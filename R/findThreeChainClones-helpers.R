#' @title Identify three TCR chain clones
#'
#' @description This internal function finds which T cell clones likely have 3
#' total TCR chains. It is a helper function for findThreeChainClones.
#'
#' @param Contigs tibble of contigs from cellranger all or
#' filtered contig files.
#' @param Chain "alpha" or "beta". This tells it whether clones with two of
#' which chain should be looked for.
#' @param Mode standard leads to modified contigs being returned.
#' If "show" is provided, two tables are returned. The first
#' shows the final set of clones with a double chain of one type, and the
#' second shows how/why some of those clones were called as such. The second
#' is pre-filtering, so some clones may not have been called as having 3 total
#' chains in the end.
#'
#' @return A list of first the clones that have 3 TCR chains, then the
#' clonotypes of select cells, then all bitypic chains. If show is provided
#' to Mode, then a list of two tables are returned showing the clones with 3
#' TCR chains.
.identifyThreeChainClones <- function(Contigs,
                                        Chain,
                                        Mode){

    if (sum(c("barcode", "high_confidence", "chain", "cdr3", "umis") %in%
        colnames(Contigs)) != 5) {
        stop(
            "Contigs should be in the all_contig_annotations or ",
            "filtered_contig_annotations format. ",
            "Please try again."
        )
    }

    if (Chain %in% c("a", "b") == FALSE) {
        stop(
            "Chain should be a or b. ",
            "Please try again."
        )
    }

    if (Mode %in% c("standard", "show") == FALSE) {
        stop(
            "Mode must be standard or show ",
            "Please try again."
        )
    }

       ### Uses scRepertoire to call clones at baseline
        tcrs <- combineTCR(Contigs)[[1]]

        if (Chain == "a") {
            tcrs <- tcrs %>% separate(
                col = "CTaa",
                into = c("bitypic", "monotypic"),
                sep = "_",
                remove = FALSE)

        } else if (Chain == "b") {

            tcrs <- tcrs %>% separate(
                col = "CTaa",
                into = c("monotypic", "bitypic"),
                sep = "_",
                remove = FALSE)

        }

        tcrs <- tcrs %>%
            relocate(all_of("monotypic"), .after = all_of("CTaa"))

        tcrs <- tcrs %>% filter(.data$monotypic != "NA" &
            .data$bitypic != "NA")

        ### Makes data frame called combos that has all cells with 3+ chains
        ### and lists the different bitypic chains with each monotypic
        combos_orig <- tcrs %>%
            filter(str_count(.data$bitypic, ";") < 2 &
                !str_detect(.data$monotypic, ";"))
        combos_full <- combos_orig %>%
            group_by(.data$monotypic, .data$bitypic) %>%
            summarise(count = n(), .groups = "drop") %>%
            ungroup() %>%
            group_by(.data$monotypic) %>%
            filter(!((.data$count / sum(.data$count)) < 0.1 &
                str_detect(.data$bitypic, ";"))) %>%
            ungroup()
        combos <- combos_full %>%
            group_by(.data$monotypic) %>%
            arrange(desc(.data$count)) %>%
            summarise(
                bitypics = paste(.data$bitypic, collapse = " "),
                counts = paste(.data$count, collapse = " "),
                total = sum(.data$count),
                chains = n()
            )
        if (nrow(combos) > 0) {
            combos <- combos %>%
                filter(.data$chains >= 3) %>%
                select(all_of(c(
                    "chains",
                    "counts",
                    "total",
                    "monotypic",
                    "bitypics"))
                ) %>%
                separate(.data$bitypics,
                    into = paste("B", seq(1, max(combos$chains)), sep = ""),
                    sep = " ", remove = FALSE, fill = "right"
                )
        } else if (nrow(combos) == 0) {
            return(Contigs)
        }

        combos <- as.data.frame(combos)
        combos <- combos %>% mutate(combo = NA)

        if (nrow(combos) > 0) {
            ### Scans if there are any clones that have multiple chains
            ### and assigns all those cells to 1 clone
            ### If there is 1+ each of M-B1B2, M-B1, and
            ### M-B2 then it is called as a double. Or also 3+ of M-B1B2.
            bitypic_cols <- colnames(combos)[6 : (ncol(combos) - 1)]


            for (i in seq_len(nrow(combos))) {

                ### Gets all chains in row
                chains <- as.character(unlist(combos[i, bitypic_cols],
                  use.names = FALSE))

                ### Gets all chains that aren't missing or combo
                valid_idx <- which(!is.na(chains) & !str_detect(chains, ";"))

                if (length(valid_idx) < 2) {
                    next
                }

                ### Gets all posible chain combos
                pairs <- do.call(
                    rbind,
                    lapply(seq_len(length(valid_idx) - 1), function(a) {
                        cbind(
                            rep(valid_idx[a], length(valid_idx) - a),
                            valid_idx[(a + 1):length(valid_idx)]
                        )
                    })
                )

                c1 <- chains[pairs[, 1]]
                c2 <- chains[pairs[, 2]]

                candidates <- as.vector(rbind(
                    paste(c1, c2, sep = ";"),
                    paste(c2, c1, sep = ";")
                ))

                ### Sees which combos are observed
                matches <- which(candidates %in% chains)

                ### Takes the first combo (i.e. the most common)
                if (length(matches) > 0) {
                    combos[i, "combo"] <- candidates[matches[1]]
                }
            }
        }

        if (Mode == "show") {
            combo_table <- combos
        }
        combos <- combos %>% select(all_of(c("monotypic", "combo")))
        combos <- combos %>% filter(!is.na(.data$combo))

        ### Uses combos_full to find more clones that may not have all 3
        ### conditions met,
        ### but still have cells with 2 bitypics that are >= 10% and
        ### also 3+ in count
        combos2 <- combos_full %>%
            filter(str_detect(.data$bitypic, ";") & .data$count >= 3)
        combos2 <- combos2 %>%
            arrange(desc(.data$count)) %>%
            group_by(.data$monotypic) %>%
            filter(row_number() == 1) %>%
            ungroup()
        combos2 <- combos2 %>%
            select(all_of("monotypic"), "combo" = all_of("bitypic")) %>%
            filter(!.data$monotypic %in% combos$monotypic)

        combos <- rbind(combos, combos2)


        ### This finds the monotypic chain that is most common
        ### for each bitypic chain.
        ### Note: this is not just 3 TCR clones, this is for all clones
        ### We want to be able to filter so that for each 2B 1M clone,
        ### each B must be most commonly seen with that M.
        pairings <- combos_orig %>%
            group_by(.data$bitypic, .data$monotypic) %>%
            summarise(count = n(), .groups = "drop")
        pairings <- pairings %>% separate(.data$bitypic,
            into = c("b1", "b2"),
            sep = ";", fill = "right"
        )
        pairings <- pairings %>%
            pivot_longer(
                cols = all_of(c("b1", "b2")),
                names_to = "original", values_to = "bitypic"
            ) %>%
            select(-all_of("original")) %>%
            filter(!is.na(.data$bitypic))
        pairings <- pairings %>%
            group_by(.data$bitypic, .data$monotypic) %>%
            summarise(
                Count = sum(.data$count),
                .groups = "drop"
            )
        pairings <- pairings %>%
            arrange(desc(.data$Count)) %>%
            group_by(.data$bitypic) %>%
            mutate(rank = min_rank(desc(.data$Count))) %>%
            group_by(.data$bitypic, .data$rank) %>%
            filter(sum(.data$rank) == 1) %>%
            ungroup() %>%
            select(-all_of("rank"))

        ### This part compares the most common monotypic for each bitypic and
        ### sees if for every 3 TCR clone both bitypics are
        ### most commonly with the monotypic
        ### If this isn't the case, they are filtered out here
        pairings2 <- pairings %>% left_join(combos, by = "monotypic")
        pairings2 <- pairings2 %>%
            filter(str_detect(
                .data$combo,
                paste("^", .data$bitypic, ";|;", .data$bitypic, "$", sep = "")
            ))
        pairings2 <- pairings2 %>%
            group_by(.data$monotypic) %>%
            mutate(combos = n_distinct(.data$combo)) %>%
            group_by(.data$monotypic, .data$combo) %>%
            summarise(
                count = n(), combos = first(.data$combos),
                .groups = "drop"
            ) %>%
            filter(.data$count == 2)

        if (nrow(pairings2) > 0) {
            if (max(pairings2$combos) > 1) {
                warning("Muliple bitypics with a monotypic remain")
            }
        }
        combos <- pairings2 %>%
            select(all_of(c("monotypic", "combo"))) %>%
            ungroup()

        ### Makes a vector special that has all of the
        ### bitypic CDR3s of chains that have 2 bitypics
        special <- as.character((combos %>%
            filter(!is.na(.data$combo)) %>%
            select(all_of("combo")))[[1]])

        if (Mode == "show") {

            combo_table <- combo_table %>%
                select(-all_of(c("chains", "bitypics"))) %>%
                rename("cell_counts_per_bitypic_pairing" = all_of("counts"),
                    "total_cells" = all_of("total"),
                    "predicted_bitypic_pairing" = all_of("combo"),
                    "monotypic_chain" = all_of("monotypic"))

            combo_table <- combo_table %>%
                select("total_cells",
                    "cell_counts_per_bitypic_pairing",
                    "monotypic_chain",
                    "predicted_bitypic_pairing",
                    everything()) %>%
                arrange(desc(.data$total_cells))

            colnames(combo_table)[5:ncol(combo_table)] <-
                str_replace(colnames(combo_table)[5:ncol(combo_table)],
                    "A", "observed_biypic_pairing")

            if (Chain == "a") {
                message("Returning a list. The first table shows all ",
                    "beta chains predicted to have a double alpha ",
                    "pairing. The second table shows how/why many ",
                    "of these clones were called as such. Please note,",
                    " the second table is pre-filtering, so not all ",
                    "clones may been in the first table since some may ",
                    "have been filtered out. In this context, the beta ",
                    "chain is the monotypic chain, and the alpha is the ",
                    "bitypic chain."
                    )
            } else if (Chain == "b"){
                message("Returning a list. The first table shows all ",
                    "alpha chains predicted to have a double beta ",
                    "pairing. The second table shows how/why many ",
                    "of these clones were called as such. Please note,",
                    " the second table is pre-filtering, so not all ",
                    "clones may been in the first table since some may ",
                    "have been filtered out. In this context, the alpha ",
                    "chain is the monotypic chain, and the beta is the ",
                    "bitypic chain."
                    )
            }
            colnames(combos) <- c("monotypic_chain",
                "predicted_bitypic_pairing")

            return(list(combos, combo_table))

        }

        return(list(combos, tcrs, special))
}







#' @title Modify three TCR chain contigs
#'
#' @description This internal function takes the clones that have 3 TCR chains
#' and modifies their contigs so that these clones can be retained as having
#' 3 TCR chains downstream. It is a helper function for findThreeChainClones.
#'
#' @param Contigs tibble of TCR contigs.
#' @param three_chain_clones a tibble of the three chain clones.
#' @param Clones the clonotypes of select cells made by
#' .identifyThreeChainClones.
#' @param Chain "alpha" or "beta". This tells it whether clones with two of
#' which chain should be looked for.
#' @param Mode modified contigs are returned with standard.
#' If "show" is provided, two tables are returned. The first
#' shows the final set of clones with a double chain of one type, and the
#' second shows how/why some of those clones were called as such. The second
#' is pre-filtering, so some clones may not have been called as having 3 total
#' chains in the end.
#'
#' @return a tibble of modified TCR contigs.
.modifyThreeChainContigs <- function(Contigs,
                                        three_chain_clones,
                                        Clones,
                                        Chain,
                                        Mode) {

    if (sum(c("barcode", "high_confidence", "chain", "cdr3", "umis") %in%
        colnames(Contigs)) != 5) {
        stop(
            "Contigs should be in the all_contig_annotations or ",
            "filtered_contig_annotations format. ",
            "Please try again."
        )
    }

    if (Chain %in% c("a", "b") == FALSE) {
        stop(
            "Chain should be a or b. ",
            "Please try again."
        )
    }

    if (Mode %in% c("standard", "show") == FALSE) {
        stop(
            "Mode must be standard or show ",
            "Please try again."
        )
    }

    if (ncol(three_chain_clones) != 2 |
        sum(colnames(three_chain_clones) %in% c("combo", "monotypic")) !=
        2) {
        stop(
             "three_chain_clones should have 2 columns named combo and ",
             "monotypic. It should be generated by ",
             ".identifyThreeChainClones. Please try again."
        )
    }

    if (ncol(Clones) != 13 |
        sum(colnames(Clones) %in% c("monotypic", "bitypic", "CTaa", "barcode",
            "CTnt")) != 5) {
      stop("Clones should have 13 columns and it should be generated by ",
           ".identifyThreeChainClones. Please try again.")
    }

        if (Mode == "show") {
            return()
        }

        ### Makes a data frame tcr_mod where combos
        ### (monotypic CDR3 + combo bitypic CDR3s) are
        ### joined to tcrs based off of shared monotypic CDR3
        ### This allows us to know which barcodes need to change.
        ### Then, if the current bitypic CDR3 is detected
        ### within the combo bitypic CDR3
        ### that has the same monotypic CDR3 then it will be assigned
        ### as then new bitypic CDR3
        ### Lastly, only the barcode and
        ### combo bitypic CDR3 are retained and cells
        ### that aren't getting assigned with a combo bitypic CDR3 are removed
        tcr_mod <- Clones %>% left_join(three_chain_clones, by = "monotypic")

        tcr_mod <- tcr_mod %>%
            mutate(bitypic_mod = case_when(
                str_detect(
                    .data$combo,
                    paste("^", .data$bitypic, ";|;",
                        .data$bitypic, "$",
                        sep = ""
                    )
                ) |
                    .data$combo == .data$bitypic ~ .data$combo,
                TRUE ~ NA
            ))

        tcr_mod <- tcr_mod %>%
            select(all_of(c("barcode", "bitypic_mod"))) %>%
            filter(!is.na(.data$bitypic_mod))

        ### A df called new_contigs is made that takes the
        ### original contigs data and left joins in tcr_mod by shared barcode
        ### Then, if the row is an bitypic chain and
        ### has a new bitypic combo it is called as being a combo
        new_contigs <- Contigs %>% left_join(tcr_mod, by = "barcode")

        if (Chain == "a"){
            new_contigs <- new_contigs %>%
                mutate(new_cdr3 = case_when(
                    !is.na(.data$bitypic_mod) & .data$chain == "TRA" ~
                        .data$bitypic_mod
                ))
        } else if (Chain == "b") {
            new_contigs <- new_contigs %>%
                mutate(new_cdr3 = case_when(
                    !is.na(.data$bitypic_mod) & .data$chain == "TRB" ~
                        .data$bitypic_mod
                ))
        }

        ### new_contigs is then filtered by reads/UMI and then,
        ### for newly assigned clones only,
        ### only the chain with the most amount of UMIs/reads for a
        ### given chain/barcode pair is retained.
        ### This filtering is necessary because essentially each row = 1 chain,
        ### so if you assign 1 row to a combo (2 chains),
        ### then you would be duplicating the data.
        ### Removing duplicates rows ensures that downstream the
        ### clone will not get duplicated.
        new_contigs <- new_contigs %>%
            group_by(.data$barcode, .data$chain) %>%
            arrange(desc(.data$umis), desc(.data$reads)) %>%
            mutate(row_num = row_number()) %>%
            ungroup()
        new_contigs <- new_contigs %>%
            filter(!(!is.na(.data$new_cdr3) & .data$row_num > 1))
        new_contigs <- new_contigs %>%
            mutate(cdr3 = case_when(
                !is.na(.data$new_cdr3) ~ .data$new_cdr3,
                is.na(.data$new_cdr3) ~ .data$cdr3
            ))

        return(new_contigs)

}








#' @title Find missing monotypic chains
#'
#' @description This internal function finds cells that may be in clones with
#' 3 TCR chains and have no monotypic chain detected. If likely, it assigns
#' them to these clones. It is a helper function for findThreeChainClones.
#'
#' @param Contigs modified TCR contigs from .modifyThreeChainContigs.
#' @param specialTCRs a vector of bitypic chains.
#' @param Chain "alpha" or "beta". This tells it whether clones with two of
#' which chain should be looked for.
#' @param Mode modified contigs are returned with standard.
#' If "show" is provided, two tables are returned. The first
#' shows the final set of clones with a double chain of one type, and the
#' second shows how/why some of those clones were called as such. The second
#' is pre-filtering, so some clones may not have been called as having 3 total
#' chains in the end.
#' @param print_stats specifies whether notes/stats
#' should be printed.
#'
#' @return a tibble of modified TCR contigs.
.findMissingMonotypicChains <- function(Contigs,
                                        specialTCRs,
                                        Chain,
                                        Mode,
                                        print_stats) {

    if (sum(c("barcode", "high_confidence", "chain", "cdr3", "umis") %in%
        colnames(Contigs)) != 5) {
        stop(
            "Contigs should be in the all_contig_annotations or ",
            "filtered_contig_annotations format. ",
            "Please try again."
        )
    }

    if (Chain %in% c("a", "b") == FALSE) {
        stop(
            "Chain should be a or b. ",
            "Please try again."
        )
    }

    if (Mode %in% c("standard", "show") == FALSE) {
        stop(
            "Mode must be standard or show ",
            "Please try again."
        )
    }

    if (!is.character(specialTCRs)) {
      stop(
          "specialTCRs must be a character vector. ",
          "Please try again."
      )
    }

    if (!is.logical(print_stats)) {
        stop(
            "printStats must be TRUE or FALSE. ",
            "Please try again."
        )
    }

        if (Mode == "show") {
            return()
        }

        special <- specialTCRs
        ab <- Contigs
        post <- combineTCR(ab)[[1]]

        ### THIS SECOND PART IS FINDING ALL THE CELLS WITH AN BITYPIC AND
        ### A MISSING MONOTYPIC THAT ARE IN CLONES WITH COMBO BITYPIC CHAINS
        ### THEN THESE CELLS ARE ASSIGNED TO THEIR COMBO CLONE SO
        ### THEY ARE CALLED DOWNSTREAM AS 1 CLONE

        if (length(special) > 0) {
            ### Makes a df vdj which summarizes the new clone sizes
            ### Then, if the bitypic isn't missing,
            ### calls the number of combo CDR3s it could be from.
            vdj <- post %>%
                group_by(.data$CTaa) %>%
                summarise(count = n()) %>%
                arrange(desc(.data$count)) %>%
                ungroup()

            if (Chain == "a") {
                vdj <- vdj %>% separate(
                    col = "CTaa",
                    into = c("bitypic", "monotypic"),
                    sep = "_",
                    remove = FALSE)

            } else if (Chain == "b") {

                vdj <- vdj %>% separate(
                    col = "CTaa",
                    into = c("monotypic", "bitypic"),
                    sep = "_",
                    remove = FALSE)

            }

            vdj <- vdj %>%
                relocate(all_of("monotypic"), .after = all_of("CTaa"))

            vdj <- vdj %>%
                rowwise() %>%
                mutate(special_count = case_when(
                    .data$bitypic != "NA"
                    ~ sum(str_detect(
                            special,
                            paste("^", .data$bitypic, ";|;",
                                .data$bitypic, "$",
                                sep = ""
                            )
                        ) | special == .data$bitypic),
                    .data$bitypic == "NA" ~ 0
                ))


            ### Makes a new column called special_tcr that has the combo
            ### bitypic that the bitypic chain belongs to.
            vdj$special_tcr <- NA
            for (i in seq_len(nrow(vdj))) {
                if (vdj$special_count[i] == 1) {
                    vdj$special_tcr[i] <-
                        special[which(str_detect(
                            special,
                            paste("^", vdj$bitypic[i], ";|;",
                                vdj$bitypic[i], "$",
                                sep = ""
                            )
                        ) | special == vdj$bitypic[i])]
                }
            }

            vdj_bad <- vdj %>% filter(.data$special_count > 1 &
                .data$monotypic == "NA")
            if (nrow(vdj_bad) > 0) {
                message(
                    "There are some bitypic chains that are in ",
                    "multiple combo clones. ",
                    "Because of this, ",
                    sum(vdj_bad$count),
                    " cells with these bitypic chains and no monotypic",
                    "won't be able to be assinged to a combo clone"
                )
            }

            ### Makes a df call special_cells that calls the cells with
            ### monotypic NA and that belong to combo bitypic clones
            ### Begins by only retaining cells that could be called
            ### and then calculates the percentage of cells with each
            ### monotypic and the same combo bitypic.
            ### You can think of these percentages as
            ### the percentage likelihood of
            ### each monotypic for cells with a NA monotypic and the same
            ### bitypic
            special_cells <- vdj %>%
                filter(!is.na(.data$special_tcr)) %>%
                mutate(mod_count = case_when(
                    .data$monotypic == "NA" ~ 0,
                    .data$monotypic != "NA" ~ .data$count
                ))
            special_cells <- special_cells %>%
                group_by(.data$special_tcr) %>%
                mutate(group_count = sum(.data$mod_count)) %>%
                ungroup()
            special_cells <- special_cells %>%
                mutate(group_pct = 100 * .data$mod_count / .data$group_count)

            ### Next we try to assign for NA monotypic cells.
            ### This starts with making a column new_bitypic which
            ### is only there for rows where the monotypic is 80%+ likely.
            ### We then add that same bitypic (that pairs with the likely
            ### monotypic) to a new column new_bitypic2 for all cells
            ### in the same bitypic group.
            ### Lastly, makes a new column new_bitypic3 that,
            ### for cells with no montoypic and a paired bitypic,
            ### has the combo bitypic chain CDR3
            special_cells <- special_cells %>%
                mutate(new_bitypic = case_when(
                    .data$group_pct >= 80 ~ .data$bitypic,
                    .data$group_pct < 80 ~ NA
                ))
            special_cells <- special_cells %>%
                arrange(desc(.data$group_pct)) %>%
                group_by(.data$special_tcr) %>%
                mutate(new_bitypic2 = dplyr::first(.data$new_bitypic)) %>%
                ungroup()
            special_cells <- special_cells %>%
                rowwise() %>%
                mutate(new_bitypic3 = case_when(
                    is.na(.data$new_bitypic2)
                    ~ NA,
                    .data$monotypic == "NA" & (str_detect(
                        .data$new_bitypic2,
                        paste("^", .data$bitypic, ";", "|", ";",
                            .data$bitypic, "$",
                            sep = ""
                        )
                    ) |
                        .data$new_bitypic2 == .data$bitypic) ~
                            .data$new_bitypic2,
                    TRUE ~ NA
                )) %>%
                ungroup()

            ### Makes a df special_cells2 that selects only the CDR3s
            ### and the new bitypic combo CDR3 and then
            ### filters rows with no new assignment
            special_cells2 <- special_cells %>%
                dplyr::select(all_of(c("CTaa", "new_bitypic3"))) %>%
                filter(!is.na(.data$new_bitypic3))

            ### Joins special_cells with post (the scRepertoire output with
            ### combo bitypic chains) by CDR3s to make vdj2
            ### This allow us to know which barcodes need to change.
            ### Barcodes that don't need to change are then removed
            vdj2 <- post %>% left_join(special_cells2, by = "CTaa")
            vdj2 <- vdj2 %>%
                select(all_of(c("barcode", "new_bitypic3"))) %>%
                filter(!is.na(.data$new_bitypic3))

            ### ab (contigs after joining bitypic chains) is
            ### joined to vdj2 to make ab2
            ### ab2 is filtered to retain only the row of that barcode/chain
            ### combo with the most UMIs/reads.
            ### This is only for bitypic chains that are changing
            ### This is necessary because, again, not doing so would result
            ### in some duplication down stream
            ### Lastly, the bitypic cdr3 for these combo bitypic - missing
            ### monotypic cells are changed
            ab2 <- ab %>% left_join(vdj2, by = "barcode")
            ab2 <- ab2 %>%
                group_by(.data$barcode, .data$chain) %>%
                arrange(desc(.data$umis), desc(.data$reads)) %>%
                mutate(row_num = row_number()) %>%
                ungroup()

            if (Chain == "a") {

                ab2 <- ab2 %>%
                    filter(!(!is.na(.data$new_bitypic3) & .data$row_num > 1 &
                        .data$chain == "TRA")) %>%
                    mutate(cdr3 = case_when(
                        .data$chain == "TRA" & !is.na(.data$new_bitypic3) ~
                            .data$new_bitypic3,
                        TRUE ~ .data$cdr3
                    ))
            } else if (Chain == "b") {

                ab2 <- ab2 %>%
                    filter(!(!is.na(.data$new_bitypic3) & .data$row_num > 1 &
                        .data$chain == "TRB")) %>%
                    mutate(cdr3 = case_when(
                        .data$chain == "TRB" & !is.na(.data$new_bitypic3) ~
                            .data$new_bitypic3,
                        TRUE ~ .data$cdr3
                    ))
            }

            ab2 <- ab2 %>% select(-all_of("new_bitypic3"))

        } else if (length(special) == 0) {
            ab2 <- ab
        }

        if (print_stats == TRUE) {

            if (Chain == "a") {

                message(
                    length(special), " clones with an extra TCR\u03B1 chain ",
                    "were detected, comprising ",
                    round(sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRA")
                    / n_distinct(ab2$barcode) * 100, 1), "% of cells (",
                    sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRA"),
                    " total cells)"
                    )
            } else if (Chain == "b") {
                message(
                    length(special), " clones with an extra TCR\u03B2 chain ",
                    "were detected, comprising ",
                    round(sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRB")
                    / n_distinct(ab2$barcode) * 100, 1), "% of cells (",
                    sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRB"),
                    " total cells)"
                )
            }
        }

        ab2 <- ab2 %>% select(-all_of(c(
            "bitypic_mod",
            "new_cdr3",
            "row_num"))
        )

        return(ab2)



}
