#' @title Read TCR contigs
#'
#' @description This internal function reads in TCR contigs and joins them.
#' It is a helper function for ECLIPSE.
#'
#' @param Contig_folder_paths vector of paths to directories that contain TCR
#' files. If you use this argument, provide either "all" or "filtered" to
#' contig_type. The order should match the numbers provided in batch.
#' @param Contig_type tells whether contig_folder_paths or contig_file_paths
#' should be used to read in contig files and how. To provide file paths,
#' this should be "manual". To provide folders, this should be "all" or
#' "filtered" to search for all or filtered contig annotations files.
#' @param Contig_file_paths vector of TCR file paths. The order should match
#' the numbers provided in batch. The files can have any name, but they must
#' be all_contig_annotations.csv or
#' filtered_contig_annotations.csv format.
#' @param Verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline.
#'
#' @return A data.frame of all contigs
.readContigs <- function(Contig_type,
                            Contig_folder_paths,
                            Contig_file_paths,
                            Verbose) {

    #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (!Contig_type %in% c(
        "ALL", "All", "all",
        "Filtered", "FILTERED", "filtered",
        "all_contig_annotations.csv",
        "filtered_contig_annotations.csv",
        "manual"
    )) {
        stop(
            "Incorrect contig type. Should be all or filtered. ",
            "Returning nothing"
        )
    }

    if (is.null(Contig_folder_paths) & Contig_type != "manual") {
        stop(
            "Please provide folders to contig_folder_paths that have the ",
            "contig annotations in them."
        )
    }

    if (is.null(Contig_file_paths) & Contig_type == "manual") {
        stop(
            "No contig file paths provided. ",
            "Please provide a vector of file path names to contig_file_paths"
        )
    }

    if (Contig_type != "manual" &
        (!is.vector(Contig_folder_paths) | is.list(Contig_folder_paths) |
        length(Contig_folder_paths) == 0)) {
        stop("contig_folder_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (Contig_type == "manual" &
        (!is.vector(Contig_file_paths) | is.list(Contig_file_paths) |
        length(Contig_file_paths) == 0)) {
        stop("contig_file_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (!is.null(Contig_file_paths) & !is.null(Contig_folder_paths)) {
        stop("Please do not use contig_file_paths and contig_folder_paths. ",
            "Please try again.")
    }

    if (!is.logical(Verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    ### Reading the data in
    if (Verbose == TRUE) message("\nReading TCR contig data in...\n")

    if (Contig_type != "manual") {
        contigs_list <- vector("list", length(Contig_folder_paths))

        for (i in seq_along(Contig_folder_paths)) {
            if (Contig_type %in% c(
                "ALL", "All", "all",
                "all_contig_annotations.csv"
            )) {
                files <- list.files(Contig_folder_paths[i],
                    pattern = "all_contig_annotations.csv",
                    full.names = TRUE
                )
            } else if (Contig_type %in%
                c(
                    "Filtered", "filtered", "FILTERED",
                    "filtered_contig_annotations.csv"
                )) {
                files <- list.files(Contig_folder_paths[i],
                    pattern = "filtered_contig_annotations.csv",
                    full.names = TRUE
                )
            }


            if (length(files) == 0) {
                stop("No contig files found, returning nothing")
            } else if (length(files) > 1) {
                stop("Multiple contig files found, returning nothing")
            } else if (length(files) == 1) {
                contigs_list[[i]] <- read.csv(files[1])
                contigs_list[[i]]$barcode <- paste(i, "_",
                    contigs_list[[i]]$barcode,
                    sep = ""
                )
            }
        }
    } else if (Contig_type == "manual") {
        contigs_list <- vector("list", length(Contig_file_paths))

        for (i in seq_along(Contig_file_paths)) {
            contigs_list[[i]] <- read.csv(Contig_file_paths[i])

            if (nrow(contigs_list[[i]]) == 0) {
                stop(
                    "File path ", i, " did not lead to any ",
                    "contigs being loaded in. ",
                    "Check the path name and that it is a .csv file"
                )
            }

            contigs_list[[i]]$barcode <-
                paste(i, "_", contigs_list[[i]]$barcode, sep = "")
        }
    }

    ### Joins together all batches and adds data to them
    ### about which sample/donor/batch each contig was from
    contigs_full <- bind_rows(contigs_list)

    return(contigs_full)
}








#' @title QC filter contigs
#'
#' @description This internal function filters contigs that do not pass QC.
#' It is a helper function for ECLIPSE.
#'
#' @param Contigs the TCR contigs data.frame made by .readContigs.
#' @param Seurat_object Seurat object that contains scRNA-seq data.
#' @param Batch name of column in Seurat meta.data. The column should contain
#' integers, not characters, that tell which position in the
#' contig_folder/file_path vectors contains the TCR contigs for the cell.
#' @param Donor name of column in Seurat meta.data. The column should contain
#' the donor, regardless of treatment/timepoint/condition.
#' @param Group name of column in Seurat meta.data. This column should specify
#' how the final statistics for each clone should be calculated. For example,
#' should they be calculated per sample, per donor, or among all cells?
#' This can and often will be the same column as donor.
#' @param Original_barcode name of column in Seurat meta.data. This must
#' contain the original 10X barcode for the cell. Typically this will end with
#' -1, but users should check their TCR contig files to be sure.
#' @param Format This tells how missing data is stored in
#' the TCR contig files. Current files likely contain "", but if older contig
#' files are used where any missing data is stored
#' as "None", then "None" should be provided here.
#' @param ReturnFilteredContigs specifies if filtered contigs should be returned
#' instead of running the whole ECLIPSE pipeline.
#' @param Verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline.
#'
#' @return A list of 3, with the filtered contigs first, the seurat object,
#' then contigs that are only seen in one cell aka orphan contigs.
#' Alternatively, if ReturnFilteredContigs == TRUE, then just the contigs are
#' returned.
.qcFilterContigs <- function(Contigs,
                                Seurat_object,
                                Batch,
                                Donor,
                                Group,
                                Original_barcode,
                                Format,
                                ReturnFilteredContigs,
                                Verbose) {

    contigs_full <- Contigs
    seurat_object <- Seurat_object
    batch <- Batch
    donor <- Donor
    group <- Group
    original_barcode <- Original_barcode
    format <- Format
    returnFilteredContigs <- ReturnFilteredContigs
    verbose <- Verbose

    rm(Seurat_object)

    ### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (sum(colnames(contigs_full) %in% c("barcode", "is_cell",
        "cdr3", "chain", "high_confidence")) != 5) {
        stop("Contigs should be in a contigs format and generated by ",
            "running .readContigs. Please try again."
        )
    }

    if (!batch %in% colnames(seurat_object@meta.data)) {
        stop(
            "Batch column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!donor %in% colnames(seurat_object@meta.data)) {
        stop(
            "Donor column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!group %in% colnames(seurat_object@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!original_barcode %in% colnames(seurat_object@meta.data)) {
        stop(
            "Original barcode column name provided is not ",
            "in the Seurat object. ",
            "Please try again."
        )
    }

    if (format %in% c("None", "none", "Blank") == FALSE) {
        stop(
            "None, none, or Blank not provided to format. ",
            "Please try again."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    if (!is.logical(returnFilteredContigs)) {
        stop(
            "TRUE or FALSE not provided to returnFilteredContigs. ",
            "Please try again."
        )
    }


    if ("sample" %in% colnames(contigs_full)) {
        contigs_full <- contigs_full %>% select(-all_of("sample"))
    }

    if (format %in% c("None", "none")) {
        contigs_full[contigs_full == "None"] <- ""
    }

    seurat_object$tcrEclipse_barcode <-
        paste(seurat_object@meta.data[[batch]],
            "_",
            seurat_object@meta.data[[original_barcode]],
            sep = ""
        )
    seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)

    cell_data <- seurat_object@meta.data %>%
        select(all_of(c(
            "tcrEclipse_barcode",
            unique(c(batch, donor, group)),
            "seurat_object_barcode"))
        )
    cell_data <- cell_data %>% filter(!is.na(.data[[batch]]))

    if (nrow(cell_data) != n_distinct(cell_data$tcrEclipse_barcode)) {
        stop(
            "There are likely multiple sequencing runs ",
            "labeled as one batch. ",
            "Please relabel the batch column such that each batch ",
            "is a unique sequncing run. ",
            "Returning nothing"
        )
    }

    all <- contigs_full %>%
        left_join(cell_data,
            by = c("barcode" = "tcrEclipse_barcode"),
            keep = FALSE
        )
    all <- all %>%
        mutate(barcode = .data$seurat_object_barcode) %>%
        select(-all_of("seurat_object_barcode"))
    if ("barcode" %in% colnames(seurat_object@meta.data)) {
        seurat_object@meta.data <- seurat_object@meta.data %>%
            mutate(originalBarcodeBeforeECLIPSE = .data$barcode) %>%
            select(-all_of("barcode"))
    }

    seurat_object@meta.data <- seurat_object@meta.data %>%
        dplyr::rename("barcode" = all_of("seurat_object_barcode"))

    if (verbose == TRUE) {
        message("QC Filtering:")
        message(" - ", nrow(all), " initial contigs found before filtering")
    }

    ### Renaming barcodes, filtering to only cells in seurat object,
    ### removing cells without a CDR3,
    ### and renaming QC metrics from cellranger so
    ### these contigs aren't thrown out
    all <- all %>% filter(.data$barcode %in% seurat_object$barcode)
    if (verbose == TRUE) {
        message(" - ", nrow(all), " contigs were from true cells ",
            "(coming from ", n_distinct(all$barcode), " cells)")
    }
    if (nrow(all) == 0) {
        stop(
            "Check that your barcodes are the original ",
            "(i.e. likely end with -1 and match the contig files)"
        )
    }

    if (verbose == TRUE) {
        message(" - ",
            sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE),
            " contigs were removed for not being high confidence as ",
            "defined by Cell Ranger"
        )
    }
    all <- all %>%
        filter(.data$high_confidence %in% c("true", "True", "TRUE"))
    all <- all %>% mutate(
        productive = "true",
        is_cell = "true",
        full_length = "true"
    )

    count_before_custom_filtering <- nrow(all)

    all <- all %>% filter(!.data$cdr3 %in% c("", "None") &
                            !.data$cdr3_nt %in% c("", "None"))
    all <- all %>% filter(.data$chain %in% c("TRA", "TRB", "TRG", "TRD"))

    ### Removing cells with abnormally long CDR3s
    all <- all %>% mutate(cdr3_length = str_length(.data$cdr3))

    all <- all %>% filter(.data$cdr3_length <= 30)

    ### There are TCRs with TRDV and then TRAJ and TRAC.
    ### Cellranger calls them as TRD,
    ### but they are actually TRA since they pair with a beta chain and
    ### have no J segment.
    all <- all %>%
        mutate(ad_hybrid_tcr = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC"))
            ~ TRUE,
            TRUE ~ FALSE
        ))

    all <- all %>%
        mutate(chain = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC")) ~ "TRA",
            TRUE ~ .data$chain
        ))

    all <- all %>% filter(.data$chain %in% c("TRA", "TRB"))

    ### Removing cells with a * at the end of their cdr3.
    ### These indicate stop codons, so removing these contigs
    all <- all %>% filter(!str_detect(.data$cdr3, "\\*"))

    ### Removing duplicated chains.
    ### Some datasets (particularly contig files made with
    ### older version of Cellranger)
    ### have lots of duplicated rows that have identical CDR3
    ### and gene segments for the same barcode.
    ### This messes up future analysis
    all <- all %>%
        arrange(desc(.data$umis), desc(.data$reads)) %>%
        group_by(.data$barcode, .data$chain, .data$cdr3) %>%
        filter(row_number() == 1) %>%
        ungroup()

    ### Noting cells with only 1 type of chain (alpha or beta) and
    ### that chain isn't seen in any other cell.
    ### We don't want these cells to be predicted.
    all <- all %>%
        group_by(.data$barcode) %>%
        mutate(chainTypes = n_distinct(.data$chain)) %>%
        ungroup()
    all <- all %>%
        group_by(.data$cdr3, .data$chain, .data[[donor]]) %>%
        mutate(maxChains = max(.data$chainTypes)) %>%
        ungroup()
    if (verbose == TRUE) {
        message(" - ",
            sum(all$maxChains == 1),
            " orphan contigs were found with no observed pairing. ",
            "Cells containing these contigs are unlikely to have ",
            "a clonotype predicted"
        )
    }
    orphanContigs <- all %>%
        filter(.data$maxChains == 1) %>%
        select(all_of(c("barcode", "chain", "cdr3")))

    if (verbose == TRUE) {
        message(" - ", count_before_custom_filtering - nrow(all),
            " contigs were removed for other reasons ",
            "(ex. irregular CDR3, TCR\u03B3/\u03B4 contig, duplicated, etc...)"
        )
        message(" - ",
            "After all filtering, ", nrow(all),
            " contigs from ", n_distinct(all$barcode),
            " true cells remained"
        )
    }

    if (returnFilteredContigs == TRUE) {

        all$orphanContig <- FALSE
        all$orphanContig[all$maxChains == 1] <- TRUE

        all <- all %>%
            select(-all_of(c(
                "cdr3_length",
                "ad_hybrid_tcr",
                "chainTypes",
                "maxChains"))
                )
        return(all)
    }


    return(list(all, seurat_object, orphanContigs))
}









#' @title Run VDJdive
#'
#' @description This internal function runs the EM algorithm of VDJdive and
#' extracts the results from the assignment matrix. It is a helper function
#' for ECLIPSE.
#'
#' @param Contigs TCR contigs made from .qcFilterContigs and
#' findThreeChainClones.
#' @param Donor name of column in Seurat meta.data. The column should contain
#' the donor, regardless of treatment/timepoint/condition.
#' @param Verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline.
#'
#' @return A data.frame describing the results from the EM algorithm.
.runVDJdive <- function(Contigs,
                        Verbose,
                        Donor) {

    donor <- Donor

    if (sum(colnames(Contigs) %in% c("barcode", "is_cell", "cdr3", "chain",
                                     "high_confidence")) != 5) {
        stop("Contigs should be in a contigs format and generated by ",
            "running .readContigs. Please try again."
        )
    }

    if (!is.logical(Verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )

      }
    ### THIS PART IS FEEDING THE EDITED CONTIG FILES INTO VDJdive
    ### WHICH CALLS AMBIGUOUS CELLS

    if (Verbose == TRUE) {
        message("Running VDJdive::clonoStats() for clone EM assignment...")
    }
    ### Making a temporary directory
    ### Writing a .csv file as filtered_contig_annotations.csv
    ### Reading that into VDJdive
    ### Using clonoStats for EM assignment
    ### Extracting output and changing format

    temp_dir <- tempfile("ECLIPSE_")
    dir.create(temp_dir)

    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    write.csv(Contigs,
        paste(temp_dir, "/filtered_contig_annotations.csv", sep = ""),
        row.names = FALSE
    )

    contigs <- readVDJcontigs(temp_dir)

    ### Please note suppressMessages is used here because
    ### there is a message stating internal coercion from
    ### Tsparse* to CsparseMatrix
    ### that I do not want to be printed.
    vdj <- suppressMessages(clonoStats(contigs,
        method = "EM",
        type = "TCR",
        assignment = TRUE,
        group = donor
    ))

    if (Verbose == TRUE) message("Extracting results from VDJdive...")
    df <- vdj@assignment
    df <- as(df, "RsparseMatrix")
    colnames(df) <- clonoNames(vdj)

    rm(vdj)

    ### Finding the most likely clonotype, its probability,
    ### and the probability of the 2nd most likely clonotype for each cell.
    p <- df@p
    j <- df@j
    x <- df@x

    all_clones <- colnames(df)
    n <- nrow(df)

    barcode <- rownames(df)
    clone <- character(n)
    maxes <- numeric(n)
    second_max <- numeric(n)


    for (i in seq_len(n)) {

        idx <- (p[i] + 1):p[i + 1]

        vals <- x[idx]
        cols <- j[idx] + 1

        best_idx <- which.max(vals)

        maxes[i] <- vals[best_idx]
        second_max[i] <- if (length(vals) == 1) {
            0
        } else {
            max(0, vals[-best_idx])
        }

        clone[i] <- all_clones[cols[best_idx]]
    }

    em <- data.frame(barcode = rownames(df), maxes, clone, second_max)

    return(em)

}








#' @title Assign final clonotypes
#'
#' @description This internal function determines and adds the final EM-based
#' clonotype to the Seurat object. It also adds the chains that were detected
#' and how the cell would have been called with standard analysis. It is a
#' helper function for ECLIPSE.
#'
#' @param em_results the results from running .runVDJdive.
#' @param Contigs post-QC filtered contigs.
#' @param Orphan_contigs orphan TCR contigs determined by .qcFilterContigs.
#' @param Seurat_object Seurat object that contains scRNA-seq data.
#' @param Donor name of column in Seurat meta.data. The column should contain
#' the donor, regardless of treatment/timepoint/condition.
#' @param Group name of column in Seurat meta.data. This column should specify
#' how the final statistics for each clone should be calculated. For example,
#' should they be calculated per sample, per donor, or among all cells?
#' This can and often will be the same column as donor.
#' @param CombineTcrArgs arguments that are provided to scRepertoire
#' combineTCR. This should be a list of 4 logical elements,
#' with the elements corresponding to the combineTCR arguments.
#' 1 is removeNA, 2 is removeMulti, 3 is filterMulti, 4 is filterNonproductive.
#' If no clonotype can be predicted with high confidence, then combineTCR
#' using these arguments will determine the final clonotype.
#' @param Verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline.
#'
#' @return A Seurat object with the clonotype information added.
.assignFinalClonotypes <- function(em_results,
                                    Contigs,
                                    Seurat_object,
                                    Orphan_contigs,
                                    Donor,
                                    Group,
                                    CombineTcrArgs,
                                    Verbose) {
    all <- Contigs
    em <- em_results
    orphanContigs <- Orphan_contigs
    seurat_object <- Seurat_object
    donor <- Donor
    group <- Group
    combineTcrArgs <- CombineTcrArgs
    verbose <- Verbose

    if (sum(colnames(all) %in%
        c("barcode", "is_cell", "cdr3", "chain", "high_confidence")) != 5) {
      stop("Contigs should be in a contigs format and generated by ",
           "running .qcFilterContigs. Please try again."
      )
    }

    if (sum(colnames(em) %in% c("barcode", "maxes", "clone", "second_max"))
        != 4 | ncol(em) != 4) {
        stop("em_results should have 4 columns and be made by .runVDJdive. ",
            "Please try again.")
    }

    if (ncol(Orphan_contigs) != 3 | sum(colnames(Orphan_contigs) %in%
        c("barcode", "chain", "cdr3")) != 3) {
        stop("Orphan_contigs should have 3 columns and be made by ",
            ".qcFilterContigs. Please try again.")
    }

    if (!donor %in% colnames(seurat_object@meta.data)) {
        stop(
            "Donor column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!group %in% colnames(seurat_object@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.list(combineTcrArgs) | length(combineTcrArgs) != 4 |
        !is.logical(combineTcrArgs[[1]]) | !is.logical(combineTcrArgs[[2]]) |
        !is.logical(combineTcrArgs[[3]]) | !is.logical(combineTcrArgs[[4]])) {
        stop(
            "Incorrect input to combineTcrArgs. ",
            "Please provide a list of 4 logical (TRUE or FALSE) elements."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    ### Summarizing the VDJdive EM data per cell
    ### and then filtering cells that don't meet QC thresholds
    if (verbose == TRUE) {
        message("Finalizing clonal assignments...\n")
    }
    em$ratio <- em$maxes / em$second_max

    if (verbose == TRUE) {
        em2 <- em %>%
            arrange(.data$maxes) %>%
            mutate(
                nm = row_number(),
                assigned_by_VDJdive = .data$maxes >= 0.8 |
                    (.data$maxes >= 0.7 & .data$ratio > 5) |
                    (.data$maxes > 0.5 & .data$ratio > 10)
            )
        plot2 <- em2 %>%
            ggplot() +
            geom_point(aes(
                x = .data$nm, y = .data$maxes,
                color = .data$assigned_by_VDJdive
            )) +
            xlab("Ranked position of cell") +
            ylab(paste("Probability of cell being\n",
                "in most likely clone", sep = "")) +
            ggtitle("Clonal Prediction QC") +
            theme_bw() +
            labs(color = "Source\nof Final\nClonotype") +
            scale_color_manual(values = c("cyan4", "darkorange2"),
                labels = c("Traditional\nAnalysis\n",
                    "High\nConfidence\nPrediction")) +
            theme(
                plot.title = element_text(
                    face = "bold", hjust = 0.5,
                    size = 20
                ),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 16, hjust = 0.5),
                legend.text = element_text(size = 12),
                legend.position.inside = c(0.7, 0.4)
            )
        plot(plot2)
    }

    final <- em %>% filter(.data$maxes >= 0.8 |
        (.data$maxes >= 0.7 & .data$ratio > 5) |
        (.data$maxes > 0.5 & .data$ratio > 10))

    ### Removing cells where an orphan TCR was assigned.
    ### There's no way for the algorithm to accurately predict these so
    ### we don't want them.
    orphanContigs <- orphanContigs %>% left_join(final, by = "barcode")
    orphanContigs <- orphanContigs %>%
        mutate(remove = case_when(
            is.na(.data$clone) ~ FALSE,
            .data$chain == "TRA" &
                str_detect(
                    .data$clone,
                    paste("^", .data$cdr3, " ", sep = "")
                ) ~ TRUE,
            .data$chain == "TRB" &
                str_detect(
                    .data$clone,
                    paste(" ", .data$cdr3, "$", sep = "")
                ) ~ TRUE,
            TRUE ~ FALSE
        ))
    orphanContigs <- orphanContigs %>%
        filter(.data$remove) %>%
        select(all_of(c("barcode", "remove")))
    orphanContigs <- unique(orphanContigs)

    final <- final %>% left_join(orphanContigs, by = "barcode")
    final <- final %>% filter(is.na(.data$remove))

    ### Adding the VDJdive and scRepertoire calls and
    ### expressed chains to the Seurat object
    seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]

    contig_list <- createHTOContigList(all,
        seurat_object,
        group.by = "tcrEclipseGroup"
    )

    no_vdjd <- combineTCR(contig_list)

    seurat_object <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )

    seurat_object$highConfChains <- seurat_object$CTaa
    seurat_object@meta.data <- seurat_object@meta.data %>%
        select(-all_of(c(
            "CTgene",
            "CTnt",
            "CTaa",
            "CTstrict",
            "clonalProportion",
            "clonalFrequency",
            "cloneSize"))
        )

    no_vdjd <- combineTCR(contig_list,
        removeNA = combineTcrArgs[[1]],
        removeMulti = combineTcrArgs[[2]],
        filterMulti = combineTcrArgs[[3]],
        filterNonproductive = combineTcrArgs[[4]]
    )
    obj <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )
    rm(seurat_object)

    columns_to_remove <- which(colnames(obj@meta.data) %in%
        c("scR_CTgene", "scR_CTnt", "scR_CTaa", "scR_CTstrict",
            "scR_clonalFrequency", "scR_clonalProportion", "scR_cloneSize"))

    if (length(columns_to_remove) > 0) {

        obj@meta.data <- obj@meta.data[, -columns_to_remove]

    }

    obj@meta.data <- obj@meta.data %>%
        rename(
            scR_CTgene = all_of("CTgene"),
            scR_CTnt = all_of("CTnt"),
            scR_CTaa = all_of("CTaa"),
            scR_CTstrict = all_of("CTstrict"),
            scR_clonalProportion = all_of("clonalProportion"),
            scR_clonalFrequency = all_of("clonalFrequency"),
            scR_cloneSize = all_of("cloneSize")
        )
    final <- final %>% mutate(clone2 = str_replace(.data$clone, " ", "_"))
    final <- final %>% mutate(clone2 = str_replace_all(
        .data$clone2,
        ";", ":::"
    ))
    final <- final %>% select(all_of("barcode"),
        "VDJdive_clone" = all_of("clone2"),
        "EM_max_prop" = all_of("maxes"),
        "EM_2ndmax_prop" = all_of("second_max"),
        "EM_max_2ndmax_ratio" = all_of("ratio")
    )

    obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")


    ### Making CTaa our final clone call which uses VDJdive
    ### if available but if not the scRepertoire call
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneCallSource = case_when(
            !is.na(.data$VDJdive_clone) ~ "VDJdive",
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                "scRepertoire",
            TRUE ~ NA
        ))

    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$VDJdive_clone) ~ .data$VDJdive_clone,
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                .data$scR_CTaa,
            TRUE ~ NA
        ))


    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = str_replace_all(.data$CTaa, ";", ":::"))

    ### Combining clones. Cells that were ambiguous may be assigned to a
    ### combo clone by VDJdive but only be listed with one chain
    ### Here those cells have their CDR3 renamed as the combo CDR3
    obj@meta.data <- obj@meta.data %>%
        ungroup() %>%
        mutate(CTaa_donor_source = case_when(
            !is.na(.data$CTaa) &
                !is.na(.data[[donor]]) &
                !is.na(.data$cloneCallSource) ~
                paste(.data$CTaa, .data[[donor]],
                    .data$cloneCallSource,
                    sep = " "
                ),
            TRUE ~ NA
        ))

    type <- c("alpha", "beta")
    opposite_type <- c("beta", "alpha")

    for (i in 1:2) {

        custom <- obj@meta.data %>%
            filter(!is.na(.data$CTaa) & !is.na(.data[[donor]]) &
                .data$cloneCallSource == "VDJdive") %>%
            group_by(.data$CTaa, .data[[donor]], .data$cloneCallSource) %>%
            summarise(count = n()) %>%
            arrange(desc(.data$count))
        custom <- custom %>%
            separate(.data$CTaa,
                into = c("alpha", "beta"),
                remove = FALSE, sep = "_"
            )

        custom <- custom %>%
            mutate(special = case_when(
                str_detect(.data[[type[i]]], ":::") &
                    .data$cloneCallSource == "VDJdive" ~ .data[[type[i]]],
                !str_detect(.data[[type[i]]], ":::") |
                    .data$cloneCallSource != "VDJdive" ~ NA))

        custom <- custom %>%
            group_by(.data[[opposite_type[i]]], .data[[donor]]) %>%
            arrange(desc(.data$count)) %>%
            mutate(special_Chain = first(.data$special)) %>%
            ungroup()

        if (type[i] == "alpha") {

            custom <- custom %>%
                mutate(new_clone = case_when(
                    str_detect(
                        .data$special_Chain,
                        paste("^", .data[[type[i]]], ":::", "|", ":::",
                              .data[[type[i]]], "$",
                            sep = ""
                        )
                    ) &
                        .data[[type[i]]] != .data$special_Chain ~
                        paste(.data$special_Chain, .data[[opposite_type[i]]],
                            sep = "_"),
                    TRUE ~ NA
                ))
        } else if (type[i] == "beta") {
            custom <- custom %>%
                mutate(new_clone = case_when(
                    str_detect(
                        .data$special_Chain,
                        paste("^", .data[[type[i]]], ":::", "|", ":::",
                              .data[[type[i]]], "$",
                            sep = ""
                        )
                    ) &
                        .data[[type[i]]] != .data$special_Chain ~
                        paste(.data[[opposite_type[i]]], .data$special_Chain,
                            sep = "_"),
                    TRUE ~ NA
                ))
        }

        custom <- custom %>%
            unite(
                col = "CTaa_donor_source",
                all_of(c("CTaa", donor, "cloneCallSource")),
                sep = " ", remove = FALSE
            )
        custom <- custom %>%
            filter(!is.na(.data$new_clone)) %>%
            select(all_of(c("CTaa_donor_source", "new_clone")))

        obj@meta.data <- obj@meta.data %>%
            left_join(custom, by = "CTaa_donor_source")

        obj@meta.data <- obj@meta.data %>%
            mutate(CTaa = case_when(
                !is.na(.data$new_clone) ~ .data$new_clone,
                is.na(.data$new_clone) ~ .data$CTaa
            )) %>%
            select(-all_of("new_clone"))
    }


    obj@meta.data <- obj@meta.data %>%
        select(-all_of("CTaa_donor_source"))

    return(obj)
}








#' @title Reformat Seurat object
#'
#' @description This internal function renames columns in the Seurat
#' meta.data, removes unneeded columns, and formats the columns for use with
#' scRepertoire downstream. It is a helper function for ECLIPSE.
#'
#' @param Seurat_object Seurat object that contains scRNA-seq data.
#' @param Group name of column in Seurat meta.data. This column should specify
#' how the final statistics for each clone should be calculated. For example,
#' should they be calculated per sample, per donor, or among all cells?
#' This can and often will be the same column as donor.
#' @param DetailedEclipseStats Provides detailed stats
#' relating to EM assignment in Seurat meta.data for each cell if TRUE.
#' @param Verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline.
#'
#' @return A Seurat object that is properly formatted and has all TCR
#' information added.
.reformatSeuratObject <- function(Seurat_object,
                                    Group,
                                    DetailedEclipseStats,
                                    Verbose) {

    obj <- Seurat_object
    group <- Group
    detailedEclipseStats <- DetailedEclipseStats
    verbose <- Verbose

    if (sum(colnames(obj@meta.data) %in% c("cloneCallSource", "CTaa",
        "scR_CTaa", "VDJdive_clone")) != 4) {
        stop("Seurat_object should be made by running ",
            ".assignFinalClonotypes. Please try again.")
    }

    if (!group %in% colnames(obj@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.logical(detailedEclipseStats)) {
        stop(
            "TRUE or FALSE not provided to detailedEclipseStats. ",
            "Please try again."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    ### Making columns so that scRepertoire visualizations can be used.
    ### Grouping by the group specified in the arguments, not the donor
    obj@meta.data <- obj@meta.data %>%
        group_by(.data[[group]]) %>%
        mutate(eclipseTcrGroupSize = sum(!is.na(.data$CTaa))) %>%
        ungroup()
    obj@meta.data <- obj@meta.data %>%
        group_by(.data$CTaa, .data[[group]]) %>%
        mutate(
            clonalProportion = case_when(
                !is.na(.data$CTaa) ~ n() / first(.data$eclipseTcrGroupSize),
                is.na(.data$CTaa) ~ NA
            ),
            clonalFrequency = case_when(
                !is.na(.data$CTaa) ~ n(),
                is.na(.data$CTaa) ~ NA
            )
        ) %>%
        ungroup() %>%
        select(-all_of("eclipseTcrGroupSize"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneSize = factor(case_when(
            .data$clonalProportion <= 1 & .data$clonalProportion > 0.1
            ~ "Hyperexpanded (0.1 < X <= 1)",
            .data$clonalProportion <= 0.1 & .data$clonalProportion > 0.01
            ~ "Large (0.01 < X <= 0.1)",
            .data$clonalProportion <= 0.01 & .data$clonalProportion > 0.001
            ~ "Medium (0.001 < X <= 0.01)",
            .data$clonalProportion <= 0.001 & .data$clonalProportion > 1e-04
            ~ "Small (1e-04 < X <= 0.001)",
            .data$clonalProportion <= 1e-04 & .data$clonalProportion > 0
            ~ "Rare (0 < X <= 1e-04)",
            .data$clonalProportion <= 0
            ~ "None ( < X <= 0)",
            is.na(.data$clonalProportion)
            ~ NA
        )))

    ### Making a column that notes whether the cells in the clone appear
    ### to be MAIT or iNKT based on TCR segment usage
    obj@meta.data <- obj@meta.data %>%
        mutate(
            scR_mait_conv = case_when(
                str_detect(.data$scR_CTgene, "TRAV1-2\\.TRAJ33") ~
                    "MAIT Conventional (TRAV1-2 TRAJ33)",
                TRUE ~
                    "Normal"
            ),
            scR_mait_unconv = case_when(
                str_detect(
                    .data$scR_CTgene,
                    "TRAV1-2\\.TRAJ12|TRAV1-2\\.TRAJ20"
                )
                ~ "MAIT Unconventional (TRAV1-2 TRAJ12 or TRAJ20)",
                TRUE ~ "Normal"
            ),
            scR_inkt = case_when(
                str_detect(.data$scR_CTgene, "TRAV10\\.TRAJ18") &
                    str_detect(.data$scR_CTgene, "TRBV25-1")
                ~ "iNKT",
                TRUE ~ "Normal"
            )
        )
    obj@meta.data <- obj@meta.data %>%
        mutate(scR_unconventional_count = 3 -
            (.data$scR_mait_conv == "Normal") -
            (.data$scR_mait_unconv == "Normal") -
            (.data$scR_inkt == "Normal"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cell_unconventional_subset = case_when(
            .data$scR_unconventional_count == 0 ~ "Conventional",
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_conv != "Normal"
            ~ .data$scR_mait_conv,
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_unconv != "Normal"
            ~ .data$scR_mait_unconv,
            .data$scR_unconventional_count == 1 & .data$scR_inkt != "Normal"
            ~ .data$scR_inkt,
            .data$scR_unconventional_count > 1
            ~ "Multiple Unconventional Types"
        ))

    df <- obj@meta.data %>%
        group_by(
            .data$CTaa,
            .data$cell_unconventional_subset,
            .data[[group]]
        ) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        filter(.data$cell_unconventional_subset != "Conventional") %>%
        arrange(desc(.data$count))

    if (nrow(df) > 0) {
        obj@meta.data <- obj@meta.data %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                clone_unconventional_subset =
                    first(.data$cell_unconventional_subset)
            )
        df <- df %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            ungroup() %>%
            select(all_of(c("CTaa_group", "clone_unconventional_subset")))

        if ("clone_unconventional_subset" %in% colnames(obj@meta.data)) {
            obj@meta.data <- obj@meta.data %>%
                select(-all_of("clone_unconventional_subset"))
        }

        obj@meta.data <- obj@meta.data %>% left_join(df, by = "CTaa_group")
        obj@meta.data <- obj@meta.data %>%
            mutate(clone_unconventional_subset = case_when(
                is.na(.data$clone_unconventional_subset)
                ~ "Conventional",
                !is.na(.data$clone_unconventional_subset)
                ~ .data$clone_unconventional_subset
            )) %>%
            select(-all_of("CTaa_group"))
        if (verbose == TRUE) {
            message(
                sum(obj$clone_unconventional_subset != "Conventional"),
                " cells were detected in clones that appear to ",
                "be MAIT or iNKT cells"
            )
        }
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of(c(
        "scR_mait_conv",
        "scR_mait_unconv",
        "scR_inkt",
        "scR_unconventional_count",
        "tcrEclipseGroup",
        "tcrEclipse_barcode"))
    )

    obj@meta.data <- as.data.frame(obj@meta.data)
    rownames(obj@meta.data) <- obj$barcode

    if (verbose == TRUE) {
        message(
            round(sum(obj$cloneCallSource == "VDJdive", na.rm = TRUE) /
                nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype using ",
            "high confidence prediction from VDJdive::clonoStats()"
        )
        message(
            round(sum(!is.na(obj$CTaa)) / nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype ",
            "in the end (either from prediction or traditional analysis)"
        )
        message("These final clonal annotations are stored in CTaa")
        message("Summary of the top 5 largest clones:")

        big_clones <- obj@meta.data %>%
            group_by(.data[[group]]) %>%
            mutate(group_count = sum(!is.na(.data$CTaa))) %>%
            ungroup() %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                Clone_Size = n(),
                Percent_of_Repertoire = 100 * n() /
                    first(.data$group_count)
            ) %>%
            select(all_of(c(
                "CTaa",
                group,
                "Clone_Size",
                "Percent_of_Repertoire"
                ))
            ) %>%
            arrange(desc(.data$Clone_Size)) %>%
            filter(!is.na(.data$CTaa)) %>%
            head(n = 5)
        print(big_clones)
    }

    ### Cleaning up object before returning to user
    if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)) {
        obj@meta.data <- obj@meta.data %>%
            mutate(barcode = .data$originalBarcodeBeforeECLIPSE) %>%
            select(-all_of("originalBarcodeBeforeECLIPSE"))
    } else if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)
    == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of("barcode"))
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of("VDJdive_clone"))

    if (detailedEclipseStats == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of(c(
            "EM_max_prop",
            "EM_2ndmax_prop",
            "EM_max_2ndmax_ratio"))
        )
    }

    return(obj)
}


