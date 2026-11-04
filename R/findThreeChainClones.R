#' @title Finds clones with 3 total TCR chains and modifies their contigs
#'
#' @description Analyzes chain pairings across a sample to find clones with
#' 2 TCR alpha or 2 TCR beta chains. Then, the contigs are modified so that
#' these cells only have 1 alpha and 1 beta contig, with ::: separating each
#' CDR3 sequence. These contigs are then returned.
#'
#'
#' @param contigs tibble of contigs from cellranger all or
#' filtered contig files.
#' @param chain "alpha" or "beta". This tells it whether clones with two of
#' which chain should be looked for.
#' @param mode "standard" is the default, where modified contigs are returned.
#' If "show" is provided, two tables are returned. The first
#' shows the final set of clones with a double chain of one type, and the
#' second shows how/why some of those clones were called as such. The second
#' is pre-filtering, so some clones may not have been called as having 3 total
#' chains in the end.
#' @param printStats specifies whether notes/stats
#' should be printed. Default is TRUE.
#'
#'
#'
#' @return a tibble of modified contigs,
#' except if mode show is selected.
#'
#' @examples
#' data(post_qc_contigs)
#'
#' modified_contigs <- findThreeChainClones(
#'                         contigs = post_qc_contigs,
#'                         chain = "a")
#'
#' @export
findThreeChainClones <- function(contigs,
                                    chain,
                                    mode = "standard",
                                    printStats = TRUE) {
    if (sum(c("barcode", "high_confidence", "chain", "cdr3", "umis") %in%
        colnames(contigs)) != 5) {
        stop(
            "Contigs should be in the all_contig_annotations or ",
            "filtered_contig_annotations format. ",
            "Please try again."
        )
    }

    if (chain %in% c(
        "alpha", "Alpha", "a", "A", "ALPHA",
        "beta", "Beta", "b", "B", "BETA"
    ) == FALSE) {
        stop(
            "Chain should be alpha or beta. ",
            "Please try again."
        )
    }

    if (mode %in% c("standard", "show") == FALSE) {
        stop(
            "Mode must be standard or show ",
            "Please try again."
        )
    }

    if (!is.logical(printStats)) {
        stop(
            "printStats must be TRUE or FALSE. ",
            "Please try again."
        )
    }

    if (chain %in% c("alpha", "Alpha", "a", "A", "ALPHA")) {
        chain <- "a"
    } else if (chain %in% c("beta", "Beta", "b", "B", "BETA")) {
        chain <- "b"
    }



        output1 <- .identifyThreeChainClones(Contigs = contigs,
                                             Chain = chain,
                                             Mode = mode)

        if (mode == "show") {
            return(output1)
        }

        output2 <- .modifyThreeChainContigs(Contigs = contigs,
                                            three_chain_clones = output1[[1]],
                                            Clones = output1[[2]],
                                            Chain = chain,
                                            Mode = mode)

        output3 <- .findMissingMonotypicChains(Contigs = output2,
                                               specialTCRs = output1[[3]],
                                               Chain = chain,
                                               Mode = mode,
                                               print_stats = printStats)

        if (mode == "standard") {
            return(output3)
        }

}

