#' @title Determines the true clonotype of each cell
#'
#' @description Analyzes scTCR-seq data and returns results in Seurat object.
#' Specifically, chain pairings across the sample are analyzed to predict
#' chains lost to sequencing dropout and to decipher whether cells with >2
#' total TCR chains are caused by technical artifacts or true biology.
#'
#' @param seurat_object Seurat object that contains scRNA-seq data.
#' @param contig_folder_paths vector of paths to directories that contain TCR
#' files. If you use this argument, provide either "all" or "filtered" to
#' contig_type. The order should match the numbers provided in batch.
#' @param contig_type tells whether contig_folder_paths or contig_file_paths
#' should be used to read in contig files and how. To provide file paths,
#' this should be "manual". To provide folders, this should be "all" or
#' "filtered" to search for all or filtered contig annotations files.
#' @param contig_file_paths vector of TCR file paths. The order should match
#' the numbers provided in batch. The files can have any name, but they must
#' be all_contig_annotations.csv or
#' filtered_contig_annotations.csv format.
#' @param batch name of column in Seurat meta.data. The column should contain
#' integers, not characters, that tell which position in the
#' contig_folder/file_path vectors contains the TCR contigs for the cell.
#' @param donor name of column in Seurat meta.data. The column should contain
#' the donor, regardless of treatment/timepoint/condition.
#' @param group name of column in Seurat meta.data. This column should specify
#' how the final statistics for each clone should be calculated. For example,
#' should they be calculated per sample, per donor, or among all cells?
#' This can and often will be the same column as donor.
#' @param original_barcode name of column in Seurat meta.data. This must
#' contain the original 10X barcode for the cell. Typically this will end with
#' -1, but users should check their TCR contig files to be sure.
#' @param format default is blank. This tells how missing data is stored in
#' the TCR contig files. Current files likely contain "" so "blank" is the
#' default, but if older contig files are used where any missing data is stored
#' as "None", then "None" should be provided here.
#' @param detailedEclipseStats default is FALSE. Provides detailed stats
#' relating to EM assignment in Seurat meta.data for each cell if TRUE.
#' @param combineTcrArgs arguments that are provided to scRepertoire
#' combineTCR. This should be a list of 4 logical elements,
#' with the elements corresponding to the combineTCR arguments.
#' 1 is removeNA, 2 is removeMulti, 3 is filterMulti, 4 is filterNonproductive.
#' If no clonotype can be predicted with high confidence, then combineTCR
#' using these arguments will determine the final clonotype. The default is
#' the combineTCR default where only filterNonproductive is TRUE.
#' @param allowDoubleChains specifies whether 2 alpha chains, 2 beta chains,
#' or either 2 alpha or 2 beta chains are allowed. Default is both.
#' @param returnFilteredContigs specifies if filtered contigs should be returned
#' instead of running the whole ECLIPSE pipeline. Default is FALSE.
#' @param verbose specifies whether to print an important
#' QC graph, a small table, and also notes through the pipeline. Default is
#' TRUE.
#'
#'
#' @return A Seurat object with TCR data provided in meta.data.
#'
#' @examples
#' data(test_seurat)
#'
#' files <- system.file("extdata/examples", package = "ECLIPSE")
#' tcr_paths <- list.files(
#'                 files,
#'                 pattern = "\\.csv.gz$",
#'                 full.names = TRUE)
#'
#' cd8_tcr <- ECLIPSE(test_seurat,
#'                     contig_file_paths = tcr_paths,
#'                     contig_type = "manual",
#'                     batch = "index",
#'                     donor = "donor",
#'                     group = "donor",
#'                     original_barcode = "Cell_Barcode")
#'
#' @export
ECLIPSE <- function(seurat_object,
                    contig_type = "all",
                    contig_folder_paths = NULL,
                    contig_file_paths = NULL,
                    batch,
                    donor,
                    group,
                    original_barcode,
                    format = "Blank",
                    detailedEclipseStats = FALSE,
                    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
                    allowDoubleChains = "both",
                    returnFilteredContigs = FALSE,
                    verbose = TRUE) {

    #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (!contig_type %in% c(
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

    if (is.null(contig_folder_paths) & contig_type != "manual") {
        stop(
            "Please provide folders to contig_folder_paths that have the ",
            "contig annotations in them."
        )
    }

    if (is.null(contig_file_paths) & contig_type == "manual") {
        stop(
            "No contig file paths provided. ",
            "Please provide a vector of file path names to contig_file_paths"
        )
    }

    if (contig_type != "manual" &
        (!is.vector(contig_folder_paths) | is.list(contig_folder_paths) |
        length(contig_folder_paths) == 0)) {
        stop("contig_folder_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (contig_type == "manual" &
        (!is.vector(contig_file_paths) | is.list(contig_file_paths) |
        length(contig_file_paths) == 0)) {
        stop("contig_file_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (!is.null(contig_file_paths) & !is.null(contig_folder_paths)) {
        stop("Please do not use contig_file_paths and contig_folder_paths. ",
            "Please try again.")
    }

    if (!batch %in% colnames(seurat_object@meta.data)) {
        stop(
            "Batch column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.numeric(seurat_object@meta.data[[batch]])) {
        stop("Batch column should be numeric. Please try again.")
    }

    if (sum((seurat_object@meta.data[[batch]] %% 1) != 0) |
       1 %in% seurat_object@meta.data[[batch]] == FALSE ) {
       stop("Batch column must contain whole numbers and start with 1. ",
           "Please try again.")
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

    if (is.numeric(seurat_object@meta.data[[donor]]) |
        is.numeric(seurat_object@meta.data[[group]])) {
        stop("Donor and group cannot be numeric. Please try again.")
    }

    if (!original_barcode %in% colnames(seurat_object@meta.data)) {
        stop(
            "Original barcode column name provided is not ",
            "in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.logical(detailedEclipseStats)) {
        stop(
            "TRUE or FALSE not provided to detailedEclipseStats. ",
            "Please try again."
        )
    }

    if (format %in% c("None", "none", "Blank") == FALSE) {
        stop(
            "None, none, or Blank not provided to format. ",
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

    if (allowDoubleChains %in% c(
        "Alpha", "alpha",
        "beta", "Beta",
        "Both", "both",
        "Neither", "neither"
    ) == FALSE) {
        stop(
            "Alpha, Beta, Both, or Neither not provided to ",
            "allowDoubleChains. ",
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

    output1 <- .readContigs(Contig_type = contig_type,
                            Contig_folder_paths = contig_folder_paths,
                            Contig_file_paths = contig_file_paths,
                            Verbose = verbose)

    output2 <- .qcFilterContigs(output1,
                                Seurat_object = seurat_object,
                                Batch = batch,
                                Donor = donor,
                                Group = group,
                                Original_barcode = original_barcode,
                                Format = format,
                                ReturnFilteredContigs = returnFilteredContigs,
                                Verbose = verbose)

    if (returnFilteredContigs == TRUE) {
      return(output2)
    }


    if (verbose == TRUE) message("\nFinding clones with extra chains...\n")

    donors <- unique(as.vector(output2[[1]][, donor])[[1]])
    new_contigs <- vector("list", length(donors))

    for (i in seq_along(donors)) {
        ### THIS PART IS FINDING THE COMBO CLONES AND
        ### REWRITING THE CONTIG FILES
        if (allowDoubleChains %in% c("Neither", "neither") == FALSE &
            verbose == TRUE) {
            message(
                "Running donor ", i, "/",
                length(donors), " (", donors[i], "):"
            )
        }

        slim <- output2[[1]] %>% filter(.data[[donor]] == donors[i])

        if (allowDoubleChains %in% c("Alpha", "alpha", "Both", "both")) {
            ### Calling clones with multiple alpha chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "alpha",
                printStats = verbose
            )
        }

        if (allowDoubleChains %in% c("Beta", "beta", "Both", "both")) {
            ### Calling clones with multiple beta chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "beta",
                printStats = verbose
            )
        }

        new_contigs[[i]] <- slim

        if (verbose == TRUE) message()
    }

    contigs_post <- bind_rows(new_contigs)


    output3 <- .runVDJdive(Contigs = contigs_post,
                           Donor = donor,
                           Verbose = verbose)

    output4 <- .assignFinalClonotypes(em_results = output3,
                                      Contigs = output2[[1]],
                                      Seurat_object = output2[[2]],
                                      Orphan_contigs = output2[[3]],
                                      Donor = donor,
                                      Group = group,
                                      CombineTcrArgs = combineTcrArgs,
                                      Verbose = verbose)

    output5 <- .reformatSeuratObject(Seurat_object = output4,
        Group = group,
        DetailedEclipseStats = detailedEclipseStats,
        Verbose = verbose)

    if (returnFilteredContigs == FALSE) {
        return(output5)
    }

}

