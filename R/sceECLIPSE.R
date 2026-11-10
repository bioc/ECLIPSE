#' @title Determines the true clonotype of each cell
#'
#' @description Analyzes scTCR-seq data and returns results as a
#' SingleCellExperiment object. Unlike regular ECLIPSE, this function
#' requires a SingleCellExperiment object as an input, not a Seurat object.
#' Specifically, chain pairings across the sample are analyzed to predict
#' chains lost to sequencing dropout and to decipher whether cells with >2
#' total TCR chains are caused by technical artifacts or true biology.
#'
#' @param sce_object SingleCellExperiment object that contains scRNA-seq data.
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
#' @param batch name of column in colData. The column should contain
#' integers, not characters, that tell which position in the
#' contig_folder/file_path vectors contains the TCR contigs for the cell.
#' @param donor name of column in colData. The column should contain
#' the donor, regardless of treatment/timepoint/condition.
#' @param group name of column in colData. This column should specify
#' how the final statistics for each clone should be calculated. For example,
#' should they be calculated per sample, per donor, or among all cells?
#' This can and often will be the same column as donor.
#' @param original_barcode name of column in colData. This must
#' contain the original 10X barcode for the cell. Typically this will end with
#' -1, but users should check their TCR contig files to be sure.
#' @param format default is blank. This tells how missing data is stored in
#' the TCR contig files. Current files likely contain "" so "blank" is the
#' default, but if older contig files are used where any missing data is stored
#' as "None", then "None" should be provided here.
#' @param detailedEclipseStats default is FALSE. Provides detailed stats
#' relating to EM assignment in colData for each cell if TRUE.
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
#' @return A SingleCellExperiment object with TCR data provided in colData.
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
#' sce_obj <- SeuratObject::as.Seurat(test_seurat)
#'
#' cd8_tcr <- sceECLIPSE(sce_obj,
#'                     contig_file_paths = tcr_paths,
#'                     contig_type = "manual",
#'                     batch = "index",
#'                     donor = "donor",
#'                     group = "donor",
#'                     original_barcode = "Cell_Barcode")
#'
#' @export
sceECLIPSE <- function(sce_object,
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

    if (inherits(sce_object, "SingleCellExperiment") == FALSE) {
        stop("sce_object is not a SingleCellExperiment object. ",
            "Please try again.")
    }

    if (requireNamespace("SeuratObject", quietly = TRUE) == FALSE) {
        stop("Please install SeuratObject to use sceECLIPSE. ",
            "Please try again.")
    }

    if (requireNamespace("SummarizedExperiment", quietly = TRUE) == FALSE) {
        stop("Please install SummarizedExperiment to use sceECLIPSE. ",
            "Please try again.")
    }

    if (requireNamespace("S4Vectors", quietly = TRUE) == FALSE) {
        stop("Please install S4Vectors to use sceECLIPSE. ",
            "Please try again.")
    }

    ### All other arguments are identical to ECLIPSE so it will check them

    seurat_obj <- SeuratObject::as.Seurat(sce_object, data = NULL)

    seurat_object <- ECLIPSE(seurat_object = seurat_obj,
        contig_type = contig_type,
        contig_folder_paths = contig_folder_paths,
        contig_file_paths = contig_file_paths,
        batch = batch,
        donor = donor,
        group = group,
        original_barcode = original_barcode,
        format = format,
        detailedEclipseStats = detailedEclipseStats,
        combineTcrArgs = combineTcrArgs,
        allowDoubleChains = allowDoubleChains,
        returnFilteredConitgs = returnFilteredContigs,
        verbose = verbose
        )

    if (returnFilteredContigs == TRUE) {
        return(seurat_object)
    }

    new_cols <- colnames(seurat_object@meta.data)[
        colnames(seurat_object@meta.data) %in% colnames(seurat_obj@meta.data)
        == FALSE]

   seurat_object$cell_barcode_for_eclipse <- rownames(seurat_object@meta.data)

    meta <- seurat_object@meta.data %>% select(
        all_of(c("cell_barcode_for_eclipse", new_cols)))

    new_meta <- as.data.frame(SummarizedExperiment::colData(sce_object)) %>%
        rownames_to_column("cell_barcode_for_eclipse") %>%
        left_join(meta, by = "cell_barcode_for_eclipse") %>%
        column_to_rownames("cell_barcode_for_eclipse")

    if (sum(rownames(new_meta) == colnames(sce_object)) !=
        length(rownames(new_meta))) {
        stop("Barcodes were changed in conversion. ",
            "Please check your object.")
    }

    SummarizedExperiment::colData(sce_object) <- S4Vectors::DataFrame(new_meta)

    return(sce_object)

}
