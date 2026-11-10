data(test_seurat)
data(post_qc_contigs)
data(unit_test_object)

files <- system.file("extdata/examples", package = "ECLIPSE")
tcr_paths <- list.files(files,
    pattern = "\\.csv.gz$",
    full.names = TRUE
)

unit_test_files <- system.file("extdata/unit_tests/D5L", package = "ECLIPSE")
unit_test_tcr_paths <- list.files(unit_test_files,
    pattern = "\\.csv",
    full.names = TRUE)

test_seurat_with_tcr <- ECLIPSE(
    contig_file_paths = tcr_paths,
    contig_type = "manual",
    seurat_object = test_seurat,
    batch = "index",
    donor = "donor",
    group = "donor",
    original_barcode = "Cell_Barcode",
    detailedEclipseStats = FALSE,
    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
    allowDoubleChains = "both",
    verbose = FALSE
)

out_contigs <- findThreeChainClones(post_qc_contigs,
    chain = "b"
)

filtered_seurat <- tcrArtifactDetect(test_seurat_with_tcr)

test_seurat$index2 <- 1.2




test_that("ECLIPSE accepts only valid inputs", {
    expect_error(
        ECLIPSE(
            contig_folder_paths = list(tcr_paths),
            contig_type = "all",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "contig_folder_paths must be a vector"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = list(tcr_paths),
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "contig_file_paths must be a vector"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = list(),
            contig_type = "all",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Please provide folders to contig_folder_paths that have"
    )

    expect_error(
        ECLIPSE(
            contig_type= "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "No contig file paths provided"
    )


    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "fake",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Batch column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "fake",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Donor column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "fake",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Group column name provided is not in the Seurat object"
    )


    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "bc",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Original barcode column name provided is not in the Seurat object"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "fake",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Incorrect contig type"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = 1,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "TRUE or FALSE not provided to detailedEclipseStats"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            format = "fake"
        ),
        "None, none, or Blank not provided to format"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(1, 2, 3, 4),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "Incorrect input to combineTcrArgs"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "fake",
            verbose = FALSE
        ),
        "Alpha, Beta, Both, or Neither not provided to allowDoubleChains"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = "fale"
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = 1
        ),
        "TRUE or FALSE not provided to returnFilteredContigs"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_folder_path = c(files, files),
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
        ),
        "Please do not use contig_file_paths and contig_folder_paths"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "donor",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
        ),
        "Batch column should be numeric"
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index2",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
        ),
        "Batch column must contain whole numbers and start with 1."
    )

    expect_error(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "index",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
        ),
        "Donor and group cannot be numeric"
    )
})




test_that("ECLIPSE returns the correct output", {
    ### Expected columns present in meta.data
    expect_equal(
        c(
            "CTaa", "scR_CTaa", "clonalFrequency",
            "orig.ident", "nCount_RNA"
        ) %in%
            colnames(test_seurat_with_tcr@meta.data),
        rep(TRUE, 5)
    )

    ### Meta data retained
    expect_true(nrow(test_seurat_with_tcr@meta.data) > 1 &
        ncol(test_seurat_with_tcr@meta.data) > 1)

    ### clonalProportion adds up to 1 for each group
    sum_freq_per_donor <- test_seurat_with_tcr@meta.data %>%
        filter(!is.na(clonalProportion)) %>%
        group_by(donor, CTaa) %>%
        summarise(prop = first(clonalProportion)) %>%
        group_by(donor) %>%
        summarise(total = sum(prop))

    expect_equal(mean(sum_freq_per_donor$total), 1)

    ### clonalFrequency is whole numbers
    expect_equal(
        sum(!is.na(test_seurat_with_tcr@meta.data$clonalFrequency) &
              test_seurat_with_tcr$clonalFrequency %% 1 != 0),
        0
    )

    ### Clone names are formatted correctly (i.e no -)
    expect_equal(
        length(grep("-", test_seurat_with_tcr$CTaa)),
        0
    )

    ### Columns didn't get removed
    expect_equal(
        sum(as.character(seq(1, nrow(test_seurat_with_tcr@meta.data)))
        == rownames(test_seurat_with_tcr@meta.data)),
        0
    )

    ### Number of cells annotated in scR_CTaa = CTaa
    expect_equal(
        sum(is.na(test_seurat_with_tcr$scR_CTaa)),
        sum(is.na(test_seurat_with_tcr$CTaa))
    )

    ### returnFilteredContigs = TRUE should return a tibble
    expect_s3_class(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = TRUE,
            returnFilteredContigs = TRUE
        ),
        "tbl_df"
    )

    expect_s4_class(
        ECLIPSE(
            contig_file_paths = unit_test_tcr_paths,
            contig_type = "manual",
            seurat_object = unit_test_object,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = TRUE,
            returnFilteredContigs = FALSE
        ),
        "Seurat"
    )

unit_test_object$barcode <- unit_test_object$Cell_Barcode

    expect_s4_class(
        ECLIPSE(
            contig_file_paths = unit_test_tcr_paths,
            contig_type = "manual",
            seurat_object = unit_test_object,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = TRUE,
            returnFilteredContigs = FALSE
        ),
        "Seurat"
    )


})


test_that("The verbose argument of ECLIPSE works correctly", {

    expect_message(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = TRUE,
        ),
        "Finding clones with extra chains"
    )

    expect_silent(
        ECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            seurat_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
        )
    )

})

test_that("findThreeChainClones accepts only valid inputs", {
    expect_error(
        findThreeChainClones(post_qc_contigs[, 1:4],
            chain = "a",
            mode = "standard",
            printStats = TRUE
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "fake",
            mode = "standard",
            printStats = TRUE
        ),
        "Chain should be alpha or beta"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "fake",
            printStats = TRUE
        ),
        "Mode must be standard or show"
    )

    expect_error(
        findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "standard",
            printStats = "fake"
        ),
        "printStats must be TRUE or FALSE"
    )
})




test_that("findThreeChainClones returns correct outputs", {

    ### Returns data.frame
    expect_true(is.data.frame(out_contigs))

    ### This should keep contigs as is or join them together
    ### leading to the same or less rows
    expect_lte(
        nrow(out_contigs),
        nrow(post_qc_contigs)
    )

    ### findThreeChainClones shouldn't add any new columns
    expect_equal(
        ncol(out_contigs),
        ncol(post_qc_contigs)
    )

    ### Expected columns present
    expect_equal(
        sum(c("barcode", "high_confidence", "chain", "cdr3", "umis")
        %in% colnames(out_contigs)),
        5
    )

    ### mode = "show" returns a list
    expect_true(
        is.list(findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "show",
            printStats = TRUE
            )
        )
    )

    expect_true(
        sum(colnames(findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "show",
            printStats = TRUE)[[1]]
        ) %in% c("monotypic_chain", "predicted_bitypic_pairing")) == 2
    )

    expect_true(
        sum(colnames(findThreeChainClones(post_qc_contigs,
            chain = "a",
            mode = "show",
            printStats = TRUE)[[2]]
        ) %in% c("total_cells", "cell_counts_per_bitypic_pairing")) == 2
    )


})




test_that("tcrArtifactDetect accepts only valid inputs", {

    expect_error(
        tcrArtifactDetect(test_seurat),
            "The seurat object must have had ECLIPSE ran on it"
    )

    expect_error(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = "fake",
            singleChainLimit = 2,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "filter, showPlots, and verbose must be TRUE or FALSE"
    )

    expect_error(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 2.5,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "singleChainLimit and totalChainLimit both must be"
    )

    expect_error(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 0,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
        ),
        "singleChainLimit must be at least 2"
    )

    expect_warning(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 1,
            totalChainLimit = -1,
            showPlots = TRUE,
            verbose = TRUE
            ),
            "Filtering cells with more than 1"
    )

    expect_error(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 0,
            totalChainLimit = "1",
            showPlots = TRUE,
            verbose = TRUE
            ),
            "singleChainLimit and totalChainLimit both must be numbers"
    )


    expect_warning(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 2,
            totalChainLimit = 2,
            showPlots = TRUE,
            verbose = TRUE
            ),
            "Filtering all cells with more than 2 total chains"
    )
})



test_that("tcrArtifactDetect returns a filtered Seurat object", {
    ### No columns are added
    expect_equal(
        ncol(test_seurat_with_tcr@meta.data),
        ncol(filtered_seurat@meta.data)
    )

    ### Rows are the same or decreased
    expect_lte(
        nrow(filtered_seurat@meta.data),
        nrow(test_seurat_with_tcr@meta.data)
    )
})

test_that("The verbose argument of tcrArtifactDetect works", {

    expect_message(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 2,
            totalChainLimit = 3,
            showPlots = TRUE,
            verbose = TRUE
            ),
        "were removed due to having >"
    )

    expect_silent(
        tcrArtifactDetect(test_seurat_with_tcr,
            filter = TRUE,
            singleChainLimit = 2,
            totalChainLimit = 3,
            showPlots = TRUE,
            verbose = FALSE
        )
    )
})




### Making all inputs
seurat_object <- test_seurat
contig_type <- "manual"
contig_folder_paths <- NULL
contig_file_paths <- tcr_paths
original_barcode <- "Cell_Barcode"
verbose <- FALSE
batch <- "index"
donor <- "donor"
group <- "donor"
format <- "Blank"
returnFilteredContigs <- FALSE
combineTcrArgs <- list(FALSE, FALSE, FALSE, TRUE)
detailedEclipseStats <- FALSE
allowDoubleChains <- TRUE

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




test_that(".readContigs works correctly", {

    expect_error(
        .readContigs(
            Contig_type = "fake",
            Contig_folder_paths = NULL,
            Contig_file_paths = tcr_paths,
            Verbose = TRUE
        ),
        "Incorrect contig type. Should be all or filtered"
    )

    expect_error(
        .readContigs(
            Contig_type = "all",
            Contig_folder_paths = NULL,
            Contig_file_paths = NULL,
            Verbose = TRUE
        ),
        "Please provide folders to contig_folder_paths that have the"
    )

    expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = NULL,
            Contig_file_paths = NULL,
            Verbose = TRUE
        ),
        "No contig file paths provided"
    )

    expect_error(
        .readContigs(
            Contig_type = "all",
            Contig_folder_paths = list("a", "b"),
            Contig_file_paths = NULL,
            Verbose = TRUE
        ),
        "contig_folder_paths must be a vector"
    )

    expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = NULL,
            Contig_file_paths = list("a", "b"),
            Verbose = TRUE
        ),
        "contig_file_paths must be a vector"
    )

    expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = c("a", "b"),
            Contig_file_paths = c("a", "b"),
            Verbose = TRUE
        ),
        "Please do not use contig_file_paths and contig_folder_paths"
    )

    expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = c("a", "b"),
            Contig_file_paths = c("a", "b"),
            Verbose = TRUE
        ),
        "Please do not use contig_file_paths and contig_folder_paths"
    )

     expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = NULL,
            Contig_file_paths = c("a", "b"),
            Verbose = 1
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_error(
        .readContigs(
            Contig_type = "manual",
            Contig_folder_paths = NULL,
            Contig_file_paths =
                list.files(system.file("extdata/",
                    package = "ECLIPSE"), full.names = TRUE)[1],
            Verbose = TRUE
        ),
        "did not lead to any contigs being loaded in"
    )

    vignette_dir <- system.file("extdata/vignette", package = "ECLIPSE")
    tcr_folder_paths <- list.dirs(path = vignette_dir, full.names = TRUE)[-1]

    expect_s3_class(
        .readContigs(
            Contig_type = "all",
            Contig_folder_paths = tcr_folder_paths,
            Contig_file_paths = NULL,
            Verbose = TRUE
        ),
        "data.frame"
    )

})



test_that(".qcFilterContigs works correctly", {

  expect_error(
      .qcFilterContigs(
          Contigs = output1[, 1:2],
          Seurat_object = test_seurat,
          Batch = "index",
          Donor = "donor",
          Group = "donor",
          Original_barcode = "original_barcode",
          Format = "Blank",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "Contigs should be in a contigs format"
  )

  expect_error(
      .qcFilterContigs(
          Contigs = output1,
          Seurat_object = test_seurat,
          Batch = "index4",
          Donor = "donor",
          Group = "donor",
          Original_barcode = "Cell_Barcode",
          Format = "Blank",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "Batch column name provided is not in the Seurat"
  )

  expect_error(
      .qcFilterContigs(
          Contigs = output1,
          Seurat_object = test_seurat,
          Batch = "index",
          Donor = "donor2",
          Group = "donor",
          Original_barcode = "Cell_Barcode",
          Format = "Blank",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "Donor column name provided is not in the Seurat object"
  )

    expect_error(
      .qcFilterContigs(
          Contigs = output1,
          Seurat_object = test_seurat,
          Batch = "index",
          Donor = "donor",
          Group = "donor2",
          Original_barcode = "Cell_Barcode",
          Format = "Blank",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "Group column name provided is not in the Seurat object"
  )

    expect_error(
      .qcFilterContigs(
          Contigs = output1,
          Seurat_object = test_seurat,
          Batch = "index",
          Donor = "donor",
          Group = "donor",
          Original_barcode = "fake",
          Format = "Blank",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "Original barcode column name provided is not"
  )

    expect_error(
      .qcFilterContigs(
          Contigs = output1,
          Seurat_object = test_seurat,
          Batch = "index",
          Donor = "donor",
          Group = "donor",
          Original_barcode = "Cell_Barcode",
          Format = "fake",
          ReturnFilteredContigs = FALSE,
          Verbose = TRUE
      ),
      "None, none, or Blank not provided to format"
  )

    expect_error(
        .qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = 2
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_error(
        .qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = 2,
            Verbose = TRUE
        ),
        "TRUE or FALSE not provided to returnFilteredContigs"
    )

    expect_true(
        is.list(.qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = TRUE
        ))
    )

    expect_true(
        sum(colnames(.qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = TRUE)[[1]])
        %in% c("barcode", "reads", "cdr3", "chainTypes", "cdr3_nt")) == 5
    )

    expect_s4_class(
        .qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = TRUE
        )[[2]],
        "Seurat"
    )

    expect_true(
        sum(colnames(.qcFilterContigs(
            Contigs = output1,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = TRUE)[[3]])
        %in% c("barcode", "chain", "cdr3")) == 3
    )

    ### This leads to no contigs being left, triggering the error
    expect_error(
        .qcFilterContigs(
            Contigs = post_qc_contigs,
            Seurat_object = test_seurat,
            Batch = "index",
            Donor = "donor",
            Group = "donor",
            Original_barcode = "Cell_Barcode",
            Format = "Blank",
            ReturnFilteredContigs = FALSE,
            Verbose = TRUE
        ),
        "Check that your barcodes are the original"
    )

})


test_that(".runVDJdive works correctly", {

    expect_error(
        .runVDJdive(
            Contigs = contigs_post[,1:2],
            Verbose = TRUE,
            Donor = "donor"
        ),
        "Contigs should be in a contigs format and generated by"
    )

    expect_error(
        .runVDJdive(
            Contigs = contigs_post,
            Verbose = 2,
            Donor = "donor"
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_true(
        is.data.frame(.runVDJdive(
            Contigs = contigs_post,
            Verbose = TRUE,
            Donor = "donor"
        ))
    )

})



test_that(".assignFinalClonotypes works correctly", {

    expect_error(
        .assignFinalClonotypes(
            em_results = output3[,1:2],
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "em_results should have 4 columns"
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]][,1:2],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "Contigs should be in a contigs format"
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output3,
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "Orphan_contigs should have 3 columns and be made by "
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor2",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "Donor column name provided is not in the Seurat object"
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor2",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "Group column name provided is not in the Seurat object"
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, 1),
            Verbose = TRUE
        ),
        "Incorrect input to combineTcrArgs"
    )

    expect_error(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = 2
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_s4_class(
        .assignFinalClonotypes(
            em_results = output3,
            Contigs = output2[[1]],
            Seurat_object = output2[[2]],
            Orphan_contigs = output2[[3]],
            Donor = "donor",
            Group = "donor",
            CombineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            Verbose = TRUE
        ),
        "Seurat"
    )
})


output42 <- output4
output42@meta.data <- output42@meta.data[, 1:2]


test_that(".reformatSeuratObject works correctly", {

    expect_error(
        .reformatSeuratObject(
            Seurat_object = output42,
            Group = "donor",
            DetailedEclipseStats = FALSE,
            Verbose = TRUE
        ),
        "Seurat_object should be made by running"
    )

    expect_error(
        .reformatSeuratObject(
            Seurat_object = output4,
            Group = "donor2",
            DetailedEclipseStats = FALSE,
            Verbose = TRUE
        ),
        "Group column name provided is not in the Seurat object"
    )

    expect_error(
        .reformatSeuratObject(
            Seurat_object = output4,
            Group = "donor",
            DetailedEclipseStats = 2,
            Verbose = TRUE
        ),
        "TRUE or FALSE not provided to detailedEclipseStats"
    )

    expect_error(
        .reformatSeuratObject(
            Seurat_object = output4,
            Group = "donor",
            DetailedEclipseStats = FALSE,
            Verbose = 2
        ),
        "TRUE or FALSE not provided to verbose"
    )

    expect_s4_class(
        .reformatSeuratObject(
            Seurat_object = output4,
            Group = "donor",
            DetailedEclipseStats = FALSE,
            Verbose = TRUE
        ),
        "Seurat"
    )

})




### Prepping inputs for findThreeChainClones helper functions
contigs <- post_qc_contigs
chain <- "a"
mode = "standard"
printStats = TRUE

output1 <- .identifyThreeChainClones(Contigs = contigs,
                                     Chain = chain,
                                     Mode = mode)


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

contigs2 <- ECLIPSE(
    contig_file_paths = unit_test_tcr_paths,
    contig_type = "manual",
    seurat_object = unit_test_object,
    batch = "index",
    donor = "donor",
    group = "donor",
    original_barcode = "Cell_Barcode",
    detailedEclipseStats = FALSE,
    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
    allowDoubleChains = "both",
    verbose = TRUE,
    returnFilteredContigs = TRUE
)

output12 <- .identifyThreeChainClones(Contigs = contigs2,
                                     Chain = chain,
                                     Mode = mode)


test_that(".identifyThreeChainClones is working correctly", {

    expect_error(
        .identifyThreeChainClones(
            Contigs = post_qc_contigs[, 1:2],
            Chain = "a",
            Mode = "standard"
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_error(
        .identifyThreeChainClones(
            Contigs = post_qc_contigs,
            Chain = "c",
            Mode = "standard"
        ),
        "Chain should be a or b"
    )

    expect_error(
        .identifyThreeChainClones(
            Contigs = post_qc_contigs,
            Chain = "a",
            Mode = "fake"
        ),
        "Mode must be standard or show"
    )

    expect_true(
        is.list(
            .identifyThreeChainClones(
                Contigs = post_qc_contigs,
                Chain = "a",
                Mode = "standard"
            )
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_true(
        sum(colnames(
            .identifyThreeChainClones(
                Contigs = post_qc_contigs,
                Chain = "a",
                Mode = "standard")[[1]])
            %in% c("monotypic", "combo")) == 2
    )

    expect_true(
        sum(colnames(
            .identifyThreeChainClones(
                Contigs = post_qc_contigs,
                Chain = "a",
                Mode = "standard")[[2]])
            %in% c("barcode", "TCR1", "CTgene", "CTaa", "CTstrict")) == 5
    )

    expect_true(
        is.vector(
            .identifyThreeChainClones(
                Contigs = post_qc_contigs,
                Chain = "a",
                Mode = "standard")[[3]]
        )
    )

    expect_message(
        .identifyThreeChainClones(
            Contigs = post_qc_contigs,
            Chain = "b",
            Mode = "show"
        ),
        "Returning a list. The first table shows all"
    )
})





test_that(".modifyThreeChainContigs works correctly", {

    expect_error(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs[, 1:2],
            three_chain_clones = output1[[1]],
            Clones = output1[[2]],
            Chain = "a",
            Mode = "standard"
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_error(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs,
            three_chain_clones = output1[[1]],
            Clones = output1[[2]],
            Chain = "c",
            Mode = "standard"
        ),
        "Chain should be a or b"
    )

    expect_error(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs,
            three_chain_clones = output1[[1]],
            Clones = output1[[2]],
            Chain = "a",
            Mode = "fake"
        ),
        "Mode must be standard or show"
    )

    expect_error(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs,
            three_chain_clones = output1[[1]][, 1],
            Clones = output1[[2]],
            Chain = "a",
            Mode = "standard"
        ),
        "three_chain_clones should have 2 columns named combo and"
    )

    expect_error(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs,
            three_chain_clones = output1[[1]],
            Clones = output1[[2]][, 1:2],
            Chain = "a",
            Mode = "standard"
        ),
        "Clones should have 13 columns and it should be generated by"
    )

    expect_true(
        is.null(
            .modifyThreeChainContigs(
                Contigs = contigs2,
                three_chain_clones = output12[[1]],
                Clones = output12[[2]],
                Chain = "a",
                Mode = "show"
            )
        )
    )

    expect_s3_class(
        .modifyThreeChainContigs(
            Contigs = post_qc_contigs,
            three_chain_clones = output1[[1]],
            Clones = output1[[2]],
            Chain = "a",
            Mode = "standard"
        ),
        "tbl_df"
    )

})



test_that(".findMissingMonotypicChains works correctly", {

    expect_error(
        .findMissingMonotypicChains(Contigs = output2[, 1:2],
            specialTCRs = output1[[3]],
            Chain = "a",
            Mode = "standard",
            print_stats = TRUE
        ),
        "Contigs should be in the all_contig_annotations or"
    )

    expect_error(
        .findMissingMonotypicChains(Contigs = output2,
            specialTCRs = output1[[3]],
            Chain = "c",
            Mode = "standard",
            print_stats = TRUE
        ),
        "Chain should be a or b"
    )

    expect_error(
        .findMissingMonotypicChains(Contigs = output2,
            specialTCRs = output1[[3]],
            Chain = "a",
            Mode = "fake",
            print_stats = TRUE
        ),
        "Mode must be standard or show"
    )

    expect_error(
        .findMissingMonotypicChains(Contigs = output2,
            specialTCRs = list(output1[[3]]),
            Chain = "a",
            Mode = "standard",
            print_stats = TRUE
        ),
        "specialTCRs must be a character vector"
    )

    expect_error(
        .findMissingMonotypicChains(Contigs = output2,
            specialTCRs = output1[[3]],
            Chain = "a",
            Mode = "standard",
            print_stats = 2
        ),
        "printStats must be TRUE or FALSE"
    )

    expect_true(
        is.null(
            .findMissingMonotypicChains(Contigs = output2,
                specialTCRs = output1[[3]],
                Chain = "a",
                Mode = "show",
                print_stats = TRUE
            )
        )
    )

    expect_s3_class(
        .findMissingMonotypicChains(Contigs = output2,
            specialTCRs = output1[[3]],
            Chain = "a",
            Mode = "standard",
            print_stats = TRUE
        ),
        "tbl_df"
    )
})



test_sce <- as.SingleCellExperiment(test_seurat)



test_that("sceECLIPSE works correctly", {

    ### scECLIPSE doesn't accept Seurat objects
    expect_error(
        sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_seurat,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "sce_object is not a SingleCellExperiment object"
    )

    ### sceECLIPSE returns SCE objects
    expect_s4_class(
        sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE
        ),
        "SingleCellExperiment"
    )

    ### sceECLIPSE returns a tibble if returnFilteredContigs == TRUE
    expect_s3_class(
        sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = TRUE
        ),
        "tbl_df"
    )

   ### sceECLIPSE returns SCE objects with expected new columns in colData
   expect_true(
        sum(colnames(colData(sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = FALSE
        ))) %in% c("CTaa", "highConfChains", "cloneCallSource",
            "scR_CTaa", "cloneSize")) == 5
    )

   ### sceECLIPSE returns SCE objects that retain all original columns
   expect_true(
        sum(colnames(colData(test_sce)) %in% colnames(colData(sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = FALSE
        )))) == ncol(colData(test_sce))
    )

   ### Ordering of cells is correct and retained with sceECLIPSE
   expect_true(
        sum(rownames(colData(test_sce)) == rownames(colData(sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "donor",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = FALSE
        )))) == nrow(colData(test_sce))
    )

   ### Normal ECLIPSE errors still work with sceECLIPSE
   expect_error(
        sceECLIPSE(
            contig_file_paths = tcr_paths,
            contig_type = "manual",
            sce_object = test_sce,
            batch = "index",
            donor = "fake",
            group = "donor",
            original_barcode = "Cell_Barcode",
            detailedEclipseStats = FALSE,
            combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
            allowDoubleChains = "both",
            verbose = FALSE,
            returnFilteredContigs = FALSE
        ),
        "Donor column name provided is not in the Seurat object"
   )
})
