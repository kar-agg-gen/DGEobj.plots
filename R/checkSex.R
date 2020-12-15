#' Plot XIST vs. highest expressed Y-chromosome gene
#'
#' Take a DGEobj as input and plot expression of XIST vs the highest expressed
#' Y-chromosome gene.
#'
#' This function uses the original unfiltered data by default, because Y-linked genes
#' are often below the low intensity threshold in tissues other than testes.
#' Nevertheless, there are usually enough reads for the plot to be interpretable.
#'
#' @param dgeObj A DGEobj (Required)
#' @param species One of "human", "mouse", "rat"
#' @param sexCol Character string name of the sex annotation column in the
#'    design table within the DGEobj (Optional)
#' @param labelCol Character string name of the design column to use to label points
#'    with ggrepel (optional if showLabels = FALSE)
#' @param showLabels Set TRUE to turn on ggrepel labels (Default = FALSE)
#' @param chrX Character string name of the X chromosome in the gene annotation
#'    within the geneData object (Default = "X")
#' @param chrY Character string name of the Y chromosome in the gene annotation
#'    within the geneData object (Default = "Y")
#' @param baseFontSize Sets the base font size for ggplot themes (Default = 14)
#' @param orig Set to FALSE to use filtered DGEobj data. Set to TRUE to use original
#'    unfiltered data (Default = TRUE)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#'    checkSex(DGEobj, species = "human")
#' }
#'
#' @import ggplot2 magrittr ggrepel
#' @importFrom assertthat assert_that
#' @importFrom edgeR calcNormFactors
#' @importFrom dplyr filter select left_join mutate arrange
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
checkSex <- function(dgeObj,
                     species,
                     sexCol,
                     labelCol,
                     showLabels = FALSE,
                     chrX = "X",
                     chrY = "Y",
                     baseFontSize = 14,
                     orig = TRUE) {

    assertthat::assert_that(!missing(species),
                            !missing(dgeObj),
                            msg = "Specify both a DGEob and a species. Both are required.")
    assertthat::assert_that(tolower(species) %in% c("human", "mouse", "rat"),
                            msg = "Species must be one of 'human', 'mouse', or 'rat'.")

    if (tolower(species) == "rat") { # No XIST gene in Rat
        x <- .getTopExpressedGene(dgeObj, chr = chrX, orig = orig)
    } else {
        x <- switch(tolower(species),
                    human = list(gene = "ENSG00000229807", genename = "XIST"),
                    mouse = list(gene = "ENSMUSG00000086503", genename = "Xist")
        )
    }

    y <- .getTopExpressedGene(dgeObj, chr = chrY, orig = orig)

    # Get data for the two x and y genes
    if (orig) {
        log2CPM <- DGEobj.utils::convertCounts(DGEobj::getItem(dgeObj, "counts_orig"), unit = "cpm", log = TRUE, normalize = "tmm")
    } else {
        log2CPM <- DGEobj.utils::convertCounts(DGEobj::getItem(dgeObj, "counts"), unit = "cpm", log = TRUE, normalize = "tmm")
    }

    idx <- rownames(log2CPM) %in% c(x$gene, y$gene)
    plotDat <- t(log2CPM[idx,]) %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname")

    # Add sample identifier data
    dtemp <- DGEobj::getItem(dgeObj, "design_orig") %>%
        tibble::rownames_to_column(var = "rowname")

    # Add labelCol and sexCol data as needed
    if (!missing(labelCol) & !missing(sexCol)) {

        dtemp %<>% dplyr::select(rowname, labelCol = labelCol, sexCol = sexCol)
        plotDat <- dplyr::left_join(plotDat, dtemp, by = "rowname")
        colnames(plotDat) <- c("SampID", x$genename, y$genename, labelCol, sexCol)
        sexPlot <- ggplot(plotDat, aes_string(x = x$genename, y = y$genename, label = labelCol, color = sexCol))

    } else if (!missing(labelCol) & missing(sexCol)) {

        dtemp %<>% dplyr::select(rowname, label = labelCol)
        plotDat <- dplyr::left_join(plotDat, dtemp, by = "rowname")
        colnames(plotDat) <- c("SampID", x$genename, y$genename, labelCol)
        sexPlot <- ggplot(plotDat, aes_string(x = x$genename, y = y$genename, label = labelCol))

    } else if (missing(labelCol) & !missing(sexCol)) {

        dtemp %<>% dplyr::select(rowname, sexCol = sexCol)
        plotDat <- dplyr::left_join(plotDat, dtemp, by = "rowname")
        colnames(plotDat) <- c("SampID", x$genename, y$genename, sexCol)
        sexPlot <- ggplot(plotDat, aes_string(x = x$genename, y = y$genename, color = sexCol))

    } else if (missing(labelCol) & missing(sexCol)) {
        colnames(plotDat) <- c("SampID", x$genename, y$genename)
        sexPlot <- ggplot(plotDat, aes_string(x = x$genename, y = y$genename))
    }

    sexPlot <- sexPlot +
        geom_point(size = 4, shape = 21) +
        xlab(paste("X (", x$genename, ")", sep = "")) +
        ylab(paste("Y (", y$genename, ")", sep = "")) +
        ggtitle("Sex Determination Plot") +
        theme_grey()
    if (showLabels == TRUE & !missing(labelCol)) {
        sexPlot <- sexPlot + ggrepel::geom_label_repel(label.size = 0.125)
    }

    return(sexPlot)
}


# Helper function
.getTopExpressedGene <- function(dgeObj,
                                 chr,
                                 orig = FALSE){
    # Get top mean expressed gene from given chromosome
    if (orig == TRUE) {
        geneData <- DGEobj::getItem(dgeObj, "geneData_orig")
    } else {
        geneData <- DGEobj::getItem(dgeObj, "geneData")
    }

    # Support function for Sex analysis
    .getChrDat <- function(dgeObj,
                           chr,
                           orig = FALSE) {
        # chr can be one or more chromosome names
        # Returns NULL if none found
        if (orig == TRUE) {
            geneData <- DGEobj::getItem(dgeObj, "geneData_orig")
            log2CPM <- DGEobj.utils::convertCounts(DGEobj::getItem(dgeObj, "counts_orig"), unit = "cpm", log = TRUE, normalize = "tmm")
        } else {
            geneData <- DGEobj::getItem(dgeObj, "geneData")
            log2CPM <- DGEobj.utils::convertCounts(DGEobj::getItem(dgeObj, "counts"), unit = "cpm", log = TRUE, normalize = "tmm")
        }

        # Filter to specified chr
        if (!is.null(chr)) {
            Chromosome <- NULL #binding variable Chromosome local to function
            geneData %<>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "geneid") %>%
                dplyr::filter(toupper(Chromosome) %in% toupper(chr)) %>%
                tibble::column_to_rownames(var = "geneid")
        }

        log2CPM_mat <- log2CPM[rownames(log2CPM) %in% rownames(geneData),]

        if (is.null(nrow(log2CPM_mat))) {
            return(NULL)
        } else {
            return(log2CPM_mat)
        }
    }

    # Filter geneData to selected chr
    if (!is.null(chr)) {
        log2CPM_mat <- .getChrDat(dgeObj, chr = chr, orig = orig)
        assertthat::assert_that(!is.null(log2CPM_mat),
                            msg = "No gene data found for specified chromosome.")
    }

    meanLogCPM <- rowMeans(log2CPM_mat)

    # Get the top expressed gene data
    # Calc mean and sort on mean descending
    log2CPM  <- log2CPM_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "geneid") %>%
        dplyr::mutate(meanLogCPM = meanLogCPM) %>%
        dplyr::arrange(dplyr::desc(meanLogCPM)) %>%
        dplyr::mutate(meanLogCPM = NULL) %>%
        tibble::column_to_rownames(var = "geneid")
    gene <- rownames(log2CPM)[1]
    genename <- geneData$GeneName[rownames(geneData) == gene]

    return(list(gene = gene, genename = genename))
}
