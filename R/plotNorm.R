#' Plot log2CPM before and after normalization
#'
#' Takes a DGEobj containing counts or a counts matrix as input. Returns a ggplot object containing
#' a faceted plot of log2CPM before and after normalization. Either a box plot or density plot
#' type can be chosen.
#'
#' Normalization is performed by edgeR::calcNormFactors. Note TMM is specifically tailored to count-based
#' data.  Thus this function is only appropriate for count-based data.
#'
#' @param DGEdata  A DGEobj or counts matrix.
#' @param plotType  One of "box" or "density." (Default = "box")
#' @param normalize Default = "TMM" and invokes TMM normalization. Other allowed
#'   values are: "RLE", "upperquartile", or "none". Invokes edgeR::calcNormFactors for
#'   normalization.
#'
#' @return A faceted ggplot plot showing before/after log2CPM normalization.
#'
#' @examples
#' \dontrun{
#'    myNormPlotBox <- plotNorm(myDGEobj, plotType = "box")
#'    myNormPlotDensity <- plotNorm(counts, plotType = "density")
#' }
#'
#' @import magrittr ggplot2
#' @importFrom stringr str_c
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom DGEobj getItem
#' @importFrom assertthat assert_that
#'
#' @export
plotNorm <- function(DGEdata,
                     plotType = "box",
                     normalize = "tmm") {

    assertthat::assert_that(any(c("matrix", "DGEobj") %in% class(DGEdata)),
                            msg = "DGEdata must be of either class 'matrix' or 'DGEobj'.")
    assertthat::assert_that(tolower(plotType) %in% c("box", "density"),
                            msg = "plotType must be one of 'box' or 'density'.")
    assertthat::assert_that(tolower(normalize) %in% c("tmm", "rle", "upperquartile", "none"),
                            msg = "normalize must be one of 'TMM', 'RLE', 'upperquartile', or 'none'.")

    if ("matrix" %in% class(DGEdata)) {
        counts <- DGEdata
    } else {
        counts <- DGEobj::getItem(DGEdata, "counts")
    }

    log2cpm <- DGEobj.utils::convertCounts(counts, unit = "cpm", log = TRUE, normalize = "none") %>%
        as.data.frame
    log2CPM_tmm <- DGEobj.utils::convertCounts(counts, unit = "cpm", log = TRUE, normalize = normalize) %>%
        as.data.frame

    log2cpm_colnames <- colnames(log2cpm)
    tall <- data.frame("GeneID" = row.names(log2cpm), log2cpm, row.names = NULL)

    tall <-  stats::reshape(data          = tall,
                            idvar         = "GeneID",
                            varying       = log2cpm_colnames,
                            v.names       = "Log2CPM",
                            direction     = "long",
                            timevar       = "SampleID",
                            times         = log2cpm_colnames,
                            new.row.names = sequence(prod(length(log2cpm_colnames), nrow(tall))))
    tall$Normalization = "none"

    log2cpm_tmm_colnames <- colnames(log2CPM_tmm)
    tall_tmm <- data.frame("GeneID" = row.names(log2CPM_tmm), log2CPM_tmm, row.names = NULL)
    tall_tmm <-  stats::reshape(data          = tall_tmm,
                                idvar         = "GeneID",
                                varying       = log2cpm_tmm_colnames,
                                v.names       = "Log2CPM",
                                direction     = "long",
                                timevar       = "SampleID",
                                times         = log2cpm_tmm_colnames,
                                new.row.names = sequence(prod(length(log2cpm_tmm_colnames), nrow(tall_tmm))))
    tall_tmm$Normalization = toupper(normalize)

    tall <- tall %>%
        rbind(tall_tmm)

    if (tolower(plotType) == "density") {
        resultPlot <- ggplot(tall, aes(x = Log2CPM, color = SampleID)) +
            geom_density() +
            facet_grid(~Normalization) +
            ggtitle(stringr::str_c("Log2CPM before/after", normalize, "normalization", sep = " "))  +
            theme_gray() +
            theme(legend.position = "none")
    }

    if (tolower(plotType == "box")) {
        resultPlot <- ggplot(tall, aes(x = SampleID, y = Log2CPM, color = SampleID)) +
            geom_boxplot(alpha = 0.5) +
            facet_grid(~Normalization) +
            ggtitle(stringr::str_c("Log2CPM before/after", normalize, "normalization", sep = " "))  +
            theme_gray() +
            theme(axis.text.x = element_blank(),
                  legend.position = "none")
    }

    return(resultPlot)
}
