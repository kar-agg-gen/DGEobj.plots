#' Plot histogram analysis of p-value distributions
#'
#' Generate a facet plot (or optionally individual plots) from a dataframe of
#' numbers. Intended to perform histogram analysis of p-value distributions,
#' but should be useful for any dataframe of numeric columns.
#'
#' @param P.Val A matrix or dataframe of numeric data; col = samples
#' @param plotType Plot type must be canvasXpress or ggplot (Default to canvasXpress).
#' @param facet Set to FALSE to print individual plots instead of a faceted plot. (Default = TRUE)
#' @param binWidth Range is always 0-1 for p-values. (Default = 0.02)
#' @param alpha Set the transparency. (Default = 0.6)
#' @param color Color for the histogram outline. (Default = "dodgerblue3")
#' @param fill Fill color for the histogram. (Default = "dodgerblue3")
#'
#' @return A ggplot2 object if facet = TRUE or a list of plots if facet = FALSE. (Default = TRUE)
#'
#' @examples
#' \dontrun{
#'    # Print to console using all defaults
#'    MyPvalMatrix <- extractCol(getType(myDGEobj, "topTable"), "P.Value")
#'    plotPvalHist(MyPvalMatrix)
#'
#'    # Use some custom arguments
#'    myplot <- plotPvalHist(MyPValMatrix)
#' }
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom canvasXpress canvasXpress
#'
#' @export
plotPvalHist <- function(P.Val,
                         plotType = "canvasXpress",
                         facet = TRUE,
                         binWidth = 0.02,
                         alpha = 0.6,
                         color = "dodgerblue3",
                         fill = "dodgerblue3") {

    assertthat::assert_that(plotType %in% c("canvasXpress", "ggplot"),
                            msg = "Plot type must be either canvasXpress or ggplot.")

    if (is.matrix(P.Val)) {
        P.Val <- P.Val %>%
            as.data.frame
    }

    NumSamples <- ncol(P.Val)
    SampNames <- colnames(P.Val)

    # Set up Tall format
    P.Val$GeneID = rownames(P.Val)
    P.Val <-  stats::reshape(data          = P.Val,
                             idvar         = "GeneID",
                             varying       = SampNames,
                             v.names       = "Pval",
                             direction     = "long",
                             timevar       = "Levels",
                             times         = SampNames,
                             new.row.names = sequence(prod(length(SampNames), nrow(P.Val))))

    title <- "P-value Histograms"
    plotlist <- list()
    if (plotType == "canvasXpress") {
        if (facet) {
            cx.data <- subset(P.Val, select = Pval)
            var.annot <- subset(P.Val, select = -c(Pval))
            plotlist <- canvasXpress::canvasXpress(data                    = cx.data,
                                                   varAnnot                = var.annot,
                                                   histogramData           = TRUE,
                                                   graphType               = "Scatter2D",
                                                   histogramBins           = binWidth,
                                                   colors                  = color,
                                                   title                   = title,
                                                   xAxisTitle              = "P-value",
                                                   yAxisTitle              = "Count",
                                                   hideHistogram           = FALSE,
                                                   showHistogramDensity    = FALSE,
                                                   showLegend              = FALSE,
                                                   segregateVariablesBy    = list("Levels"))
        } else {
            for (i in 1:NumSamples) {
                s <- SampNames[i]
                MyPVal <- dplyr::filter(P.Val, grepl(s, Levels))
                cx.data <- subset(MyPVal, select = Pval)
                var.annot <- subset(MyPVal, select = -c(Pval))
                Hist_Pval <- canvasXpress::canvasXpress(data                    = cx.data,
                                                        varAnnot                = var.annot,
                                                        histogramData           = TRUE,
                                                        graphType               = "Scatter2D",
                                                        histogramBins           = binWidth,
                                                        colors                  = color,
                                                        title                   = paste(title, "\n", s),
                                                        xAxisTitle              = "P-value",
                                                        yAxisTitle              = "Count",
                                                        hideHistogram           = FALSE,
                                                        showHistogramDensity    = FALSE,
                                                        showLegend              = FALSE,
                                                        segregateVariablesBy    = list("Levels"))

                plotlist[[i]] = Hist_Pval
            }
        }
    } else {
        if (facet) {
            numcol <- 3
            numrow <- (NumSamples / numcol) %>% ceiling

            plotlist <- ggplot2::ggplot(data = P.Val, aes(x = Pval)) +
                ggplot2::geom_histogram(alpha = alpha,
                                        fill = fill,
                                        color = color,
                                        binwidth = binWidth) +
                ggplot2::xlab("P-value") +
                ggplot2::ylab("Count") +
                ggtitle(title) +
                ggplot2::scale_fill_brewer(palette = "Set1") +
                ggplot2::facet_wrap(~Levels, nrow = numrow, scales = "free")
        } else {
            for (i in 1:NumSamples) {
                s <- SampNames[i]
                MyPVal <- dplyr::filter(P.Val, grepl(s, Levels))

                Hist_Pval <- ggplot2::ggplot(data = MyPVal, aes(x = Pval)) +
                    ggplot2::geom_histogram(alpha = alpha,
                                            fill = fill,
                                            color = color,
                                            binwidth = binWidth) +
                    ggplot2::xlab("P-value") +
                    ggplot2::ylab("Count") +
                    ggplot2::ggtitle(paste(title, "\n", s))

                plotlist[[i]] = Hist_Pval
            }
        }
    }
    return(plotlist)
}
