#' Plot edgeR dispersion from a DGEobj
#'
#' Creates an edgeR dispersion plot for RNA-Seq QC purposes.  Takes a counts matrix
#' or DGEList for input.  Dispersion is plotted against AveLogCPM.  Optionally,
#' the plot can instead be Biological Coefficient of Variation (BCV is the square root of
#' dispersion) against AveLogCPM.
#'
#' @param DGEdata Counts matrix or DGEList (Required)
#' @param designMatrix A design matrix created by stats::model.matrix (Required)
#' @param plotType Plot type must be canvasXpress or ggplot (Default to canvasXpress).
#' @param plotCategory One of "dispersion" or "BCV" (Default = "dispersion")
#' @param symbolSize (Default = 1)
#' @param symbolShape see
#'   http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ (Default = 1)
#' @param symbolColor Default = "darkblue"
#' @param symbolFill  Default = "darkblue"
#' @param symbolAlpha Transparency for the points. Value from 0 to 1. Smaller
#'   indicate more transparency (Default = 0.3)
#' @param linefitSize (Default = 1)
#' @param linefitColor (Default = "yellow")
#' @param lineFit (Default = NULL) Any type supported by geom_smooth (Loess is
#'   recommended).
#' @param lineType Any ggplot2 lineType supported by geom_smooth (Default = solid)
#' @param lineAlpha Alpha desired for geom_smooth line (Default = 1)
#' @param ... Extra parameters to pass to edgeR::estimateDisp
#'
#' @return canvasxpress or ggplot object based on plotType selection
#'
#' @examples
#' \dontrun{
#'    # canvasxpress plot
#'    myCxplot <- plotDispersion(myDGElist)
#'    myCxplot <- plotDispersion(myDGEobj)
#'
#'    # ggplot
#'    myGgplot <- plotDispersion(myDGElist, plotType = "ggplot")
#'    myGgplot <- plotDispersion(myDGEobj, plotType = "ggplot")
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom edgeR calcNormFactors estimateDisp DGEList
#' @importFrom canvasXpress canvasXpress
#'
#' @export
plotDispersion <- function(DGEdata,
                           designMatrix,
                           plotType = "canvasXpress",
                           plotCategory = "dispersion",
                           symbolSize = 1,
                           symbolShape = 1,
                           symbolColor = "darkblue",
                           symbolFill = "darkblue",
                           symbolAlpha = 0.3,
                           linefitSize = 1,
                           linefitColor = "yellow",
                           lineFit = NULL,
                           lineType = "solid",
                           lineAlpha = 1,
                           ...) {

    assertthat::assert_that(!missing(DGEdata),
                            msg = "DGEdata must be specified.")
    assertthat::assert_that(!missing(designMatrix),
                            "matrix" %in% class(designMatrix),
                            msg = "designMatrix must be specified and should be of class 'matrix'.")
    assertthat::assert_that(plotType %in% c("canvasXpress", "ggplot"),
                            msg = "Plot type must be either canvasXpress or ggplot.")
    assertthat::assert_that(tolower(plotCategory) %in% c("dispersion", "bcv"),
                            msg = "Plot value must be either dispersion or BCV.")

    if (class(DGEdata)[[1]] == "DGEList") {
        dgelist <- DGEdata %>%
            edgeR::calcNormFactors() %>%
            edgeR::estimateDisp(design = designMatrix, robust = TRUE, ...)
    } else {
        dgelist <- DGEdata %>%  # Process a counts matrix
            as.matrix %>%
            edgeR::DGEList() %>%
            edgeR::calcNormFactors() %>%
            edgeR::estimateDisp(design = designMatrix, robust = TRUE, ...)
    }

    if (tolower(plotCategory) == "dispersion") {
        plotdata <- data.frame(AveLogCPM = dgelist$AveLogCPM, Dispersion = dgelist$tagwise.dispersion)
        title <- "EdgeR Dispersion Plot"
        ylab  <- "Dispersion"
    } else {
        plotdata <- data.frame(AveLogCPM = dgelist$AveLogCPM, Dispersion = sqrt(dgelist$tagwise.dispersion))
        title <- "EdgeR BCV Plot"
        ylab  <- "BCV"
    }
    rownames(plotdata) <- rownames(dgelist$counts)

    if (plotType == "canvasXpress") {
        showLoessFit <- FALSE
        afterRender <- NULL
        if (!is.null(lineFit)) {
            lineFit <- cxSupportedLineFit(lineFit)
            if (lineFit == "lm") {
                afterRender <- list(list("addRegressionLine"))
            } else if (lineFit == "loess") {
                showLoessFit <- TRUE
            }
        }

        if (tolower(plotCategory) == "bcv") {
            colnames(plotdata) <- c("AveLogCPM", "BCV")
        }

        MyDispPlot <- canvasXpress::canvasXpress(data                    = plotdata,
                                                 graphType               = "Scatter2D",
                                                 colors                  = symbolColor,
                                                 sizes                   = symbolSize,
                                                 title                   = title,
                                                 yAxisTitle              = ylab,
                                                 showLoessFit            = showLoessFit,
                                                 fitLineColor            = linefitColor,
                                                 fitLineStyle            = lineType,
                                                 afterRender             = afterRender)
    } else {
        MyDispPlot <- ggplot(plotdata, aes(x = AveLogCPM, y = Dispersion)) +
            geom_point(size = symbolSize,
                       shape = symbolShape,
                       fill = symbolFill,
                       color = symbolColor,
                       alpha = symbolAlpha)

        if (!is.null(lineFit)) {
            MyDispPlot <- MyDispPlot +
                geom_smooth(method = lineFit,
                            size = linefitSize,
                            color = linefitColor,
                            linetype = lineType,
                            alpha = lineAlpha)
        }

        MyDispPlot <- MyDispPlot +
            ylab(ylab) +
            ggtitle(title) +
            expand_limits(x = 0, y = 0)
    }

    return(MyDispPlot)
}
