#' Create QC metric plots
#'
#' Takes a dataframe with metric names in the first column and
#' samples in the second through nth column.  Each row is a different QC metric.
#' Returns a list of QC plots as defined by other arguments.  If only one metric is plotted,
#' a ggplot object is returned. By default, horizontal reference lines are drawn
#' at the median and +/- n SDs based on the hlineSD argument. These are
#' statistical reference points, NOT pass/fail limits.
#'
#' @param qcdata A dataframe or tibble with metric names in the first column and
#'   samples in columns 2-n.  Each row is a different QC metric. (Required)
#' @param metricNames A list of metrics to plot.  Values must exist in column 1
#'   of the data frame. (Required)
#' @param sampleNames By default will use the sample names in qcdata (colnames).
#'   Optionally use this argument to provide different sample names.  Supply the
#'   alternative sample names here in the order they appear in the qcdata
#'   data.frame (columns 2-n). The plot order is based on an alphabetical sort
#'   of the original column names so might be different from the order in the
#'   supplied data.frame.
#' @param plotType One of "bar", "point", "pointline".  For a different
#'   plot type for each metric, pass a list of plotTypes with length equal to
#'   length(metricNames). (Default = "bar")
#' @param barColor Color for the bar outline (Default = "dodgerblue4")
#' @param barFill Color for the bar area (Default = "dodgerblue3")
#' @param barSize Set the bar size (thickness of each bar perimeter; Default =
#'   0.1)
#' @param barWidth Set the bar width (Default = 0.8)
#' @param barAlpha Transparency for the bar layer (Range = 0-1) (Default = 1)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Range = 0-1) (Default = 1)
#' @param pointSize Size of the points (Default = 4)
#' @param lineColor Color of the line (Default = "dodgerblue4")
#' @param lineSize Size of the line fit (Default = 1)
#' @param lineAlpha Transparency for the line layer (Range = 0-1) (Default=1)
#' @param lineType One of c("solid", "dashed", "dotted", "dotdash", "longdash",
#'   "twodash"). (Default = "solid")
#' @param histColor Outline color for the histogram (Default = "dodgerblue4")
#' @param histFill Fill color for the histogram (Default = "dodgerblue3")
#' @param histSize Thickness of the bar borders (Default = 1)
#' @param histAlpha Transparency of the histogram (Default = 1)
#' @param xAngle Angle to set the sample labels on the X axis (Default = 90;
#'   Range = 0 - 90)
#' @param baseTextSize Default = 14
#' @param hlineSD Draw two reference lines 1) at the median value 2) the number of
#'   SDs defined by the value of hlineSD. (Default = 3; 0 to disable the reference lines.)
#' @param winsorize This implements a robust method to calculate standard
#'   deviations.  It is used to calculate the standard deviation for the
#'   placement of horizontal reference lines (hlineSD argument).  The adaptive
#'   winsorization used here only trims extreme values when normality is
#'   violated. see https://www.r-bloggers.com/winsorization/ for details.
#'   (Default = TRUE).
#'
#' @return ggplot object if one plot is specified.  A list of ggplot objects if 2 or more metrics specified.
#'
#' @examples
#' \dontrun{
#'   someFaveMetrics <- c("Alignment_MappedRate", "Alignment_PairedRate",
#'                               "Source_rRNA", "Strand_Read1AntiSense",
#'                               "Strand_ReadPairAntiSense", "Profile_ExonRate",
#'                               "Profile_InterGene_FPK")
#'
#'   MyQCplots <- QCplots(qcdata, metricNames = someFaveMetrics) # All defaults
#'
#'   # Draw the first plot
#'   print(MyQCplots[[1]])
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#' @importFrom stats median sd
#'
#' @export
QCplots <- function(qcdata,
                    metricNames,
                    sampleNames,
                    plotType = "bar",
                    barColor = "dodgerblue4",
                    barFill = "dodgerblue3",
                    barSize = 0.1,
                    barAlpha = 1,
                    barWidth = 0.9,
                    pointColor = "dodgerblue4",
                    pointFill = "dodgerblue3",
                    pointShape = 21,
                    pointAlpha = 1,
                    pointSize = 4,
                    lineColor = "dodgerblue4",
                    lineSize = 1,
                    lineType = "solid",
                    lineAlpha = 1,
                    histColor = "dodgerblue4",
                    histFill = "dodgerblue3",
                    histSize = 1,
                    histAlpha = 1,
                    xAngle = 90,
                    baseTextSize = 14,
                    hlineSD = 3,
                    winsorize = TRUE) {

    assertthat::assert_that("data.frame" %in% class(qcdata),
                            msg = "qcdata must be of class 'data.frame'.")
    assertthat::assert_that(tolower(plotType) %in% c("bar", "point", "pointline", "histogram"),
                            msg = "plotType must be one of 'bar', 'point', 'pointline', or 'histogram'.")
    assertthat::assert_that(xAngle >= 0  && xAngle <= 90,
                            msg = "xAngle must be greater than or equal to 0 but also less than or equal to 90.")

    # Replace "." in column 2-n with NA and convert columns 2-n to numeric
    dotIdx <- qcdata == "."
    if (sum(dotIdx, na.rm = TRUE) > 0) {
        qcdata[dotIdx] <- NA
        for (columnName in colnames(qcdata)[2:ncol(qcdata)])
            qcdata[columnName] <- as.numeric(qcdata[[columnName]])
    }

    # Convert first col to rownames and transpose
    qcdata %<>% column_to_rownames(var = colnames(qcdata)[1]) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "Sample")

    assertthat::assert_that(all(metricNames %in% colnames(qcdata)),
                            msg = "All of the specified metricNames must be included in the colnames of qcdata.")

    # If only one plotType, apply to all plots
    if (length(plotType) == 1) {
        plotType <- rep(plotType, length(metricNames))
    }

    plots <- list()

    # Helper function
    winsorize <- function(x, multiple = 3, na.rm = FALSE, ...){

        med <- median(x, na.rm = na.rm)
        y <- x - med # Median centered
        sc <- stats::mad(y, center = 0, na.rm = na.rm, ...) * multiple
        y[ y > sc ] <- sc
        y[ y < -sc ] <- -sc
        y + med
    }

    for (metric in metricNames) {
        idx <- metric == metricNames
        plot_type <- plotType[idx]

        # Calculate mean and sd for hline yintercepts
        if (winsorize == TRUE) {
            thisMetric <- winsorize(qcdata[,metric], na.rm = TRUE)
        } else {
            thisMetric <- qcdata[, metric, drop = TRUE]
        }
        metricMedian <- median(thisMetric, na.rm = TRUE)
        metricMean   <- mean(thisMetric, na.rm = TRUE)
        metricSD     <- sd(thisMetric, na.rm = TRUE)

        if (!tolower(plot_type) == "histogram") {
            p <- ggplot(qcdata, aes_string(x = "Sample", y = metric, group = 1))
        }

        p <- switch(tolower(plot_type),
                    bar = {p + geom_bar(stat = "identity", color = barColor, fill = barFill, alpha = barAlpha, width = barWidth)},
                    point = {p + geom_point(color = pointColor, fill = pointFill, shape = pointShape, alpha = pointAlpha, size = pointSize)},
                    pointline = {p +
                            geom_point(color = pointColor, fill = pointFill, shape = pointShape, alpha = pointAlpha, size = pointSize) +
                            geom_line(color = lineColor, size = lineSize, linetype = lineType, alpha = lineAlpha)},
                    histogram = {p <- ggplot(qcdata, aes_string(x = metric)) +
                        geom_histogram(colour = histColor, fill = histFill, size = histSize, alpha = histAlpha)}
        )

        # Draw hline xSD above or below the mean
        SD <- metricSD * hlineSD
        if (hlineSD > 0 & tolower(plot_type) == "histogram") {

            # Plot vlines for the histogram
            p <- p +
                geom_vline(xintercept = metricMedian, color = "grey70") +
                geom_vline(xintercept = metricMean + SD, color = "firebrick3", linetype = "longdash") +
                geom_vline(xintercept = metricMean - SD, color = "firebrick3", linetype = "longdash")
        } else if (hlineSD > 0) {  # Use hlines for other plot types
            p <- p +
                geom_hline(yintercept = metricMedian, color = "grey70") +
                geom_hline(yintercept = metricMean + SD, color = "firebrick3", linetype = "longdash")
            # Most qc plots floor to 0 so no point plotting -SD if it goes below zero.
            if (metricMean - SD > 0) {
                p <- p + geom_hline(yintercept = metricMean - SD, color = "firebrick3", linetype = "longdash")
            }
        }

        # Set x axis text angle
        hjust = 1
        vjust = 1
        if (xAngle == 0) {
            hjust = 0.5
        }
        if (xAngle == 90) {
            vjust = 0.5
        }
        p <- p +
            ggtitle(metric) +
            theme_gray(baseTextSize) +
            theme(axis.text.x = element_text(angle = xAngle, hjust = hjust, vjust = vjust))

        # Set user supplied axis tick labels
        if (!missing(sampleNames)) {
            if (length(sampleNames) == nrow(qcdata)) {
                p <- p + scale_x_discrete(breaks = qcdata$Sample,
                                          labels = sampleNames)
            }
        }

        plots[[metric]] <- p
    }

    # Put the warning outside the loop so it only gets called once.
    if (!missing(sampleNames)) {
        if (length(sampleNames) != nrow(qcdata)) {
            warning("sampleNames not the right length, thus ignored")
        }
    }

    if (length(plots) == 1) {
        plots <- plots[[1]]
    }

    return(plots)
}
