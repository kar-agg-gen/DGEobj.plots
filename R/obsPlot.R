#' Function obsPlot
#'
#' Provides a summary plot for each observation (gene), showing data for each
#' experiment group. The plot can optionally include one or more of the
#' following layers: boxplot, violin plot, individual points, and/or mean of all
#' points.  The layers are built up in the order listed with user settable
#' transparency, colors etc.  By default, the Boxplot, Point, and Mean layers are
#' active. Also, by default, the plots are faceted.  Faceted plot can be turned
#' off to return a list of individual ggplot2 graphics for each gene.
#'
#' Input is a dataframe or matrix of observations (rows; usually genes) by
#' Samples (columns) and requires rownames to identify the observations (genes).
#' A Block vector is required to define which Samples (columns) belong to the
#' same group and generate a summary plot by group with a separate plot for
#' each observation (gene).
#'
#' Rownames are required and all the data columns should be numeric.  Each
#' observation (gene) generates a separate plot, so you should pass a smallish
#' list of genes unless you want a lot of output.  By default the plot is
#' faceted. You should turn facet off if you have more than ~25 genes to plot.
#' Colnames and Rownames will be used to label the plots by default but custom
#' col/row names can also be supplied as separate vectors.
#'
#' @param data Matrix or dataframe of whatever data you want to plot with
#'   samples in columns and observations (genes) in rows (with rownames and colnames).  Bug: Currently requires at least 2 rows of data.
#' @param block Character or numeric vector. Determines which samples belong to the same group.  Must be same
#'   length as ncol(data).  Assign the same value to each member of a group.
#'   Note that the block values are used to label the X Axis to identify the
#'   groups.  Thus short mnemonic labels are useful here.
#' @param blockOrder Character vector. Use this to explicitly set the display order of the blocks
#'   in the plots.  blockOrder should contain unique(block) in the order you
#'   wish to items to be arranged along the X-axis.  (Optional: Default = unique(block))
#' @param obsNames Character vector of row (observation) names to use instead of actual
#'   rownames. (Default = NULL)
#' @param sampNames Character vector of column (sample) names to use instead of actual
#'   colnames. (Default = NULL)
#' @param plotBy Name for the type of observation being plotted (Default = "Gene")
#' @param valType Name for the Value being plotted.  (Default = "Log2CPM")
#' @param boxLayer Adds a boxplot layer (Default = TRUE)
#' @param violinLayer Adds a violin layer (Default = FALSE)
#' @param pointLayer Adds a point layer (Default = TRUE)
#' @param meanLayer Adds a mean layer (Default = TRUE)
#' @param xlab X axis label (defaults to "Block")
#' @param ylab Y axis label (defaults to valType)
#' @param title Plot title (optional; Defaults = NULL)
#' @param boxColor Color for the boxplot layer (Default = "grey30")
#' @param boxFill Fill Color for the boxplot layer (Default = "deepskyblue3")
#' @param boxAlpha Transparency for the box layer (Default = 0.5)
#' @param violinColor Color for the violin layer (Default = "grey30")
#' @param violinFill Fill Color for the violin layer (Default = "yellow")
#' @param violinAlpha Alpha value for the violin layer (Default = 0.5)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Default = 1)
#' @param boxNotch Turn on/off confidence interval notches on boxplots (Default = FALSE)
#' @param boxNotchWidth Set the width of box notches (0-1) (Default = 0.8)
#' @param pointSize Size of the points (Default = 4)
#' @param pointJitter Amount to jitter the points (Default = 0) Try 0.2 if you
#'   have a lot of points.
#' @param meanColor Color for the mean layer (Default = "red2")
#' @param meanFill Fill color for the mean layer (Default = "goldenrod1")
#' @param meanShape Shape for the mean layer (Default = 21; fillable circle)
#' @param meanAlpha Transparency for the mean layer (Default = 0.7)
#' @param meanSize Size of the mean points (Default = 3)
#' @param baseFontSize The smallest size font in the figure in points. (Default =
#'   12)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. Default = bw"
#' @param facet Specifies whether to facet (TRUE) or print individual plots
#'   (FALSE)  (Default = TRUE)
#' @param facetCol If facet = TRUE, used to specify a number of columns of faceted graphs.
#'   (Default = NULL)
#' @param xAngle Angle to set the sample labels on the X axis. (Default =  30; Range = 0-90)
#' @param scales Specify same scales or independent scales for each subplot (Default = "free_y";
#'   Allowed values: "fixed", "free_x", "free_y", "free")
#' @param returnPlotDat Returns the dataframe used for the plot as a list member (Default = FALSE)
#'
#' @return ggplot object. If Facet = TRUE (default) returns a faceted ggplot object. If
#'   facet = FALSE or returnPlotDat = TRUE, returns a list of ggplot objects indexed
#'   by observation (gene) names. If returnPlotDat = TRUE, the last element
#'   of the list is the dataframe used to generate the plot.
#'
#' @examples
#' \dontrun{
#'    # Simple faceted plot with custom title
#'    # Get Log2CPM from a DGEobj object
#'    dgelist <- getItem(DGEobj, "DGEList")
#'    Log2CPM <- edgeR::cpm(dgelist, log = TRUE)
#'    genes <- Log2CPM[1:12,] #first dozen genes
#'    # Define six treatment groups in triplicate
#'    MyBlock <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6)
#'    MyPlot <- obsPlot(genes, MyBlock, title = "Plot Title")
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom assertthat assert_that
#'
#' @export
obsPlot <- function(data,
                    block,
                    blockOrder,
                    obsNames,
                    sampNames,
                    plotBy = "Gene",
                    valType = "Log2CPM",
                    boxLayer = TRUE,
                    violinLayer = FALSE,
                    pointLayer = TRUE,
                    meanLayer = TRUE,
                    xlab = NULL, ylab = NULL, title = NULL,
                    boxColor = "grey30",
                    boxFill = "deepskyblue3",
                    boxAlpha = 0.5,
                    boxNotch = FALSE,
                    boxNotchWidth = 0.8,
                    violinColor = "grey30",
                    violinFill = "goldenrod1",
                    violinAlpha = 0.5,
                    pointColor = "grey30",
                    pointFill = "dodgerblue4",
                    pointShape = 21,
                    pointAlpha = 1,
                    pointSize = 2,
                    pointJitter = 0,
                    meanColor = "red2",
                    meanFill = "goldenrod1",
                    meanShape = 22,
                    meanAlpha = 0.7,
                    meanSize = 3,
                    baseFontSize = 12,
                    themeStyle = "grey",
                    facet = TRUE,
                    facetCol = NULL,
                    xAngle = 30,
                    scales = "free_y",
                    returnPlotDat = FALSE) {

    .addGeoms <- function(MyPlot) {
        if (boxLayer == TRUE) {
            MyPlot <- MyPlot + geom_boxplot(alpha = boxAlpha,
                                            color = boxColor,
                                            fill = boxFill,
                                            notch = boxNotch,
                                            notchwidth = boxNotchWidth,
                                            outlier.shape = outlier.shape,
                                            outlier.size = outlier.size)
        }

        if (violinLayer == TRUE) {
            MyPlot <- MyPlot + geom_violin(alpha = violinAlpha,
                                           color = violinColor,
                                           fill = violinFill)
        }

        if (pointLayer == TRUE) {
            if (pointJitter > 0) {
                MyPlot <- MyPlot + geom_point(position = position_jitter(width = pointJitter),
                                              alpha = pointAlpha,
                                              color = pointColor,
                                              fill = pointFill,
                                              size = pointSize,
                                              shape = pointShape)
            } else {
                MyPlot <- MyPlot + geom_point(alpha = pointAlpha,
                                              color = pointColor,
                                              fill = pointFill,
                                              size = pointSize,
                                              shape = pointShape)
            }
        }

        if (meanLayer == TRUE) {
            MyPlot <- MyPlot +
                stat_summary(fun.y = mean,
                             geom = "point",
                             shape = meanShape,
                             size = meanSize,
                             color = "red",
                             fill = "goldenrod1",
                             alpha = meanAlpha)
        }

        return(MyPlot)
    }

    # Argument checks
    assertthat::assert_that(!missing(block),
                            !missing(data),
                            msg = "block and data are both required and must be specified.")
    if (!class(block) == "character") {
        block <- as.character(block)
    }

    if (is.matrix(data)) {
        data <- as.data.frame(data, stringsAsFactors = FALSE)
    }

    assertthat::assert_that(length(block) == ncol(data),
                            msg = "block should be the same length as the number of columns in data.")

    if (missing(sampNames)) {
        sampNames <- colnames(data)
    } else {
        assertthat::assert_that(length(sampNames) == ncol(data),
                                msg = "sampNames should be the same length as the number of columns in data.")
    }

    if (missing(obsNames)) {
        obsNames <- rownames(data)
    } else {
        assertthat::assert_that(length(obsNames) == nrow(data),
                                msg = "obsNames should be the same length as the number of rows in data.")
    }

    if (missing(blockOrder)) {
        blockOrder <- unique(block)
    } else {
        assertthat::assert_that(length(blockOrder) == length(unique(block)),
                                msg = "blockOrder should be the same length as the number of unique entries in block.")
    }

    # Reduce box outliers to a dot if geom_points turned on.
    outlier.size <- 1.5
    outlier.shape <- 19
    if (pointLayer) {
        outlier.size <- 1
        outlier.shape <- "."
    }
    # Get axis labels if not provided
    if (is.null(xlab)) {
        xlab <- "Group"
    }

    if (is.null(ylab)) {
        ylab <- valType
    }

    # Build tall data
    groupdf <- data.frame(cbind(Samples = sampNames, Block = block), stringsAsFactors = FALSE)
    data[[plotBy]] = as.character(obsNames)
    data %<>% tidyr::gather(key = "Samples", value = "valType", -Gene, na.rm = TRUE)
    # Attach the group info
    data %<>% dplyr::left_join(groupdf)
    data$Block <- factor(data$Block, levels = blockOrder)

    # Plot code
    if (facet) {
        # Set facet columns to sqrt of unique observations (rounded up)
        if (is.null(facetCol)) {
            numcol <- obsNames %>% unique %>% length %>% sqrt %>% ceiling
        } else {
            numcol <- facetCol
        }

        if (numcol > 6) {
            warning("You're putting a lot of plots into a Facet Plot")
        }

        MyPlot <- ggplot2::ggplot(data, aes_string(x = "Block", y = valType))
        MyPlot <- .addGeoms(MyPlot)
        MyPlot <- MyPlot + ggplot2::facet_wrap(~Gene, nrow = numcol, scales = scales)

        MyPlot <- MyPlot + ggplot2::xlab(xlab)
        MyPlot <- MyPlot + ggplot2::ylab(ylab)
        MyPlot <- MyPlot + ggplot2::ggtitle(title)
        if (tolower(themeStyle) == "bw" ) {
            MyPlot <- MyPlot + theme_bw() + baseTheme(baseFontSize)
        } else {
            MyPlot <- MyPlot + theme_grey() + baseTheme(baseFontSize)
        }

        # Rotate xaxis group labels
        if (xAngle > 0) {
            MyPlot <- MyPlot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
        }

    } else {

        plotlist <- list()

        for (obs in data[[plotBy]]) {
            dat <- data[data[[plotBy]] == obs, ]
            aplot <- ggplot(dat, aes_string(x = "Block", y = valType)) + # Samples vs Log2CPM
                xlab(xlab) +
                ylab(ylab) +
                ggtitle(obs) +
                theme_grey() +
                facetTheme(baseFontSize)
            aplot <- .addGeoms(aplot)
            if (meanLayer == TRUE) {
                aplot <- aplot +
                    stat_summary(fun.y = mean,
                                 geom = "point",
                                 shape = meanShape,
                                 size = meanSize,
                                 color = "red",
                                 fill = "goldenrod1",
                                 alpha = 1.0)
            }
            # Rotate xaxis group labels
            if (xAngle > 0) {
                aplot <- aplot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
            }
            plotlist[[obs]] <- aplot
        }

        MyPlot = plotlist
    }

    # Add dataframe as specified by argument returnPlotDat
    if (returnPlotDat == TRUE) {
        if ("list" %in% class(MyPlot)) {
            MyPlot$PlotDat <- data
        } else {
            MyPlot <- list(MyPlot, data)
        }
    }

    return(MyPlot)
}
