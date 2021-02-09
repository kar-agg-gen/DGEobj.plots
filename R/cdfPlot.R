#' Create deluxe CDF Plots
#'
#' CDF plots are a good complement to p-value histograms as a way to evaluate
#' model performance and examine support for differential expression. Results
#' are ranked by p-value on the x-axis and the p-value plotted on the y-axis.
#' Since p-value distributions should be flat, this type of plot should produce a
#' straight line.  Any observations that fail to meet the null hypothesis will
#' appear as a break in the line at the low end of the curve.
#'
#' This function is designed to take topTable dataframes and display the
#' corresponding CDF plot. Data for the p-values below 0.1 (configurable via
#' pvalMax argument) are show in a full size plot. An inset figure shows the
#' whole p-value scale and highlights the portion shown in the full plot. Points
#' below 0.01 are a different color by default (threshold set by pThreshold
#' argument; shape/color attributes customizable through other arguments).
#'
#' The output can be controlled using printPlot = TRUE, which outputs the compound plot
#' to the console/knitr.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable.  The column
#' "P.Value" is used by default to accommodate the column names used in topTable
#' dataframes.  Any other dataframe can be used with by explicitly defining the
#' p-value column name with the appropriate argument.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color, and Fill).
#' There are optional arguments that allow these to be adjusted. A length of 2
#' is required for these arguments which applies the attributes in
#' this order: Significant, Not Significant.
#'
#' @param contrastDF A dataframe with LogRatio and LogIntensity columns and optionally a p-value or FDR column.
#' @param plotType Plot type must be canvasXpress or ggplot (Default to canvasXpress).
#' @param pvalCol Name of the p-value or FDR column (Default = "P.Value")
#' @param pvalMax Limit the range of the main plot (Default = 0.10)
#' @param pThreshold Used to color points (default = 0.01)
#' @param xlab X axis label (default is "LogIntensity column name "Rank")
#' @param ylab Y axis label (default is p-value column name)
#' @param title Plot title (Optional)
#' @param insetTitle Title for the inset plot (Optional)
#' @param symbolSize Size of symbols for up, no change, and down. Default = c(4, 1, 4);
#'        Note: All three cannot be the same size. Decimal values are acceptable to help offset that
#'        (e.g. 4, 4.1, 4.2).
#' @param symbolShape Shape of the symbols for up, no change, and down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle); Note: The same symbol shape cannot
#'        be selected for all three symbols.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: Colors cannot be duplicated.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: Colors cannot be duplicated.
#' @param alpha Controls the transparency of the plotted points (0-1; Default =
#'   0.5)
#' @param referenceLine Color for an horizontal line drawn at the p-threshold
#'   (Default = NULL; NULL disables, set to desired color to enable)
#' @param refLineThickness Set thickness of the reference line (Default = 1)
#' @param legendPosition (Default = "right")
#' @param baseFontSize The smallest size font in the figure in points. (Default
#'   = 12)
#' @param viewportX (Default = 0.15)
#' @param viewportY (Default = 0.93)
#' @param viewportWidth (Default = 0.35)
#' @param printPlot Specify printing the combined plot to the console/knitr
#'   (Default = TRUE)
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0 - 1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return A list containing main plot, inset plot for both plotType. For plotType ='ggplot' list contains viewport, The plot can be
#'    reconstructed with print(cdfPlot$main); print(inset, vp = cdfPlot$viewport)
#'
#' @examples
#' \dontrun{
#'    # Plot to console (contrastDF is a topTable dataframe)
#'    cdfPlot(contrastDF, title = "My CDF Plot")
#' }
#' @import ggplot2 magrittr
#' @importFrom grid viewport
#' @importFrom dplyr arrange left_join
#' @importFrom assertthat assert_that
#' @importFrom canvasXpress canvasXpress
#'
#' @export
cdfPlot <- function(contrastDF,
                    plotType = "canvasXpress",
                    pvalCol = "P.Value",
                    pThreshold = 0.01,
                    xlab = NULL,
                    ylab = NULL,
                    title = NULL,
                    insetTitle = NULL,
                    symbolSize = c(2, 1),
                    symbolShape = c(20, 20),
                    symbolColor = c("red3", "deepskyblue4"),
                    symbolFill = c("red3", "deepskyblue4"),
                    alpha = 1,
                    referenceLine = NULL,
                    refLineThickness = 3,
                    legendPosition = "right",
                    baseFontSize = 12,
                    viewportX = 0.15,
                    viewportY = 0.93,
                    viewportWidth = 0.35,
                    pvalMax = 0.10,
                    printPlot = TRUE,
                    footnote,
                    footnoteSize = 3,
                    footnoteColor = "black",
                    footnoteJust = 1) {

    assertthat::assert_that(pvalCol %in% colnames(contrastDF),
                            msg = "Specified pvalCol not found in the supplied dataframe (contrastDF).")
    assertthat::assert_that(plotType %in% c("canvasXpress", "ggplot"),
                            msg = "Plot type must be either canvasXpress or ggplot.")

    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 2,
                                !length(symbolShape) == 2,
                                !length(symbolColor) == 2,
                                !length(symbolFill) == 2,
                                msg = "All specified symbol arguments must be of length 2, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    groupNames <- c("Significant", "Not Significant")
    # Storing column names in x and y variable
    x <- "Rank"
    y <- pvalCol

    if (is.null(xlab)) {
        xlab <- x
    }

    if (is.null(ylab)) {
        ylab <- y
    }

    if (is.null(title)) {
        title = ""
    }

    if (is.null(insetTitle)) {
        insetTitle = ""
    }

    # Combo PLOT: full data inset, most significant data in main plot
    # Rank by p-value
    contrastDF <- contrastDF %>%
        dplyr::arrange(!!sym(y))
    contrastDF$Rank <- c(1:nrow(contrastDF))

    # Let's plot the p-value subsets
    Sig <- contrastDF[[y]] <= pThreshold
    NotSig <- contrastDF[[y]] > pThreshold

    # Create group factor column in contrastDF
    contrastDF$group <- NA
    contrastDF$group[Sig] <- "Significant"
    contrastDF$group[NotSig] <- "Not Significant"

    contrastDF$group <- contrastDF$group %>%
        factor(levels = groupNames)

    # Set an order field to force plotting of 'not significant' first
    contrastDF$order <- NA
    contrastDF$order[NotSig] <- 1
    contrastDF$order[Sig] <- 2

    # Rows to include in the zoomed in plot
    subsetRows <- sum(contrastDF[[y]] <= pvalMax)
    contrastDFsubset <- contrastDF[1:subsetRows,]

    if (plotType == "canvasXpress") {
        ## Create the canvasXpress cx.data and var.annot
        # Main plot
        cx.data <- data.frame(a = contrastDF[colnames(contrastDF) == x],
                              b = contrastDF[colnames(contrastDF) == y])
        colnames(cx.data) <- c(x, y)
        var.annot <- data.frame(Group = contrastDF$group)
        rownames(var.annot) <- rownames(cx.data)

        # Inst plot
        cx.data.subset <- data.frame(a = contrastDFsubset[colnames(contrastDFsubset) == x],
                              b = contrastDFsubset[colnames(contrastDFsubset) == y])
        colnames(cx.data.subset) <- c(x, y)
        var.annot.subset <- data.frame(Group = contrastDFsubset$group)
        rownames(var.annot.subset) <- rownames(cx.data.subset)

        decorations <- list()
        if (!is.null(referenceLine)) {
            referenceLine <- paste(c("rgba(", paste(c(paste(col2rgb(referenceLine, alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")
            decorations <- list(line = list(list(color = referenceLine, width = refLineThickness, y = pThreshold)))
        }

        # Footnote
        if (missing(footnote)) {
            footnote <- NULL
        }
        maxY <- pThreshold + pThreshold*0.2

        cdfMain <- canvasXpress::canvasXpress(data              = cx.data,
                                              varAnnot          = var.annot,
                                              decorations       = decorations,
                                              graphType         = "Scatter2D",
                                              colorBy           = "Group",
                                              colors            = symbolFill,
                                              title             = title,
                                              xAxisTitle        = xlab,
                                              yAxisTitle        = ylab,
                                              citation          = footnote,
                                              citationFontSize  = footnoteSize,
                                              citationColor     = footnoteColor,
                                              setMaxY           = maxY)

        cdfInset <- canvasXpress::canvasXpress(data             = cx.data.subset,
                                               varAnnot         = var.annot.subset,
                                               graphType        = "Scatter2D",
                                               colorBy          = "Group",
                                               colors           = symbolFill,
                                               title            = insetTitle,
                                               xAxisTitle       = xlab,
                                               yAxisTitle       = ylab,
                                               setMaxY          = max(contrastDFsubset[[y]]))
        cdfPlot <- list("main" = cdfMain, "inset" = cdfInset)
    } else {
        names(symbolShape) <- groupNames
        names(symbolSize)  <- groupNames
        names(symbolColor) <- groupNames
        names(symbolFill)  <- groupNames

        ssc <- data.frame(group = groupNames,
                          symbolShape = symbolShape,
                          symbolSize = symbolSize,
                          symbolColor = symbolColor,
                          symbolFill = symbolFill,
                          stringsAsFactors = FALSE)

        contrastDF <- contrastDF %>%
            dplyr::left_join(ssc)

        # Plot subset percent of the data for the main plot
        cdfMain <- ggplot(contrastDFsubset, aes_string(x = x, y = y)) +
            aes(shape = group, size = group, color = group, fill = group, order = order) +
            # Scale lines tell it to use the actual values, not treat them as factors
            scale_shape_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolShape) +
            scale_size_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolSize) +
            scale_color_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolColor) +
            scale_fill_manual(name = "Group", guide = "legend", labels = ssc$group, values = ssc$symbolFill) +
            geom_point(alpha = alpha)

        # Optional Decorations
        if (!is.null(referenceLine)) {
            cdfMain <- cdfMain +
                geom_hline(yintercept = pThreshold, color = referenceLine,
                           size = refLineThickness, alpha = 0.5)
        }

        # Add Labels
        cdfMain <- cdfMain +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title)

        cdfMain <- setLegendPosition(cdfMain, legendPosition)

        if (!missing(footnote)) {
            cdfMain <- addFootnote(cdfMain,
                                   footnoteText = footnote,
                                   footnoteSize = footnoteSize,
                                   footnoteColor = footnoteColor,
                                   footnoteJust = footnoteJust)
        }

        # Set up the inset plot with All Data
        cdfInset <- ggplot(contrastDF, aes_string(x = x, y = y)) +
            aes(shape = group, size = group,
                color = group, fill = group,
                order = order) +
            # Scale lines tell it to use the actual values, not treat them as factors
            scale_shape_manual(name = "Group", guide = "none", labels = ssc$group,
                               values = ssc$symbolShape) +
            scale_size_manual(name = "Group", guide = "none", labels = ssc$group,
                              values = ssc$symbolSize) +
            scale_color_manual(name = "Group", guide = "none", labels = ssc$group,
                               values = ssc$symbolColor) +
            scale_fill_manual(name = "Group", guide = "none", labels = ssc$group,
                              values = ssc$symbolFill) +
            geom_rect(xmin = 0, xmax = subsetRows,
                      ymin = 0, ymax = max(contrastDFsubset[[y]]), color = "lightblue",
                      fill = "lightblue", alpha = 0.2) +
            geom_point(alpha = alpha)

        # Add Labels and title
        cdfInset <- cdfInset +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(insetTitle)

        # Adjust font size for inset.
        factor <- 4/baseFontSize
        cdfInset <- cdfInset + baseTheme(baseFontSize*factor)

        # Adjust viewport Y if main Title present
        vy <- viewportY
        if (!is.null(title)) {
            adjust <- (baseFontSize/150)
            vy <- viewportY - adjust
        }

        # A viewport taking up a fraction of the plot area (upper left)
        vp <- grid::viewport(width = viewportWidth, height = viewportWidth,
                             x = viewportX, y = vy,
                             just = c("left", "top"))

        if (printPlot == TRUE) {
            print(cdfMain)
            print(cdfInset, vp = vp)
        }

        cdfPlot <- list(main = cdfMain, inset = cdfInset, viewport = vp)
    }
    return(cdfPlot)
}
