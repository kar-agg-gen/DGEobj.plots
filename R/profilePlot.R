#' Plot log intensity versus log ratio
#'
#' A profile plot shows Log Intensity on the X axis and Log Ratio on the Y axis.
#' This function is intended to show the profile plot from a dataframe
#' created by topTable or topTreat.  Properly normalized data will generally be
#' centered around LogRatio = 0.
#'
#' By default, the plot places "logFC" on the Y axis and "AveExpr" on the X
#' axis. By default, a reference horizontal line is shown at LogRatio = 0 on the
#' Y axis. Optionally, additional reference lines will be drawn at +/- a user
#' supplied LogRatio threshold. A Loess line fit is drawn through the actual
#' data. The points are color coded using the significance measure (p-value or
#' FDR threshold) supplied by the user. By default, the P.Value field is used
#' with a threshold of 0.01 to color code the points.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable and topTreat.  The columns named "logFC",
#' "AveExpr", and "P.Value" are used by default to accommodate the column
#' names used in topTable/topTreat dataframes.  Any other dataframe
#' can be used with fold-change, intensity, and significance measures, with appropriate
#' arguments to define the column names to use provided. By default, the
#' column names will be used for the axis labels, but can be overridden with xlab and ylab arguments.
#'
#' A significance measure (which defaults to P.Value <= 0.01) is used to color code genes that
#' are significantly increased or decreased.  Use the appropriate arguments to use an FDR measure instead.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color, and Fill), but they can be
#' adjusted through the use of optional arguments. A length of 3 is
#' required for these arguments which applies the attributes in this order:
#' Increased, NoChange, Decreased.
#'
#' @param contrastDF Dataframe with LogIntensity and LogRatio columns and optionally a p-value or FDR column.
#' @param plotType Plot type must be canvasXpress or ggplot (Default to canvasXpress).
#' @param logRatioCol Name of the LogRatio column (Default = "logFC")
#' @param logIntCol Name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol Name of the p-value or FDR column (Default = "P.Value")
#' @param xlab X axis label (Defaults to the LogIntensity column name)
#' @param ylab Y axis label (Defaults to the LogRatio column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (Default = 0.01)
#' @param geneSymLabels Labels for gene symbols in contrastDF.
#' @param geneSymCol Name of the gene symbol column in contrastDF.  The gene symbol is
#'    not in topTable output by default so the user has to bind this column
#'    to the dataframe in advance.  Then this column will be used to label
#'    significantly changed points.
#' @param symbolSize Size of symbols for Up, no change, and Down. default = c(4, 1, 4);
#'        Note: All three cannot be the same size. Decimal values are acceptable to help offset that
#'        (e.g. 4, 4.1, 4.2).
#' @param symbolShape Shape of the symbols for Up, no change, and Down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle) Note: The same symbol shape cannot
#'        be selected for all three symbols.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: Colors cannot be duplicated.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: Colors cannot be duplicated.
#' @param alpha Controls the transparency of the plotted points (0-1; Default = 0.5)
#' @param sizeBySignificance Set to TRUE to size points by the negative Log10 of the
#'        Significance measure (Default = FALSE)
#' @param referenceLine Color for an intercept = 0 horizontal reference line
#'        (Default = "blue"; NULL disables)
#' @param foldChangeLines Position of reference horizontal lines for fold change
#'        (Default = log2(1.5); NULL disables)
#' @param refLineThickness Controls size of the horizontal reference line through geom_hline.
#'   (Default = 1)
#' @param lineFitType Enable a line fit through the data (Default = "loess";
#'        "lm" produces a linear fit. NULL disables)
#' @param lineFitColor Color for the fit line (Default = "goldenrod1")
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se", "nw", "sw", NULL.
#'        top/bottom/left/right place the legend outside the figure.  ne/se/nw/sw place the figure
#'        inside the figure. NULL disables the legend. Default = "right"
#' @param footnote Optional string placed right justified at bottom of plot.
#' @param footnoteSize Applies to footnote. (Default = 3)
#' @param footnoteColor Applies to footnote. (Default = "black")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (Default = 1)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'    # Simple plot with custom title(contrastDF is a topTable dataframe)
#'    myPlot <- profilePlot(contrastDF, title = "Plot Title")
#'
#'    # Some options with a custom datafile
#'    myPlot <- profilePlot(contrastDF,
#'                          pthreshold = 0.1,
#'                          logRatioCol = "Log2ratio",
#'                          logIntCol = "AverageIntensity",
#'                          pvalCol = "BHFDR",
#'                          xlab = "Log2Intensity", ylab = "Log2Ratio",
#'                          title = "Profile Plot Title",
#'                          referenceLine = "blue",
#'                          legendPosition = "ne")
#' }
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join filter arrange
#' @importFrom assertthat assert_that
#' @importFrom ggrepel geom_text_repel
#'
#' @export
profilePlot <- function(contrastDF,
                        plotType = "canvasXpress",
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold = 0.01,
                        geneSymLabels,
                        geneSymCol,
                        xlab = NULL, ylab = NULL, title = NULL,
                        symbolSize = c(4, 1, 4),
                        symbolShape = c(21, 20, 21),
                        symbolColor = c("black", "grey0", "grey25"),
                        symbolFill = c("red3", "deepskyblue4", "grey25"),
                        alpha = 0.5,
                        sizeBySignificance = FALSE,
                        referenceLine = "grey25",
                        foldChangeLines = log2(1.5),
                        refLineThickness = 1,
                        lineFitType = "loess",
                        lineFitColor = "goldenrod1",
                        legendPosition = "right",
                        footnote,
                        footnoteSize = 3,
                        footnoteColor = "black",
                        footnoteJust = 1) {

    # Make sure specified columns exist
    assertthat::assert_that(logRatioCol %in% colnames(contrastDF),
                            msg = "logRatioCol column not found in contrastDF.")
    assertthat::assert_that(logIntCol %in% colnames(contrastDF),
                            msg = "logIntCol column not found in contrastDF.")
    assertthat::assert_that(pvalCol %in% colnames(contrastDF),
                            msg = "pvalCol column not found in contrastDF.")
    assertthat::assert_that(plotType %in% c("ggplot", "canvasXpress"),
                            msg = "Plot type must be either ggplot or canvasXpress.")
    if (!missing(geneSymCol)) {
        assertthat::assert_that(geneSymCol %in% colnames(contrastDF),
                                msg = "geneSymCol column not found in contrastDF.")
    }
    if (!missing(symbolSize) || !missing(symbolShape) || !missing(symbolColor) || !missing(symbolFill)) {
        assertthat::assert_that(!length(symbolSize) == 3,
                                !length(symbolShape) == 3,
                                !length(symbolColor) == 3,
                                !length(symbolFill) == 3,
                                msg = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    }

    groupNames <- c("Increased", "Decreased", "No Change")

    # Columns to plot
    if (is.null(xlab)) { # Use colname unless supplied as argument
        xlab <- logIntCol
    }

    if (is.null(ylab)) {
        ylab <- logRatioCol
    }

    if (is.null(title)) {
        title <- ""
    }

    x <- make.names(colnames(contrastDF)[colnames(contrastDF) == logIntCol])
    y <- make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])
    colnames(contrastDF)[colnames(contrastDF) == logIntCol] <- make.names(colnames(contrastDF)[colnames(contrastDF) == logIntCol])
    colnames(contrastDF)[colnames(contrastDF) == logRatioCol] <- make.names(colnames(contrastDF)[colnames(contrastDF) == logRatioCol])

    # Need a NLP column for sizing
    contrastDF$negLog10P <- -log10(contrastDF[[pvalCol]])

    contrastDF$group <- NA
    for (i in seq(nrow(contrastDF))) {
        if (contrastDF[i, pvalCol] <= pthreshold) {
            if (contrastDF[i, logRatioCol] > 0) {
                contrastDF$group[i] <- "Increased"
            } else if (contrastDF[i, logRatioCol] < 0) {
                contrastDF$group[i] <- "Decreased"
            }
        } else {
            contrastDF$group[i] <- "No Change"
        }
    }

    contrastDF$group <- contrastDF$group %>%
        factor(levels = groupNames)

    # Set an order field to force plotting of NoChange first
    contrastDF$order <- 0
    contrastDF$order[contrastDF$group %in% c("Increased", "Decreased")] <- 1


    # plotType
    if (plotType == "canvasXpress") {
        symbolFill[1] <- paste(c("rgba(", paste(c(paste(col2rgb(symbolFill[1], alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")
        symbolFill[2] <- paste(c("rgba(", paste(c(paste(col2rgb(symbolFill[2], alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")
        symbolFill[3] <- paste(c("rgba(", paste(c(paste(col2rgb(symbolFill[3], alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")

        ## Create the canvasXpress df and var annotation
        cx.data <- data.frame(a = contrastDF[colnames(contrastDF) == x],
                              b = contrastDF[colnames(contrastDF) == y])
        colnames(cx.data) <- c(x, y)
        var.annot <- data.frame(Group = contrastDF$group, nLog10pVal = contrastDF$negLog10P)
        rownames(var.annot) <- rownames(cx.data)
        events <- NULL

        if (!missing(geneSymCol)) {
            var.annot <- cbind(var.annot, GeneLabel = contrastDF[[geneSymCol]])
            events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
                                                if (o != null && o != false) {
                                                    if (o.objectType == null) {
                                                        t.showInfoSpan(e, '<b>' + o.y.vars + '</b> <br/>' +
                                                        '<b>' + 'GeneLabel'  + '</b>' + ': ' + o.z.GeneLabel[0] + '<br/>' +
                                                        '<b>' + o.y.smps[0]  + '</b>' + ': ' + o.y.data[0][0] + '<br/>' +
                                                        '<b>' + o.y.smps[1]  + '</b>' + ': ' + o.y.data[0][1]);
                                                    } else {
                                                        t.showInfoSpan(e, o.display);
                                                    };
                                                }; }}")
        }

        # Optional Decorations
        sizeBy <- NULL
        sizeByShowLegend <- FALSE
        if (sizeBySignificance == TRUE) {
            sizeBy <- "nLog10pVal"
            sizeByShowLegend <- TRUE
        }

        if (!is.null(referenceLine)) {
            referenceLine <- paste(c("rgba(", paste(c(paste(col2rgb(referenceLine, alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")
            decorations <- list(line = list(list(color = referenceLine, width = refLineThickness, y = 0)))
        }

        if (!is.null(foldChangeLines)) {
            decorations <- list(line = append(decorations$line, list(list(color = symbolFill[which(groupNames == "Increased")],
                                                                          width = refLineThickness,
                                                                          y     = foldChangeLines),
                                                                     list(color = symbolFill[which(groupNames == "Decreased")],
                                                                          width = refLineThickness,
                                                                          y     = -foldChangeLines)

            )))
        }

        showLoessFit <- FALSE
        if (!is.null(lineFitType)) {
            lineFitColor <- paste(c("rgba(", paste(c(paste(col2rgb(lineFitColor, alpha = FALSE), collapse = ","), 0.5), collapse = ","), ")"), collapse = "")
            showLoessFit <- TRUE
        }

        # Footnote
        if (missing(footnote)) {
            footnote <- NULL
        }

        foldChangeMargin <- (foldChangeLines + (foldChangeLines * 0.2))
        profilePlot <- canvasXpress(data                    = cx.data,
                                    varAnnot                = var.annot,
                                    decorations             = decorations,
                                    graphType               = "Scatter2D",
                                    colorBy                 = "Group",
                                    colors                  = symbolFill,
                                    legendPosition          = legendPosition,
                                    showDecorations         = TRUE,
                                    showLoessFit            = showLoessFit,
                                    loessColor              = lineFitColor,
                                    sizes                   = c(4, 10, 12, 14, 16, 18, 20, 22, 24, 26),
                                    sizeByShowLegend        = sizeByShowLegend,
                                    title                   = title,
                                    xAxisTitle              = list(xlab),
                                    yAxisTitle              = list(ylab),
                                    sizeBy                  = sizeBy,
                                    setMaxY                 = foldChangeMargin,
                                    setMinY                 = -1*foldChangeMargin,
                                    citation                = footnote,
                                    citationFontSize        = footnoteSize,
                                    citationColor           = footnoteColor,
                                    events                  = events)

    } else {
        names(symbolShape) = groupNames
        names(symbolSize)  = groupNames
        names(symbolColor) = groupNames
        names(symbolFill)  = groupNames

        ssc = data.frame(group = groupNames,
                         symbolShape = symbolShape,
                         symbolSize = symbolSize,
                         symbolColor = symbolColor,
                         symbolFill = symbolFill,
                         stringsAsFactors = FALSE)


        contrastDF <- contrastDF %>%
            dplyr::left_join(ssc)

        profilePlot <- ggplot(contrastDF, aes_string(x = x, y = y)) +
            aes(shape = group,
                size = group,
                color = group,
                fill = group,
                order = order) +
            # Scale lines tell it to use the actual values, not treat them as factors
            scale_shape_manual(name = "Group", guide = "legend", labels = ssc$group,
                               values = ssc$symbolShape) +
            scale_size_manual(name = "Group", guide = "legend", labels = ssc$group,
                              values = ssc$symbolSize) +
            scale_color_manual(name = "Group", guide = "legend", labels = ssc$group,
                               values = ssc$symbolColor) +
            scale_fill_manual(name = "Group", guide = "legend", labels = ssc$group,
                              values = ssc$symbolFill) +
            geom_point(alpha = alpha)

        # Optional Decorations
        if (sizeBySignificance == TRUE) {
            profilePlot <- profilePlot + aes(size = negLog10P) +
                scale_size_continuous()
        }

        if (!is.null(referenceLine)) {
            profilePlot <- profilePlot +
                geom_hline(yintercept = 0,
                           color = referenceLine,
                           size = refLineThickness,
                           alpha = 0.5)
        }

        if (!is.null(foldChangeLines)) {
            profilePlot <- profilePlot +
                geom_hline(yintercept = foldChangeLines,
                           color = symbolFill["Increased"],
                           size = refLineThickness,
                           alpha = 0.5) +
                geom_hline(yintercept = -foldChangeLines,
                           color = symbolFill["Decreased"],
                           size = refLineThickness,
                           alpha = 0.5)
        }

        if (!is.null(lineFitType)) {
            profilePlot <- profilePlot +
                geom_smooth(aes(group = NULL, shape = NULL, size = NULL, color = NULL, fill = NULL),
                            method = tolower(lineFitType),
                            size = refLineThickness,
                            color = lineFitColor,
                            alpha = alpha,
                            se = FALSE,
                            show.legend = FALSE)
        }

        # Add genesym labels to increased, decreased genes
        if (!missing(geneSymLabels) & !missing(geneSymCol)) {
            idx <- contrastDF[[geneSymCol]] %in% geneSymLabels
            contrastDFsubset <- contrastDF[idx,]
            profilePlot <- profilePlot +
                geom_text_repel(data = contrastDFsubset,
                                aes_string(x = x, y = y, label = geneSymCol),
                                show.legend = FALSE)
        }


        if (!is.null(title)) {
            profilePlot <- profilePlot +
                ggtitle(title)
        }

        # Footnote
        if (!missing(footnote)) {
            profilePlot <- addFootnote(profilePlot,
                                       footnoteText = footnote,
                                       footnoteSize = footnoteSize,
                                       footnoteColor = "black",
                                       footnoteJust = footnoteJust)
        }

        profilePlot <- setLegendPosition(profilePlot, legendPosition)
    }
    return(profilePlot)
}
