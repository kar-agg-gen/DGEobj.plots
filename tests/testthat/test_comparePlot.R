context("DGEobj.plots - tests for comparePlot.R functions")


test_that("comparePlot.R: comparePlot()", {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    contrastList <- getType(t_obj1, "topTable")[1:2]
    # Capture the default logFC and P.Value
    compareDat <- comparePrep(contrastList)

    ## canvasXpress plot testing
    cPlot <- comparePlot(compareDat)
    expect_s3_class(cPlot , c("canvasXpress", "htmlwidget"))

    # testing function with compareDat without significance measures supplied
    cPlot <- comparePlot(compareDat[, 1:2])
    expect_s3_class(cPlot , c("canvasXpress", "htmlwidget"))

    # testing aesthetics of plots
    cPlot <- comparePlot(compareDat,
                         pThreshold = 0.001,
                         title    = "MyPlot",
                         xlab     = "xaxis-title",
                         ylab     = "yaxis-title",
                         symbolSize = c(5, 5, 2, 2),
                         symbolShape = c(21, 21, 21, 20),
                         symbolFill = c("red", "blue", "yellow", "grey25"),
                         alpha = 0.5,
                         crosshair = "grey50",
                         referenceLine = "darkgoldenrod1",
                         refLineThickness = 1,
                         legendPosition = "right",
                         footnote = "This is my footnote")
    expect_s3_class(cPlot , c("canvasXpress", "htmlwidget"))

    ## ggplot plot testing
    cPlot <- comparePlot(compareDat, plotType = "ggplot")
    expect_s3_class(cPlot , c("gg", "ggplot"))



    cPlot <- comparePlot(compareDat[,1:2], plotType = "ggplot")
    expect_s3_class(cPlot , c("gg", "ggplot"))

    # testing aesthetics of plots
    cPlot <- comparePlot(compareDat,
                         plotType = "ggplot",
                         pThreshold = 0.001,
                         title    = "MyPlot",
                         xlab     = "xaxis-title",
                         ylab     = "yaxis-title",
                         symbolSize = c(5, 5, 2, 2),
                         symbolShape = c(21, 21, 21, 20),
                         symbolFill = c("red", "blue", "yellow", "grey25"),
                         alpha = 0.5,
                         crosshair = "grey50",
                         referenceLine = "darkgoldenrod1",
                         refLineThickness = 1,
                         legendPosition = "right",
                         footnote = "This is my footnote")
    expect_setequal(unlist(cPlot$labels[c("title","y", "x")]), c("MyPlot", "yaxis-title", "xaxis-title"))
    expect_setequal(cPlot$layers[[2]]$aes_params$colour, "grey50")

    # testing assert statement
    expect_error(comparePlot(compareDat[, 1, drop = FALSE]),
                 regexp = "Need at least two numeric columns in compareDF.")
    expect_error(comparePlot(t_obj1$design),
                 regexp = "Need at least two numeric columns in compareDF.")
    expect_error(comparePlot(compareDat, plotType = "cx"),
                 regexp = "Plot type must be either canvasXpress or ggplot.")

    # failing function with Symbol argument length < 4
    expect_error(comparePlot(compareDat, symbolSize = 1),
                 regexp = "All specified symbol arguments must be of length 4, including symbolSize, symbolShape, symbolColor, and symbolFill.")
})

test_that("comparePlot.R: comparePrep()", {
    suppressWarnings(skip_if(is.null(getType(t_obj1, "topTable"))))

    contrastList <- getType(t_obj1, "topTable")[1:2]
    # Capture the default logFC and P.Value
    compareDat <- comparePrep(contrastList)
    expect_s3_class(compareDat,"data.frame")

    expect_error(comparePrep(contrastList[[1]]),
                 regexp = "contrastList must be a named list of length 2 where both items are of class 'data.frame'.")
    expect_error(comparePrep(contrastList, valueCol = "P.val"),
                 regexp = "The valueCol must be included in the colnames of both items of contrastList.")
    expect_error(comparePrep(contrastList, significanceCol = "P.val"),
                 regexp = "The significanceCol must be included in the colnames of both items of contrastList.")
    contrastList_uncommon_ids <- list("BDL_vs_Sham" = contrastList$BDL_vs_Sham[1:10,], "EXT1024_vs_BDL" = contrastList$EXT1024_vs_BDL[21:30,])
    expect_error(comparePrep(contrastList_uncommon_ids),
                 regexp = "No common gene IDs were found between the two dataframes in contrastList.")
})
