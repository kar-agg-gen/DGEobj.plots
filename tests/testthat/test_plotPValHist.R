context("DGEobj.plots - tests for plotPValHist.R functions")


test_that("plotPValHist.R: plotPvalHist()", {
    skip_if(is.null(t_obj1$DGEList))

    # testing plotPvalHist with savePlot and facet = TRUE
    pvalMatrix <- extractCol(getType(t_obj1, "topTable"), colName = "P.Value", robust = FALSE)
    pval_plot <- plotPvalHist(pvalMatrix)
    expect_s3_class(pval_plot, c("canvasXpress", "htmlwidget"))

    pval_plot <- plotPvalHist(pvalMatrix, plotType = "ggplot")
    expect_s3_class(pval_plot, c("gg","ggplot"))

    # testing plotPvalHist with savePlot and facet = FALSE
    pval_plot <- plotPvalHist(as.matrix(pvalMatrix),
                              facet     = FALSE)
    expect_length(pval_plot, 4)
    expect_s3_class(pval_plot[[1]], c("canvasXpress", "htmlwidget"))

    pval_plot <- plotPvalHist(as.matrix(pvalMatrix),
                              plotType = "ggplot",
                              facet     = FALSE)
    expect_length(pval_plot, 4)
    expect_s3_class(pval_plot[[1]], c("gg","ggplot"))

    expect_error(plotPvalHist(pvalMatrix, plotType = "cx"),
                 regexp = "Plot type must be either canvasXpress or ggplot.")
})
