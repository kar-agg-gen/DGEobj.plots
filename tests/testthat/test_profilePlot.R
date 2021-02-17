context("DGEobj.plots - tests for profilePlot.R functions")


test_that("profilePlot.R: profilePlot()", {
    skip_if(!("ReplicateGroupDesign_fit" %in% names(t_obj1)))

    contrastDF <- t_obj1$BDL_vs_Sham
    profile_plot <- profilePlot(contrastDF)
    expect_s3_class(profile_plot, c("canvasXpress", "htmlwidget"))

    profile_plot <- profilePlot(contrastDF, plotType = "ggplot")
    expect_s3_class(profile_plot, c("gg", "ggplot"))

    contrastDF$GeneSym <- rep(c("sym1", "sym2","sym3","sym4"),nrow(contrastDF)/4)
    profile_plot <- profilePlot(contrastDF         = contrastDF,
                                title              = "Plot Title",
                                sizeBySignificance = TRUE,
                                geneSymCol         = "GeneSym",
                                geneSymLabels      = c("sym1", "sym2"),
                                symbolSize         = c(10, 10, 4),
                                xlab               = "XLabel",
                                ylab               = "YLabel",
                                footnote           = "This is footnote")
    expect_s3_class(profile_plot, c("canvasXpress", "htmlwidget"))

    profile_plot <- profilePlot(contrastDF         = contrastDF,
                                plotType           = "ggplot",
                                title              = "Plot Title",
                                sizeBySignificance = TRUE,
                                geneSymCol         = "GeneSym",
                                geneSymLabels      = c("sym1", "sym2"),
                                symbolSize         = c(10, 10, 4),
                                xlab               = "XLabel",
                                ylab               = "YLabel",
                                footnote           = "This is footnote")
    expect_s3_class(profile_plot, c("gg","ggplot"))

    expect_error(profilePlot(contrastDF, logRatioCol = "xyz"),
                 regexp =  "logRatioCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", logIntCol = "xyz"),
                 regexp =  "logIntCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", pvalCol = "xyz"),
                 regexp = "pvalCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", geneSymCol = "xyz"),
                 regexp = "geneSymCol column not found in contrastDF.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", symbolSize = 2),
                 regexp = "All specified symbol arguments must be of length 3, including symbolSize, symbolShape, symbolColor, and symbolFill.")
    expect_error(profilePlot(contrastDF, logRatioCol = "adj.P.Val", plotType = "cx"),
                 regexp = "Plot type must be either ggplot or canvasXpress.")
})
