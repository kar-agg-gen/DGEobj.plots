context("DGEobj.plots - tests for plotDispersion.R functions")


test_that("plotDispersion.R: plotDispersion()", {
    skip_if(is.null(t_obj1$DGEList))

    # creating designMatrix and designlist
    dgelist <- t_obj1$DGEList
    designMatrix <- stats::model.matrix(~ 0 + ReplicateGroup, getItem(t_obj1, "design"))

    # testing with input DGEList
    plot_disp <- plotDispersion(DGEdata = dgelist, designMatrix = designMatrix)
    expect_s3_class(plot_disp, c("canvasXpress", "htmlwidget"))

    # testing with input DGEList
    plot_disp <- plotDispersion(DGEdata = dgelist, designMatrix = designMatrix, plotType = "ggplot")
    expect_s3_class(plot_disp, c("gg", "ggplot"))

    # testing with input countMatrix
    plot_disp <- plotDispersion(DGEdata       = t_obj1$counts,
                                designMatrix  = designMatrix,
                                lineFit       = "loess",
                                plotValue     =  "BCV")
    expect_s3_class(plot_disp, c("canvasXpress", "htmlwidget"))

    plot_disp <- plotDispersion(DGEdata       = t_obj1$counts,
                                plotType      = "ggplot",
                                designMatrix  = designMatrix,
                                lineFit       = "loess",
                                plotValue     =  "BCV")
    expect_s3_class(plot_disp, c("gg", "ggplot"))

    # testing assert statements
    expect_error(plotDispersion(),
                 regexp = "DGEdata must be specified.")
    expect_error(plotDispersion(DGEdata = t_obj1$counts),
                 regexp = "argument \"designMatrix\" is missing, with no default")
    expect_error(plotDispersion(DGEdata = dgelist, designMatrix = designMatrix, plotValue = "cx"),
                 regexp = "Plot value must be either dispersion or BCV.")
    expect_error(plotDispersion(DGEdata = dgelist, designMatrix = designMatrix, plotType = "cx"),
                 regexp = "Plot type must be either ggplot or canvasXpress.")
})
