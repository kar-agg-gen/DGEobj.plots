context("DGEobj.plots - tests for ggplotMDS.R functions")


test_that("ggplotMDS.R: ggplotMDS()", {
    skip_if(is.null(t_obj1$DGEList))
    skip_if(is.null(t_obj1$design$organism))

    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$organism)
    expect_length(mds_plot, 2)
    expect_named(mds_plot, c("plot", "mdsobj"))
    expect_type(mds_plot, "list")
    expect_s3_class(mds_plot$plot, c("canvasXpress", "htmlwidget"))

    mds_plot <- ggplotMDS(DGEdata  = t_obj1,
                          plotType = "ggplot",
                          colorBy  = t_obj1$design$organism)
    expect_length(mds_plot, 2)
    expect_named(mds_plot, c("plot", "mdsobj"))
    expect_type(mds_plot, "list")
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing parameter sizeBy and shapeBy
    mds_plot <- ggplotMDS(DGEdata        = t_obj1,
                          colorBy        = t_obj1$design$organism,
                          shapeBy        = t_obj1$design$ReplicateGroup,
                          sizeBy         = rep(1:4,12),
                          hlineIntercept = 0.25,
                          vlineIntercept = 0.25)
    expect_s3_class(mds_plot$plot, c("canvasXpress", "htmlwidget"))

    mds_plot <- ggplotMDS(DGEdata        = t_obj1,
                          plotType       = "ggplot",
                          colorBy        = t_obj1$design$organism,
                          shapeBy        = t_obj1$design$ReplicateGroup,
                          sizeBy         = rep(1:4,12),
                          hlineIntercept = 0.25,
                          vlineIntercept = 0.25)
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing assert statements
    expect_error(ggplotMDS(DGEdata = 1:10),
                 regexp = "DGEdata must be of class 'DGEList', 'DGEobj', or 'matrix'.")
    expect_error(ggplotMDS(DGEdata = t_obj1),
                 regexp = "colorBy must be specified and should be the length of the number of columns in DGEdata.")
    expect_error(ggplotMDS(DGEdata = t_obj1,
                           colorBy = t_obj1$design$organism,
                           shapeBy = 1:10),
                 regexp = "shapeBy should be the length of the number of columns in DGEdata.")
    expect_error(ggplotMDS(DGEdata = t_obj1,
                           colorBy = t_obj1$design$organism,
                           sizeBy  = 1:10),
                 regexp = "sizeBy should be the length of the number of columns in DGEdata.")
    expect_error(ggplotMDS(DGEdata  = t_obj1,
                           plotType = "cx",
                           colorBy  = t_obj1$design$organism),
                 regexp = "Plot type must be either canvasXpress or ggplot.")
})

test_that("ggplotMDS.R: MDS_var_explained()", {
    skip_if(is.null(t_obj1$DGEList))
    skip_if(is.null(t_obj1$design$organism))

    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$organism)
    var_result <- MDS_var_explained(mds_plot$mdsobj)
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("canvasXpress", "htmlwidget"))
    expect_s3_class(var_result$cumvar, c("canvasXpress", "htmlwidget"))

    var_result <- MDS_var_explained(mds_plot$mdsobj,
                                    plotType = "ggplot")
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("gg", "ggplot"))
    expect_s3_class(var_result$cumvar, c("gg", "ggplot"))

    # testing mds matrix
    var_result <- MDS_var_explained(t_obj1$counts)
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("canvasXpress", "htmlwidget"))
    expect_s3_class(var_result$cumvar, c("canvasXpress", "htmlwidget"))

    var_result <- MDS_var_explained(t_obj1$counts,
                                    plotType = "ggplot")
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("gg", "ggplot"))
    expect_s3_class(var_result$cumvar, c("gg", "ggplot"))

    # testing assert statements
    expect_error(MDS_var_explained(), regexp = "mds is required and must be specified.")
    expect_error(MDS_var_explained(t_obj1$counts, plotType = "cx"),
                 regexp = "Plot type must be either canvasXpress or ggplot.")

})
