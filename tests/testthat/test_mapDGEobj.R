context("DGEobj.plots - tests for mapDGEobj.R functions")


test_that('mapDGEobj.R: mapDGEobj()', {
    map_DGEobj <- mapDGEobj(t_obj1)
    expect_s3_class(map_DGEobj, c("canvasXpress", "htmlwidget"))

    map_DGEobj <- mapDGEobj(t_obj1, plotType = "ggplot")
    expect_s3_class(map_DGEobj, "igraph")

    expect_error(mapDGEobj(obj),
                 regexp = "object 'obj' not found")
    expect_error(mapDGEobj(t_obj1, plotType = "cx"),
                 regexp = "Plot type must be either canvasXpress or ggplot.")
})
