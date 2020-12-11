require(testthat)
require(stats)

require(DGEobj)
require(DGEobj.plots)
require(DGEobj.utils)

t_obj1 <- readRDS(system.file("exampleObj.RDS", package = "DGEobj", mustWork = TRUE))
