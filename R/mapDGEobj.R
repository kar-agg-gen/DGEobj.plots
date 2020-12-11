#' Plot network of DGEobj relationships
#'
#' Reads a DGEobj and produces a node pair file defining parent/child relationships between
#' data items in the DGEobj.
#'
#' @param dgeObj DGEobj to find the parent/child relationships between data items.
#' @param directed Passed to igraph::graph_from_data_frame. Indicates if the graph should
#'     be directed or not. Default = TRUE.
#'
#' @return A class igraph network object.
#'
#' @examples
#' \dontrun{
#'   library(igraph)
#'   library(RColorBrewer)
#'
#'   # Prepare an iGraph object for plotting
#'   mynet <- mapDGEobj(dgeObj)
#'
#'   # Define a layout
#'   lay.sug <- layout_with_sugiyama(mynet)
#'
#'   # 2D Plot
#'   plot(mynet,
#'        edge.arrow.size = .3,
#'        vertex.color = "dodgerblue2",
#'        vertex.frame.color = "dodgerblue4",
#'        layout = lay.sug$layout)
#'
#'   # 2D Plot; colorby the basetype attribute
#'   pal <- brewer.pal(length(unique(V(mynet)$basetype)), "Set1")
#'   myPallet <- pal[as.numeric(as.factor(vertex_attr(mynet, "basetype")))]
#'
#'   plot(mynet,
#'        edge.arrow.size = .3,
#'        vertex.color = myPallet,
#'        vertex.label.family = "Helvetica",
#'        layout=lay.sug$layout)
#'
#'   # 2D Interactive plot
#'   plotHandle <- tkplot(net,
#'                        vertex.color = myPallet,
#'                        vertex.frame.color = "dodgerblue4",
#'                        canvas.width = 800,
#'                        canvas.height = 800,
#'                        layout = lay.sug$layout)
#'   tk_close(plotHandle)
#' }
#'
#' @import magrittr
#' @importFrom assertthat assert_that
#' @importFrom igraph graph_from_data_frame
#'
#' @export
mapDGEobj <- function(dgeObj, directed = TRUE) {

    assertthat::assert_that("DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be of class 'DGEobj'.")

    child <- names(dgeObj)
    parent <- attr(dgeObj, "parent") %>% as.character()
    type <- attr(dgeObj, "type") %>% as.character()
    basetype <- attr(dgeObj, "basetype") %>% as.character()
    edges <- as.data.frame(cbind(parent = parent,
                                 child = child),
                           stringsAsFactors = FALSE)
    # Remove incomplete edges (children without parents)
    idx <- lapply(parent, nchar) != 0
    edges <- edges[idx,]

    nodes <- data.frame(node = names(dgeObj),
                        type = as.character(type),
                        basetype = as.character(basetype))

    igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
}
