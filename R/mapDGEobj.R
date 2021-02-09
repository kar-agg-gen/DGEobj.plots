#' Plot network of DGEobj relationships
#'
#' Reads a DGEobj and produces a node pair file defining parent/child relationships between
#' data items in the DGEobj.
#'
#' @param dgeObj DGEobj to find the parent/child relationships between data items.
#' @param plotType Plot type must be canvasXpress or ggplot (Default to canvasXpress).
#' @param directed Passed to igraph::graph_from_data_frame. Indicates if the graph should
#'     be directed or not. Default = TRUE.
#'
#' @return A class igraph network object for plotType ggplot and canvasxpress network plot for plotType canvasxpress.
#'
#' @examples
#' \dontrun{
#'   library(igraph)
#'   library(RColorBrewer)
#'
#'   # Prepare canvasxpress network plot
#'   mynet <- mapDGEobj(dgeObj)
#'
#'   # Prepare an iGraph object for plotting
#'   mynet <- mapDGEobj(dgeObj, plotType = "ggplot")
#'
#'   # Define a layout
#'   lay.sug <- layout_with_sugiyama(mynet)
#'
#'   # 2D Plot
#'   plot(mynet,
#'        plotType = "ggplot",
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
#'        plotType = "ggplot",
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
#' @importFrom canvasXpress canvasXpress
#' @importFrom htmlwidgets JS
#'
#' @export
mapDGEobj <- function(dgeObj, plotType = "canvasXpress", directed = TRUE) {

    assertthat::assert_that("DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be of class 'DGEobj'.")
    assertthat::assert_that(plotType %in% c("canvasXpress", "ggplot"),
                            msg = "Plot type must be either canvasXpress or ggplot.")

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

    if (plotType == "canvasXpress") {

        colnames(nodes) <- c("id", "Type", "BaseType")
        colnames(edges) <- c("id1", "id2")

        events <- htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
                                                if (o != null && o != false) {
                                                    if (o.objectType == null) {
                                                        if (o.nodes != null) {
                                                            t.showInfoSpan(e, '<b>' + 'Node' + ': ' + o.nodes[0].id + '</b> <br/>' +
                                                             '<b>' + 'Type'  + '</b>' + ': ' + o.nodes[0].Type + '<br/>' +
                                                             '<b>' + 'Base type'  + '</b>' + ': ' + o.nodes[0].BaseType + '<br/>');

                                                        } else if (o.edges != null) {
                                                            t.showInfoSpan(e, '<b>' + o.edges[0].id1 + '&#10230;' + o.edges[0].id2 + '</b>');
                                                        }
                                                    } else {
                                                        t.showInfoSpan(e, o.display);
                                                    };
                                                }; }}")

        mapDGEplot <- canvasXpress::canvasXpress(data              = list(nodeData = nodes, edgeData = edges),
                                                 colorNodeBy       = "Type",
                                                 edgeWidth         = 2,
                                                 graphType         = "Network",
                                                 nodeSize          = 30,
                                                 networkLayoutType = "forceDirected",
                                                 events            = events)

    } else {
        mapDGEplot <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
    }

    return(mapDGEplot)
}
