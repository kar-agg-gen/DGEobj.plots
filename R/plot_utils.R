baseTheme <- function(base_size = 18,
                      base_family = "",
                      scale_legend = TRUE) {
    if (scale_legend == TRUE && base_size > 17) {
        Legend.ScaledSize <- 14/base_size
        Legend.ScaledFont <- 12/base_size
    } else {
        Legend.ScaledSize = 1.2 # theme_grey defaults
        Legend.ScaledFont = 0.8
    }

    baseTheme <- theme(  # Base size is the axis tick mark labels; other elements are scaled
        axis.text.x = element_text(size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(size = rel(1.25), vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = rel(1.25), vjust = 0.5, hjust = 0.5, color = "black"),
        plot.title = element_text(face = "bold", size = rel(1.5)),
        legend.text = element_text(colour = "Black", size = rel(Legend.ScaledFont)),
        legend.title = element_text(colour = "Black", size = rel(1.2)),
        legend.title.align = 0.5,
        legend.key.size = unit(Legend.ScaledSize, "lines"),
        strip.text.x = element_text(size = rel(1.6)),
        strip.text.y = element_text(size = rel(1.6)),
        text = element_text(size = base_size, family = base_family)
    )
}


facetTheme <- function(base_size = 18,
                       base_family = "") {
    facetTheme = theme(
        axis.text.x = element_text(size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(face = "bold", colour = "Black", size = rel(2.0)),
        axis.title.y = element_text(face = "bold", colour = "Black", size = rel(2.0)),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(lineheight = .8, face = "bold", size = rel(2.4)),
        legend.text = element_text(colour = "Black", size = rel(1.6)),
        legend.title = element_text(colour = "Black", size = rel(1.4)),
        legend.title.align = 0.5,
        strip.text.x = element_text(size = rel(1.6)),
        strip.text.y = element_text(size = rel(1.6)),
        text = element_text(size = base_size, family = base_family)
    )
}


baseFont <- function(base_size = 12,
                     base_family = ""){
    baseFont <- theme(text = element_text(size = base_size, family = base_family))
}


setLegendPosition <- function(ggplot,
                              legendPosition = "right",
                              themeStyle = "grey") {

    # 8 legend positions, with/without grey background
    if (tolower(legendPosition) %in% c("top", "bottom", "left", "right")) {

        ggplot <- ggplot +
            theme(legend.position = tolower(legendPosition))

    } else { # Place legend inside figure; with or without a grey background
        #Four corners
        if (tolower(legendPosition) == "ne") {
            ggplot <- ggplot +
                theme(legend.justification = c(1,1), legend.position = c(1,1))
        } else if (tolower(legendPosition) == "se") {
            ggplot <- ggplot +
                theme(legend.justification = c(1,0), legend.position = c(1,0))
        } else if (tolower(legendPosition) == "sw") {
            ggplot <- ggplot +
                theme(legend.justification = c(0,0), legend.position = c(0,0))
        } else if (tolower(legendPosition) == "nw") {
            ggplot <- ggplot +
                theme(legend.justification = c(0,1), legend.position = c(0,1))
        }

    }

    # Set Legend background
    if (tolower(themeStyle) == "bw") {
        ggplot <- ggplot +
            theme(legend.background = element_rect(color = "grey30", fill = "white"))
    } else { # Grey style
        ggplot <- ggplot +
            theme(legend.background = element_rect(color = "grey30", fill = "grey90"))
    }

    return(ggplot)

}



#' @importFrom utils packageVersion
yrange <- function(my.ggp) {
    assertthat::assert_that(class(my.ggp)[[2]] == "ggplot",
                            msg = "my.ggp must be of class 'ggplot'.")
    # Method used is ggplot2 version-dependent
    ggplot_version <- stringr::str_sub(as.character(utils::packageVersion("ggplot2")),1,1)
    if (ggplot_version == 2) {
        # ggplot2 v2 solution:
        range <- ggplot2::ggplot_build(my.ggp)$layout$panel_ranges[[1]]$y.range
    } else {
        # ggplot2 v3 solution:
        range <- ggplot2::ggplot_build(my.ggp)$layout$panel_params[[1]]$y.range
    }
    return(range)
}

xrange <- function(my.ggp) {
    assertthat::assert_that(class(my.ggp)[[2]] == "ggplot",
                            msg = "my.ggp must be of class 'ggplot'.")
    # Method used is ggplot2 version-dependent
    ggplot_version <- stringr::str_sub(as.character(utils::packageVersion("ggplot2")),1,1)
    if (ggplot_version == 2) {
        # ggplot2 v2:
        range <- ggplot2::ggplot_build(my.ggp)$layout$panel_ranges[[1]]$x.range
    } else {
        # ggplot2 v3 solution:
        range <-  ggplot2::ggplot_build(my.ggp)$layout$panel_params[[1]]$x.range
    }
    return(range)
}


addFootnote <- function(my.ggp,
                        footnoteText,
                        footnoteSize,
                        footnoteColor,
                        footnoteJust = 1,
                        yoffset = 0) {

    yr <- yrange(my.ggp)
    xr <- xrange(my.ggp)
    yoffset <- yoffset * (yr[2] - yr[1])
    xcoord <- ifelse(footnoteJust < 0.50, xr[1], xr[2])
    if (footnoteJust == 0.5) {# Special case = center
        xcoord <- xr[1] + ((xr[2] - xr[1]) / 2)
    }
    my.ggp <- my.ggp +
        annotate("text", label = footnoteText, x = xcoord, y = yr[1] + yoffset,
                 size = footnoteSize,
                 colour = footnoteColor,
                 hjust = footnoteJust,
                 vjust = 1)
}
