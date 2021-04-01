#' @export
#'
#' @title Plot the results of radsex "distrib" analysis
#'
#' @description Plot a distribution of markers between two groups obtained with
#' radsex "distrib". In the resulting plot, the color of a tile at coordinates
#' (x, y) indicates the number of markers present in x individuals from group1
#' and y individuals from group2.
#'
#' @param input_file Path to a table of distribution of markers between groups
#' obtained with radsex "distrib".
#'
#' @param groups Vector of length 2 with the groups to plot on the x and y axes,
#' should match column headers in the input data. If NA, groups will be infered
#' from the input data (default: NA).
#'
#' @param group_labels Vector of length 2 with the labels to associate to each
#' group in the plot axis titles, or NA to use the group names (default: NA).
#'
#' @param output_file Path to an output file for the generated distrib plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 8).
#'
#' @param height Plot height when plotting to an output file, in inches
#' (default: 7).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param autoscale If TRUE, the width and height of the output file will be
#' adjusted so that horizontal and vertical scales are the same. The resulting
#' plot will not be square if there are different number of individuals in each
#' group (default: TRUE).
#'
#' @param title Plot title (default: NA, i.e. no title).
#'
#' @param show_significance If TRUE, tiles for which association with group is
#' significant are highlighted using the color defined with the
#' "significance_color" parameter (default: TRUE).
#'
#' @param significance_color Color of the border for tiles in which association
#' with group is significant (default: "red3").
#'
#' @param significance_thickness Thickness of the border for tiles in which
#' association with group is significant (default: 0.5).

#' @param bins A vector of values to use as bins in the color palette
#' (default: c(0, 1, 5, 25, 100, 1000)).
#'
#' @param colors A vector of two colors used to create the color palette
#' gradient (default: c("white", "navyblue")).
#'
#' @return A ggplot object for the tile plot
#'
#' @examples
#' distrib_plot <- radsex_distrib("distrib.tsv",
#'                                groups = c("M", "F"),
#'                                group_labels = c("Males", "Females"),
#'                                title = "Distribution of markers",
#'                                significance_color = "green",
#'                                bins = c(0, 10, 100, 1000),
#'                                colors = c("white", "red3"))

radsex_distrib <- function(input_file,
                           groups = NA,
                           group_labels = NA,
                           output_file = NA,
                           width = 8,
                           height = 7,
                           res = 300,
                           autoscale = TRUE,
                           title = NA,
                           show_significance = TRUE,
                           significance_color = "red3",
                           significance_thickness = 0.5,
                           bins = c(0, 1, 5, 25, 100, 1000),
                           colors = c("white", "navyblue")) {

    # Load the data
    data <- load_marker_distribution(input_file,
                                     groups = groups,
                                     group_labels = group_labels)

    t <- draw_distrib(data,
                      output_file = output_file,
                      width = width,
                      height = height,
                      res = res,
                      autoscale = autoscale,
                      title = title,
                      show_significance = show_significance,
                      significance_color = significance_color,
                      significance_thickness = significance_thickness,
                      bins = bins,
                      colors = colors)

    return(t)
}





#' @export
#'
#' @title Draw a tile plot for radsex "distrib" results
#'
#' @description Draw a tile plot of the distribution of markers between two
#' groups. In the resulting plot, the color of a tile at coordinates (x, y)
#' indicates the number of markers present in x individuals from group1 and
#' y individuals from group2.
#'
#' @param data Distribution of markers between groups obtained with the
#' \code{\link{load_marker_distribution}} function.
#'
#' @param output_file Path to an output file for the generated distrib plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param height Plot height when plotting to an output file, in inches
#' (default: 4).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param autoscale If TRUE, the width and height of the output file will be
#' adjusted so that horizontal and vertical scales are the same. The resulting
#' plot will not be square if there are different number of individuals in each
#' group (default: TRUE).
#'
#' @param title Plot title (default: NA, i.e. no title).
#'
#' @param show_significance If TRUE, tiles for which association with group is
#' significant are highlighted using the color defined with the
#' "significance_color" parameter (default: TRUE).
#'
#' @param significance_color Color of the border for tiles in which association
#' with group is significant (default: "red3").
#'
#' @param significance_thickness Thickness of the border for tiles in which
#' association with group is significant (default: 0.5).

#' @param bins A vector of values to use as bins in the color palette
#' (default: c(0, 1, 5, 25, 100, 1000)).
#'
#' @param colors A vector of two colors used to create the color palette
#' gradient (default: c("white", "navyblue")).
#'
#' @return A ggplot object for the tile plot
#'
#' @examples
#' distrib_plot <- draw_distrib(data,
#'                              title = "Distribution of markers",
#'                              significance_color = "green",
#'                              bins = c(0, 10, 100, 1000),
#'                              colors = c("white", "red3"))

draw_distrib <- function(data,
                         output_file = NA,
                         width = 8,
                         height = 7,
                         res = 300,
                         autoscale = TRUE,
                         title = NA,
                         show_significance = TRUE,
                         significance_color = "red3",
                         significance_thickness = 0.5,
                         bins = c(0, 1, 5, 25, 100, 1000),
                         colors = c("white", "navyblue")) {

    group1 <- data$groups[1]
    group2 <- data$groups[2]

    # Check that there are at least two bins in user-defined bins
    if (length(bins) < 2) {

        stop(paste0("Incorrect number of values (",
                    length(bins), ") in parameter \"bins\" from ",
                    "\"draw_distrib\": at least two values are required"))

    }

    # Check that there are exactly two colors in user-defined colors
    if (length(colors) != 2) {

        stop(paste0("Incorrect number of values (",
                    length(bins), ") in parameter \"colors\" from ",
                    "\"draw_distrib\": exactly two values are required"))

    }

    # Generate color palette from the user-defined parameters
    color_palette <- generate_color_palette(bins, colors)

    # Associate the corresponding bin to each row of the data frame
    data_bins <- lapply(data$data$Markers,
                        function(x) names(color_palette)[tail(which(x >= bins),
                                                              n=1)])

    data$data$Bin <- factor(unlist(data_bins), levels = names(color_palette))

    # Remove significance for tiles without markers
    data$data$Signif <- as.factor(as.logical(data$data$Signif) &
                                      (data$data$Markers > 0))

    # Interval for labels on the x and y axis (to avoid crowded axes)
    x_label_interval <- round(data$counts[group1] / 5)
    y_label_interval <- round(data$counts[group2] / 5)

    # Generate the base plot
    tile_plot <- ggplot2::ggplot(data$data, ggplot2::aes(x = data$data[[group1]],
                                                         y = data$data[[group2]])) +
        ggplot2::geom_tile(ggplot2::aes(fill = Bin),
                           color = "grey50",
                           size = 0.1) +
        ggplot2::theme_bw() +
        # Font sizes work well for plot widths between 6 and 12. In the future,
        # try to implement a system to scale font size to plot width.
        ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5),
                       panel.border = ggplot2::element_rect(size = 0.5,
                                                            color = "black"),
                       panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_text(size = 16,
                                                         color = "black",
                                                         face = "bold"),
                       axis.title.x = ggplot2::element_text(size = 18,
                                                            face = "bold",
                                                            margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(size = 18,
                                                            face = "bold",
                                                            margin = ggplot2::margin(0, 10, 0, 0)),
                       legend.margin = ggplot2::margin(0, 0, 10, 0),
                       legend.title = ggplot2::element_text(size = 14,
                                                            face = "bold"),
                       legend.text = ggplot2::element_text(size = 11),
                       legend.key.height = ggplot2::unit(0.05, "npc"),
                       legend.key.width = ggplot2::unit(0.05, "npc"),
                       legend.key = ggplot2::element_rect(size = 0.5,
                                                          color = "grey80"),
                       legend.position = "right",
                       legend.text.align = 0) +
        ggplot2::scale_fill_manual(name = "Markers",
                                   breaks = names(color_palette),
                                   values = color_palette,
                                   labels = names(color_palette),
                                   drop = FALSE) +
        ggplot2::scale_x_continuous(name = paste0("Number of ", data$group_labels[group1]),
                                    breaks = seq(0,
                                                 data$counts[group1],
                                                 x_label_interval),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(name = paste0("Number of ", data$group_labels[group2]),
                                    breaks = seq(0,
                                                 data$counts[group2],
                                                 y_label_interval),
                                    expand = c(0, 0))

    # Add highlight to significant tiles if specified
    if (show_significance) {

            tile_plot <- tile_plot +
            ggplot2::geom_tile(data = data$data,
                               ggplot2::aes(x = data$data[[group1]],
                                            y = data$data[[group2]],
                                            color = Signif),
                               fill = "NA",
                               size = significance_thickness) +
            ggplot2::scale_color_manual(name = ggplot2::element_blank(),
                                        values = c("TRUE"=significance_color,
                                                   "FALSE"="NA",
                                                   color_palette),
                                        breaks = c("TRUE"),
                                        labels = c("Signif."))

    }

    # Add title if specified
    if (!is.na(title)) {

        tile_plot <- tile_plot +
            ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = 20,
                                                              face = "bold",
                                                              margin = ggplot2::margin(0, 0, 10, 0)))

    } else {

        tile_plot <- tile_plot +
            ggplot2::theme(plot.title = ggplot2::element_blank())

    }

    # Adjust plot dimensions if autoscale was set
    if (autoscale) {

        ratio <- data$counts[group2] / data$counts[group1]

        if (is.na(title)) {

            height <- width * (ratio - 0.05)

        } else {

            height <- width * ratio

        }

        # Adjust height to account for the legend (~ 0.1 total plot width)
        height <- 0.9 * height

    }

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = tile_plot,
                        width = width,
                        height = height,
                        dpi = res)

    } else {

        print(tile_plot)

    }

    return(tile_plot)

}



#' @title Generate a color palette for "distrib" plot
#'
#' @description Generate a color palette from the specified bins and colors to
#' use in radsex "distrib" results plot.
#'
#' @param bins A vector of values to use as bins in the color palette.
#'
#' @param colors A vector of two colors used to create the color palette
#' gradient.
#'
#' @return A color palette vector.


generate_color_palette <- function(bins, colors) {

    bin_labels <- c()

    # Create names for all bins except the last.
    # Names are the bin value if the bin has size 1, or an interval if the
    # bin has size > 1
    for (i in 1:(length(bins) - 1)) {

        if (bins[i] == bins[i + 1] - 1) {

            temp <- as.character(bins[i])

        } else {

            temp <- paste(as.character(bins[i]),
                          as.character(bins[i + 1] - 1), sep = "-")

        }

        bin_labels <- c(bin_labels, temp)

    }

    # Last bin is special ("> last bin")
    bin_labels <- c(bin_labels,
                    paste0(">", as.character(bins[length(bins)])))

    # Generate color palette from the bins and the color values
    palette <- setNames(grDevices::colorRampPalette(colors)(length(bin_labels)),
                        bin_labels)

    return(palette)
}
