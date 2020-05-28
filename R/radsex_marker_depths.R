#' @export
#'
#' @title Plot the results of radsex "subset", "signif", or "process"
#'
#' @description Plot a heatmap of individual depths for a subset of markers from
#' radsex "process", "subset", or "signif". In the resulting heatmap, the color
#' of a tile at coordinates (x, y) indicates the depth of the y-th marker in
#' the x-th individual.
#'
#' @param input_file Path to a table of individual marker depths obtained with
#' radsex "subset", "signif", or "process".
#'
#' @param output_file Path to an output file for the generated marker depths
#' plot, or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 10).
#'
#' @param height Plot height when plotting to an output file, in inches
#' in inches (default: 8).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param title Plot title (default: NA, i.e. no title).
#'
#' @param group_info_file Path to a tabulated file with columns
#' "Individual | Group", or NA to not use group information (default: NA).
#'
#' @param group_labels Vector of length 2 with the labels to associate to each
#' group in the plot legend, or NA to use the group names (default: NA).
#'
#' @param label_colors Vector of length 2 specifying label color for individuals
#' from each group if individual group information is used
#' (default: c("dodgerblue3", "firebrick2)).
#'
#' @param min_depth Minimum depth to consider a marker present in an individual:
#' depths lower than this value will be set to 0 (default: NA, i.e. get value
#' from input data).
#'
#' @param max_depth Maximum marker depth allowed in an individual: depths higher
#' than this value will be set to this value (default: 100).
#'
#' @param distance_method Method to use to compute the distance matrix, directly
#' passed to \code{\link[stats]{dist}}. Possible values: "euclidean", "maximum",
#' "manhattan", "canberra", "binary", and "minkowski" (default: "euclidean")
#'
#' @param clustering_method Method to use for clustering, directly passed to
#' \code{\link[stats]{hclust}}. Possible values: "ward.D", "ward.D2", "single",
#' "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC),
#' and "centroid" (= UPGMC) (default: "ward.D").
#'
#' @param depth_colors Vector of length 5 specifying colors for the following
#' values: [0, (1:mean), (mean:3rd quartile), (3rd quartile:(max - 1)), max]
#' (default: c("white", "royalblue2", "black", "gold2", "red3"))
#'
#' @param color_by_presence If TRUE, color markers by presence / absence using
#' the parameters "presence_color" and "presence_min_depth" (default: FALSE).
#'
#' @param presence_color Color of present markers when coloring by presence /
#' absence of markers (default: "grey30").
#'
#' @param presence_min_depth Minimum depth to consider a marker present in an
#' individual when coloring by presence / absence (default: 5).
#'
#' @param show_individuals If TRUE, display individual names on the x-axis
#' (default: TRUE).
#'
#' @param show_markers If TRUE, display marker names on the y-axis
#' (defautlt: FALSE).
#'
#' @param show_individuals_dendrogram If TRUE, display the individual clustering
#' dendrogram on the x-axis (default: TRUE).
#'
#' @param show_markers_dendrogram If TRUE, display the marker clustering
#' dendrogram on the y-axis (default: FALSE).
#'
#' @examples
#' data <- load_table("subset.tsv")
#' group_info <- load_group_info("popmap.tsv")
#' marker_depths <- draw_marker_depths(data,
#'                                     group_info = group_info,
#'                                     title = "Marker depths clustered",
#'                                     label_colors = c("green", "purple"),
#'                                     show_individuals = TRUE,
#'                                     show_markers = TRUE,
#'                                     show_individuals_dendrogram = TRUE,
#'                                     show_markers_dendrogram = TRUE)

radsex_marker_depths <- function(input_file,
                                 output_file = NA,
                                 width = 10,
                                 height = 8,
                                 res = 300,
                                 title = NA,
                                 group_info_file = NA,
                                 group_labels = NA,
                                 label_colors = c("dodgerblue3", "red3"),
                                 min_depth = NA,
                                 max_depth = 150,
                                 distance_method = "euclidean",
                                 clustering_method = "ward.D",
                                 depth_colors = c("white", "royalblue2", "black",
                                                  "gold2", "red3"),
                                 color_by_presence = FALSE,
                                 presence_color = "grey30",
                                 presence_min_depth = 5,
                                 show_individuals = TRUE,
                                 show_markers = FALSE,
                                 show_individuals_dendrogram = TRUE,
                                 show_markers_dendrogram = FALSE) {

    if (!is.na(group_info_file)) {

        group_info <- load_group_info(group_info_file,
                                      group_labels = group_labels)

    } else {

        group_info <- NA

    }

    data <- load_table(input_file)

    data <- cluster_marker_depths(data,
                                  min_depth = min_depth,
                                  max_depth = max_depth,
                                  distance_method = distance_method,
                                  clustering_method = clustering_method)

    d <- draw_marker_depths(data = data,
                            group_info = group_info,
                            output_file = output_file,
                            width = width,
                            height = height,
                            res = res,
                            title = title,
                            label_colors = label_colors,
                            depth_colors = depth_colors,
                            color_by_presence = color_by_presence,
                            presence_color = presence_color,
                            presence_min_depth = presence_min_depth,
                            show_individuals = show_individuals,
                            show_markers = show_markers,
                            show_individuals_dendrogram = show_individuals_dendrogram,
                            show_markers_dendrogram = show_markers_dendrogram)

    return(d)
}





#' @export
#'
#' @title RADSex marker depths plot
#'
#' @description Draw a heatmap of individual depths for a subset of markers from
#' radsex "process", "subset", or "signif". In the resulting heatmap, the color
#' of a tile at coordinates (x, y) indicates the depth of the y-th marker in
#' the x-th individual.
#'
#' @param data Table of marker depths obtained with the \code{\link{load_table}}
#' function.
#'
#' @param group_info A table of individual group information obtained with the
#' \code{\link{load_group_info}} function, or NA to not use group information
#' (default: NA).
#'
#' @param output_file Path to an output file for the generated marker depths
#' plot, or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 10).
#'
#' @param height Plot height when plotting to an output file, in inches
#' in inches (default: 8).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param title Plot title (default: NA, i.e. no title).
#'
#' @param label_colors Vector of length 2 specifying label color for individuals
#' from each group if individual group information is used
#' (default: c("dodgerblue3", "firebrick2)).
#'
#' @param depth_colors Vector of length 5 specifying colors for the following
#' values: [0, (1:mean), (mean:3rd quartile), (3rd quartile:(max - 1)), max]
#' (default: c("white", "royalblue2", "black", "gold2", "red3"))
#'
#' @param color_by_presence If TRUE, color markers by presence / absence using
#' the parameters "presence_color" and "presence_min_depth" (default: FALSE).
#'
#' @param presence_color Color of present markers when coloring by presence /
#' absence of markers (default: "grey30").
#'
#' @param presence_min_depth Minimum depth to consider a marker present in an
#' individual when coloring by presence / absence (default: 5).
#'
#' @param show_individuals If TRUE, display individual names on the x-axis
#' (default: TRUE).
#'
#' @param show_markers If TRUE, display marker names on the y-axis
#' (defautlt: FALSE).
#'
#' @param show_individuals_dendrogram If TRUE, display the individual clustering
#' dendrogram on the x-axis (default: TRUE).
#'
#' @param show_markers_dendrogram If TRUE, display the marker clustering
#' dendrogram on the y-axis (default: FALSE).
#'
#' @examples
#' data <- load_table("subset.tsv")
#' group_info <- load_group_info("popmap.tsv")
#' marker_depths <- draw_marker_depths(data,
#'                                     group_info = group_info,
#'                                     title = "Marker depths clustered",
#'                                     label_colors = c("green", "purple"),
#'                                     show_individuals = TRUE,
#'                                     show_markers = TRUE,
#'                                     show_individuals_dendrogram = TRUE,
#'                                     show_markers_dendrogram = TRUE)

draw_marker_depths <- function(data,
                               group_info = NA,
                               output_file = NA,
                               width = 10,
                               height = 8,
                               res = 300,
                               title = NA,
                               label_colors = c("dodgerblue3", "red3"),
                               depth_colors = c("white", "royalblue2", "black",
                                                "gold2", "red3"),
                               color_by_presence = FALSE,
                               presence_color = "grey30",
                               presence_min_depth = 5,
                               show_individuals = TRUE,
                               show_markers = FALSE,
                               show_individuals_dendrogram = TRUE,
                               show_markers_dendrogram = FALSE) {

    # Define individual label colors
    individual_names <- data$individuals$labels[data$individuals$order]

    if (is.na(c(group_info)[1])) {

        # Set all labels to black
        labels_palette <- rep("black", length(individual_names))
        names(labels_palette) <- individual_names

    } else {

        # Color labels based on groups according to label_colors
        temp <- stats::setNames(label_colors,
                                c(group_info$groups[1], group_info$groups[2]))
        labels_palette <- temp[group_info$individual_groups[individual_names]]

    }

    if (color_by_presence) {

        # Create a binary presence variable (depth >= min_depth)
        data$data$presence <- data$data$depth >= presence_min_depth

        # Generate basic heatmap for presence of markers
        heatmap <- ggplot2::ggplot(data$data,
                                   ggplot2::aes(x = individual,
                                                y = id,
                                                fill = presence)) +
            ggplot2::theme_bw() +
            ggplot2::scale_fill_manual(name = paste0("Min. depth ", presence_min_depth),
                                       breaks = c("TRUE", "FALSE"),
                                       labels = c("Present", "Absent"),
                                       values=c("TRUE" = presence_color,
                                                "FALSE"="white")) +
            ggplot2::geom_tile(color = "grey10", size = 0.05)

    } else {

        heatmap <- ggplot2::ggplot(data$data,
                                   ggplot2::aes(x = individual,
                                                y = id,
                                                fill = depth)) +
            ggplot2::theme_bw() +
            # Create color scale from a normalized distribution of depth values
            # Distribution vector:
            # Min, 1st Q., Median, Mean, 3rd Q., Max
            ggplot2::scale_fill_gradientn(name = "Depth",
                                          colours = depth_colors,
                                          values = c(0,
                                                     0.000001,
                                                     data$distribution[3] / data$distribution[6],
                                                     data$distribution[5] / data$distribution[6],
                                                     1)) +
            ggplot2::geom_tile(color = "grey30", size = 0.02) +
            ggplot2::theme(legend.title = ggplot2::element_text(size = 14, face = "bold"),
                           legend.key.height = ggplot2::unit(0.1, "npc"),
                           legend.key.width = ggplot2::unit(0.06, "npc"))

    }


    # Compute the main heatmap object
    heatmap <- heatmap +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(5, 15, 15, 30),
                       panel.border = ggplot2::element_rect(size = 0.75,
                                                            color = "black"),
                       legend.position = "right",
                       legend.title = ggplot2::element_text(size = 12, face = "bold"),
                       legend.text = ggplot2::element_text(size = 11)) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0))

    # Add plot title if specified in settings
    if (!is.na(title)) {

        heatmap <- heatmap +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = 20,
                                                              face = "bold",
                                                              margin = ggplot2::margin(0, 0, 10, 0))) +
            ggplot2::ggtitle(title)

    } else {

        heatmap <- heatmap +
            ggplot2::theme(plot.title = ggplot2::element_blank())

    }

    # Add individual labels if specified in settings
    if (show_individuals) {

        # ggplot issues a warning when passing a vector of colors, this is not
        # officially supported but works with ggplot2 v <= 3.3.0 so far.
        # Monitor ggplot changes but I couldn't find another proper way to do it
        heatmap <- heatmap +
            suppressWarnings(ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                                                vjust = 0.5,
                                                                                color = labels_palette,
                                                                                size = 10,
                                                                                margin = ggplot2::margin(2.5,0,0,0))))

    }

    # Add marker labels if specified in settings
    if (show_markers) {

        heatmap <- heatmap +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))

    }

    # Create heatmap "grob" for the combined gtable
    combined <- ggplot2::ggplotGrob(heatmap)

    # Base position of dendrograms in the gtable
    individual_dendrogram_top <- combined$layout$b[which(combined$layout$name == "axis-b")]
    individual_dendrogram_left <- combined$layout$l[which(combined$layout$name == "axis-b")]
    marker_dendrogram_top <- combined$layout$t[which(combined$layout$name == "axis-l")]
    marker_dendrogram_left <- combined$layout$l[which(combined$layout$name == "ylab-l")] - 1

    # Increment top position of individual dendrogram in the gtable if
    # individual labels are displayed.
    if (show_individuals) {

        individual_dendrogram_top <- individual_dendrogram_top + 1

    }

    # Decrement left position of marker dendrogram in the gtable if marker
    # labels are displayed.
    if (show_markers) {

        marker_dendrogram_left <- marker_dendrogram_left - 1

    }

    if (show_individuals_dendrogram) {

        # Compute individual dendrogram object
        individual_dendrogram <- suppressMessages(ggdendro::ggdendrogram(data$individuals,
                                                                         labels = FALSE,
                                                                         leaf_labels = FALSE,
                                                                         theme_dendro = TRUE,
                                                                         rotate = FALSE) +
                                                      ggplot2::theme(plot.margin = grid::unit(c(0.1, 0.01, 0, 0), 'npc'),
                                                                     axis.text.x = ggplot2::element_blank(),
                                                                     axis.text.y = ggplot2::element_blank(),
                                                                     axis.title.x = ggplot2::element_blank()) +
                                                      ggplot2::scale_y_reverse(expand = c(0, 0.5)) +
                                                      ggplot2::scale_x_continuous(expand = c(0, 0)))

        # Add a row to the combined gtable for the dendrogram
        combined <- gtable::gtable_add_rows(combined,
                                            grid::unit(0.1, "npc"),
                                            pos = individual_dendrogram_top)

        # Add the dendrogram to the combined gtable
        combined <- gtable::gtable_add_grob(combined,
                                            ggplot2::ggplotGrob(individual_dendrogram),
                                            t = individual_dendrogram_top,
                                            l = individual_dendrogram_left,
                                            b = individual_dendrogram_top + 1,
                                            r = individual_dendrogram_left + 1)
    }

    if (show_markers_dendrogram) {

        marker_dendrogram <- suppressMessages(ggdendro::ggdendrogram(data$markers,
                                                                     labels = FALSE,
                                                                     leaf_labels = FALSE,
                                                                     theme_dendro = TRUE,
                                                                     rotate = FALSE) +
                                                    ggplot2::theme(plot.margin = grid::unit(c(0.0, 0.1, 0.0, 0), 'npc'),
                                                                   axis.text.x = ggplot2::element_blank(),
                                                                   axis.text.y = ggplot2::element_blank(),
                                                                   axis.title.x = ggplot2::element_blank()) +
                                                    ggplot2::coord_flip() +
                                                    ggplot2::scale_y_reverse(expand = c(0, 0)) +
                                                    ggplot2::scale_x_continuous(expand = c(0, 0.5)))

        # Add a row to the combined gtable for the dendrogram
        combined <- gtable::gtable_add_cols(combined,
                                            grid::unit(0.04, "npc"),
                                            pos = 0)

        # Add the dendrogram to the combined gtable
        combined <- gtable::gtable_add_grob(combined,
                                            ggplot2::ggplotGrob(marker_dendrogram),
                                            t = marker_dendrogram_top,
                                            l = marker_dendrogram_left,
                                            b = marker_dendrogram_top,
                                            r = marker_dendrogram_left + 1)
    }

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = combined,
                        width = width,
                        height = height,
                        dpi = res)

    } else {

        plot(combined)

    }

    return(combined)

}





#' @title Cluster marker depths data
#'
#' @description Cluster individuals and markers by depth values in a table of
#' marker depths.
#'
#' @param data A table of marker depths obtained with the
#' \code{\link{load_table}} function.
#'
#' @param min_depth Minimum depth to consider a marker present in an individual:
#' depths lower than this value will be set to 0 (default: NA, i.e. get value
#' from input data).
#'
#' @param max_depth Maximum marker depth allowed in an individual: depths higher
#' than this value will be set to this value (default: 100).
#'
#' @param distance_method Method to use to compute the distance matrix, directly
#' passed to \code{\link[stats]{dist}}. Possible values: "euclidean", "maximum",
#' "manhattan", "canberra", "binary", and "minkowski" (default: "euclidean")
#'
#' @param clustering_method Method to use for clustering, directly passed to
#' \code{\link[stats]{hclust}}. Possible values: "ward.D", "ward.D2", "single",
#' "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC),
#' and "centroid" (= UPGMC) (default: "ward.D").
#'
#' @return A list with the following elements:
#' \item{data}{A data frame of marker depths with individuals and markers
#' ordered based on the clustering results}
#' \item{individuals}{Individuals clustering results}
#' \item{markers}{Markers clustering results}
#' \item{distribution}{Distribution of depth values}
#'
#' @examples
#' clustering_data <- cluster_marker_depths(data,
#'                                          min_depth = 1,
#'                                          max_depth = 100,
#'                                          distance_method = "binary",
#'                                          clustering_method = "complete")

cluster_marker_depths <- function(data,
                                  min_depth = NA,
                                  max_depth = 150,
                                  distance_method = "euclidean",
                                  clustering_method = "ward.D") {

    properties <- data$properties
    data <- data$data

    # If min. depth was not specified, get it from data properties or set it
    # to 1 if there was no min_depth value in the properties
    if (is.na(min_depth)) {

        if ("min_depth" %in% names(properties)) {

            min_depth <- as.numeric(properties$min_depth)

        } else {

            warning("Did not find \"min_depth\" in input file properties,
                    setting minimum depth to 1")

            min_depth <- 1

        }

    }

    # Remove the marker column which is not useful here
    data <- drop_columns(data, "sequence")

    # Count numbers of individuals and markers
    n_individuals <- ncol(data) - 1  # First column is marker id
    n_markers <- nrow(data)

    # Extract depth values in a matrix for clustering
    depth <- as.matrix(data[, -1], rownames.force = TRUE)
    rownames(depth) <- data$id

    # Set depth value boundaries according to min_depth and max_depth
    depth[which(depth > max_depth)] <- max_depth
    depth[which(depth < min_depth)] <- 0
    data[, -1] <- depth

    #### INDIVIDUAL CLUSTERING
    # Compute distances
    distances <- dist(t(depth), method = distance_method)
    # Cluster individuals based on distances
    individual_clusters <- hclust(distances, method = clustering_method)
    # Reorder marker depth columns in the data based on clustering results
    data <- data[, c(1, individual_clusters$order + 1)]
    # Order individuals based on clustering results
    individuals <- individual_clusters$labels[individual_clusters$order]

    #### MARKERS CLUSTERING
    # Compute distances
    distances <- dist(depth, method = distance_method)
    # Cluster markers based on distances
    marker_clusters <- hclust(distances, method = clustering_method)
    # Reorder marker depth rows in the data based on clustering results
    data <- data[marker_clusters$order,]
    # Order markers based on clustering results
    markers <- marker_clusters$labels[marker_clusters$order]

    # Prepare data for ggplot
    melted <- suppressMessages(reshape2::melt(data,
                                              id.vars = c("id"),
                                              variable.name = "individual",
                                              value.name = "depth"))
    melted$id <- factor(melted$id, levels = markers)
    melted$individual <- factor(as.character(melted$individual),
                                levels = individuals)

    # Compute quantiles for color scale
    distribution <- summary(replace(melted$depth,
                                    which(melted$depth == 0),
                                    NA),
                            na.rm=TRUE)

    # Generate output list
    output <- list(data = melted,
                   individuals = individual_clusters,
                   markers = marker_clusters,
                   distribution = distribution)

    return(output)
}

