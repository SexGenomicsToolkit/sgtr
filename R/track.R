#' @export
#'
#' @title Create a track object
#'
#' @description Generate an object storing all properties for a plot track
#'
#' @param metrics Metrics to represent in the track; metric names are given by
#' column names in the input data frame (i.e. the output of the
#' \code{\link{load_genome_metrics}} function). Possible values: a string if
#' the track represents a single metric (e.g. "Fst"), or a vector if the track
#' represents several metrics (e.g. c("Snps_females", "Snps_males")).
#'
#' @param label Label to use for the track, will be displayed on the y axis. Can
#' be set to NA for single metric tracks to set the label to the name of the
#' metric, but label has to be specified for multi-metrics tracks.
#'
#' @param bg_colors Background colors for circos and manhattan plots. Values can
#' be a string (e.g. "white") and will be applied to all sectors / chromosomes,
#' or a vector of colors (e.g. c("white", "grey80")) to assign
#' alternating values to sectors / chromosomes. Not used in region plots.
#'
#' @param ylim Y-axis limits, a vector with values c(min, max).
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis.
#' The value for this parameter is directly passed to "major.lines.y" in
#' \code{\link{ggplot2::theme}} for manhattan and region plots.
#'
#' @param major_lines_x If TRUE, reference lines will be plotted for the x axis.
#' The value for this parameter is directly passed to "major.lines.x" in
#' \code{\link{ggplot2::theme}} for region plots. Not used in manhattan plots.
#'
#' @param legend_position Position of the legend, directly passed to
#' "legend.position" in \code{\link{ggplot2::theme}} for region plots. Not used
#' in manhattan and in circos plots.
#'
#' @param h_lines List of objects defined with \code{\link{h_line}} adding
#' horizontal lines to the track (default: NA).
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Single metric
#' track_data <- track("Fst", color = "grey70", point_size = 0.75)
#'
#' # Multiple metrics
#' track_data <- track(c("Snp_females", "Snp_males"), label = "SNPs",
#'                     colors = c("firebrick2", "dodgerblue3"))

track <- function(metrics,
                  label = NA,
                  bg_colors = NA,
                  ylim = NA,
                  major_lines_x = NA,
                  major_lines_y = NA,
                  legend_position = NA,
                  h_lines = NA) {

    n_metrics <- length(metrics)

    # Assign label if necessary
    if (is.na(list(label))) {

        if (n_metrics == 1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    # Create track object from specified values
    track <- list(metrics = metrics,
                  label = label,
                  bg_colors = bg_colors,
                  ylim = ylim,
                  major_lines_x = major_lines_x,
                  major_lines_y = major_lines_y,
                  legend_position = legend_position,
                  h_lines = h_lines)

    return(track)

}





#' @export
#'
#' @title Create a track object with a single metric
#'
#' @description Simplified interface to generate a track representing a
#' single metric. All parameters have default values as NA, meaning their value
#' will be taken from default values specified in high-level functions.
#'
#' @param metric Metric to represent in the track, should correspond to a
#' column name in the input data frame (i.e. the output of the
#' \code{\link{load_genome_metrics}} function).
#'
#' @param colors Either a single color (e.g. "grey20") which will be
#' applied to all values, a vector of colors (e.g. c("blue", "yellow")) which
#' will be alternatingly applied to linkage groups / chromosomes, or a function
#' returning a color for a given value.
#'
#' @param point_size Point size for plots of type "points", a float.
#'
#' @param alpha Alpha value, a float.
#'
#' @param type Plot type for circos and region tracks, either "points" or
#' "ribbon"; not used in manhattan plots.
#'
#' @param label Label to use for the track, will be displayed on the y axis. If
#' NA, will use the metric name as label.
#'
#' @param bg_colors Background colors for circos and manhattan plots, not used
#' in region plots. Values can be a string (e.g. "white") and will be applied
#' to all sectors / chromosomes, or a vector of strings
#' (e.g. c("white", "grey80")) to assign alternating values to
#' sectors / chromosomes.
#'
#' @param ylim Y-axis limits, a vector with values c(min, max) or NA to infer
#' y limits from the data.
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis,
#' equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}.
#'
#' @param major_lines_x If TRUE, reference lines will be plotted for the x axis,
#' equivalent to panel.grid.major.x in \code{\link{ggplot2::theme}}.
#'
#' @param legend_position Position of the legend, directly passed to
#' "legend.position" in \code{\link{ggplot2::theme}} in region plots, not used
#' in manhattan and in circos plots.
#'
#' @param h_lines List of objects defined with \code{\link{h_line}} adding
#' horizontal lines to the track (default: NA).
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Create a track object for the metric "Fst" in grey with a point size of
#' # 0.75 and y-axis limits of (0, 1).
#' fst_track <- single_metric_track("Fst",
#'                                  colors = "grey70",
#'                                  point_size = 0.75,
#'                                  ylim = c(0, 1))

single_metric_track <- function(metric,
                                colors = NA,
                                point_size = NA,
                                alpha = NA,
                                type = NA,
                                label = NA,
                                bg_colors = NA,
                                ylim = NA,
                                major_lines_x = NA,
                                major_lines_y = NA,
                                legend_position = NA,
                                h_lines = NA) {

    # Assign label if necessary
    if (is.na(list(label))) { label <- metric }

    metrics <- list(metric(metric,
                           colors = colors,
                           point_size = point_size,
                           alpha = alpha,
                           type = type))

    # Create track object from specified values
    track <- track(metrics = metrics,
                   label = label,
                   bg_colors = bg_colors,
                   ylim = ylim,
                   major_lines_x = major_lines_x,
                   major_lines_y = major_lines_y,
                   legend_position = legend_position,
                   h_lines = h_lines)

    return(track)

}





#' @export
#'
#' @title Create a track object with multiple metrics
#'
#' @description Simplified interface to generate a track representing multiple
#' metric.s All parameters have default values as NA, meaning their value
#' will be taken from default values specified in high-level functions.
#'
#' @param metrics Vector of metrics to represent in the track, should correspond
#' to column names in the input data frame (i.e. the output of the
#' \code{\link{load_genome_metrics}} function).
#'
#' @param metric_labels Vector of metric labels, will be used in plot legends.
#' If NA, metric names will be used (default: NA).
#'
#' @param colors Either a single color (e.g. "grey20") which will be
#' applied to all values, a vector of colors (e.g. c("blue", "yellow")) which
#' will be alternatingly applied to linkage groups / chromosomes, or a function
#' returning a color for a given value.
#'
#' @param point_size Point size for plots of type "points", can be a single
#' float value (e.g. 0.5) and will be applied to all metrics, or a vector of
#' floats of length equal to the number of metrics (e.g. c(0.5, 1, 0.75) for
#' three metrics).
#'
#' @param alpha Alpha value, can be a single float value (e.g. 0.5) and will be
#' applied to all metrics, or a vector of floats of length equal to the number
#' of metrics (e.g. c(0.5, 0.8) for two metrics).
#'
#' @param type Plot type for circos and region tracks, either "points" or
#' "ribbon"; can be a single value and will be applied to all metrics, or a
#' vector of strings of length equal to the number of metrics
#' (e.g. c("points", "points", "ribbon") for three metrics).
#'
#' @param label Label to use for the track, will be displayed on the y axis.
#'
#' @param bg_colors Background colors for circos and manhattan plots, not used
#' in region plots. Values can be a string (e.g. "white") and will be applied
#' to all sectors / chromosomes, or a vector of strings
#' (e.g. c("white", "grey80")) to assign alternating values to
#' sectors / chromosomes.
#'
#' @param ylim Y-axis limits, a vector with values c(min, max) or NA to infer
#' y limits from the data.
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis,
#' equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}.
#'
#' @param major_lines_x If TRUE, reference lines will be plotted for the x axis,
#' equivalent to panel.grid.major.x in \code{\link{ggplot2::theme}}.
#'
#' @param legend_position Position of the legend, directly passed to
#' "legend.position" in \code{\link{ggplot2::theme}} in region plots, not used
#' in manhattan and in circos plots.
#'
#' @param h_lines List of objects defined with \code{\link{h_line}} adding
#' horizontal lines to the track (default: NA).
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Create a track object labeled "SNPs" for the metrics "Snp_females" in blue
#' # and "Snp_males" in red plotted as ribbons with alpha value of 0.5.
#' snps_track <- track(c("Snp_females", "Snp_males"),
#'                     label = "SNPs",
#'                     colors = c("firebrick2", "dodgerblue3"),
#'                     type = "ribbon",
#'                     alpha = 0.5)

multi_metrics_track <- function(metrics,
                                metric_labels = NA,
                                colors = NA,
                                point_size = NA,
                                alpha = NA,
                                type = NA,
                                label = NA,
                                bg_colors = NA,
                                ylim = NA,
                                major_lines_x = NA,
                                major_lines_y = NA,
                                legend_position = NA,
                                h_lines = NA) {

    n_metrics <- length(metrics)
    metric_names <- metrics

    # Create an empty list of the right length
    metrics <- vector(mode = "list", length = n_metrics)

    if (is.na(c(metric_labels)[1])) {

        metric_labels <- rep(NA, n_metrics)

    }

    for (i in 1:n_metrics) {

        metrics[[i]] <- metric(metric_names[i],
                               metric_labels[i],
                               NA, NA, NA, NA)

    }

    # Assign property values for all metrics
    metrics <- assign_values(metrics, "colors", colors)
    metrics <- assign_values(metrics, "point_size", point_size)
    metrics <- assign_values(metrics, "alpha", alpha)
    metrics <- assign_values(metrics, "type", type)

    # Create track object from specified values
    track <- track(metrics = metrics,
                   label = label,
                   bg_colors = bg_colors,
                   ylim = ylim,
                   major_lines_x = major_lines_x,
                   major_lines_y = major_lines_y,
                   legend_position = legend_position,
                   h_lines = h_lines)

    return(track)

}





#' @export
#'
#' @title Create a metric object
#'
#' @description Generate an object storing all properties for a metric
#'
#' @param metric Metric name, should match a column in the input data frame
#' (i.e. the output of the \code{\link{load_genome_metrics}} function).
#'
#' @param label Metric label, will be used in plot legends. If NA, the metric
#' name will be used (default: NA).
#'
#' @param colors Colors, either be a string (e.g. "grey20") or a vector of
#' strings (e.g. c("blue", "yellow") (default: "grey60").
#'
#' @param point_size Point size for plots of type "points", a float
#' (default: 0.5).
#'
#' @param alpha Alpha value, a float (default: 1).
#'
#' @param type Plot type for circos and region tracks, either "points" or
#' "ribbon"; not used in manhattan plots (default: "points").
#'
#' @return A named list with the value of each metric property
#'
#' @examples
#' fst <- metric("Fst", color = "grey70", point_size = 0.75)

metric <- function(metric,
                   label = NA,
                   colors = "grey60",
                   point_size = 0.5,
                   alpha = 1,
                   type = "points") {

    if (is.na(label)) { label <- metric }

    metric <- list(name = metric,
                   label = label,
                   colors = colors,
                   point_size = point_size,
                   alpha = alpha,
                   type = type)

    return(metric)

}





#' @export
#'
#' @title Create a h_line object
#'
#' @description Generate an object storing all properties for a horizontal
#' line
#'
#' @param y Y-axis coordinate for the line.
#'
#' @param label Label to add to the line, or NA for no label (default: NA).
#'
#' @param label_x Normalized x coordinate for the label (default: 0.1).
#'
#' @param label_font_size Font size for label (default: 7).
#'
#' @param color Line color (default: "grey60").
#'
#' @param type Line type as defined in R "lty" parameter (default: 1).
#'
#' @param size Line thickness (default: 1).
#'
#' @return A named list with the value of each line property
#'
#' @examples
#' signif_line <- h_line("p>0.05", color = "red")

h_line <- function(y,
                   label = NA,
                   label_x = 0.1,
                   label_font_size = 7,
                   color = "grey60",
                   type = 1,
                   size = 1) {

    h_line <- list(y = y,
                   label = label,
                   label_x = label_x,
                   label_font_size = label_font_size,
                   color = color,
                   type = type,
                   size = size)

    return(h_line)

}




#' @title Configure track data from default values
#'
#' @description Wrapper around several track value assignment functions to
#' assign default values and generate track data.
#'
#' @param track A track object generated with the \code{\link{track}} function).
#'
#' @param defaults A named list with default values for all properties.
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_metrics}} function).
#'
#' @param region Region to plot, output of the \code{\link{parse_region}}
#' function, or NA to retain the entire genome (default: NA)
#'
#' @return A track object with values for all properties and track data stored
#' in track$data.
#'
#' @examples
#' track <- track("SNP_males", label = "Male-specific SNPs")
#' defaults <- list(colors = "red", alpha = 0.5)
#' data <- load_genome_metrics("window.tsv")
#' track <- configure_track(track, defaults, data)

configure_track <- function(track, defaults, data, region = NA) {

    # Assign default values to track properties when not set by user
    track <- assign_defaults(track, defaults)

    # Create input data for track
    track <- create_track_data(track, data, region = region)

    return(track)

}





#' @title Assign default values to a track object
#'
#' @description Assign default values to all properties for which the value
#' was not specified by the user (e.g. value is NA).
#'
#' @param track A track object generated with the \code{\link{track}} function)
#'
#' @param defaults A named list with default values for all properties.
#'
#' @return A track object with values for all properties.
#'
#' @examples
#' track <- track("Fst", label = "Window FST")
#' defaults <- list(colors = "red", alpha = 0.5)
#' track <- assign_defaults(track, defaults)

assign_defaults <- function(track, defaults) {

    properties <- names(track)
    default_properties <- names(defaults)
    n_metrics <- length(track$metrics)

    # Iterate through all properties in a track
    for (i in 1:length(properties)) {

        property <- properties[i]
        # Transform to vector and take first element because NA behave
        # strangely in containers
        current_value <- c(track[[property]])[1]

        # Assign a value to properties not defined in track creation (NA value)
        # Use 'list' because checking if an expression is NA triggers a warning,
        # which does not happen when the expression is part of a list
        if (is.na(list(current_value)) & property %in% default_properties) {

            track[[property]] <- defaults[[property]]

        }

    }

    track$metrics <- assign_values(track$metrics, "colors", defaults$colors)
    track$metrics <- assign_values(track$metrics, "point_size",
                                   defaults$point_size)
    track$metrics <- assign_values(track$metrics, "alpha", defaults$alpha)
    track$metrics <- assign_values(track$metrics, "type", defaults$type)

    return(track)

}





#' @title Create track data
#'
#' @description Create an input data frame for the track plotting functions
#'
#' @param track Track object for the current plot, generated with the
#' \code{\link{track}} function.
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_metrics}} function).
#'
#' @param region Region to plot, output of the \code{\link{parse_region}}
#' function, or NA to retain the entire genome (default: NA)
#'
#' @return A data frame with columns:
#' Contig | Position | Metric 1 | Metric 1 colors |  Metric N | Metric N colors
#'
#' @examples
#' genomic_data <- load_genome_metrics("psass_window.tsv")
#' track <- track("Fst")
#' track <- create_track_data(genomic_data, track)

create_track_data <- function(track, data, region = NA) {

    metrics <- c()
    for (i in 1:length(track$metrics)) {
        metrics <- c(metrics, track$metrics[[i]]$name)
    }

    # Extract required columns and create color columns
    data <- data[, c("Contig_plot", "Position_plot", metrics)]

    names(data) <- c("Contig", "Position", metrics)

    # Combine data for multiple metrics
    data <- reshape2::melt(data, measure.vars = metrics,
                           variable.name="Metric", value.name = "Value")

    # Assign color to data points
    for (i in 1:length(metrics)) {

        metric <- track$metrics[[i]]

        colors <- metric$colors

        if (is.function(colors)) {

            indices <- which(data$Metric == metric$name)
            color_values <- metric$colors(data$Value[indices])
            data$Color[indices] <- color_values

        } else if (length(colors) == 1) {

            # Only one color, apply it to all values in metric
            data$Color[which(data$Metric == metric$name)] <- colors

        } else {

            # Multiple colors, create alternating values for LG / chromosomes
            n_contigs <- length(unique(data$Contig))
            color_values <- stats::setNames(rep(colors, n_contigs)[1:n_contigs],
                                            unique(data$Contig))
            data$Color <- color_values[data$Contig]

        }

    }

    # Sort data by Contig then Position
    data <- data[order(data$Contig, data$Position), ]

    # Remove row names
    rownames(data) <- c()

    # Extract region if specified
    if (is.list(region)) {

        data <- subset(data,
                       data$Contig == region$contig &
                           data$Position >= region$start &
                           data$Position <= region$end)

    }

    track$data <- data

    return(track)

}






#' @title Assign value to track property
#'
#' @description Assign value to a single track property, handling single vs
#' multi metrics value definitions.
#'
#' @param metrics Named list of metrics properties/
#'
#' @param colors Colors, can be a string (e.g. "grey20") and will be
#' applied to all metrics, a vector of size equal to the number of metrics
#' (e.g. c("red", "blue") for two metrics), or a color scale object <TBD>.
#'
#' @param ylim Y-axis limits, a vector with values c(min, max).
#'
#' @param point_size Point size for plots of type "points". Values can be a
#' float (e.g. 0.5) and will be applied to all metrics, or a vector of
#' size equal to the number of metrics (e.g. c(1, 1.5, 3) for three metrics).
#'
#' @param alpha Alpha value. Values can be a float (e.g. 0.5) and will
#' be applied to all metrics, or a vector of size equal to the number of metrics
#' (e.g. c(1, 0.5, 0.75) for three metrics).
#'
#' @param type Plot type for circos and region plots, either "ribbon" or
#' "points". Values can be a string (e.g. "ribbon") and will be applied to all
#' metrics, or a vector of size equal to the number of metrics
#' (e.g. c("points", "ribbon") for two metrics). Not used in manhattan plots.
#'
#' @return The correct value to assign to the property.
#'
#' @examples
#' track <- single_metric_track("Fst")
#' assign_value(track$metrics, "alpha", 0.6)

assign_values <- function(metrics, property, values) {

    n_metrics <- length(metrics)

    if (n_metrics == 1) {

        if (is.na(c(metrics[[1]][[property]])[1])) {

            metrics[[1]][[property]] <- values

        }

    } else {

        if (length(values) == 1) {

            for (i in 1:n_metrics) {

                if (is.na(metrics[[i]][[property]])) {

                    metrics[[i]][[property]] <- values

                }

            }

        } else if (length(values) == n_metrics) {

            for (i in 1:n_metrics) {

                if (is.na(metrics[[i]][[property]])) {

                    metrics[[i]][[property]] <- values[i]

                }

            }

        } else {

            stop(paste0("Incorrect value \"", values, "\" for property \"",
                        property, "\" in track definition."))

        }

    }

    return(metrics)

}