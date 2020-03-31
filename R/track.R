#' @title Create a track object
#'
#' @description Generate an object storing all properties for a plot track
#'
#' @param metrics Metrics to represent in the track; metric names are given by
#' column names in the input data frame (i.e. the output of the
#' \code{\link{load_genome_input}} function). Possible values: a string if
#' the track represents a single metric (e.g. "Fst"), or a vector if the track
#' represents several metrics (e.g. c("Snps_females", "Snps_males")).
#'
#' @param label Label to use for the track, will be displayed on the y axis. Can
#' be set to NA for single metric tracks to set the label to the name of the
#' metric, but label has to be specified for multi-metrics tracks.
#'
#' @param colors Colors, can be a string (e.g. "grey20") and will be
#' applied to all metrics, a vector of size equal to the number of metrics
#' (e.g. c("red", "blue") for two metrics), or a color scale object <TBD>.
#'
#' @param bg_colors Background colors for circos and manhattan plots. Values can
#' be a string (e.g. "white") and will be applied to all sectors / chromosomes,
#' or a vector of colors (e.g. c("white", "grey80")) to assign
#' alternating values to sectors / chromosomes. Not used in region plots.
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
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis.
#' The value for this parameter is directly passed to "major.lines.y" in
#' \code{\link{ggplot2::theme}} for manhattan and region plots.
#'
#' @param major_lines_x If TRUE, reference lines will be plotted for the x axis.
#' The value for this parameter is directly passed to "major.lines.x" in
#' \code{\link{ggplot2::theme}} for manhattan and region plots.
#'
#' @param legend_position Position of the legend, directly passed to
#' "legend.position" in \code{\link{ggplot2::theme}} for region plots. Not used
#' in manhattan and in circos plots.
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

track <- function(metrics, label = NA, colors = NA, bg_colors = NA,
                  ylim = NA, point_size = NA, type = NA,
                  alpha = NA, major_lines_x = NA, major_lines_y = NA,
                  legend_position = NA) {

    n_metrics <- length(metrics)

    # Assign label if necessary
    if (is.na(label)) {

        if (n_metrics == 1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    # Create track object from specified values
    track <- list(metrics = metrics,
                  label = label,
                  colors = assign_value(colors, "colors", n_metrics),
                  bg_colors = assign_value(bg_colors, "bg_colors", n_metrics),
                  point_size = assign_value(point_size, "point_size",
                                            n_metrics),
                  ylim = assign_value(ylim, "ylim", n_metrics),
                  alpha = assign_value(alpha, "alpha", n_metrics),
                  type = assign_value(type, "type", n_metrics),
                  major_lines_x = assign_value(major_lines_x, "major_lines_x",
                                               n_metrics),
                  major_lines_y = assign_value(major_lines_y, "major_lines_y",
                                               n_metrics),
                  legend_position = assign_value(legend_position,
                                                 "legend_position", n_metrics))

    return(track)

}





#' @title Configure track data from default values
#'
#' @description Wrapper around several track value assignment functions to
#' assign default values and generate track data.
#'
#' @param track A track object generated with the \code{\link{track}} function)
#'
#' @param defaults A named list with default values for all properties.
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_input}} function).
#'
#' @return A track object with values for all properties and track data stored
#' in track$data.
#'
#' @examples
#' track <- track("SNP_males", label = "Male-specific SNPs")
#' defaults <- list(colors = "red", alpha = 0.5)
#' data <- load_genome_metrics("window.tsv")
#' track <- configure_track(track, defaults, data)

configure_track <- function(track, defaults, data) {

    # Assign default values to track properties when not set by user
    track <- assign_defaults(track, defaults)

    # Create input data for track
    track <- create_track_data(data, track)

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
    for (i in 1:length(track)) {

        property <- properties[i]
        current_value <- c(track[[property]])[1]

        # Assign a value to properties not defined in track creation (NA value)
        if (is.na(current_value) & property %in% default_properties) {

            default_value <- defaults[[property]]
            track[[property]] <- assign_value(default_value, property,
                                              n_metrics)

        }

    }

    return(track)

}





#' @title Create track data
#'
#' @description Create an input data frame for the track plotting functions
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_input}} function).
#'
#' @param track Track object for the current plot, generated with the
#' \code{\link{track}} function.
#'
#' @return A data frame with columns:
#' Contig | Position | Metric 1 | Metric 1 colors |  Metric N | Metric N colors
#'
#' @examples
#' genomic_data <- load_genome_metrics("psass_window.tsv")
#' track <- track("Fst")
#' track <- create_track_data(genomic_data, track)

create_track_data <- function(data, track) {

    # Extract required columns and create color columns
    track$data <- data$data[, c("Contig_plot", "Position_plot", track$metrics)]

    names(track$data) <- c("Contig", "Position", track$metrics)

    # Combine data for multiple metrics
    track$data <- reshape2::melt(track$data, measure.vars = track$metrics,
                                 variable.name="Metric")

    # Assign color to data points
    track$data$Color <- stats::setNames(track$color,
                                        track$metrics)[track$data$Metric]

    # Sort data by Contig then Position
    track$data <- track$data[order(track$data$Contig, track$data$Position), ]

    return(track)

}






#' @title Assign value to track property
#'
#' @description Assign value to a single track property, handling single vs
#' multi metrics value definitions.
#'
#' @param value Value to assign, either a single value or a vector of values.
#'
#' @param proprety Name of the property to assign the value to, a string.
#'
#' @param n_metrics Number of metrics in the track, an integer.
#'
#' @return The correct value to assign to the property.
#'
#' @examples
#' n_metrics <- 2
#' colors <- assign_value("red", "colors", n_metrics)
#' bg_colors <- assign_value(c("white", "grey80"), "bg_colors", n_metrics)

assign_value <- function(value, property, n_metrics) {

    # Define globality for each track property: if TRUE, the property is global,
    # meaning only one value is requried (e.g. ylim); if FALSE, the property
    # applies to each metric and requires multiple values (e.g. colors)
    global <- c("colors" = FALSE, "bg_colors" = TRUE, "point_size" = FALSE,
                "ylim" = TRUE, "alpha" = FALSE, "type" = FALSE,
                "major_lines_x" = TRUE, "major_lines_y" = TRUE,
                "legend_position" = TRUE)

    # Assign multiple values when required
    if (n_metrics > 1 & global[property] == FALSE & length(value) < n_metrics) {

        value <- rep(value, n_metrics)

    } else if (!(length(value) %in% c(1, n_metrics))) {

        stop(paste0("Incorrect value \"", value, "\" for property \"",
                    property, "\" in track definition."))

    }

    return(value)

}