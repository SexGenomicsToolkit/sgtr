#' @title Plot Manhattan
#'
#' @description Generate a Manhattan plot with multiple track for the entire genome
#'
#' @param input.file.path Path to a genomic data input file (e.g. result of PSASS or RADSex)
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{manhattan_track}} function
#'
#' @param chromosomes.file.path Path to a tabulated file specifying chromosome names (default: NULL)
#'
#' @param detect.chromosomes If TRUE, will consider contigs starting with LG, CH, or NC as chromosomes if no chromosomes were specified (default: TRUE)
#'
#' @param unplaced.label Label for unplaced contigs (default: "U.")
#'
#' @param output.file Path to an output file for the generated Manhattan plot, or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in inches (default: 12)
#'
#' @param track.height Height of a single track when plotting to an output file, in inches (default: 6)
#'
#' @param res Image resolution when plotting to an output file, in dpi (default: 300)
#'
#' @param x.title Title of the x-axis
#'
#' @param show.chromosomes.names If TRUE, display chromosome names on the x axis
#'
#' @param chromosomes.as.numbers If TRUE, display chromosome numbers instead of names for readability
#'
#' @param default.point.color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default.bg.color Default background color for a track when not specified in the track (default: "white")
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.25)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @examples
#' a

plot_manhattan <- function(input.file.path, tracks,
                           chromosomes.file.path = NULL, detect.chromosomes = TRUE, unplaced.label = "U.",
                           output.file = NULL, width = 14, track.height = 6, res = 300,
                           x.title = "Chromosomes",
                           show.chromosomes.names = TRUE, chromosomes.as.numbers = FALSE,
                           default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                           default.bg.color = c("grey85", "white"),
                           default.point.size = 0.5, default.ylim = NULL) {

    # Load chromosome names (return NULL when chromosomes file path is NULL)
    chromosomes <- load_chromosome_names(chromosomes.file.path)

    # Load genomic data
    data <- load_genome_input(input.file.path, chromosomes = chromosomes, detect.chromosomes = detect.chromosomes, unplaced.label = unplaced.label)

    # Draw the plot
    draw_manhattan_plot(data$data,
                        data$contig.lengths,
                        tracks,
                        output.file = output.file,
                        width = width,
                        track.height = track.height,
                        res = res,
                        x.title = x.title,
                        show.chromosomes.names = show.chromosomes.names,
                        chromosomes.as.numbers = chromosomes.as.numbers,
                        default.point.color = default.point.color,
                        default.bg.color = default.bg.color,
                        default.point.size = default.point.size,
                        default.ylim = default.ylim
    )

}






#' @title Draw manhattan plot
#'
#' @description Generate a Manhattan plot with multiple track for the entire genome
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig.lengths Contig lengths from the output of the \code{\link{load_genome_input}} function
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{manhattan_track}} function
#'
#' @param output.file Path to an output file for the generated Manhattan plot, or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in inches (default: 12)
#'
#' @param track.height Height of a single track when plotting to an output file, in inches (default: 6)
#'
#' @param res Image resolution when plotting to an output file, in dpi (default: 300)
#'
#' @param x.title Title of the x-axis
#'
#' @param show.chromosomes.names If TRUE, display chromosome names on the x axis
#'
#' @param chromosomes.as.numbers If TRUE, display chromosome numbers instead of names for readability
#'
#' @param default.point.color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default.bg.color Default background color for a track when not specified in the track (default: "white")
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.25)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @examples
#' a

draw_manhattan_plot <- function(data, contig.lengths, tracks,
                                output.file = NULL, width = 14, track.height = 6, res = 300,
                                x.title = "Chromosomes",
                                show.chromosomes.names = TRUE, chromosomes.as.numbers = FALSE,
                                default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                                default.bg.color = c("grey85", "white"),
                                default.point.size = 0.5, default.ylim = NULL) {

    # Compute increment to add to each contig (i.e. cumulative length of each contig before this one)
    increments <- setNames(c(0, cumsum(head(contig.lengths, -1))), names(contig.lengths))

    # Adjust x-axis position for each point in the data based on increments
    data$Position_plot = data$Position_plot + increments[data$Contig_plot]

    # Create data for background rectangles (for alternating background color)
    backgrounds <- data.frame(contig = names(increments), start = increments, end = increments + contig.lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE  # For x-axis labels and title

        # Assign default values to track properties that were not specified by the user
        tracks[[i]] <- assign_manhattan_track_default(tracks[[i]], default.point.color = default.point.color, default.bg.color = default.bg.color,
                                                      default.point.size = default.point.size, default.ylim = default.ylim)

        # Generate track data
        track_data <- create_manhattan_track_data(data, tracks[[i]])

        # Generate background data
        track_background_data <- assign_manhattan_background_colors(backgrounds, tracks[[i]])

        # Generate track plot
        plots[[i]] <- plot_track_manhattan(track_data, track_background_data, tracks[[i]], bottom.track = bottom_track)

    }

    # Combine all tracks in a single plot
    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    # Output to file if specified or print in current R device otherwise
    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = track.height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    return(combined)

}





#' @title Create manhattan track object
#'
#' @description Generate an object storing all properties for a manhattan track
#'
#' @param metric Metric represented by the track, which should correspond to a column name in the data frame used as input
#' data in the plot (output of the \code{\link{load_genome_input}} function)
#'
#' @param label Track label, NULL to set the label to the metric name (default: NULL)
#'
#' @param point.color Point color, either a string (e.g. "grey20") or a vector of alternating colors (e.g. c("red", "blue") for two colors)
#'
#' @param bg.color Background color, either a string (e.g. "white") or a vector of alternating colors (e.g. c("grey85", "white") for two colors)
#'
#' @param point.size Point size, directly passed to ggplot
#'
#' @param ylim Vector of y-axis limits for the track; if NULL, infer directly from data
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Single metric
#' track_data <- manhattan_track("Fst", point.color = c("blue", "yellow"), bg.color = c("white", "grey50"), point.size = 0.75)

manhattan_track <- function(metric, label = NULL, point.color = NULL, bg.color = NULL, point.size = NULL, ylim = NULL) {

    # Set a default label if not specified
    if (is.null(label)) { label <- metric }

    # Create track object
    track_info <- list(metric=metric, label=label, point.color=point.color, bg.color=bg.color, point.size=point.size, ylim=ylim)

    return(track_info)
}





#' @title Assign default values to a manhattan track object
#'
#' @description Assign default values to all properties for which the value was not specified by the user (i.e. value is NULL)
#'
#' @param track A track object generated with the \code{\link{manhattan_track}} function)
#'
#' @param default.point.color Default point color when not specified in track data (default: c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param default.bg.color Default background color when not specified in track data (default: c("white", "grey85"))
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.5)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @return A track object with default values for properties not specified by the user
#'
#' @examples
#' track_data <- assign_manhattan_track_default(track_data, default.point.size = 0.75)
#'

assign_manhattan_track_default <- function(track, default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                                           default.bg.color = c("white", "grey85"),
                                           default.point.size = 0.5, default.ylim = NULL) {

    if (is.null(track$point.color)) { track$point.color <- default.point.color }
    if (is.null(track$bg.color)) { track$bg.color <- default.bg.color }
    if (is.null(track$point.size)) { track$point.size <- default.point.size }
    if (is.null(track$ylim)) { track$ylim <- default.ylim }

    return(track)

}





#' @title Create manhattan track data
#'
#' @description Create input data frame for the \code{\link{plot_track_manhattan}} function from the genomic data and the track information
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param track Track object for the current plot, generated with the \code{\link{manhattan_track}} function
#'
#' @return A data frame with columns:
#' Contig_plot | Position_plot | Value | Color
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' track <- manhattan_track("Fst")
#'
#' track_data <- create_manhattan_track_data(genomic_data, track)
#'

create_manhattan_track_data <- function(data, track) {

    # Extract required columns and create color columns
    track_data <- data[, c("Contig_plot", "Position_plot", track$metric)]
    names(track_data)[3] <- "Value"
    # Sort data by Contig then Position
    track_data <- track_data[order(track_data$Contig_plot, track_data$Position_plot), ]
    # Assign point color
    contigs <- unique(track_data$Contig_plot)
    n_contigs <- length(contigs)
    points_palette <- setNames(rep(track$point.color, n_contigs)[1:n_contigs], contigs)
    track_data$Color <- points_palette[track_data$Contig_plot]

    return(track_data)

}





#' @title Assign color to manhattan background data
#'
#' @description Assign color to each background rectangle from a backgrounds data frame for the \code{\link{plot_track_manhattan}} function
#'
#' @param data Data frame with contig, start, and end for each background rectangle
#'
#' @param track Track object for the current plot, generated with the \code{\link{manhattan_track}} function
#'
#' @return A data frame with columns:
#' contig | start | end | color
#'
#' @examples
#' track <- manhattan_track("Fst")
#' background_data <- assign_manhattan_background_colors(background_data, track)
#'

assign_manhattan_background_colors <- function(data, track) {

    background_data <- data

    # Assign point and background color
    n_contigs <- nrow(data)
    bg_palette <- setNames(rep(track$bg.color, n_contigs)[1:n_contigs], data$contig)
    background_data$color <- bg_palette[background_data$contig]

    return(background_data)
}





#' @title Plot a manhattan track
#'
#' @description Plot a single track for a manhattan plot
#'
#' @param data Data to draw in the manhattan track, generated with the \code{\link{create_manhattan_track_data}} function)
#'
#' @param data Backgrounds data generated with the \code{\link{assign_manhattan_background_colors}} function)
#'
#' @param track Track object storing properties for the current track, generated with the \code{\link{manhattan_track}} function
#'
#' @param bottom.track If TRUE, x-axis labels and title will be added to the plot
#'
#' @param x.title Title of the x-axis
#'
#' @param show.chromosomes.names If TRUE, display chromosome names on the x axis
#'
#' @param chromosomes.as.numbers If TRUE, display chromosome numbers instead of names for readability
#'
#' @return A ggplot object for the manhattan plot
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' fst_track <- manhattan_track("Fst")
#' manhattan_data <- create_manhattan_track_data(genomic_data, fst_track)
#'
#' fst_plot <- plot_track_manhattan(manhattan_data, fst_track, bottom.track=TRUE)
#'

plot_track_manhattan <- function(data, backgrounds, track, bottom.track = FALSE,
                                 x.title = "Chromosomes",
                                 show.chromosomes.names = TRUE, chromosomes.as.numbers = FALSE) {

    # Maximum / minimum Y value
    if (is.null(track$ylim)) { track$ylim = c(min(data[, 3]), 1.025 * max(data[, 3]) + 0.01) }
    ymin <- track$ylim[1]
    ymax <- track$ylim[2]

    # Assign labels to chromosomes
    if (chromosomes.as.numbers) {

        chromosomes_names <- c(seq(1, nrow(backgrounds) - 1), "U")
        x.labels.angle = 0

    } else {

        chromosomes_names <- backgrounds$contig
        x.labels.angle = 90

    }

    # Create a fake color palette for both backgrounds and points
    merged_color_palette <- setNames(c(unique(data$Color), unique(backgrounds$color)), c(unique(data$Color), unique(backgrounds$color)))

    manhattan_plot <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        # Backgrounds with alternating colors
        ggplot2::geom_rect(data = backgrounds,
                           ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color),
                           alpha = 0.5) +
        # Data points
        ggplot2::geom_point(data = data,
                            ggplot2::aes(x = Position_plot, y = Value, color = Color),
                            size = track$point.size,
                            alpha = 1) +
        # Attribute color values from merged color scale for points and backgrounds
        ggplot2::scale_color_manual(values = merged_color_palette) +
        ggplot2::scale_fill_manual(values = merged_color_palette) +
        # Generate y-axis
        ggplot2::scale_y_continuous(name = track$label,
                                    limits = c(ymin, ymax),
                                    expand = c(0, 0)) +
        # Adjust theme elements
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(face = "bold", margin = ggplot2::margin(0, 10, 0, 0)),
                       axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                       axis.text.x = ggplot2::element_text(color = "black", face = "bold", vjust = 0.5, angle = x.labels.angle, hjust = 1))

    if (show.chromosomes.names) {

        # Generate x-axis, use background start and end to place chromosome labels
        manhattan_plot <- manhattan_plot + ggplot2::scale_x_continuous(name = x.title,
                                                                       breaks = backgrounds$start + (backgrounds$end - backgrounds$start) / 2,
                                                                       labels = chromosomes_names,
                                                                       expand = c(0, 0))

    } else {

        # Generate x-axis
        manhattan_plot <- manhattan_plot +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

    }

    if (!(bottom.track)) {

        manhattan_plot <- manhattan_plot +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())

    }

    return(manhattan_plot)

}
