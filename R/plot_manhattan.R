#' @export
#'
#' @title Plot genome metrics in a manhattan plot
#'
#' @description Generate a manhattan plot with one or multiple tracks from a
#' genome metrics file
#'
#' @param input_file Path to a genome metrics input file (e.g. result of
#' PSASS or RADSex).
#' Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
#' corresponding position on the corresponding contig.
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects for manhattan plots are generated with the
#' \code{\link{track}} function. Tracks can represent one or mutiple
#' metrics from the genome metrics input file.
#'
#' @param chromosomes_file Path to the chromosome names file
#' (i.e. tab-separated file without header and with columns
#' <Contig ID> | <Chromosome name>). If NA, all contigs will be considered
#' unplaced except if detect_chromosomes is set to TRUE, in which case
#' chromosomes will be detected automatically from contig identifiers
#' (default: NA).
#'
#' @param detect_chromosomes If TRUE, will consider contigs starting with
#' "LG", "CHR", or "NC" as chromosomes if no chromosomes were specified
#' (default: TRUE).
#'
#' @param unplaced_label Label for unplaced contigs superscaffold
#' (default: "Unplaced").
#'
#' @param comment_char Character indicating a comment line in the input file
#' (default: "#").
#'
#' @param comment_sep Character separating two fields in a comment line
#' (default: ";").
#'
#' @param comment_internal_sep Character separating property and value in a
#' field from a comment line (default: ":").#'
#'
#' @param output_file Path to an output file for the generated manhattan plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param track_height Height of a single track when plotting to an output file,
#' in inches (default: 6).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param x_title Title of the x-axis (default: "Chromosomes").
#'
#' @param show_chromosomes_names If TRUE, display chromosome names on the x axis
#' (default: TRUE)
#'
#' @param chromosomes_as_numbers If TRUE, display chromosome numbers instead of
#' names for readability (default: FALSE)
#'
#' @param default_colors Default colors when not specified in track data
#' (default: "dodgerblue3", "darkgoldenrod2").
#'
#' @param default_bg_colors Default background colors for standard sectors
#' when not specified in track data (default: c("white")).
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 0.5).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 1).
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @examples
#' plot_manhattan("data/psass_window.tsv",
#'                tracks = list(track("Fst",
#'                                    label = expression("F"["ST"])),
#'                              track("Snps_females",
#'                                    label = "Pool-specific SNPs",
#'                                    colors = "firebrick2"),
#'                              track("Snps_males",
#'                                    label = "Pool-specific SNPs",
#'                                    colors = "dodgerblue3"))
#'                default_bg_colors = c("white", "grey60"),
#'                output_file = "manhattan.png")

plot_manhattan <- function(input_file, tracks,
                           chromosomes_file = NA,
                           detect_chromosomes = TRUE,
                           unplaced_label = "U.",
                           comment_char = "#",
                           comment_sep = ";",
                           comment_internal_sep = ":",
                           output_file = NA,
                           width = 12,
                           track_height = 6,
                           res = 300,
                           x_title = "Chromosomes",
                           show_chromosomes_names = TRUE,
                           chromosomes_as_numbers = FALSE,
                           default_colors = c("dodgerblue3", "darkgoldenrod2"),
                           default_bg_colors = c("grey85", "white"),
                           default_point_size = 0.5,
                           default_ylim = NA,
                           default_alpha = 1,
                           default_major_lines_y = TRUE) {

    # Load chromosome names (return NA if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes,
                                unplaced_label = unplaced_label,
                                comment_char = comment_char,
                                comment_sep = comment_sep,
                                comment_internal_sep = comment_internal_sep)
    # Draw the plot
    m <- draw_manhattan_plot(data$data,
                             data$lengths,
                             tracks,
                             output_file = output_file,
                             width = width,
                             track_height = track_height,
                             res = res,
                             x_title = x_title,
                             show_chromosomes_names = show_chromosomes_names,
                             chromosomes_as_numbers = chromosomes_as_numbers,
                             default_colors = default_colors,
                             default_bg_colors = default_bg_colors,
                             default_point_size = default_point_size,
                             default_ylim = default_ylim,
                             default_alpha = default_alpha,
                             default_major_lines_y = default_major_lines_y)

    return(m)
}






#' @export
#'
#' @title Draw genome metrics in a manhattan plot
#'
#' @description Generate a manhattan plot with one or multiple tracks from a
#' genome metrics data frame and a contig lengths vector
#'
#' @param data Genome metrics data frame (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_input}} function).
#' Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
#' corresponding position on the corresponding contig.
#'
#' @param contig_lengths Contig lengths data frame from the output of the
#' \code{\link{load_genome_input}} function. Format: <Contig ID> | <Length>
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects for manhattan plots are generated with the
#' \code{\link{track}} function. Tracks can represent one or mutiple
#' metrics from the genome metrics input file.
#'
#' @param output_file Path to an output file for the generated manhattan plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param track_height Height of a single track when plotting to an output file,
#' in inches (default: 6).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param x_title Title of the x-axis (default: "Chromosomes").
#'
#' @param show_chromosomes_names If TRUE, display chromosome names on the x axis
#' (default: TRUE)
#'
#' @param chromosomes_as_numbers If TRUE, display chromosome numbers instead of
#' names for readability (default: FALSE)
#'
#' @param default_colors Default colors when not specified in track data
#' (default: "dodgerblue3", "darkgoldenrod2").
#'
#' @param default_bg_colors Default background colors for standard sectors
#' when not specified in track data (default: c("white")).
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 0.5).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 1).
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @examples
#' metrics <- load_genome_metrics("psass_window.tsv")
#' draw_manhattan(metrics$data,
#'                metrics$lengths,
#'                tracks = list(track("Fst",
#'                                    label = expression("F"["ST"])),
#'                              track("Snps_females",
#'                                    label = "Pool-specific SNPs",
#'                                    colors = "firebrick2"),
#'                              track("Snps_males",
#'                                    label = "Pool-specific SNPs",
#'                                    colors = "dodgerblue3"))
#'                output_file = "manhattan.png")

draw_manhattan_plot <- function(data,
                                contig_lengths,
                                tracks,
                                output_file = NA,
                                width = 12,
                                track_height = 6,
                                res = 300,
                                x_title = "Chromosomes",
                                show_chromosomes_names = TRUE,
                                chromosomes_as_numbers = FALSE,
                                default_colors = c("dodgerblue3", "darkgoldenrod2"),
                                default_bg_colors = c("grey85", "white"),
                                default_point_size = 0.5,
                                default_ylim = NA,
                                default_alpha = 1,
                                default_major_lines_y = TRUE) {

    defaults <- list(colors = default_colors,
                     bg_colors = default_bg_colors,
                     point_size = default_point_size,
                     ylim = default_ylim,
                     alpha = default_alpha,
                     major_lines_y = default_major_lines_y,
                     h_lines = NA)

    # Compute increment to add to each contig (i.e. cumulative length of each
    # contig before this one)
    increments <- setNames(c(0, cumsum(head(contig_lengths, -1))),
                           names(contig_lengths))

    # Adjust x-axis position for each point in the data based on increments
    data$Position_plot = data$Position_plot + increments[data$Contig_plot]

    # Create data for background rectangles (for alternating background color)
    backgrounds <- data.frame(contig = names(increments),
                              start = increments,
                              end = increments + contig_lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE  # For x-axis labels and title

        # Configure track values
        tracks[[i]] <- configure_track(tracks[[i]], defaults, data)

        # Add backgrounds to track
        tracks[[i]] <- manhattan_bg_colors(tracks[[i]], backgrounds)

        # Generate track plot
        plots[[i]] <- draw_manhattan_track(tracks[[i]],
                                           bottom_track = bottom_track,
                                           x_title = x_title,
                                           show_chromosomes_names = show_chromosomes_names,
                                           chromosomes_as_numbers = chromosomes_as_numbers)

    }

    # Combine all tracks in a single plot
    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = combined,
                        width = width,
                        height = track_height * n_tracks,
                        dpi = res)

    } else {

        print(combined)

    }

    # Return plot list
    return(plots)

}





#' @export
#'
#' @title Plot a manhattan track
#'
#' @description Plot a single track for a manhattan plot
#'
#' @param track Track object storing properties for the current track,
#' generated with the \code{\link{track}} function.
#'
#' @param bottom_track If TRUE, x-axis labels and title will be added to
#' the plot (default: FALSE).
#'
#' @param x_title Title of the x-axis (default: "Chromosomes")
#'
#' @param show_chromosomes_names If TRUE, display chromosome names on the x axis
#' (default: TRUE)
#'
#' @param chromosomes_as_numbers If TRUE, display chromosome numbers instead of
#' names for readability (default: FALSE)
#'
#' @return A ggplot object for the track
#'
#' @examples
#' data <- load_genome_metrics("window.tsv")
#' defaults <- list(colors = c("yellow", "blue"))
#' track <- track("SNP_males", label = "Male-specific SNPs")
#' track <- configure_track(track, defaults, data)
#' track_plot <- draw_manhattan_track(track, bottom_track=TRUE)

draw_manhattan_track <- function(track,
                                 bottom_track = FALSE,
                                 x_title = "Chromosomes",
                                 show_chromosomes_names = TRUE,
                                 chromosomes_as_numbers = FALSE) {

    data <- track$data
    backgrounds <- track$backgrounds
    metric <- track$metrics[[1]]

    # Assign values for y-axis limits
    if (is.na(track$ylim)) {

        ymin <- min(data$Value)
        if (ymin < 0) { ymin <- 1.025 * ymin } else { ymin <- 0.976 * ymin }
        ymax <- max(data$Value)
        if (ymax < 0) { ymax <- 0.976 * ymax } else { ymax <- 1.025 * ymax }
        track$ylim <- c(ymin, ymax)

    }

    # Assign labels to chromosomes
    if (chromosomes_as_numbers) {

        chromosomes_names <- c(seq(1, nrow(backgrounds) - 1), "U")
        x_labels_angle <- 0
        x_labels_offset <- 0.5

    } else {

        chromosomes_names <- backgrounds$contig
        x_labels_angle <- 90
        x_labels_offset <- 1

    }

    # Create a fake color palette for both backgrounds and points
    merged_color_palette <- setNames(c(unique(data$Color),
                                       unique(backgrounds$color)),
                                     c(unique(data$Color),
                                       unique(backgrounds$color)))

    manhattan_plot <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        # Backgrounds with alternating colors
        ggplot2::geom_rect(data = backgrounds,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = track$ylim[1],
                                        ymax = track$ylim[2],
                                        fill = color),
                           alpha = 0.5) +
        # Data points
        ggplot2::geom_point(data = data,
                            ggplot2::aes(x = Position,
                                         y = Value,
                                         color = Color),
                            size = metric$point_size,
                            alpha = metric$alpha) +
        # Set color values from merged color scale for points and backgrounds
        ggplot2::scale_color_manual(values = merged_color_palette) +
        ggplot2::scale_fill_manual(values = merged_color_palette) +
        # Generate y-axis
        ggplot2::scale_y_continuous(name = track$label,
                                    limits = c(track$ylim[1], track$ylim[2]),
                                    expand = c(0, 0)) +
        # Adjust theme elements
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       axis.title.x = ggplot2::element_text(face = "bold",
                                                            margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(face = "bold",
                                                            margin = ggplot2::margin(0, 10, 0, 0)),
                       axis.text.y = ggplot2::element_text(color = "black",
                                                           face = "bold"),
                       axis.text.x = ggplot2::element_text(color = "black",
                                                           face = "bold",
                                                           vjust = 0.5,
                                                           angle = x_labels_angle,
                                                           hjust = x_labels_offset))

    if (show_chromosomes_names) {

        # Generate x-axis, use background start and end to place chromosome labels
        breaks <- backgrounds$start + (backgrounds$end - backgrounds$start) / 2
        manhattan_plot <- manhattan_plot +
            ggplot2::scale_x_continuous(name = x_title,
                                        breaks = breaks,
                                        labels = chromosomes_names,
                                        expand = c(0, 0))

    } else {

        # Generate x-axis
        manhattan_plot <- manhattan_plot +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank())

    }

    if (!(bottom_track)) {

        manhattan_plot <- manhattan_plot +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank())

    }

    if (is.na(x_title)) {

        manhattan_plot <- manhattan_plot +
            ggplot2::theme(axis.title.x = ggplot2::element_blank())

    }

    if (!(track$major_lines_y)) {

        manhattan_plot <- manhattan_plot +
            ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

    } else {

        major_lines <-  ggplot2::element_line(linetype=3, color = "grey30")
        manhattan_plot <- manhattan_plot +
            ggplot2::theme(panel.grid.major.y = major_lines)

    }

    # Add horizontal lines if defined
    if (!is.na(c(track$h_lines)[1])) {

        for (i in 1:length(track$h_lines)) {

            l <- track$h_lines[[i]]

            manhattan_plot <- manhattan_plot + ggplot2::geom_hline(yintercept = l$y,
                                                                   color = l$color,
                                                                   linetype = l$type,
                                                                   size = l$size)

            xrange <- ggplot2::ggplot_build(manhattan_plot)$layout$panel_scales_x[[1]]$range$range
            label_x <- l$label_x * (xrange[2] - xrange[1])
            label_y <- (track$ylim[2] - track$ylim[1]) / 20 + l$y
            manhattan_plot <- manhattan_plot + ggplot2::annotate("text",
                                                                 x = label_x,
                                                                 y = label_y,
                                                                 label = l$label,
                                                                 color = l$color,
                                                                 size = l$label_font_size)

        }

    }

    return(manhattan_plot)

}





#' @title Assign color to manhattan background data
#'
#' @description Assign color to each background rectangle from a backgrounds
#' data frame for the \code{\link{plot_track_manhattan}} function.
#'
#' @param track Track object for the current plot, generated with the
#' \code{\link{track}} function.
#'
#' @param data Data frame with contig, start, and end for each
#' background rectangle.
#'
#' @return A track object containing a backgrounds data frame with columns:
#' contig | start | end | color.
#'
#' @examples
#' track <- single_metric_track("Fst")
#' background_data <- data.frame(contig = c("Chr1", "Chr2"),
#'                               start = c(0, 1000)
#'                               end = c(999, 2635))
#' background_data <- manhattan_bg_colors(track, background_data)

manhattan_bg_colors <- function(track, data) {

    n_contigs <- nrow(data)

    # Create palette of alternating background colors
    bg_palette <- setNames(rep(track$bg_colors, n_contigs)[1:n_contigs],
                           data$contig)

    # Assign color to each background
    data$color <- bg_palette[data$contig]

    # Store the backgrounds data frame in the track object
    track$backgrounds <- data

    return(track)
}
