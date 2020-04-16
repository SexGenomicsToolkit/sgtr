#' @export
#'
#' @title Plot genome metrics for a specificied genomic region
#'
#' @description Generate a linear plot for a specified genomic region with one
#' or multiple tracks from a genome metrics file
#'
#' @param input_file Path to a genome metrics input file (e.g. result of
#' PSASS or RADSex).
#' Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
#' corresponding position on the corresponding contig.
#'
#' @param region Region to plot, defined with the syntax "Contig" or
#' "Contig:start-end"
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects for region plots are generated with the
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
#' @param output_file Path to an output file for the generated region plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param track_height Height of a single track when plotting to an output file,
#' in inches (default: 4).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param default_colors Default colors when not specified in track data
#' (default: "grey20").
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 1).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 0.6).
#'
#' @param default_type Default plot type when not specified in track data,
#' either "ribbon" or "points" (default: "ribbon").
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @param default_major_lines_x Default value for drawing major lines on the
#' x axis. If TRUE, reference lines will be plotted (default: FALSE).
#'
#' @param default_legend_position Default value for position of the legend,
#' directly passed to "legend.position" in \code{\link{ggplot2::theme}}
#' (default: "right").
#'
#' @return
#'
#' @examples
#' region_plot <- plot_region("psass_window.tsv",
#'                            "Chr24:0-6000000",
#'                            tracks = list(single_metric_track("Fst",
#'                                                              label = expression("F"["ST"])),
#'                                          multi_metrics_track(c("Snp_females", "Snp_males"),
#'                                                              label = "Pool-specific SNPs",
#'                                                              colors = c("firebrick2",
#'                                                                         "dodgerblue3"))),
#'                            output_file = "region_plot.png")

plot_region <- function(input_file,
                        region,
                        tracks,
                        chromosomes_file = NA,
                        detect_chromosomes = TRUE,
                        output_file = NA,
                        width = 12,
                        track_height = 4,
                        res = 300,
                        default_colors = "grey20",
                        default_point_size = 1,
                        default_ylim = NA,
                        default_alpha = 0.6,
                        default_type = "ribbon",
                        default_major_lines_y = TRUE,
                        default_major_lines_x = FALSE,
                        default_legend_position = "right") {


    # Load chromosome names (return NA if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes)
    # Draw the plot
    r <- draw_region(data$data,
                     data$lengths,
                     region,
                     tracks,
                     output_file = output_file,
                     width = width,
                     track_height = track_height,
                     res = res,
                     default_colors = default_colors,
                     default_point_size = default_point_size,
                     default_ylim = default_ylim,
                     default_alpha = default_alpha,
                     default_major_lines_y = default_major_lines_y,
                     default_major_lines_x = default_major_lines_x,
                     default_legend_position = default_legend_position)

    return(r)

}





#' @export
#'
#' @title Draw genome metrics for a specified genomic region
#'
#' @description Generate a linear plot for a specified genomic region with one
#' or multiple tracks from a genome metrics data frame and a contig lengths
#' vector
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
#' @param region Region to plot, defined with the syntax "Contig" or
#' "Contig:start-end"
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects for region plots are generated with the
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
#' (default: "U.").
#'
#' @param output_file Path to an output file for the generated manhattan plot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param track_height Height of a single track when plotting to an output file,
#' in inches (default: 4).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param default_colors Default colors when not specified in track data
#' (default: "grey20").
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 1).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 0.6).
#'
#' @param default_type Default plot type when not specified in track data,
#' either "ribbon" or "points" (default: "ribbon").
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @param default_major_lines_x Default value for drawing major lines on the
#' x axis. If TRUE, reference lines will be plotted (default: FALSE).
#'
#' @param default_legend_position Default value for position of the legend,
#' directly passed to "legend.position" in \code{\link{ggplot2::theme}}
#' (default: "right").
#'
#' @return
#'
#' @examples
#' metrics <- load_genome_metrics("psass_window.tsv")
#' region_plot <- draw_region(metrics$data,
#'                            metrics$lengths,
#'                            "Chr24:0-6000000",
#'                            tracks = list(single_metric_track("Fst",
#'                                                              label = expression("F"["ST"])),
#'                                          multi_metrics_track(c("Snp_females", "Snp_males"),
#'                                                              label = "Pool-specific SNPs",
#'                                                              colors = c("firebrick2",
#'                                                                         "dodgerblue3"))),
#'                            output_file = "region_plot.png")

draw_region <- function(data,
                        contig_lengths,
                        region,
                        tracks,
                        chromosomes_file = NA,
                        detect_chromosomes = TRUE,
                        unplaced_label = "U.",
                        output_file = NA,
                        width = 12,
                        track_height = 4,
                        res = 300,
                        default_colors = "grey20",
                        default_point_size = 1,
                        default_ylim = NA,
                        default_alpha = 0.6,
                        default_type = "ribbon",
                        default_major_lines_y = TRUE,
                        default_major_lines_x = FALSE,
                        default_legend_position = "right") {

    defaults <- list(colors = default_colors,
                     point_size = default_point_size,
                     ylim = default_ylim,
                     alpha = default_alpha,
                     type = default_type,
                     major_lines_y = default_major_lines_y,
                     major_lines_x = default_major_lines_x,
                     legend_position = default_legend_position)

    # Add original contig names to contig lengths so that the user can
    # specify both chromosome names or contig names in region
    tmp <- unique(data[, c("Contig", "Length")])
    contig_lengths <- c(contig_lengths,
                        setNames(tmp$Length, tmp$Contig))

    # Get contig, min, and max from the region string
    region_info <- parse_region(region, contig_lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE  # For x-axis labels and title

        # Configure track values
        tracks[[i]] <- configure_track(tracks[[i]], defaults, data, region_info)

        # Generate track plot
        plots[[i]] <- draw_region_track(tracks[[i]], region_info,
                                        bottom_track = bottom_track)

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
#' @title Plot a region track
#'
#' @description Plot a single track for a genomic region
#'
#' @param track Track object storing properties for the current track,
#' generated with the \code{\link{track}} function.
#'
#' @param region Region to plot, output of the \code{\link{parse_region}}
#' function, or NA to retain the entire genome.
#'
#' @param bottom_track If TRUE, x-axis labels and title will be added to
#' the plot (default: FALSE).
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' data <- load_genome_metrics("psass_window.tsv")
#' defaults <- list(colors = "blue", alpha = 0.5)
#' track <- single_metric_track("SNP_males", label = "Male-specific SNPs")
#' track <- configure_track(track, defaults, data)
#' draw_region_track(track, "Chr01", bottom_track=TRUE)

draw_region_track <- function(track, region, bottom_track = FALSE) {

    data <- track$data
    n_metrics <- length(track$metrics)

    if (n_metrics > 1 |
        (n_metrics == 1 & !is.function(track$metrics[[1]]$colors))) {

        # Color palettes are not allowed in multi-metrics tracks, expect at most
        # one color per metric. Create a merged color scale for all metrics and
        # assign proper labels.

        colors <- c()
        metrics <- c()

        for (i in 1:length(track$metrics)) {

            metric <- track$metrics[[i]]
            colors <- c(colors, metric$colors)
            metrics <- c(metrics, metric$label)

        }

        data$Color <- factor(data$Color, levels = colors)

        color_scales <- list(ggplot2::scale_color_manual(name = "",
                                                         values = colors,
                                                         labels = metrics),
                             ggplot2::scale_fill_manual(name = "",
                                                        values = colors,
                                                        labels = metrics))

    } else {

        # Colors were defined using a palette function in a single metric track,
        # generate the proper color scale legend.

        colors <- stats::setNames(unique(data$Color), unique(data$Color))
        color_scales <- list(ggplot2::scale_color_manual(name = "",
                                                         values = colors),
                             ggplot2::scale_fill_manual(name = "",
                                                        values = colors))

    }

    # Create major grid lines for y axis if specified
    if (track$major_lines_y) {

        major_lines_y <- ggplot2::element_line(color = "grey95", linetype = 1)

    } else {

        major_lines_y <- ggplot2::element_blank()

    }

    # Create major grid lines for x axis if specified
    if (track$major_lines_x) {

        major_lines_x <- ggplot2::element_line(color = "grey95", linetype = 1)

    } else {

        major_lines_x <- ggplot2::element_blank()

    }

    # Add x axis if bottom track
    if (!bottom_track) {

        axis_title_x <- ggplot2::element_blank()

    } else {

        axis_title_x <- ggplot2::element_text(face = "bold",
                                              margin = ggplot2::margin(10, 5, 5, 5))

    }

    # Assign values for y-axis limits
    if (is.na(c(track$ylim)[1])) {

        ymin <- min(data$Value)
        if (ymin < 0) { ymin <- 1.025 * ymin } else { ymin <- 0.976 * ymin }
        ymax <- max(data$Value)
        if (ymax < 0) { ymax <- 0.976 * ymax } else { ymax <- 1.025 * ymax }
        track$ylim <- c(ymin, ymax)

    }

    # Initialize the plot
    g <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        ggplot2::scale_y_continuous(name = track$label,
                                    expand = ggplot2::expansion(c(0, 0.01), 0),
                                    limits = track$ylim) +
        generate_x_scale(region) +
        ggplot2::theme(legend.position = track$legend_position,
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x,
                       axis.title.x = axis_title_x)

    # Draw data for each metric
    for (i in 1:length(track$metrics)) {

        metric <- track$metrics[[i]]
        plot_data <- subset(data, data$Metric == metric$name)

        if (metric$type == "ribbon") {

            g <- g + ggplot2::geom_ribbon(data = plot_data,
                                          ggplot2::aes(x = Position,
                                                       ymin = track$ylim[1],
                                                       ymax = Value,
                                                       fill = Color,
                                                       color = Color),
                                          size = 0.4,
                                          alpha = metric$alpha)

        } else if (metric$type == "points") {

            g <- g + ggplot2::geom_point(data = plot_data,
                                         ggplot2::aes(x = Position,
                                                      y = Value,
                                                      fill = Color,
                                                      color = Color),
                                         size = metric$point_size,
                                         alpha = metric$alpha)

        }

    }

    g <- g + color_scales

    if (n_metrics == 1) { g <- g + ggplot2::guides(color = FALSE,
                                                   fill = FALSE) }

    return(g)

}
