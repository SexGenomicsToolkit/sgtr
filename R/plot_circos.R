#' @title Plot genome metrics in a circular layout
#'
#' @description Generate a circos plot with multiple tracks from a
#' genome metrics file
#'
#' @param input_file_path Path to a genome metrics input file (e.g. result of
#' PSASS or RADSex).
#' Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
#' corresponding position on the corresponding contig.
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects for circos plots are generated with the
#' \code{\link{circos_track}} function. Tracks can represent one or mutiple
#' metrics from the genome metrics input file.
#'
#' @param chromosomes_file_path Path to the chromosome names file
#' (i.e. tab-separated file without header and with columns
#' <Contig ID> | <Chromosome name>). If NULL, all contigs will be considered
#' unplaced except if detect_chromosomes is set to TRUE, in which case
#' chromosomes will be detected automatically from contig identifiers
#' (default: NULL)
#'
#' @param detect_chromosomes If TRUE, will consider contigs starting with
#' "LG", "CHR", or "NC" as chromosomes if no chromosomes were specified
#' (default: TRUE)
#'
#' @param unplaced_label Label for unplaced contigs superscaffold
#' (default: "Unplaced")
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot (default NULL)
#'
#' @param output_file Path to an output file for the generated circos plot,
#' or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in pixel
#' (default: 2400)
#'
#' @param height Plot height when plotting to an output file, in pixel
#' (default: 2400)
#'
#' @param res Image resolution when plotting to an output file, in ppi
#' (default: 120)
#'
#' @param default_color Default color for a track when not specified
#' in track data (default: "grey20")
#'
#' @param default_point_size Default point size for a track when not specified
#' in track data (default: 0.25)
#'
#' @param default_ylim Default y-axis limits for a track when not specified
#' in track data (default: NULL, i.e. infer from data)
#'
#' @param default_bg_color Default background color for standard sectors
#' in the track when not specified in track data (default: "white")
#'
#' @param default_bg_highlight_color Default background color for highlighted
#' sectors in the track when not specified in track data (default: "grey80")
#'
#' @param sector_titles_expand Manually set the space between sector titles
#' and x-axis as a multiple of ymax (default: NULL)
#'
#' @examples
#' plot_circos("data/psass_window.tsv",
#'             tracks = list(track("Fst",
#'                                 label = expression("F"["ST"])),
#'                           track(c("Snps_females", "Snps_males"),
#'                                 label = "Pool-specific SNPs",
#'                                 color = c("firebrick2", "dodgerblue3")),
#'                           track(c("Depth_ratio"),
#'                                 label = "Depth ratio",
#'                                 color = "grey50")),
#'             output_file = "circos.png")

plot_circos <- function(input_file_path, tracks,
                        chromosomes_file_path = NULL, detect_chromosomes = TRUE,
                        unplaced_label = "Unplaced",
                        highlight = NULL,
                        output_file = NULL, width = 2400, height = 2400,
                        res = 120,
                        default_color = "grey20", default_point_size = 0.25,
                        default_ylim = NULL, default_bg_color = "white",
                        default_bg_highlight_color = "grey80",
                        sector_titles_expand = NULL) {

    # Load chromosome names (return NULL if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file_path)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file_path,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes,
                                unplaced_label = unplaced_label)

    # Draw the plot
    draw_circos(data$data,
                data$lengths,
                tracks,
                highlight = highlight,
                output_file = output_file,
                width = width,
                height = height,
                res = res,
                default_color = default_color,
                default_point_size = default_point_size,
                default_ylim = default_ylim,
                default_bg_color = default_bg_color,
                default_bg_highlight_color = default_bg_highlight_color,
                sector_titles_expand = sector_titles_expand
    )

}





#' @title Draw genome metrics in a circular layout
#'
#' @description Generate a circos plot with multiple tracks from a
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
#' plot. Track data objects for circos plots are generated with the
#' \code{\link{circos_track}} function. Tracks can represent one or mutiple
#' metrics from the genome metrics input file.
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot (default NULL)
#'
#' @param output_file Path to an output file for the generated circos plot,
#' or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in pixel
#' (default: 2400)
#'
#' @param height Plot height when plotting to an output file, in pixel
#' (default: 2400)
#'
#' @param res Image resolution when plotting to an output file, in ppi
#' (default: 120)
#'
#' @param default_color Default color for a track when not specified
#' in track data (default: "grey20")
#'
#' @param default_point_size Default point size for a track when not specified
#' in track data (default: 0.25)
#'
#' @param default_ylim Default y-axis limits for a track when not specified
#' in track data (default: NULL, i.e. infer from data)
#'
#' @param default_bg_color Default background color for standard sectors
#' in the track when not specified in track data (default: "white")
#'
#' @param default_bg_highlight_color Default background color for highlighted
#' sectors in the track when not specified in track data (default: "grey80")
#'
#' @param sector_titles_expand Manually set the space between sector titles
#' and x-axis as a multiple of ymax (default: NULL)
#'
#' @examples
#' metrics <- load_genome_metrics("psass_window.tsv")
#' plot_circos(metrics$data,
#'             metrics$lengths,
#'             tracks = list(track("Fst",
#'                                 label = expression("F"["ST"])),
#'                           track(c("Snps_females", "Snps_males"),
#'                                 label = "Pool-specific SNPs",
#'                                 color = c("firebrick2", "dodgerblue3")),
#'                           track(c("Depth_ratio"),
#'                                 label = "Depth ratio",
#'                                 color = "grey50")),
#'             output_file = "circos.png")

draw_circos <- function(data, contig_lengths, tracks,
                        highlight = NULL,
                        output_file = NULL, width = 2400, height = 2400,
                        res = 120,
                        default_color = "grey20", default_point_size = 0.25,
                        default_ylim = NULL, default_bg_color = "white",
                        default_bg_highlight_color = "grey80",
                        sector_titles_expand = NULL) {

    # Check that sectors to highlight exist
    highlight <- check_highlight_sectors(data, contig_lengths, highlight)

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, length(contig_lengths)), contig_lengths),
                             ncol=2)

    # Setup gaps between sectors
    gaps <- rep(1, nrow(sector_lengths))
    # Bigger gap after the last sector to leave some space for track titles
    gaps[length(gaps)] <- 12

    # Set a small offset to starting angle to make room for axes labels
    starting_angle_offset <- 7.5

    # Set track width sector_title_expand factor based on number of tracks
    n_tracks <- length(tracks)
    # Width of 0.25 if single track is plotted
    track_height <- 0.25
    # For multiple tracks, total width 0.5 (half the circos' width)
    if (n_tracks > 1) { track_height <- 0.5 / n_tracks }
    if (is.null(sector_titles_expand)) {
        sector_titles_expand <- 1.1 + n_tracks / 10
    }

    # Reset and setup circos parameters
    circlize::circos.clear()
    circlize::circos.par("track.height" = track_height,
                         "start.degree" = 90 - starting_angle_offset,
                         "gap.degree" = gaps,
                         "cell.padding" = c(0.01, 0.1, 0.01, 0.1),
                         "points.overflow.warning" = FALSE,
                         "track.margin" = c(0, 0.015),
                         "unit.circle.segments" = 100)

    # Open output file if specified
    if (!is.null(output_file)) {
        png(output_file, width=width, height=height,  res=res)
    }

    # Initialize circos plot
    circlize::circos.initialize(factors = names(contig_lengths),
                                xlim = sector_lengths,
                                sector.width = contig_lengths)

    defaults <- set_track_defaults(track,
                                   color = default_color,
                                   point_size = default_point_size,
                                   ylim = default_ylim,
                                   bg_color = default_bg_color,
                                   bg_highlight_color = default_bg_highlight_color)

    # Draw specified tracks
    top_track <- TRUE
    for (i in c(1:length(tracks))) {

        configure_circos_track_data(tracks[[i]])

        # Plot a single track
        plot_circos_track(tracks[[i]],
                          top.track = top_track, sector_titles_expand = sector_titles_expand,
                          first.sector = names(contig_lengths)[1])

        top_track <- FALSE
    }

    # Close output file
    if (!is.null(output_file)) { dev.off() }

}





#' @title Check highlight sectors
#'
#' @description Check that sectors to be highlighted exist, allowing both contig names and chromosomes names
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig_lengths Contig lengths from the output of the \code{\link{load_genome_input}} function
#'
#' @param highlight Vector of names of sectors to highlight
#'
#' @return A clean vector of names for sectors to highlight
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#'
#' highlight <- check_highlight_sectors(genomic_data$data, genomic_data$lengths, c("Chr01", "NC_002364.1"))
#'

check_highlight_sectors <- function(data, contig_lengths, highlight) {

    if (is.null(highlight)) {return(c())}

    for (i in 1:length(highlight)) {

        if (!(highlight[i] %in% names(contig_lengths))) {  # Sector to highlight not found in list of sector lengths

            contigs <- stats::setNames(unique(data$Contig_plot), unique(data$Contig))

            if (highlight[i] %in% names(contigs)) {  # Look for specified sector to highlight in original contig names

                highlight[i] <- contigs[highlight[i]]

            } else {

                print(paste0("Could not find sector to highlight \"", highlight[i], "\"."))
                exit(1)

            }
        }
    }

    return(highlight)

}





#' @title Create circos track data
#'
#' @description Create input data frame for the \code{\link{plot_track_circos}} function from the genomic data and the track information
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param track Track object for the current plot, generated with the \code{\link{circos_track}} function
#'
#' @return A data frame with columns:
#' Contig_plot | Position_plot | Metric 1 | Metric 1 colors | ... | Metric N | Metric N colors
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' track <- circos_track("Fst")
#'
#' region_data <- create_circos_track_data(genomic_data, track)
#'

create_circos_track_data <- function(data, track) {

    # Extract required columns and create color columns
    track_data <- data[, c("Contig_plot", "Position_plot", track$metrics)]
    # Combine data for multiple metrics
    track_data <- reshape2::melt(track_data, measure.vars = track$metrics, variable.name="Color")
    # Assign color to data points
    track_data$Color <- stats::setNames(track$color, track$metrics)[track_data$Color]
    # Sort data by Contig then Position
    track_data <- track_data[order(track_data$Contig_plot, track_data$Position_plot), ]

    return(track_data)

}





#' @title Create circos track object
#'
#' @description Generate an object storing all properties for a circos track
#'
#' @param metrics Metrics included in the track. Metrics should correspond to column names in the data frame used as input
#' data in the plot (output of the \code{\link{load_genome_input}} function). Possible values: a string if the track includes
#' a single metric (e.g. "Fst"), or a vector if the track includes several metrics (e.g. c("Snps_females", "Snps_males"))
#'
#' @param label Track label, NULL to set the label to the metric name for single-metric tracks. Label has to be specified for multi-metrics tracks (default: NULL)
#'
#' @param color Track color. Values can be a string (e.g. "grey20") and will then be applied to all metrics, or a vector of size equal to the number of metrics
#' (e.g. c("red", "blue") for two metrics)
#'
#' @param point.size Point size for plots of type "points". Values can be a float (e.g. 0.5) and will then be applied to all metrics,
#' or a vector of size equal to the number of metrics (e.g. c(1, 1.5, 3) for three metrics)
#'
#' @param ylim Vector of y-axis limits for the track; if NULL, infer directly from data
#'
#' @param bg.color Background color for standard sectors in the track (single value or vector of length equal to the number of sectors in the plot)
#'
#' @param bg.highlight.color Background color for highlighted sectors in the track (single value)
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Single metric
#' track_data <- circos_track("Fst", color = "grey70", point.size = 0.75)
#'
#' # Multiple metrics
#' track_data <- circos_track(c("Snp_females", "Snp_males"), color = c("firebrick2", "dodgerblue3"))
#'

circos_track <- function(metrics, label = NULL, color = NULL, point.size = NULL, ylim = NULL,
                         bg.color = NULL, bg.highlight.color = NULL) {

    n_metrics <- length(metrics)

    if (is.null(label)) {

        if (n_metrics == 1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    # Handle multi-values properties for multiple metrics (i.e. metric-specific properties)
    if (n_metrics > 1) {

        # For each option, assign value to each metric if single value was defined in multi-metrics track
        if (length(color) == 1) { color <- rep(color, n_metrics) }
        if (length(point.size) == 1) { point.size <- rep(point.size, n_metrics) }

    }

    # Single-value properties (not metric-specific)
    ylim <- ylim
    bg.color <- bg.color
    bg.highlight.color <- bg.highlight.color

    track_info <- list(metrics=metrics, label=label, color=color, point.size=point.size, ylim = ylim,
                       bg.color = bg.color, bg.highlight.color = bg.highlight.color)

    return(track_info)

}





#' @title Assign default values to a circos track object
#'
#' @description Assign default values to all properties for which the value was not specified by the user (e.g. value is NULL)
#'
#' @param track A track object generated with the \code{\link{circos_track}} function)
#'
#' @param default_color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default_point_size Default point size for a track of type "points" when not specified in track data (default: 0.5)
#'
#' @param default_ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @param default_bg_color Default background color for standard sectors in the track (default: "white")
#'
#' @param default_bg_highlight_color Default background color for highlighted sectors in the track (default: "grey80")
#'
#' @return A track object with default values for properties not specified by the user
#'
#' @examples
#' track_data <- assign_circos_track_default(track_data, default.alpha = 0.75)
#'

assign_circos_track_default <- function(track, default_color = "grey20", default_point_size = 0.5, default_ylim = NULL,
                                        default_bg_color = "white", default_bg_highlight_color = "grey80") {

    n_metrics <- length(track$metrics)

    # Metrics-specific options (create vector)
    if (is.null(track$color)) { track$color <- rep(default_color, n_metrics) }
    if (is.null(track$point.size)) { track$point.size <- rep(default_point_size, n_metrics) }

    # Track-specific options (single value)
    if (is.null(track$ylim)) { track$ylim <- default_ylim }
    if (is.null(track$bg.color)) { track$bg.color <- default_bg_color }
    if (is.null(track$bg.highlight.color)) { track$bg.highlight.color <- default_bg_highlight_color }

    return(track)

}





configure_circos_track_data <- function(track, defaults) {

    # Assign default values to track properties when not set by user
    track <- assign_circos_track_default(defaults)

    # Setup sector background colors for the track
    if (length(tracks[[i]]$bg.color) == 1) {

        # Single background color value for normal and highlighted sectors
        tracks[[i]]$bg_colors <- rep(tracks[[i]]$bg.color,
                                     length(contig_lengths))
        tracks[[i]]$bg_colors[which(names(contig_lengths) %in% highlight)] <- tracks[[i]]$bg.highlight.color

    } else {  # Background colors defined as vector of values for "bg.color"

        tracks[[i]]$bg_colors <- tracks[[i]]$bg.color

    }

    # Create input data for track
    create_circos_track_data(data, tracks[[i]])

}




#' @title Plot a circos track
#'
#' @description Plot a single track for a circos plot
#'
#' @param data Input data generated with the \code{\link{create_region_track_data}} function)
#'
#' @param track Track object storing properties for the current track, generated with the \code{\link{circos_track}} function
#'
#' @param top.track If TRUE, x-axis labels and sector names will be added to the track
#'
#' @param sector_titles_expand Manually set the space between sector titles and x-axis as a multiple of ymax (default: 1.3)
#'
#' @param first.sector Name of the first sector to draw y-axis on (default: NULL)
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' fst_track <- circos_track("Fst")
#' circos_data <- create_circos_track_data(genomic_data, fst_track)
#'
#' plot_track_region(region_data, fst_track, top.track=TRUE)
#'

plot_circos_track <- function(data, track, top.track = FALSE, sector_titles_expand = 1.3, first.sector=NULL) {

    # Assign values for y-axis parameters
    if (is.null(track$ylim)) { track$ylim = c(0.975 * min(data[, 4]) - 0.01, 1.025 * max(data[, 4]) + 0.01) }
    ylim <- track$ylim
    ylabel <- track$label

    # Draw the top track of the plot
    circlize::circos.track(factors = data$Contig_plot,
                           x = data$Position_plot,
                           y = as.vector(unlist(data[,4])),
                           ylim = track$ylim,
                           bg.col = track$bg_colors,
                           panel.fun = function(x, y) {  # panel.fun is the function drawing the track

                               # Get useful sector information
                               sector.index <- circlize::get.cell.meta.data("sector.index")
                               xcenter <- circlize::get.cell.meta.data("xcenter")
                               ymin <- circlize::get.cell.meta.data("ylim")[1]
                               ymax <- circlize::get.cell.meta.data("ylim")[2]
                               xmin <- circlize::get.cell.meta.data("xlim")[1]
                               xmax <- circlize::get.cell.meta.data("xlim")[2]
                               xplot <- circlize::get.cell.meta.data("xplot")

                               # Add top axis and titles to sectors
                               if (top.track) {

                                   # Create x axis on top of sectors
                                   circlize::circos.axis(h = "top",
                                                         major.at = c(0, xmax / 3, 2 * xmax / 3, xmax),  # Label every 1/3 of the axis
                                                         labels.cex = 1.2,
                                                         labels.facing = "outside",
                                                         direction="outside",
                                                         labels = convert_to_mb(c(0, xmax / 3, 2 * xmax / 3, xmax)),  # Conversion to Mb
                                                         minor.ticks = 4,
                                                         labels.pos.adjust = TRUE)

                                   # Add sector names
                                   circlize::circos.text(xcenter,
                                                         sector_titles_expand * ymax,
                                                         sector.index,
                                                         cex = 1.5,
                                                         facing = "bending.inside",
                                                         niceFacing = TRUE)
                               }

                               # Plot the data
                               circlize::circos.points(x, y, cex = track$point.size,
                                                       col = data$Color,
                                                       bg = data$Color,
                                                       pch = 21)

                               # Add Y axis on the first sector only
                               if (sector.index == first.sector) {

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          at = c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]),  # 3 labels
                                                          sector.index = first.sector,
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = FALSE,
                                                          labels = round(c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]), 0))

                                   #Add y axis labels
                                   label_offset <- - 7.5 * (xmax - xmin) / (xplot[1] - xplot[2])  # Axis title will be plotted 7.5° on the left of the axis
                                   circlize::circos.text(label_offset,
                                                         0.5 * (ymax - ymin) + ymin,
                                                         ylabel,
                                                         sector.index = first.sector,
                                                         facing = "clockwise",
                                                         cex = 1.3,
                                                         font = 2)
                               }
                           }
    )

}
