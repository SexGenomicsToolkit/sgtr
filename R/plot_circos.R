#' @export
#'
#' @title Plot genome metrics in a circular layout
#'
#' @description Generate a circos plot with multiple tracks from a
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
#' plot. Track data objects are generated with the \code{\link{track}},
#' \code{\link{single_metric_track}}, or \code{\link{multi_metrics_track}}
#' functions.
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
#' field from a comment line (default: ":").
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot (default NA).
#'
#' @param highlight_bg_color Background color for highlighted sectors
#' (default: "grey80").
#'
#' @param output_file Path to an output file for the generated circos plot,
#' or NA to plot in the current R device (default: NA)
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
#' @param default_colors Default colors when not specified in track data
#' (default: "grey20").
#'
#' @param default_bg_colors Default background colors for standard sectors
#' when not specified in track data (default: "white").
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 0.25).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 1).
#'
#' @param default_type Default plot type when not specified in track data,
#' either "ribbon" or "points" (default: "points").
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @param default_major_lines_x Default value for drawing major lines on the
#' x axis. If TRUE, reference lines will be plotted (default: FALSE).
#'
#' @param sector_titles_expand Manually set the space between sector titles
#' and x-axis as a multiple of ymax (default: NA)
#'
#' @examples
#' plot_circos("data/psass_window.tsv",
#'             tracks = list(single_metric_track("Fst",
#'                                               label = expression("F"["ST"])),
#'                           multi_metrics_track(c("Snp_females", "Snp_males"),
#'                                               label = "Pool-specific SNPs",
#'                                               colors = c("firebrick2",
#'                                                          "dodgerblue3")),
#'             default_type = "ribbon",
#'             output_file = "circos.png")

plot_circos <- function(input_file, tracks,
                        chromosomes_file = NA,
                        detect_chromosomes = TRUE,
                        unplaced_label = "U.",
                        comment_char = "#",
                        comment_sep = ";",
                        comment_internal_sep = ":",
                        highlight = NA,
                        highlight_bg_color = "grey80",
                        output_file = NA,
                        width = 2400,
                        height = 2400,
                        res = 120,
                        default_colors = "grey20",
                        default_bg_colors = "white",
                        default_point_size = 0.25,
                        default_ylim = NA,
                        default_alpha = 1,
                        default_type = "points",
                        default_major_lines_y = TRUE,
                        default_major_lines_x = FALSE,
                        sector_titles_expand = NA) {

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
    draw_circos(data$data,
                data$lengths,
                tracks,
                highlight = highlight,
                highlight_bg_color = highlight_bg_color,
                output_file = output_file,
                width = width,
                height = height,
                res = res,
                default_colors = default_colors,
                default_bg_colors = default_bg_colors,
                default_point_size = default_point_size,
                default_ylim = default_ylim,
                default_alpha = default_alpha,
                default_type = default_type,
                default_major_lines_y = default_major_lines_y,
                default_major_lines_x = default_major_lines_x,
                sector_titles_expand = sector_titles_expand)

}





#' @export
#'
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
#' plot. Track data objects are generated with the \code{\link{track}},
#' \code{\link{single_metric_track}}, or \code{\link{multi_metrics_track}}
#' functions.
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot (default NA).
#'
#' @param highlight_bg_color Background color for highlighted sectors
#' (default: "grey80").
#'
#' @param output_file Path to an output file for the generated circos plot,
#' or NA to plot in the current R device (default: NA)
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
#' @param default_colors Default colors when not specified in track data
#' (default: "grey20").
#'
#' @param default_bg_colors Default background colors for standard sectors
#' when not specified in track data (default: "white").
#'
#' @param default_point_size Default point size when not specified in track data
#' (default: 0.25).
#'
#' @param default_ylim Default y-axis limits when not specified in track data
#' (default: NA, i.e. infer from data).
#'
#' @param default_alpha Default alpha value when not specified in track data
#' (default: 1).
#'
#' @param default_type Default plot type when not specified in track data,
#' either "ribbon" or "points" (default: "points").
#'
#' @param default_major_lines_y Default value for drawing major lines on the
#' y axis. If TRUE, reference lines will be plotted (default: TRUE).
#'
#' @param default_major_lines_x Default value for drawing major lines on the
#' x axis. If TRUE, reference lines will be plotted (default: FALSE).
#'
#' @param sector_titles_expand Manually set the space between sector titles
#' and x-axis as a multiple of ymax (default: NA)
#'
#' @examples
#' metrics <- load_genome_metrics("psass_window.tsv")
#' draw_circos(metrics$data,
#'             metrics$lengths,
#'             tracks = list(single_metric_track("Fst",
#'                                               label = expression("F"["ST"])),
#'                           multi_metrics_track(c("Snp_females", "Snp_males"),
#'                                               label = "Pool-specific SNPs",
#'                                               colors = c("firebrick2",
#'                                                          "dodgerblue3")),
#'             output_file = "circos.png")

draw_circos <- function(data,
                        contig_lengths,
                        tracks,
                        highlight = NA,
                        highlight_bg_color = "grey80",
                        output_file = NA,
                        width = 2400,
                        height = 2400,
                        res = 120,
                        default_colors = "grey20",
                        default_bg_colors = "white",
                        default_point_size = 0.25,
                        default_ylim = NA,
                        default_alpha = 1,
                        default_type = "points",
                        default_major_lines_y = TRUE,
                        default_major_lines_x = FALSE,
                        sector_titles_expand = NA) {

    defaults <- list(colors = default_colors,
                     bg_colors = default_bg_colors,
                     point_size = default_point_size,
                     ylim = default_ylim,
                     alpha = default_alpha,
                     type = default_type,
                     major_lines_x = default_major_lines_x,
                     major_lines_y = default_major_lines_y,
                     h_lines = NA)

    # Check that sectors to highlight exist
    highlight <- check_highlight_sectors(data, contig_lengths, highlight)

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, length(contig_lengths)), contig_lengths),
                             ncol = 2)

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
    if (is.na(sector_titles_expand)) {
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

    # Open output file in the correct device
    open_output_device(output_file, width, height, res)

    # Initialize circos plot
    circlize::circos.initialize(factors = names(contig_lengths),
                                xlim = sector_lengths,
                                sector.width = contig_lengths)

    # Draw specified tracks
    top_track <- TRUE
    for (i in c(1:length(tracks))) {

        tracks[[i]] <- configure_track(tracks[[i]], defaults, data)

        tracks[[i]] <- configure_backgrounds(tracks[[i]], highlight,
                                             highlight_bg_color)

        # Plot a single track
        draw_circos_track(tracks[[i]],
                          top_track = top_track,
                          sector_titles_expand = sector_titles_expand,
                          first_sector = names(contig_lengths)[1])

        top_track <- FALSE
    }

    # Close output file
    if (!is.na(output_file)) { dev.off() }

}





#' @export
#'
#' @title Plot a circos track
#'
#' @description Plot a single track for a circos plot
#'
#' @param track Track object storing properties for the current track,
#' generated with the \code{\link{track}}, \code{\link{single_metric_track}},
#' or \code{\link{multi_metrics_track}} functions.
#'
#' @param top_track If TRUE, x-axis labels and sector names will be added to
#' the track (default: FALSE).
#'
#' @param sector_titles_expand Manually set the space between sector titles and
#' x-axis as a multiple of ymax (default: 1.3).
#'
#' @param first_sector Name of the first sector to draw y-axis on
#' (default: NA).
#'
#' @examples
#' data <- load_genome_metrics("window.tsv")
#' defaults <- list(colors = "blue", alpha = 0.5)
#' track <- single_metric_track("SNP_males", label = "Male-specific SNPs")
#' track <- configure_track(track, defaults, data)
#' circlize::circos.par("track.height" = 12)
#' draw_circos_track(track, top_track=TRUE)

draw_circos_track <- function(track,
                              top_track = FALSE,
                              sector_titles_expand = 1.3,
                              first_sector=NA) {

    # Get data from track
    data <- track$data

    # Assign values for y-axis limits
    if (is.na(c(track$ylim)[1])) {

        ymin <- min(data[, 4])
        if (ymin < 0) { ymin <- 1.025 * ymin } else { ymin <- 0.976 * ymin }
        ymax <- max(data[, 4])
        if (ymax < 0) { ymax <- 0.976 * ymax } else { ymax <- 1.025 * ymax }
        track$ylim <- c(ymin, ymax)

    }

    ylim <- track$ylim
    ylabel <- track$label

    # Draw the top track of the plot
    circlize::circos.track(factors = data$Contig,
                           x = data$Position,
                           y = data$Value,
                           ylim = ylim,
                           bg.col = track$bg_colors,
                           panel.fun = function(x, y) {

                               # Get useful sector information
                               sector_index <- circlize::CELL_META$sector.index
                               xcenter <- circlize::CELL_META$xcenter
                               ymin <- circlize::CELL_META$ylim[1]
                               ymax <- circlize::CELL_META$ylim[2]
                               xmin <- circlize::CELL_META$xlim[1]
                               xmax <- circlize::CELL_META$xlim[2]
                               xplot <- circlize::CELL_META$xplot

                               x_ticks <- c(0, xmax / 3, 2 * xmax / 3, xmax)
                               y_ticks <- c(ylim[1],
                                            (ylim[2] - ylim[1]) / 2 + ylim[1],
                                            ylim[2])
                               major_y <- c(ylim[1],
                                            (ylim[2] - ylim[1]) / 4 + ylim[1],
                                            (ylim[2] - ylim[1]) / 2 + ylim[1],
                                            3 * (ylim[2] - ylim[1]) / 4 + ylim[1],
                                            ylim[2])

                               # Add top axis and titles to sectors
                               if (top_track) {

                                   x_ticks_mb <- convert_to_mb(x_ticks)

                                   # Create x axis on top of sectors
                                   circlize::circos.axis(h = "top",
                                                         # Label every 1/3 of the axis
                                                         major.at = x_ticks,
                                                         labels.cex = 1.2,
                                                         labels.facing = "outside",
                                                         direction="outside",
                                                         # Conversion to Mb
                                                         labels = x_ticks_mb,
                                                         minor.ticks = 4,
                                                         labels.pos.adjust = TRUE)

                                   # Add sector names
                                   circlize::circos.text(xcenter,
                                                         sector_titles_expand * ymax,
                                                         sector_index,
                                                         cex = 1.5,
                                                         facing = "bending.inside",
                                                         niceFacing = TRUE)
                               }

                               for (i in 1:length(track$metrics)) {

                                   metric <- track$metrics[[i]]

                                   colors <- subset(data$Color,
                                                    data$Metric == metric$name &
                                                    data$Contig == sector_index)

                                   if (metric$type == "points") {

                                       # Plot the data as points
                                       circlize::circos.points(x, y,
                                                               cex = metric$point_size,
                                                               col = colors,
                                                               bg = colors,
                                                               pch = 21)

                                   } else if (metric$type == "ribbon") {

                                       # Plot the data as ribbon
                                       circlize::circos.lines(x, y,
                                                              col = colors[1],
                                                              border = colors[1],
                                                              lwd = 1,
                                                              area = TRUE,
                                                              baseline = ymin)

                                   }
                               }

                               # Add Y axis on the first sector only
                               if (sector_index == first_sector) {

                                   y_ticks_lab <- y_ticks
                                   # Round y-axis labels if values are big
                                   # enough
                                   if (y_ticks[3] - y_ticks[1] >= 10) {

                                       y_ticks_lab <- round(y_ticks, 0)

                                   } else {

                                       y_ticks_lab <- round(y_ticks, 2)

                                   }

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          at = y_ticks,  # 3 labels
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = FALSE,
                                                          labels = y_ticks_lab)

                                   # Add y axis labels
                                   # Axis title will be plotted 7.5% on the
                                   # left of the axis
                                   label_offset <- - 7.5 * (xmax - xmin) / (xplot[1] - xplot[2])
                                   circlize::circos.text(label_offset,
                                                         0.5 * (ymax - ymin) + ymin,
                                                         ylabel,
                                                         sector.index = first_sector,
                                                         facing = "clockwise",
                                                         cex = 1.3,
                                                         font = 2)
                               }

                               if (track$major_lines_x == TRUE) {

                                   n_major <- length(x_ticks)

                                   circlize::circos.segments(x_ticks,
                                                             rep(ymin, n_major),
                                                             x_ticks,
                                                             rep(ymax, n_major),
                                                             col = "grey50",
                                                             lty = 1)

                               }

                               if (track$major_lines_y == TRUE) {

                                   n_major <- length(major_y)
                                   circlize::circos.segments(rep(xmin, n_major),
                                                             major_y,
                                                             rep(xmax, n_major),
                                                             major_y,
                                                             col = "grey70",
                                                             lty = 3)

                               }

                           })

}





#' @title Check highlighted sectors
#'
#' @description Check that sectors to be highlighted exist,
#' allowing both contig names and chromosomes names
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded
#' with the \code{\link{load_genome_input}} function).
#'
#' @param contig_lengths Contig lengths from the output of the
#' \code{\link{load_genome_input}} function
#'
#' @param highlight Vector of names of sectors to highlight
#'
#' @return A clean vector of names for sectors to highlight
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#'
#' highlight <- check_highlight_sectors(genomic_data$data,
#'                                      genomic_data$lengths,
#'                                      c("Chr01", "NC_002364.1"))

check_highlight_sectors <- function(data, contig_lengths, highlight) {

    if (is.na(highlight)) { return(c()) }

    # Contig correspondence table
    contigs <- unique(data$Contig)

    for (i in 1:length(highlight)) {

        # Sector to highlight not found in list of sector lengths
        if (!(highlight[i] %in% names(contig_lengths))) {

            # Look for specified sector to highlight in original contig names
            if (highlight[i] %in% contigs) {

                tmp <- data$Contig_plot[which(data$Contig == highlight[i])]
                highlight[i] <- unique(tmp)[1]

            } else {

                stop(paste0("Could not find sector to highlight \"",
                             highlight[i], "\"."))

            }

        }

    }

    return(highlight)

}



#' @title Open output device
#'
#' @description Open an R output device based on an output file extension.
#'
#' @param output_file Path to the output file.
#'
#' @param width Plot width in pixel
#'
#' @param height Plot height in pixel
#'
#' @param res Image resolution in ppi
#'
#' @examples
#' open_output_device('circos.png', width = 1200, height = 1200, res = 120)
#' open_output_device('circos.pdf', width = 1400, height = 800, res = 160)

open_output_device <- function(output_file, width, height, res) {

    if (!is.na(output_file)) {

        extension <- strsplit(output_file, "\\.")[[1]]

        if (length(extension) < 2) {

            stop(paste0("Invalid output file name \"", output_file,
                        "\": could not detect extension"))

        } else {

            # Get last element in file name split by "."
            extension <- extension[length(extension)]

            if (extension == "pdf") {

                pdf(output_file, width = width, height = height)

            } else if (extension == "svg") {

                svg(output_file, width = width, height = height)

            } else if (extension == "png") {

                png(output_file, width = width, height = height, res = res)

            } else if (extension %in% c("jpeg", "jpg")) {

                jpeg(output_file, width = width, height = height, res = res)

            } else if (extension == "bmp") {

                bmp(output_file, width = width, height = height, res = res)

            } else if (extension == "tiff") {

                tiff(output_file, width = width, height = height, res = res)

            } else {

                stop(paste0("Invalid extension \"", extension,
                            "\" in output file\"", output_file,
                            ". Valid extensions: ",
                            ".pdf, .png, .svg, .bmp, .tiff, .jpeg, .jpg"))

            }

        }

    }

}



#' @title Configure sector background colors
#'
#' @description Configure background colors for all sectors in circos plot
#'
#' @param tracks List of track data objects for each track to include in the
#' plot. Track data objects are generated with the \code{\link{track}},
#' \code{\link{single_metric_track}}, or \code{\link{multi_metrics_track}}
#' functions.
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot.
#'
#' @param highlight_bg_color Background color for highlighted sectors.
#'
#' @return
#'
#' @examples
#' track <- single_metric_track("Fst", colors = "grey70", ylim = c(0, 1))
#' highlight <- c("LG07")
#' track <- configure_backgrounds(track, highlight, "grey70")

configure_backgrounds <- function(track, highlight, highlight_bg_color) {

    sectors <- unique(track$data$Contig)
    n_sectors <- length(sectors)
    track$bg_colors <- rep(track$bg_colors, n_sectors)[1:n_sectors]
    names(track$bg_colors) <- sectors
    track$bg_colors[highlight] <- highlight_bg_color

    return(track)

}
