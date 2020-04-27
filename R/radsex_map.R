#' @export
#'
#' @title RADSex circos plot
#'
#' @description Generates a circular plot of radsex "map" results in which each
#' sector represents a linkage group / chromosome and the x-axis represents the
#' position on the linkage group. The y-axis on the first track shows the bias
#' of a marker between groups and the second track shows the probability of
#' association with group for a marker.
#'
#' @param input_file Path to the result of radsex "map".
#' Format: Contig | Position | Length | Marker ID | Bias | P | Signif
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N at the
#' corresponding position on the corresponding contig.
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
#' @param output_file Path to an output file for the generated circos plot,
#' or NA to plot in the current R device (default: NA)
#'
#' @param width Plot width when plotting to an output file, in pixel
#' (default: 2400).
#'
#' @param height Plot height when plotting to an output file, in pixel
#' (default: 2400).
#'
#' @param res Image resolution when plotting to an output file, in ppi
#' (default: 120).
#'
#' @param tracks Vector of tracks to plot, in order from top to bottom; possible
#' tracks: "bias" for bias between groups, "p" for probability of association
#' with group (default: c("p", "bias")).
#'
#' @param bias_colors A vector of length 3 defining colors for extreme values in
#' the bias track, i.e. -1, 0, and 1
#' (default: c("firebrick1", "black", "dodgerblue2")).
#'
#' @param p_color Color of the track showing -log(p of association with group)
#' (default: "grey20").
#'
#' @param point_size Point size, a float (default: 0.3).
#'
#' @param highlight Vector containing the names or identifiers of contigs or
#' chromosomes to highlight in the circos plot (default NA).
#'
#' @param highlight_bg_color Background color for highlighted sectors
#' (default: "grey80").
#'
#' @param sector_titles_expand Manually set the space between sector titles
#' and x-axis as a multiple of ymax (default: NA).
#'
#' @examples
#' radsex_map_circos("radsex_map.tsv",
#'                   chromosomes_file = "chromosomes.tsv")

radsex_map_circos <- function(input_file,
                              chromosomes_file = NA,
                              detect_chromosomes = TRUE,
                              unplaced_label = "U.",
                              output_file = NA,
                              width = 2200,
                              height = 2200,
                              res = 120,
                              tracks = c("p", "bias"),
                              bias_colors = c("firebrick1",
                                              "black",
                                              "dodgerblue2"),
                              p_color = "grey40",
                              point_size = 0.3,
                              highlight = NA,
                              highlight_bg_color = "grey80",
                              sector_titles_expand = NA) {

    # Load chromosome names (return NA if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes,
                                unplaced_label = unplaced_label,
                                comment_char = "#",
                                comment_sep = ";",
                                comment_internal_sep = "#")

    defined_tracks <- tracks
    tracks <- list()
    track_n = 1

    for (i in 1:length(defined_tracks)) {

        if (defined_tracks[[i]] == "p") {

            # Transform the P-value column into -log10(p-value)
            data$data$P <- -log(data$data$P, 10)

            # Generate y-axis label for the probability of association with
            # group track
            p_label <- expression(paste("-log"[10], "(p)"))

            track <- single_metric_track("P",
                                         colors = p_color,
                                         point_size = point_size,
                                         alpha = 1,
                                         type = "points",
                                         label = p_label,
                                         bg_colors = "white",
                                         ylim = NA,
                                         major_lines_x = FALSE,
                                         major_lines_y = TRUE,
                                         legend_position = NA)

            tracks[[track_n]] <- track
            track_n <- track_n + 1

        } else if (defined_tracks[[i]] == "bias") {

            # Generate a color palette function for bias between groups
            bias_palette <- create_palette(bias_colors)

            track <- single_metric_track("Bias",
                                         colors = bias_palette,
                                         point_size = point_size,
                                         alpha = 1,
                                         type = "points",
                                         label =  "Bias",
                                         bg_colors = "white",
                                         ylim = c(-1, 1),
                                         major_lines_x = FALSE,
                                         major_lines_y = TRUE,
                                         legend_position = NA)

            tracks[[track_n]] <- track
            track_n <- track_n + 1

        } else {

            warning(paste0("Invalid track \"", defined_tracks[i],
                           "\" in radsex circos plot definition"))

        }

    }

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
                sector_titles_expand = sector_titles_expand)

}





#' @export
#'
#' @title RADSex manhattan plot
#'
#' @description Generates a Manhattan plot of radsex "map" results in which
#' -log10(probability of association with group) is plotted for all markers.
#'
#' @param input_file Path to the result of radsex "map".
#' Format: Contig | Position | Length | Marker ID | Bias | P | Signif
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N at the
#' corresponding position on the corresponding contig.
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
#' @param output_file Path to an output file for the generated circos plot,
#' or NA to plot in the current R device (default: NA)
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 12).
#'
#' @param height Plot height when plotting to an output file in inches,
#' (default: 6).
#'
#' @param res Image resolution when plotting to an output file, in dpi
#' (default: 300).
#'
#' @param colors A single color value for points that will be applied to all
#' chromosomes or a vector of alternating colors
#' (default: c("dodgerblue3", "darkgoldenrod2")).
#'
#' @param bg_colors A single color value for background that will be applied to
#' all chromosomes or a vector of alternating colors
#' (default: c("grey85", "white")).
#'
#' @param point_size Point size, a float (default: 1).
#'
#' @param x_title Title of the x-axis (default: "Chromosomes").
#'
#' @param show_chromosomes_names If TRUE, display chromosome names on the x axis
#' (default: TRUE).
#'
#' @param chromosomes_as_numbers If TRUE, display chromosome numbers instead of
#' names for readability (default: FALSE).
#'
#' @param show_signif_line If TRUE, show a horizontal line at the p-value
#' threshold for significance with Bonferroni correction (default: TRUE).
#'
#' @param signif_line_color Color of the significance line (default: "red").
#'
#' @return A ggplot object for the plot.
#'
#' @examples
#' radsex_map_manhattan("radsex_map.tsv",
#'                      output_file = "manhattan.svg",
#'                      colors = c("green", "purple"),
#'                      chromosomes_as_numbers = TRUE)

radsex_map_manhattan <- function(input_file,
                                 chromosomes_file = NA,
                                 detect_chromosomes = TRUE,
                                 unplaced_label = "U.",
                                 output_file = NA,
                                 width = 12,
                                 height = 6,
                                 res = 300,
                                 colors = c("dodgerblue3", "darkgoldenrod2"),
                                 bg_colors = c("grey85", "white"),
                                 point_size = 0.5,
                                 x_title = NA,
                                 show_chromosomes_names = TRUE,
                                 chromosomes_as_numbers = FALSE,
                                 show_signif_line = TRUE,
                                 signif_line_color = "red") {

    # Load chromosome names (return NA if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes,
                                unplaced_label = unplaced_label,
                                comment_char = "#",
                                comment_sep = ";",
                                comment_internal_sep = ":")

    # Transform the P-value column into -log10(p-value)
    data$data$P <- -log(data$data$P, 10)

    # Generate y-axis label for the probability of association with
    # group track
    y_label <- expression(bold(paste("-log"[10],
                                     "(p)"["association with group"])))

    tracks <- list(single_metric_track("P",
                                       colors = colors,
                                       point_size = point_size,
                                       alpha = 1,
                                       type = "points",
                                       label = y_label,
                                       bg_colors = bg_colors,
                                       ylim = NA,
                                       major_lines_x = FALSE,
                                       major_lines_y = TRUE,
                                       legend_position = "none"))

    if (show_signif_line) {

        s <- as.numeric(data$properties$signif_threshold)
        n_markers <- as.numeric(data$properties$n_markers)
        signif_threshold <- -log(s / n_markers, 10)
        tracks[[1]]$h_lines <- list(h_line(signif_threshold,
                                           label = paste0("p<", s),
                                           color = signif_line_color,
                                           type = 2,
                                           size = 0.75,
                                           label_font_size = 5))

    }

    # Draw the plot
    m <- draw_manhattan_plot(data$data,
                             data$lengths,
                             tracks,
                             output_file = output_file,
                             width = width,
                             track_height = height,
                             res = res,
                             x_title = x_title,
                             show_chromosomes_names = show_chromosomes_names,
                             chromosomes_as_numbers = chromosomes_as_numbers)

    return(m)

}





#' @export
#'
#' @title RADSex region plot
#'
#' @description Generates a linear plot of radsex "map" results for a given
#' genomic region, where the x-axis represents the position on the contig,
#' the y-axis on the first track shows the probability of association with group
#' for a marker and the second track shows the bias of a marker between groups,
#'
#' @param input_file Path to the result of radsex "map".
#' Format: Contig | Position | Length | Marker ID | Bias | P | Signif
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N at the
#' corresponding position on the corresponding contig.
#'
#' @param region Region to plot, defined with the syntax "Contig" or
#' "Contig:start-end"
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
#' or NA to plot in the current R device (default: NA)
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
#' @param tracks Vector of tracks to plot, in order from top to bottom; possible
#' tracks: "bias" for bias between groups, "p" for probability of association
#' with group (default: c("p", "bias")).
#'
#' @param bias_colors A vector of length 3 defining colors for extreme values in
#' the bias track, i.e. -1, 0, and 1
#' (default: c("firebrick1", "black", "dodgerblue2")).
#'
#' @param p_color Color of the track showing -log(p of association with group)
#' (default: "grey20").
#'
#' @param point_size Point size, a float (default: 1).
#'
#' @param alpha Alpha value, a float (default: 1).
#'
#' @param show_signif_line If TRUE, show a horizontal line at the p-value
#' threshold for significance with Bonferroni correction (default: TRUE).
#'
#' @param signif_line_color Color of the significance line (default: "red").
#'
#' @return A ggplot object for the plot.
#'
#' @examples
#' region_plot <- radsex_map_region("radsex_map.tsv",
#'                                  region = "Contig0001:1500-3500",
#'                                  tracks = c("bias"))

radsex_map_region <- function(input_file,
                              region,
                              chromosomes_file = NA,
                              detect_chromosomes = TRUE,
                              output_file = NA,
                              width = 12,
                              track_height = 4,
                              res = 300,
                              tracks = c("p", "bias"),
                              bias_colors = c("firebrick1",
                                              "black",
                                              "dodgerblue2"),
                              p_color = "grey40",
                              type = "points",
                              point_size = 1,
                              alpha = 1,
                              show_signif_line = TRUE,
                              signif_line_color = "red") {

    # Load chromosome names (return NA if no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes_file)

    # Load genomic metrics data
    data <- load_genome_metrics(input_file,
                                chromosomes = chromosomes,
                                detect_chromosomes = detect_chromosomes,
                                comment_char = "#",
                                comment_sep = ";",
                                comment_internal_sep = ":")

    defined_tracks <- tracks
    tracks <- list()
    track_n = 1

    for (i in 1:length(defined_tracks)) {

        if (defined_tracks[[i]] == "p") {

            # Transform the P-value column into -log10(p-value)
            data$data$P <- -log(data$data$P, 10)

            # Generate y-axis label for the probability of association with
            # group track
            p_label <- expression(bold(paste("-log"[10],
                                             "(p)"["association with group"])))

            track <- single_metric_track("P",
                                         colors = p_color,
                                         point_size = point_size,
                                         alpha = alpha,
                                         type = type,
                                         label = p_label,
                                         bg_colors = "white",
                                         ylim = NA,
                                         major_lines_x = FALSE,
                                         major_lines_y = TRUE,
                                         legend_position = "none")

            if (show_signif_line) {

                s <- as.numeric(data$properties$signif_threshold)
                n_markers <- as.numeric(data$properties$n_markers)
                signif_threshold <- -log(s / n_markers, 10)
                track$h_lines <- list(h_line(signif_threshold,
                                             label = paste0("p<", s),
                                             color = signif_line_color,
                                             type = 2,
                                             size = 0.75,
                                             label_font_size = 5))

            }

            tracks[[track_n]] <- track
            track_n <- track_n + 1



        } else if (defined_tracks[[i]] == "bias") {

            # Generate a color palette function for bias between groups
            bias_palette <- create_palette(bias_colors)

            track <- single_metric_track("Bias",
                                         colors = bias_palette,
                                         point_size = point_size,
                                         alpha = alpha,
                                         type = type,
                                         label = expression(bold("Group Bias")),
                                         bg_colors = "white",
                                         ylim = c(-1, 1),
                                         major_lines_x = FALSE,
                                         major_lines_y = TRUE,
                                         legend_position = "none")

            tracks[[track_n]] <- track
            track_n <- track_n + 1

        } else {

            warning(paste0("Invalid track \"", defined_tracks[i],
                           "\" in radsex circos plot definition"))

        }

    }

    r <- draw_region(data$data,
                     data$lengths,
                     region = region,
                     tracks = tracks,
                     output_file = output_file,
                     width = width,
                     track_height = track_height,
                     res = res)

    return(r)

}





#' @title Create bias color palette
#'
#' @description Generate a function to assign a color to a bias value based on
#' colors defined in bias_colors.
#'
#' @param bias_colors A vector of length 3 defining colors for the extremes in
#' the bias track, i.e. -1, 0, and 1.
#'
#' @return A function returning the color associated to a bias value.
#'
#' @examples
#' palette <- create_palette(c("blue", "black", "yellow"))
#' palette(0.25)

create_palette <- function(bias_colors) {

    # Create a palette made of two color ramping palettes meeting at 0
    palette <- c(grDevices::colorRampPalette(bias_colors[1:2])(10),
                 grDevices::colorRampPalette(bias_colors[2:3])(10)[-1])

    function(x) {

        # Convert a value between -1 and 1 into a value between 1 and 19 with
        # -1 --> 1, 0 --> 10, 1 --> 19
        x <- round((x + 1) * 9 + 1)
        x <- palette[x]

    }

}
