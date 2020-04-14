#' @title Draw a boxplot for a metric
#'
#' @description Draw a boxplot of the distribution of a metric for each group
#' in a data frame with columns: "Group | <Metric>".
#'
#' @param data A data frame with columns "Group | <Metric>". Groups can be
#' factors and will be plotted in order given by the factor levels or
#' strings and will be plotted in order of appearance in the data. The metric
#' is assumed to have continuous values.
#'
#' @param group_labels Vector specifying the label to associate to each group
#' in the plot x-axis, or NA to use the group names (default: NA).
#'
#' @param group_colors Vector specyfing the color to use for each group in the
#' plot (default: c("dodgerblue3", "firebrick2")).
#'
#' @param alpha Alpha value for boxplots, a float (default: 0.8).
#'
#' @param output_file Path to an output file for the generated boxplot,
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
#' @param title Plot title, or NA for no title (default: NA).
#'
#' @param ylab Title of the y axis, or NA to use the metric name (default: NA).
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis,
#' equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}
#' (default: TRUE).
#'
#' @param pretty_y_scale If TRUE, generate pretty labels for the y-axis scale
#' and add the unit (M. or K.) to the axis title (default: TRUE).
#'
#' @param boxplot_jitter If TRUE, superpose the distribution of values as points
#' in the boxplot (default: TRUE).
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' data <- load_table("depth.tsv")
#' data <- data[, c("Group", "Reads")]
#' reads_plot <- draw_boxplot(data,
#'                            group_labels = c("Males", "Females"),
#'                            group_colors = c("green", "purple"),
#'                            ylab = "Number of reads")

draw_boxplot <- function(data,
                         group_labels = NA,
                         group_colors = c("dodgerblue3", "firebrick2"),
                         alpha = 0.8,
                         output_file = NA,
                         width = 8,
                         height = 8,
                         res = 300,
                         title = NA,
                         ylab = NA,
                         major_lines_y = TRUE,
                         pretty_y_scale = TRUE,
                         jitter = TRUE) {

    # Get metric label and rename data columns to generic names
    metric <- names(data)[2]
    names(data) <- c("Group", "Metric")

    # Get ordered groups
    if (is.factor(data$Group)) {

        groups <- levels(data$Group)

    } else {

        groups <- unique(data$Group)
    }

    # Create named vector of group labels
    if (is.na(c(group_labels)[1])) {

        group_labels <- stats::setNames(groups, groups)

    } else {

        group_labels <- stats::setNames(group_labels, groups)

    }

    # Create named vector of group colors
    group_colors <- stats::setNames(group_colors, groups)

    # Assign y-axis label if not defined by user
    if (is.na(ylab)) { ylab <- metric }

    # Create major grid lines for y axis if specified
    if (major_lines_y) {

        major_lines_y <- ggplot2::element_line(color = "grey75", linetype = 3)

    } else {

        major_lines_y <- ggplot2::element_blank()

    }

    # Create a nice y-axis scale if specified
    if (pretty_y_scale) {

        y_axis <- generate_y_scale(data$Metric, ylab)

    } else {

        y_axis <- ggplot2::scale_y_continuous(name = metric)

    }

    # Remove outliers if jitter (to prevent double points)
    if (jitter) { outlier_alpha <- 0 } else { outlier_alpha <- 1 }

    # Generate the plot
    boxplot <- ggplot2::ggplot(data, ggplot2::aes(x = Group,
                                                  y = Metric,
                                                  fill = Group)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_boxplot(alpha = alpha, width = 0.75,
                              outlier.alpha = outlier_alpha) +
        ggplot2::scale_fill_manual(name = "", values = group_colors) +
        ggplot2::scale_x_discrete(name = "", labels = group_labels) +
        y_axis +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_text(face = "bold",
                                                            size = 22,
                                                            margin = ggplot2::margin(0, 10, 0, 5)),
                       axis.text.x = ggplot2::element_text(face = "bold",
                                                           size = 20,
                                                           margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.text.y = ggplot2::element_text(face = "bold",
                                                           size = 18),
                       panel.grid.major.y = major_lines_y,
                       legend.position = "none")

    if (jitter) {

        boxplot <- boxplot +
            ggplot2::geom_jitter(ggplot2::aes(fill = Group),
                                 fill = "grey10",
                                 color = "grey10",
                                 shape = 21,
                                 alpha = 1,
                                 width = 0.15, size = 3) +
            ggplot2::scale_color_manual(name = "", values = group_colors)
    }

    # Add title if specified
    if (!is.na(title)) {

        boxplot <- boxplot +
            ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = 20,
                                                              face = "bold",
                                                              margin = ggplot2::margin(0, 0, 10, 0)))

    } else {

        boxplot <- boxplot +
            ggplot2::theme(plot.title = ggplot2::element_blank())

    }

    if (is.na(c(groups)[1])) {

        boxplot <- boxplot + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                            axis.ticks.length.x = ggplot2::element_blank())

    }

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = boxplot,
                        width = width,
                        height = height,
                        dpi = res)

    } else {

        print(boxplot)

    }

    return(boxplot)

}





#' @title Draw a barplot for a metric
#'
#' @description Draw a barplot of the value of a metric for each Sample
#' in a data frame with columns: "Sample | Group | <Metric>".
#'
#' @param data A data frame with columns "Sample | Group | <Metric>". The
#' metric is assumed to have continuous values.
#'
#' @param group_labels Vector specifying the label to associate to each group
#' in the plot legend, or NA to use the group names (default: NA).
#'
#' @param group_colors Vector specyfing the color to use for each group in the
#' plot (default: c("dodgerblue3", "firebrick2")).
#'
#' @param alpha Alpha value for the bars, a float (default: 0.8).
#'
#' @param output_file Path to an output file for the generated barplot,
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
#' @param title Plot title, or NA for no title (default: NA).
#'
#' @param ylab Title of the y axis, or NA to use the metric name (default: NA).
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis,
#' equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}
#' (default: TRUE).
#'
#' @param pretty_y_scale If TRUE, generate pretty labels for the y-axis scale
#' and add the unit (M. or K.) to the axis title (default: TRUE).
#'
#' @param order_by Value to use to order bars in the plot, either "Group",
#' "Sample", or "Metric" (default: "Metric").
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' data <- load_table("depth.tsv")
#' data <- data[, c("Sample, "Group", "Markers")]
#' markers_plot <- draw_barplot(data,
#'                              group_labels = c("Males", "Females"),
#'                              group_colors = c("blue", "red"),
#'                              ylab = "Number of reads")

draw_barplot <- function(data,
                         group_labels = NA,
                         group_colors = c("dodgerblue3", "firebrick2"),
                         alpha = 0.8,
                         output_file = NA,
                         width = 8,
                         height = 8,
                         res = 300,
                         title = NA,
                         ylab = NA,
                         major_lines_y = TRUE,
                         pretty_y_scale = TRUE,
                         order_by = "Metric") {


    # Get metric label and rename data columns to generic names
    metric <- names(data)[3]
    names(data) <- c("Sample", "Group", "Metric")

    # Get ordered groups
    if (is.factor(data$Group)) {

        groups <- levels(data$Group)

    } else {

        groups <- unique(data$Group)

    }

    # Create named vector of group labels
    if (is.na(c(group_labels)[1])) {

        group_labels <- stats::setNames(groups, groups)

    } else {

        group_labels <- stats::setNames(group_labels, groups)

    }

    group_colors <- stats::setNames(group_colors, groups)

    if (order_by == "Group") {

        order <- order(data$Group, data$Metric)

    } else {

        order <- order(unlist(data[, order_by]))

    }

    # Assign y-axis label if not defined by user
    if (is.na(ylab)) { ylab <- metric }

    # Create major grid lines for y axis if specified
    if (major_lines_y) {

        major_lines_y <- ggplot2::element_line(color = "grey75", linetype = 3)

    } else {

        major_lines_y <- ggplot2::element_blank()

    }

    # Create a nice y-axis scale if specified
    if (pretty_y_scale) {

        y_axis <- generate_y_scale(data$Metric, ylab, min = 0)

    } else {

        y_axis <- ggplot2::scale_y_continuous(name = metric)

    }

    data <- data[order, ]
    data$Sample <- factor(data$Sample, levels = data$Sample)

    # Generate the plot
    barplot <- ggplot2::ggplot(data, ggplot2::aes(x = Sample,
                                                  y = Metric,
                                                  fill = Group)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_col(alpha = alpha) +
        ggplot2::scale_fill_manual(name = "", values = group_colors) +
        y_axis +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_text(face = "bold",
                                                            size = 20,
                                                            margin = ggplot2::margin(0, 10, 0, 5)),
                       axis.text.x = ggplot2::element_text(size = 16,
                                                           angle = 90,
                                                           vjust = 0.5,
                                                           margin = ggplot2::margin(0, 0, 0, 0)),
                       axis.text.y = ggplot2::element_text(face = "bold",
                                                           size = 18),
                       panel.grid.major.y = major_lines_y,
                       legend.position = "right",
                       legend.key.size = ggplot2::unit(0.06, "npc"),
                       legend.text = ggplot2::element_text(face = "bold",
                                                           size = 18))

    # Add title if specified
    if (!is.na(title)) {

        barplot <- barplot +
            ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = 20,
                                                              face = "bold",
                                                              margin = ggplot2::margin(0, 0, 10, 0)))

    } else {

        barplot <- barplot +
            ggplot2::theme(plot.title = ggplot2::element_blank())

    }

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = barplot,
                        width = width,
                        height = height,
                        dpi = res)

    } else {

        print(barplot)

    }

    return(barplot)

}





#' @title Draw density for a metric
#'
#' @description Draw a density plot for a vector of values.
#'
#' @param data A vector of values.
#'
#' @param color Color of the plot (default: "grey60").
#'
#' @param alpha Alpha value for the bars, a float (default: 0.8).
#'
#' @param output_file Path to an output file for the generated barplot,
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
#' @param title Plot title, or NA for no title (default: NA).
#'
#' @param ylab Title of the y axis, or NA to use the metric name (default: NA).
#'
#' @param xlab Title of the x axis, or NA for no title (default: NA).
#'
#' @param major_lines_y If TRUE, reference lines will be plotted for the y axis,
#' equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}
#' (default: TRUE).
#'
#' @param pretty_y_scale If TRUE, generate pretty labels for the y-axis scale
#' and add the unit (M. or K.) to the axis title (default: TRUE).
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' data <- load_table("freq.tsv")
#' data <- rep(data$Frequency, data$Count)
#' freq_plot <- draw_density(data)

draw_density <- function(data,
                         color = "grey60",
                         alpha = 0.8,
                         output_file = NA,
                         width = 12,
                         height = 6,
                         res = 300,
                         title = NA,
                         ylab = NA,
                         xlab = NA,
                         major_lines_y = TRUE) {

    # Create major grid lines for y axis if specified
    if (major_lines_y) {

        major_lines_y <- ggplot2::element_line(color = "grey75", linetype = 3)

    } else {

        major_lines_y <- ggplot2::element_blank()

    }

    # Convert data into single column data frame
    data <- data.frame(Value = data)

    # Generate the plot
    density <- ggplot2::ggplot(data, ggplot2::aes(x = Value)) +
        cowplot::theme_cowplot() +
        ggplot2::geom_density(color = color, fill = color, alpha = alpha) +
        ggplot2::scale_x_continuous(expand = c(0.02, 0),
                                    breaks = ggplot2::waiver()) +
        ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
        ggplot2::theme(plot.title = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(face = "bold",
                                                           size = 18),
                       axis.text.y = ggplot2::element_text(face = "bold",
                                                           size = 18),
                       panel.grid.major.y = major_lines_y,
                       legend.position = "none")

    # Add ylab if specified
    if (!is.na(ylab)) {

        density <- density +
            ggplot2::ylab(ylab) +
            ggplot2::theme(axis.title.y = ggplot2::element_text(face = "bold",
                                                                size = 20,
                                                                angle = 90,
                                                                margin = ggplot2::margin(0, 10, 0, 5)))

    }

    # Add xlab if specified
    if (!is.na(xlab)) {

        density <- density +
            ggplot2::xlab(xlab) +
            ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold",
                                                                size = 20,
                                                                margin = ggplot2::margin(5, 0, 5, 0)))

    }

    # Add title if specified
    if (!is.na(title)) {

        density <- density +
            ggplot2::ggtitle(title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = 20,
                                                              face = "bold",
                                                              margin = ggplot2::margin(0, 0, 10, 0)))

    }

    # Output to file if specified or print in current R device otherwise
    if (!is.na(output_file)) {

        ggplot2::ggsave(output_file,
                        plot = density,
                        width = width,
                        height = height,
                        dpi = res)

    } else {

        print(density)

    }

    return(density)

}