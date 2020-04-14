#' @export
#'
#' @title Plot the results of radsex "depth" analysis
#'
#' @description Plot any metric computed with radsex "depth", either as a
#' boxplot showing the distribution of the metric for each group or as a barplot
#' showing the metric value for each individual.
#'
#' @param input_file Path to a table of marker depth statistics per individual
#' obtained with radsex "depth".
#'
#' @param metric Metric to display on the y axis, one of the columns of the
#' output of radsex "depth": "Reads", "Markers", "Min_depth", "Max_depth",
#' "Median_depth", "Average_depth" (default: "Markers").
#'
#' @param type Plot type, either "boxplot" or "barplot" (default: "boxplot").
#'
#' @param groups Vector specifying which groups to include in the plot, ordered.
#' Should match group names in the input data. If NA, groups will be infered
#' from the input data (default: NA).
#'
#' @param group_labels Vector specifying the label to associate to each group
#' in the plot axis / legend, or NA to use the group names (default: NA).
#'
#' @param group_colors Vector specyfing the color to use for each group in the
#' plot (default: c("dodgerblue3", "firebrick2")).
#'
#' @param alpha Alpha value for plot elements, a float (default: 0.8).
#'
#' @param output_file Path to an output file for the generated boxplot,
#' or NA to plot in the current R device (default: NA).
#'
#' @param width Plot width when plotting to an output file, in inches
#' (default: 8).
#'
#' @param height Plot height when plotting to an output file, in inches
#' (default: 8).
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
#' @param barplot_order_by Value to use to order bars in barplots, either
#' "Group", "Individual", or "Metric" (default: "Metric").
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' reads_plot <- radsex_depth("depth.tsv",
#'                            groups = c("M", "F"),
#'                            group_labels = c("Males", "Females"),
#'                            group_colors = c("green", "purple"),
#'                            metric = "Reads",
#'                            ylab = "Number of reads")
#'
#' depth_plot <- radsex_depth("depth.tsv",
#'                            by_group = FALSE,
#'                            group_colors = c("dodgerblue3", "red2"),
#'                            metric = "Median_depth",
#'                            ylab = "Number of reads",
#'                            output_file = "median_depth.png")

radsex_depth <- function(input_file,
                         metric = "Markers",
                         type = "boxplot",
                         groups = NA,
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
                         boxplot_jitter = TRUE,
                         barplot_order_by = "Metric") {

    # Load the data
    data <- load_table(input_file, get_properties = FALSE)$data

    # Infer groups from the data if not specified
    if (is.na(c(groups)[1])) { groups <- unique(data$Group) }

    # Order groups for plotting later
    data$Group <- factor(data$Group, levels = groups)

    # Generate a boxplot of distribution for each group
    if (type == "boxplot") {

        # Retain only useful columns
        data <- data[, c("Group", metric)]

        b <- draw_boxplot(data,
                          group_labels = group_labels,
                          group_colors = group_colors,
                          alpha = alpha,
                          output_file = output_file,
                          width = width,
                          height = height,
                          res = res,
                          title = title,
                          ylab = ylab,
                          major_lines_y = major_lines_y,
                          pretty_y_scale = pretty_y_scale,
                          jitter = boxplot_jitter)

    # Generate a barplot of the metric for each individual
    } else if (type == "barplot") {

        # Retain only useful columns
        data <- data[, c("Sample", "Group", metric)]

        b <- draw_barplot(data,
                          group_labels = group_labels,
                          group_colors = group_colors,
                          alpha = alpha,
                          output_file = output_file,
                          width = width,
                          height = height,
                          res = res,
                          title = title,
                          ylab = ylab,
                          major_lines_y = major_lines_y,
                          pretty_y_scale = pretty_y_scale,
                          order_by = barplot_order_by)

    } else {

        stop(paste0("Incorrect plot type \"", type,
                    "\" in radsex_depth",
                    "(availble plot types: \"boxplot\", \"barplot\")"))

    }

    return(b)
}

