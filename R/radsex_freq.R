#' @title Plot the results of radsex "freq" analysis
#'
#' @description Plot any metric computed with radsex "freq", either as a single
#' boxplot showing the distribution of markers in all individuals or as a
#' barplot showing the metric value for each number of individuals.
#'
#' @param input_file Path to a table of marker presence in all individuals
#' obtained with radsex "freq".
#'
#' @param type Plot type, either "boxplot" or "barplot" (default: "barplot").
#'
#' @param color Plot main color (default: "grey60").
#'
#' @param alpha Alpha value for plot elements, a float (default: 1).
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
#' @param barplot_frequencies If TRUE, display values as frequency between 0
#' and 1 instead of absolute values (default: FALSE).
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' freq_plot <- radsex_freq("depth.tsv",
#'                          type = "boxplot",
#'                          color = "purple",
#'                          ylab = "Number of markers",
#'                          output_file = "markers_freq.png")

radsex_freq <- function(input_file,
                        color = "grey60",
                        alpha = 0.8,
                        output_file = NA,
                        width = 10,
                        height = 6,
                        res = 300,
                        title = NA,
                        ylab = "Marker presence (density)",
                        xlab = "Number of samples",
                        major_lines_y = TRUE) {

    # Load the data
    data <- load_table(input_file, get_properties = FALSE)$data

    # Convert the data into a vector of values
    data <- rep(data$Frequency, data$Count)

    # Draw the density plot
    b <- draw_density(data,
                      color = color,
                      alpha = alpha,
                      output_file = output_file,
                      width = width,
                      height = height,
                      res = res,
                      title = title,
                      ylab = ylab,
                      xlab = xlab,
                      major_lines_y = major_lines_y)

    return(b)
}

