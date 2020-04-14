#' @title Parse region
#'
#' @description Parse a region string with format "Contig:start-end"
#'
#' @param region String defining the region, either "Contig" or
#' "Contig:start-end"
#'
#' @param contig_lengths Named vector with contig identifers as names and
#' contig lengths as values
#' (e.g. output of the \code{\link{load_genome_input}} function)
#'
#' @return A list of three elements: contig name, region start, and region end
#'
#' @examples
#' data <- load_genome_input("psass_window.tsv")
#' region_info <- parse_region("Chr01:150000-250000", data$lengths)

parse_region <- function(region, contig_lengths) {

    # Get contig
    tmp = strsplit(region, ":")
    contig = tmp[[1]][1]

    # Check that the specified contig is in list of contig lengths
    # Note: in the draw_region function, unplaced contigs are added to
    # contig_lengths before parsing the region string, so that unplaced
    # contigs can be plotted (by default only chromosomes are in the vector)
    if (!(contig %in% names(contig_lengths))) {
        stop(paste0("Error: invalid contig in region string <", region, ">"))
    }

    if (length(tmp[[1]]) == 2) {  # Case "Contig:start-end"

        tmp = strsplit(tmp[[1]][2], "-")
        start = as.numeric(tmp[[1]][1])
        end = as.numeric(tmp[[1]][2])

    } else {  # Case "Contig": get end from vector of contig lengths

        start = 0
        end = unname(contig_lengths[contig])

    }

    return(list(contig = contig, start = start, end = end))

}





#' @title Convert position to Mb
#'
#' @description Convert a genomic position to Mega bp (round to n digits).
#' This function can be called directly to generate labels in plot scales.
#' e.g. ggplot2::scale_y_continuous(labels = convert_to_mb)
#'
#' @param x Genomic position to convert, integer
#'
#' @param n Number of digits to round the value, integer (default: 0)
#'
#' @return The converted value, floating point number
#'
#' @examples
#' x_mb <- convert_to_mb(1267356, 2)
#' x_mb
#' > 1.27

convert_to_mb <- function(x, n = 0) {

    round(x / 10^6, n)

}





#' @title Convert position to Kb
#'
#' @description Convert a genomic position to Kilo bp (round to n digits).
#' This function can be called directly to generate labels in plot scales.
#' e.g. ggplot2::scale_y_continuous(labels = convert_to_kb)
#'
#' @param x Genomic position to convert, integer
#'
#' @param n Number of digits to round the value, integer (default: 0)
#'
#' @return The converted value, floating point number
#'
#' @examples
#' x_mb <- convert_to_kb(1267356, 1)
#' x_mb
#' > 1267.4

convert_to_kb <- function(x, n = 0) {

    round(x / 10^3, n)

}





#' @title Generate x-axis scale
#'
#' @description Generate an x-axis scale with a given number of ticks,
#' automatically converting genomic positions to the best unit (Mbp, Kbp ...).
#'
#' @param region_info A list with three elements: contig name, region start,
#' region end (e.g. output of the \code{\link{parse_region}} function).
#'
#' @param n_ticks Number of ticks to include in the scale (default: 10).
#'
#' @return A ggplot2::scale_x_continuous object encoding the scale.
#'
#' @examples
#' data <- load_genome_input("psass_window.tsv")
#' region_info <- parse_region("Chr01:150000-250000", data$lengths)
#' nice_x_scale <- generate_x_scale(region_info, n_ticks = 5)

generate_x_scale <- function(region_info, n_ticks = 10) {

    # Compute the size of an inter-break interval
    S <- (region_info$end - region_info$start) / n_ticks
    if (S <= 0) {
        stop(paste0("Error: the region size has to be > 0 (value: ",
                    S * n_ticks, ")"))
    }

    # Find the order of magnitude of an inter-break interval
    N <- floor(log(S, 10))
    N10 <- 10 ^ N  # N10 is an 'exponentially-rounded' inter-break interval

    # Generate a scale of n_ticks values within the region
    # Start, end, and break size values are rounded for nicer display
    scale <- seq(N10 * floor(region_info[[2]] / N10),
                 N10 * ceiling(region_info[[3]] / N10),
                 N10 * round(S / N10, 1))

    # Adjust the labels based on the size of the values (Mbp or Kbp)
    if (region_info[[2]] < 10 ^ 6 & region_info[[3]] < 10 ^ 6) {

        scale_labels <- round(scale / 10 ^ 3, 3 - N)
        bp_unit <- "K"

    } else {

        scale_labels <- round(scale / 10 ^ 6, 6 - N)
        bp_unit <- "M"

    }

    # Generate the scale object
    output <- ggplot2::scale_x_continuous(name = paste0("Position on ",
                                                        region_info[[1]],
                                                        " (", bp_unit, "bp)"),
                                          expand = c(0.01, 0.01),
                                          breaks = scale,
                                          labels = scale_labels,
                                          limits = c(min(scale[1],
                                                         region_info[[2]]),
                                                     max(tail(scale, 1),
                                                         region_info[[3]])))

    return(output)

}





#' @title Generate y-axis scale
#'
#' @description Generate an x-axis scale with a given number of ticks,
#' automatically converting values to the best unit (M., K., ...).
#'
#' @param values A vector of values that will be plotted on the y axis.
#'
#' @param metric The metric name.
#'
#' @param n_ticks Number of ticks to include in the scale (default: 5)
#'
#' @param min Minimum value for the scale, or NA to infer from the data
#' (default: NA)
#'
#' @param max Maximum value for the scale, or NA to infer from the data
#' (default: NA)
#'
#' @return A ggplot2::scale_y_continuous object encoding the scale.
#'
#' @examples
#' nice_y_scale <- generate_y_scale(seq(1000000, 10000000, 1000000),
#'                                  metric = "Reads", n_ticks = 4)

generate_y_scale <- function(values, metric, n_ticks = 5, min = NA, max = NA) {

    start = min(values)
    end = max(values)

    # Compute the size of an inter-break interval
    S <- (end - start) / n_ticks
    if (S <= 0) {
        stop(paste0("Values range has to be > 0 (range: ", S * n_ticks, ")"))
    }

    # Find the order of magnitude of an inter-break interval
    N <- floor(log(S, 10))
    N10 <- 10 ^ N  # N10 is an 'exponentially-rounded' inter-break interval

    # Generate a scale of n_ticks values within the region
    # Start, end, and break size values are rounded for nicer display
    scale <- seq(N10 * floor(start / N10),
                 N10 * ceiling(end / N10),
                 N10 * round(S / N10, 1))

    # Adjust the labels based on the size of the values (Mbp or Kbp)
    if (start < 10 ^ 3 & end < 10 ^ 3) {

        N <- max(0, -N)
        scale_labels <- format(round(scale, N), nsmall = N)
        unit <- ""

    } else if (start < 10 ^ 6 & end < 10 ^ 6) {

        N <- max(0, 3 - N)
        scale_labels <- format(round(scale / 10 ^ 3, N), nsmall = N)
        unit <- "(K.)"

    } else {

        N <- max(0, 6 - N)
        scale_labels <- format(round(scale / 10 ^ 6, N), nsmall = N)
        unit <- "(M.)"

    }

    # Assign values for y-axis limits
    ymin <- min(scale[1], start)
    if (ymin < 0) { ymin <- 1.025 * ymin } else { ymin <- 0.976 * ymin }
    ymax <- max(tail(scale, 1), end)
    if (ymax < 0) { ymax <- 0.976 * ymax } else { ymax <- 1.025 * ymax }
    ylim <- c(ymin, ymax)

    # Override min / max if specified in function call
    if (!is.na(min)) { ylim[1] <- min}
    if (!is.na(max)) { ylim[2] <- max}

    # Generate the scale object
    output <- ggplot2::scale_y_continuous(name = paste0(metric, " ", unit),
                                          expand = c(0.01, 0.01),
                                          breaks = scale,
                                          labels = scale_labels,
                                          limits = ylim)

    return(output)

}





#' @title Drop columns from a data frame
#'
#' @description Remove columns from a data frame by name. Generate a warning
#' if columns to remove do not exist in the data frame.
#'
#' @param data A data frame.
#'
#' @param columns A string (i.e. "value") or vector of strings (
#' i.e. c("X1", "X2")) of columns to drop from the data frame.
#'
#' @return A data frame in which the specified columns were removed.
#'
#' @examples
#' data <- data.frame(x = c(1), y = c(2), z = c(3))
#' data <- drop_columns(data, c("x", "z"))

drop_columns <- function(data, columns) {

    columns <- c(columns)

    for (i in 1:length(columns)) {

        if (!(columns[i] %in% names(data))) {

            warning(paste0("Ccolumn \"", columns[i], "\" was not found in the ",
                           "data frame while trying to remove the column"))

        }

    }

    return(data[, -which(names(data) %in% columns)])

}