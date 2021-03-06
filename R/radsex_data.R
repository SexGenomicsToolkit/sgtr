#' @export
#'
#' @title Load a marker distribution table
#'
#' @description Load a distribution of markers between groups obtained with
#' radsex "distrib"
#'
#' @param input_file Path to a table of distribution of markers between groups
#' obtained with radsex "distrib".
#'
#' @param groups Vector of groups to include in the plots, in order; group names
#' should match column headers in the input data. If NA, groups will be infered
#' from the input data (default: NA).
#'
#' @param group_labels Vector of length 2 with the labels to associate to each
#' group in the plot axis titles, or NA to use the group names (default: NA).
#'
#' @return A list with the following elements:
#' \item{data}{A data frame with the distribution of markers between groups}
#' \item{groups}{Vector of group names}
#' \item{group_labels}{Vector with groups as names and group labels as values}
#' \item{counts}{Vector with groups as names and individual counts as values}
#'
#' @examples
#' data <- load_marker_distribution("distrib.tsv",
#'                                  groups = c("Pop1", "Pop2"),
#'                                  group_labels = c("France", "Spain"))

load_marker_distribution <- function(input_file,
                                     groups = NA,
                                     group_labels = NA) {

    data <- load_table(input_file)
    properties <- data$properties
    data <- data$data

    # Convert Signif to boolean for old versions of readr
    if (!is.logical(data$Signif)) {

        data$Signif <- c("True" = TRUE, "False" = FALSE)[data$Signif]

    }

    # Assign group names if not specified by user
    if (is.na(c(groups)[1])) { groups <- names(data)[1:2] }

    # Assign group labels if not specified by the user
    if (is.na(c(group_labels)[1])) { group_labels = groups }

    # Create group -> label correspondence table
    group_labels <- stats::setNames(group_labels, groups)

    # Compute number of individuals in each group
    counts <- stats::setNames(c(max(data[[groups[1]]]), max(data[[groups[2]]])),
                              groups[1:2])

    if (counts[groups[1]] == 0) {

        stop("No individuals were found in specified group \"", groups[1], "\"")

    }

    if (counts[groups[2]] == 0) {

        stop("No individuals were found in specified group \"", groups[2], "\"")

    }

    return(list(data = data,
                properties = properties,
                groups = groups,
                group_labels = group_labels,
                counts = counts))

}





#' @export
#'
#' @title Load a standard table file
#'
#' @description Loads a tabulated file with readr
#'
#' @param input_file Path to a tabulated file
#'
#' @return A data frame containing the tabulated file data.
#'
#' @examples
#' data <- load_table("subset.tsv")

load_marker_depths <- function(input_file) {

    data <- suppressMessages(readr::read_delim(input_file,
                                               "\t",
                                               comment = "#",
                                               escape_double = FALSE,
                                               col_names = TRUE,
                                               trim_ws = TRUE))

    return(data)

}





#' @export
#'
#' @title Load group info
#'
#' @description Load group information from a popmap file.
#'
#' @param input_file Path to a tabulated file with columns "Individual | Group"
#'
#' @param groups Vector of groups to include in the plots, in order; group names
#' should match column headers in the input data. If NA, groups will be infered
#' from the input data (default: NA).
#'
#' @param group_labels Vector of length 2 with the labels to associate to each
#' group in the plot axis titles, or NA to use the group names (default: NA).
#'
#' @return A list with the following elements:
#' \item{individual_groups}{Vector with individuals as names and groups as values}
#' \item{groups}{Vector of group names}
#' \item{group_labels}{Vector with groups as names and group labels as values}
#' \item{counts}{Vector with groups as names and individual counts as values}
#'
#' @examples
#' popmap <- load_group_info("popmap.tsv",
#'                           groups = c("M", "F"),
#'                           group_labels = c("Males", "Females"))

load_group_info <- function(input_file,
                            groups = NA,
                            group_labels = NA) {

    data <- suppressMessages(readr::read_delim(input_file,
                                               "\t",
                                               escape_double = FALSE,
                                               col_names = FALSE,
                                               trim_ws = TRUE))

    individual_groups <- stats::setNames(data$X2, data$X1)

    counts <- table(individual_groups)

    # Assign group names if not specified by user
    if (is.na(c(groups)[1])) { groups <- names(counts) }

    # Assign group labels if not specified by the user
    if (is.na(c(group_labels)[1])) { group_labels = groups }

    # Create group -> label correspondence table
    group_labels <- stats::setNames(group_labels, groups)

    # Generate output
    output <- list(individual_groups = individual_groups,
                   groups = groups,
                   group_labels = group_labels,
                   counts = counts)

    return(output)

}
