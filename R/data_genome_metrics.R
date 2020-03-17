#' @title Load chromosome names
#'
#' @description Loads chromosome names from a tabulated file
#'
#' @param input_file_path Path to the chromosome names file (i.e. tab-separated
#' file without header and with columns <Contig ID> | <Chromosome name>).
#'
#' @return A named vector with contig identifiers as names and chromosome names
#' as values, or NULL if the input file path is NULL.
#'
#' @examples
#' chromosomes <- load_chromosome_names("chromosomes_names.tsv")
#'

load_chromosome_names <- function(input_file_path) {

    if (is.null(input_file_path)) return(NULL)

    raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t",
                                                   escape_double = FALSE,
                                                   col_names = FALSE,
                                                   trim_ws = TRUE))

    data <- gtools::mixedsort(setNames(raw_data$X2, raw_data$X1))

    return(data)

}





#' @title Load a genome metrics data file
#'
#' @description Loads a file containing metric values along a genome.
#' Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
#' Contig = contig identifier, Position = position on the contig, Length =
#' length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
#' corresponding position on the corresponding contig.
#'
#' @param input_file_path Path to a genome metrics file (format described above)
#'
#' @param chromosomes Named vector of chromosome names from
#' \code{\link{load_chromosome_names}} (default: NULL)
#'
#' @param detect_chromosomes If TRUE, will consider contigs starting with
#' "chromosome", "CHR", or "NC" as chromosomes if no chromosomes were specified
#' (default: TRUE)
#'
#' @param unplaced_label Label for unplaced contigs superscaffold
#' (default: "Unplaced")
#'
#' @return A list with two elements: "data" = parsed genome metrics data,
#' "lengths" = lengths of contigs included in the plot
#'
#' @examples
#' chromosomes <- load_chromosome_names("chromosomes_names.tsv")
#' genome_data <- load_genome_metrics("psass_window.tsv",
#'                                    chromosomes = chromosomes)

load_genome_metrics <- function(input_file_path,
                                chromosomes = NULL,
                                detect_chromosomes = TRUE,
                                unplaced_label = "Unplaced") {

    data <- suppressMessages(readr::read_delim(input_file_path,
                                               "\t",
                                               escape_double = FALSE,
                                               trim_ws = TRUE))

    # Create copies of Contig and Position columns before transformation
    data$Contig_plot <- data$Contig
    data$Position_plot <- data$Position

    # Detect chromosomes automatically if specified
    if (is.null(chromosomes) & detect_chromosomes) {

        chromosomes <- detect_chromosomes(data)

    }

    if (!is.null(chromosomes)) {

        # Separate data between chromosomes and unplaced contigs
        chromosome_data <- subset(data, data$Contig %in% names(chromosomes))
        chromosome_data$Contig_plot <- chromosomes[chromosome_data$Contig_plot]
        unplaced_data <- subset(data, !(data$Contig %in% names(chromosomes)))

        # Set chromosome color index to 2 for plotting
        chromosome_data$Color <- rep(2, nrow(chromosome_data))

        # Create a sorted data frame of contig lengths
        contig_lengths <- get_contig_lengths_from_data(data,
                                                       chromosomes,
                                                       sortby = "Contig")

        contig_lengths$Contig <- chromosomes[contig_lengths$Contig]

    } else {

        # No chromosomes: entire data is unplaced contigs
        unplaced_data <- data
        chromosome_data <- NULL

    }

    # If there are unplaced contigs, concatenate them in a superscaffold
    if (nrow(unplaced_data) > 0) {

        # Order unplaced data by contig length and then by position
        unplaced_lengths <- data.frame(unique(unplaced_data[, c(1, 3)]))
        unplaced_order <- unique(data.frame(unplaced_data[,3]))
        unplaced_lengths <- unplaced_lengths[order(unplaced_order,
                                                   decreasing = TRUE),]
        unplaced_data <- unplaced_data[order(match(unplaced_data$Contig,
                                                   unplaced_lengths$Contig),
                                             unplaced_data$Position), ]

        # Attribute an alternating color index to each unplaced contig (0/1)
        unplaced_color_index <- setNames(seq(1, nrow(unplaced_lengths)),
                                         unplaced_lengths$Contig)
        unplaced_data$Color <- unplaced_color_index[unplaced_data$Contig] %% 2

        # Transform position on each contig into position on cumulated contig
        increments <- c(0, head(cumsum(unplaced_lengths$Length), -1))
        cumulated_lengths <- setNames(increments, unplaced_lengths$Contig)
        unplaced_data$Position_plot <- (unplaced_data$Position +
                                        cumulated_lengths[unplaced_data$Contig])
        unplaced_data$Contig_plot <- unplaced_label

        # Compute total length of unplaced contigs
        total_unplaced_length <- sum(unplaced_lengths$Length)

    }

    # Create final data frame
    if (!is.null(chromosome_data) & nrow(unplaced_data) > 0) {

        data <- rbind(chromosome_data, unplaced_data)

        contig_lengths <- rbind(contig_lengths,
                                data.frame(Contig = unplaced_label,
                                           Length = total_unplaced_length))

    } else if (!is.null(chromosome_data)) {

        # Entire dataset is in chromosomes
        data <- chromosome_data

    } else if (nrow(unplaced_data) > 0) {

        # Entire dataset is unplaced
        data <- unplaced_data
        contig_lengths <- data.frame(Contig = c(unplaced_label),
                                     Length = c(total_unplaced_length))

    } else {

        stop(paste0("Error loading dataset <", input_file_path, ">",
                    " (impossible conditions)"))

    }

    return(list(data = data, lengths = setNames(contig_lengths$Length,
                                                contig_lengths$Contig)))

}





#' @title Detect chromosomes
#'
#' @description Automatically detect chromosomes in a genome metrics data frame.
#' Contigs are identified as chromosomes if their name starts with
#' "CHR", "LG", or "NC" (case unsensitive).
#'
#' @param data Genome metrics data frame: data frame with columns
#' Contig | Position | Length | <Metric1> | <Metric2> ...
#'
#' @return A named vector with contig identifiers as names and
#' chromosome names as values
#'
#' @examples
#' chromosome_names <- detect_chromosomes(data)

detect_chromosomes <- function(data) {

    valid_chromosome_start <- c("CHR", "LG", "NC")

    # Order contigs by length in a data frame
    contig_lengths <- get_contig_lengths_from_data(data)

    # Identify chromosomes based on name:
    # start with CHR, LG, or NC (case-unsensitive)
    contig_identifier_start <- toupper(substr(contig_lengths$Contig, 1, 2))
    detected_chr <- subset(contig_lengths,
                           contig_identifier_start %in% valid_chromosome_start)

    if (nrow(detected_chr) > 0) {

        # Sometimes mitochondria are also called NC_xxx.
        # If one chromosome is > 50 times smaller than the average of all other
        # chromosomes, or if it is smaller than 50000 bp,
        # it is considered to be the mitochondrion and is removed
        potential_mt <- tail(detected_chr, 1)
        if (50 * potential_mt$Length < median(detected_chr$Length) |
            potential_mt$Length < 50000) {

            detected_chr <- detected_chr[detected_chr$Contig != potential_mt,]

        }

    }

    # Create the named vector of chromosome names with detected chromosomes
    return(setNames(detected_chr$Contig, detected_chr$Contig))

}




#' @title Get contig lengths from genome metrics data
#'
#' @description Extract contig identifier and contig length columns from a
#' genome metrics data frame and return a sorted data frame.
#'
#' @param data Genome metrics data frame: data frame with columns
#' Contig | Position | Length | <Metric1> | <Metric2> ...
#'
#' @param chromosomes Named vector of chromosome names from
#' \code{\link{load_chromosome_names}}. If not NULL, only returns lengths for
#' chromosomes (default: NULL)
#'
#' @param sortby Column to sort the data frame by, either "Length" or "Contig"
#' (default: "Length")
#'
#' @return A data frame with columns <Contig ID> | <Contig Length> sorted by
#' contig length
#'
#' @examples
#' contig_lengths <- get_contig_lengths_from_data(data)

get_contig_lengths_from_data <- function(data, chromosomes = NULL, sortby = "Length") {

    contig_lengths <- data.frame(unique(data[, c("Contig", "Length")]))

    if (sortby == "Length") {

        contig_lengths <- contig_lengths[order(contig_lengths$Length,
                                               decreasing = TRUE),]

    } else if (sortby == "Contig") {

        contig_order <- gtools::mixedorder(contig_lengths$Contig)
        contig_lengths <- contig_lengths[contig_order,]

    } else {

        stop("Invalid value for parameter \"sortby\" in
             get_contig_lengths_from_data
             (accepted values: \"Length\", \"Contig\"")

    }

    if (!is.null(chromosomes)) {
        contig_lengths <- subset(contig_lengths,
                                 contig_lengths$Contig %in% names(chromosomes))
    }

    return(contig_lengths)

}
