% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_manhattan.R
\name{plot_manhattan}
\alias{plot_manhattan}
\title{Plot genome metrics in a manhattan plot}
\usage{
plot_manhattan(
  input_file,
  tracks,
  chromosomes_file = NA,
  detect_chromosomes = TRUE,
  unplaced_label = "U.",
  comment_char = "#",
  comment_sep = ";",
  comment_internal_sep = ":",
  output_file = NA,
  width = 12,
  track_height = 6,
  res = 300,
  x_title = "Chromosomes",
  show_chromosomes_names = TRUE,
  chromosomes_as_numbers = FALSE,
  default_colors = c("dodgerblue3", "darkgoldenrod2"),
  default_bg_colors = c("grey85", "white"),
  default_point_size = 0.5,
  default_ylim = NA,
  default_alpha = 1,
  default_major_lines_y = TRUE
)
}
\arguments{
\item{input_file}{Path to a genome metrics input file (e.g. result of
PSASS or RADSex).
Format: Contig | Position | Length | <Metric1> | <Metric2> ... with
Contig = contig identifier, Position = position on the contig, Length =
length of the contig, <MetricN> = value of metric N (e.g. Fst) at the
corresponding position on the corresponding contig.}

\item{tracks}{List of track data objects for each track to include in the
plot. Track data objects for manhattan plots are generated with the
\code{\link{track}} function. Tracks can represent one or mutiple
metrics from the genome metrics input file.}

\item{chromosomes_file}{Path to the chromosome names file
(i.e. tab-separated file without header and with columns
<Contig ID> | <Chromosome name>). If NA, all contigs will be considered
unplaced except if detect_chromosomes is set to TRUE, in which case
chromosomes will be detected automatically from contig identifiers
(default: NA).}

\item{detect_chromosomes}{If TRUE, will consider contigs starting with
"LG", "CHR", or "NC" as chromosomes if no chromosomes were specified
(default: TRUE).}

\item{unplaced_label}{Label for unplaced contigs superscaffold
(default: "Unplaced").}

\item{comment_char}{Character indicating a comment line in the input file
(default: "#").}

\item{comment_sep}{Character separating two fields in a comment line
(default: ";").}

\item{comment_internal_sep}{Character separating property and value in a
field from a comment line (default: ":").#'}

\item{output_file}{Path to an output file for the generated manhattan plot,
or NA to plot in the current R device (default: NA).}

\item{width}{Plot width when plotting to an output file, in inches
(default: 12).}

\item{track_height}{Height of a single track when plotting to an output file,
in inches (default: 6).}

\item{res}{Image resolution when plotting to an output file, in dpi
(default: 300).}

\item{x_title}{Title of the x-axis (default: "Chromosomes").}

\item{show_chromosomes_names}{If TRUE, display chromosome names on the x axis
(default: TRUE)}

\item{chromosomes_as_numbers}{If TRUE, display chromosome numbers instead of
names for readability (default: FALSE)}

\item{default_colors}{Default colors when not specified in track data
(default: "dodgerblue3", "darkgoldenrod2").}

\item{default_bg_colors}{Default background colors for standard sectors
when not specified in track data (default: c("white")).}

\item{default_point_size}{Default point size when not specified in track data
(default: 0.5).}

\item{default_ylim}{Default y-axis limits when not specified in track data
(default: NA, i.e. infer from data).}

\item{default_alpha}{Default alpha value when not specified in track data
(default: 1).}

\item{default_major_lines_y}{Default value for drawing major lines on the
y axis. If TRUE, reference lines will be plotted (default: TRUE).}
}
\description{
Generate a manhattan plot with one or multiple tracks from a
genome metrics file
}
\examples{
plot_manhattan("data/psass_window.tsv",
               tracks = list(track("Fst",
                                   label = expression("F"["ST"])),
                             track("Snps_females",
                                   label = "Pool-specific SNPs",
                                   colors = "firebrick2"),
                             track("Snps_males",
                                   label = "Pool-specific SNPs",
                                   colors = "dodgerblue3"))
               default_bg_colors = c("white", "grey60"),
               output_file = "manhattan.png")
}
