% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radsex_marker_depths.R
\name{draw_marker_depths}
\alias{draw_marker_depths}
\title{RADSex marker depths plot}
\usage{
draw_marker_depths(
  data,
  group_info = NA,
  output_file = NA,
  width = 10,
  height = 8,
  res = 300,
  title = NA,
  label_colors = c("dodgerblue3", "red3"),
  depth_colors = c("white", "royalblue2", "black", "gold2", "red3"),
  color_by_presence = FALSE,
  presence_color = "grey30",
  presence_min_depth = 5,
  show_individuals = TRUE,
  show_markers = FALSE,
  show_individuals_dendrogram = TRUE,
  show_markers_dendrogram = FALSE
)
}
\arguments{
\item{data}{Table of marker depths obtained with the \code{\link{load_table}}
function.}

\item{group_info}{A table of individual group information obtained with the
\code{\link{load_group_info}} function, or NA to not use group information
(default: NA).}

\item{output_file}{Path to an output file for the generated marker depths
plot, or NA to plot in the current R device (default: NA).}

\item{width}{Plot width when plotting to an output file, in inches
(default: 10).}

\item{height}{Plot height when plotting to an output file, in inches
in inches (default: 8).}

\item{res}{Image resolution when plotting to an output file, in dpi
(default: 300).}

\item{title}{Plot title (default: NA, i.e. no title).}

\item{label_colors}{Vector of length 2 specifying label color for individuals
from each group if individual group information is used
(default: c("dodgerblue3", "firebrick2)).}

\item{depth_colors}{Vector of length 5 specifying colors for the following
values: [0, (1:mean), (mean:3rd quartile), (3rd quartile:(max - 1)), max]
(default: c("white", "royalblue2", "black", "gold2", "red3"))}

\item{color_by_presence}{If TRUE, color markers by presence / absence using
the parameters "presence_color" and "presence_min_depth" (default: FALSE).}

\item{presence_color}{Color of present markers when coloring by presence /
absence of markers (default: "grey30").}

\item{presence_min_depth}{Minimum depth to consider a marker present in an
individual when coloring by presence / absence (default: 5).}

\item{show_individuals}{If TRUE, display individual names on the x-axis
(default: TRUE).}

\item{show_markers}{If TRUE, display marker names on the y-axis
(defautlt: FALSE).}

\item{show_individuals_dendrogram}{If TRUE, display the individual clustering
dendrogram on the x-axis (default: TRUE).}

\item{show_markers_dendrogram}{If TRUE, display the marker clustering
dendrogram on the y-axis (default: FALSE).}
}
\description{
Draw a heatmap of individual depths for a subset of markers from
radsex "process", "subset", or "signif". In the resulting heatmap, the color
of a tile at coordinates (x, y) indicates the depth of the y-th marker in
the x-th individual.
}
\examples{
data <- load_table("subset.tsv")
group_info <- load_group_info("popmap.tsv")
marker_depths <- draw_marker_depths(data,
                                    group_info = group_info,
                                    title = "Marker depths clustered",
                                    label_colors = c("green", "purple"),
                                    show_individuals = TRUE,
                                    show_markers = TRUE,
                                    show_individuals_dendrogram = TRUE,
                                    show_markers_dendrogram = TRUE)
}
