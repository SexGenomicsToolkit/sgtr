% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track.R
\name{single_metric_track}
\alias{single_metric_track}
\title{Create a track object with a single metric}
\usage{
single_metric_track(
  metric,
  colors = NA,
  point_size = NA,
  alpha = NA,
  type = NA,
  label = NA,
  bg_colors = NA,
  ylim = NA,
  major_lines_x = NA,
  major_lines_y = NA,
  legend_position = NA,
  h_lines = NA
)
}
\arguments{
\item{metric}{Metric to represent in the track, should correspond to a
column name in the input data frame (i.e. the output of the
\code{\link{load_genome_metrics}} function).}

\item{colors}{Either a single color (e.g. "grey20") which will be
applied to all values, a vector of colors (e.g. c("blue", "yellow")) which
will be alternatingly applied to linkage groups / chromosomes, or a function
returning a color for a given value.}

\item{point_size}{Point size for plots of type "points", a float.}

\item{alpha}{Alpha value, a float.}

\item{type}{Plot type for circos and region tracks, either "points" or
"ribbon"; not used in manhattan plots.}

\item{label}{Label to use for the track, will be displayed on the y axis. If
NA, will use the metric name as label.}

\item{bg_colors}{Background colors for circos and manhattan plots, not used
in region plots. Values can be a string (e.g. "white") and will be applied
to all sectors / chromosomes, or a vector of strings
(e.g. c("white", "grey80")) to assign alternating values to
sectors / chromosomes.}

\item{ylim}{Y-axis limits, a vector with values c(min, max) or NA to infer
y limits from the data.}

\item{major_lines_x}{If TRUE, reference lines will be plotted for the x axis,
equivalent to panel.grid.major.x in \code{\link{ggplot2::theme}}.}

\item{major_lines_y}{If TRUE, reference lines will be plotted for the y axis,
equivalent to panel.grid.major.y in \code{\link{ggplot2::theme}}.}

\item{legend_position}{Position of the legend, directly passed to
"legend.position" in \code{\link{ggplot2::theme}} in region plots, not used
in manhattan and in circos plots.}

\item{h_lines}{List of objects defined with \code{\link{h_line}} adding
horizontal lines to the track (default: NA).}
}
\value{
A named list with the value of each track property
}
\description{
Simplified interface to generate a track representing a
single metric. All parameters have default values as NA, meaning their value
will be taken from default values specified in high-level functions.
}
\examples{
# Create a track object for the metric "Fst" in grey with a point size of
# 0.75 and y-axis limits of (0, 1).
fst_track <- single_metric_track("Fst",
                                 colors = "grey70",
                                 point_size = 0.75,
                                 ylim = c(0, 1))
}
