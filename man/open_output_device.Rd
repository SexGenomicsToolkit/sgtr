% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_circos.R
\name{open_output_device}
\alias{open_output_device}
\title{Open output device}
\usage{
open_output_device(output_file, width, height, res)
}
\arguments{
\item{output_file}{Path to the output file.}

\item{width}{Plot width in pixel}

\item{height}{Plot height in pixel}

\item{res}{Image resolution in ppi}
}
\description{
Open an R output device based on an output file extension.
}
\examples{
open_output_device('circos.png', width = 1200, height = 1200, res = 120)
open_output_device('circos.pdf', width = 1400, height = 800, res = 160)
}
