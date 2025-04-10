% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npranov.R
\name{main_effect_boxplot}
\alias{main_effect_boxplot}
\title{Generate an Interactive Boxplot}
\usage{
main_effect_boxplot(data, factor, response, title)
}
\arguments{
\item{data}{Data frame containing the factors and response variable.}

\item{factor}{Factor variable for grouping.}

\item{response}{Numeric response variable.}

\item{title}{Title for the plot.}
}
\value{
An interactive boxplot using ggplot2 and plotly.
}
\description{
Generate an Interactive Boxplot
}
\examples{
main_effect_boxplot(data, "Location", "Midichlorians", "Main Effect of Location")
}
