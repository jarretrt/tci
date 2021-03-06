% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkpd_assist.R
\name{create_intvl}
\alias{create_intvl}
\title{dosing schedule
Create dosing schedule}
\usage{
create_intvl(dose, inittm = 0)
}
\arguments{
\item{dose}{Data frame with columns "time" and "infrt".}

\item{inittm}{Starting time of initial infusion}
}
\description{
Create a dosing schedule object with columns "infrt", "begin", "end" from
vectors of infusions and infusion end times. The argument "inittm" is used
to specify the starting time of the first infusion.
}
