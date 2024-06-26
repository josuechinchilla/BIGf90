% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_blup.R
\name{run_blup}
\alias{run_blup}
\title{Run blupf90+}
\usage{
run_blup(path_2_execs)
}
\arguments{
\item{path_2_execs}{path to a folder that holds the renumf90 executable. This field should be in quotes "".}
}
\description{
This function runs blupf90+ using a pre-processed parameter file called renf90.par.
}
\details{
This function runs blupf90+ using a parameter file named renf90.par.
Since this function depends only on the renf90.par parameter file, the only input needed from the user is a path where the blupf90+ executable is located. A log file called run_blup.log is also produced.
}
\examples{
## Example

# run_blup(path_2_execs = "/Users/johndoe/Desktop/bf90_execs/")

}
