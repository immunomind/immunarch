% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmers.R
\name{getKmers}
\alias{getKmers}
\alias{get.kmers}
\alias{makeKmerTable}
\title{Calculate the k-mer statistics of immune repertoires}
\usage{
getKmers(.data, .k, .col = c("aa", "nt"), .coding = TRUE)
}
\arguments{
\item{.data}{The data to be processed. Can be \link{data.frame},
\link{data.table}, or a list of these objects.

Every object must have columns in the immunarch compatible format.
\link{immunarch_data_format}

Competent users may provide advanced data representations:
DBI database connections, Apache Spark DataFrame from \link{copy_to} or a list
of these objects. They are supported with the same limitations as basic objects.

Note: each connection must represent a separate repertoire.}

\item{.k}{Integer. Length of k-mers.}

\item{.col}{Character. Which column to use, pass "aa" (by default) for CDR3 amino acid sequence,
pass "nt" for CDR3 nucleotide sequences.}

\item{.coding}{Logical. If TRUE (by default) then removes all non-coding sequences from input data first.}
}
\value{
Data frame with two columns (k-mers and their counts).
}
\description{
Calculate the k-mer statistics of immune repertoires
}
\examples{
data(immdata)
kmers <- getKmers(immdata$data[[1]], 5)
kmers \%>\% vis()
}
\concept{k-mers}
