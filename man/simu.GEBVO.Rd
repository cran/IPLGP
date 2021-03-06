% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.GEBVO.R
\name{simu.GEBVO}
\alias{simu.GEBVO}
\title{GEBV-O Strategy}
\usage{
simu.GEBVO(
  phe.t,
  geno.t,
  marker,
  geno.c = NULL,
  npl = NULL,
  weight = NULL,
  direction = NULL,
  nprog = 50,
  nsele = NULL,
  ngen = 10,
  nrep = 30,
  console = TRUE
)
}
\arguments{
\item{phe.t}{matrix. An n*t matrix denotes the phenotypic values of the
training population with n individuals and t target traits.}

\item{geno.t}{matrix. An n*p matrix denotes the marker score matrix of the
training population. The markers must be coded as 1, 0, or -1 for alleles
AA, Aa, or aa. The missing value must have been already imputed.}

\item{marker}{matrix. An p*2 matrix whose first column indicates the chromosome
number to which a marker belongs; and second column indicates the position of
the marker in centi-Morgan (cM).}

\item{geno.c}{matrix. An nc*p matrix denotes the marker score matrix of the
candidate population with nc individuals and p markers. The markers must be
coded as 1, 0, or -1 for alleles AA, Aa, or aa. The missing value must have
been already imputed. If geno.c is set to be NULL, the candidate population
is exactly the training population.}

\item{npl}{vector. A vector denotes how many parental lines with the top GEBVs
will be chosen from each trait. If npl is set to be NULL, there will be 4
parental lines chosen from each trait.}

\item{weight}{vector. A vector with length t indicates the weights of target
traits in selection index. If weight is set to be NULL, the equal weight will
be assigned to all the target traits.}

\item{direction}{vector. A vector with length t indicates the selecting
directions for target traits. The elements of direction are 1, or -1
representing the rule that the larger the better; or the smaller the better.
If direction is set to be NULL, the selecting direction will be the same as
weight.}

\item{nprog}{integer. An integer indicates the number of progenies which
will be produced for each of the best individuals at every generation.}

\item{nsele}{integer. An integer indicates the number of the best individuals
which will be selected at each generation. If nsele is set to be NULL, the
number will be the same as the number of F1 individuals.}

\item{ngen}{integer. An integer indicates the number of generations in the
simulation process.}

\item{nrep}{integer. An integer indicates the number of repetitions in the
simulation process.}

\item{console}{logical. A logical variable, if console is set to be TRUE,
the simulation process will be shown in the R console.}
}
\value{
\item{method}{The GEBV-O strategy.}
\item{weight}{The weights of target traits in selection index.}
\item{direction}{The selecting directions of target traits in selection index.}
\item{mu}{The mean vector of target traits.}
\item{sd}{The standard deviation vector of target traits.}
\item{GEBV.value}{The GEBVs of target traits in each generation and each
repetition.}
\item{parental.lines}{The IDs and D-score of parental lines selected in
each repetition.}
\item{suggested.subset}{The most frequently selected parental lines by this
strategy.}
}
\description{
Identify parental lines based on GEBV-O strategy and simulate their offsprings.
}
\note{
The function output.best and output.gain can be used to summarize the result.
}
\examples{
# generate simulated data
set.seed(2000)
phe.test <- data.frame(trait1 = rnorm(10,30,10), trait2 = rnorm(10,10,5))
geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
geno.candidate <- matrix(sample(c(1,-1), 300, replace = TRUE), 15, 20)

# run and output
result <- simu.GEBVO(phe.test, geno.test, marker.test, geno.candidate,
nprog = 5, nsele = 10, ngen = 5, nrep = 5)
result$suggested.subset
}
\references{
Chung PY, Liao CT. 2020. Identification of superior parental lines for
biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
}
\seealso{
\code{\link[IPLGP]{phe.sd}}
\code{\link[IPLGP]{GBLUP.fit}}
\code{\link[IPLGP]{GA.Dscore}}
\code{\link[IPLGP]{simu.gamete}}
\code{\link[IPLGP]{simu.GDO}}
\code{\link[IPLGP]{simu.GEBVGD}}
\code{\link[IPLGP]{output.best}}
\code{\link[IPLGP]{output.gain}}
}
