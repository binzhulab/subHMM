
\name{subHMM-package}
\alias{subHMM}
\docType{package}
\title{
 The SubHMM Model for Identifying Tumor Subclones
}
\description{
 A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
}
\details{
Allele-specific copy number alteration (ASCNA) analysis is for identifying copy
number abnormalities in tumor cells. Unlike normal cells, tumor cells are 
heterogeneous with a combination of dominant and minor subclones with distinct copy number
profiles. Estimating the clonal proportion and identifying main and subclone 
genotypes across the genome is important for understanding tumor progression. Several
ASCNA tools have recently been developed, but they have been limited to the 
identification of subclone regions, and not the genotype of the subclones. 
This package uses a hidden Markov model-based approach that estimates both 
sub-clone region as well as region-specific subclone genotype and clonal proportion. 
A hidden state variable is specified which represents the conglomeration of clonal genotypes
and subclone status. A two-step algorithm for parameter estimation is implemented,
where in the first step, a standard hidden Markov model with this conglomerated
state space is fit. 
Then, in the second step, region-specific estimates of the clonal
proportions are obtained by maximizing region-specific pseudo-likelihoods.
}
\references{
  Choo-Wosoba H., Albert P.S., Zhu B.
  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
}
\author{Hyoyoung Choo-Wosoba, Paul S. Albert, and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}