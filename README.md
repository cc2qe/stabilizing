# Stabilizing selection forward simulations

Scripts, data, and figures for forward evolutionary simulations
from Zhang et al. manuscript.

## Forward simulations
The `r06.1` directories contain the simulations included in the paper.
These apply a stabilizing selection fitness model and measure LD as well
as snp-pair effect correlation. Scripts adapted from an early implementation
credited to Arun Durvasula.

## emeraLD build
I used Corbin Quick's [emeraLD](https://github.com/statgen/emeraLD) software to
compute LD. I am using a [custom build](https://github.com/cc2qe/emeraLD/tree/dfc414e4f3026aab42e3e7b212e2d502a45de4fc)
of that includes a bugfix and the ability to receive a list of SNPs to
interrogate.

## Citation
Zhang et al. (2023).
Pervasive correlations between causal disease effects of proximal SNPs vary with functional annotations and implicate stabilizing selection.
In medRxiv (p. 2023.12.04.23299391). https://doi.org/10.1101/2023.12.04.23299391