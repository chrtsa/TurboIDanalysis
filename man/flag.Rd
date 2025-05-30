% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{flag}
\alias{flag}
\title{True Positive and False Positive Rates}
\usage{
flag(data, sample, TP, FP = "FP", gene_name = "gene_name", ratio, TP_flag)
}
\arguments{
\item{data}{The dataset that will be used.}

\item{sample}{The name of the sample (e.g. SYP or NES)}

\item{TP}{A list or dataframe containing the true positive proteins, used for flagging. Contains proteins that should be biotinylated by TurboID. Ensure that the name of the column in the TP list matches the dataset.}

\item{FP}{A list or dataframe containing the false positive proteins, used for flagging. Contains proteins that do not belong in the cellular compartment being analysed. Ensure that the name of the column in the FP list matches the dataset.}

\item{gene_name}{The name of the column in data that contains the gene or protein names. Make sure that the name of the column in the TP and FP list matches the dataset}

\item{ratio}{The name of the column in data that contains the ratio of interest}

\item{TP_flag}{The name of the column to be created or modified in data for the true positive flag, so whether this protein is a TP for SYP or NES respectively}
}
\value{
The initial dataset, including TPR and FPR
}
\description{
This function calculates the True Positive Rate (TPR) and False Positive Rate (FPR) for the provided dataset. First, the proteomics hits are flagged as TP or FP, based on the user-provided TP and FP lists, generated based on prior annotation and literature evidence. The proteins are then normalised against the distribution of FP proteins, and ranked in descending order. Then, the TPR and FPR are calculated at each potential cutoff. The optimal cutoff will be calculated based on the ROC at a later stage.
}
