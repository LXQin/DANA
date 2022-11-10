#' Chromosomal coordinate for miRNA (hsa) from MiRBase v22
#'
#' Chromosomal coordinates of Homo sapiens microRNAs. Note that all
#' lower-case notation is used for miRNA names.
#' \itemize{
#'   \item microRNAs: miRBase v22
#'   \item genome-build-id: GRCh38
#'   \item genome-build-accession: NCBI_Assembly:GCA_000001405.15
#' }
#' Hairpin precursor sequences have type "miRNA_primary_transcript".
#' Note, these sequences do not represent the full primary transcript,
#' rather a predicted stem-loop portion that includes the precursor
#' miRNA. Mature sequences have type "miRNA".
#'
#' @format 'miRNA.coordinates'
#' A data frame with 3,619 rows and 7 columns:
#' \describe{
#'   \item{name}{miRNA name in \emph{all lower-case}
#'     (e.g. "hsa-let-7b", "hsa-mir-100", ...)}
#'   \item{chromosome}{The chromosome the gene is located on.}
#'   \item{position}{Nucleotide coordinate on the chromosome.}
#'   \item{type}{Either "miRNA_primary_transcript" for heirpin precursor
#'     sequences or "miRNA" for mature sequences.}
#'   \item{start.position, end.position}{Starting and ending nucleotide
#'     coordinate on the chromosome.}
#'   \item{id}{Corresponding id (e.g. "MI0000063", "MI0000102")}
#'   ...
#' }
"miRNA.coordinates"
