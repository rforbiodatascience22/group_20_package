# Create a polypeptide sequence from the codons based on the dataset.
#' The creation of a polypeptide sequence from the codon.
#'
#' Creates a polypeptide sequence of amino acids based on the codon.
#' @param codons
#'
#' @return Return character strings representing the polypeptide
#' @author Oriade Latifah Simpson, \email{s172084@@dtu.dk}
#' @export
#'
#' @examples
#' polypeptide_sequence <- make_polypeptide(c("AUG", "CAA", "CGA", "UCA"))
#' polypeptide_sequence
make_polypeptide <- function(codons){
  amino_acid <- paste0(the_codon_table[codons], collapse = "")
  return(amino_acid)
}
