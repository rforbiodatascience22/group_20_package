#' Title Creation of a sequence of RNA nucleotides from DNA nucleotides
#'
#' Substitutes the T for a U when converting from a DNA sequence to an RNA sequence.
#' @param the_dna_sequence indicating DNA strand as input to the function.
#'
#' @return Return character strings representing the RNA
#' @author Oriade Latifah Simpson, \email{s172084@@dtu.dk}
#' @export
#'
#' @examples
#' my_rna <- make_rna(my_dna)
#' my_rna
make_rna <- function(the_dna_sequence){
  the_rna_sequence <- gsub("T", "U", the_dna_sequence)
  return(the_rna_sequence)
}
