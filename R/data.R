#' The Standard Code
#'
#' Translation Table 1 from Genetic Codes
#'
#' @docType data
#'
#' @usage data(the_codon_table)
#'
#' @format A dataframe
#'
#' @keywords datasets
#'
#' @references (Takahashi et al. (2005). Genomics 2005 Mar;85(3):360-71.doi: 10.1016/j.ygeno.2004.11.012.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/15718103/}{PubMed})
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1}{The Genetic Codes}
#'
#' @examples
#' data(the_codon_table)
#' amino_acid_sq <- paste0(the_codon_table[codons], collapse = "")
the_codon_table <- data.frame("UUU" = "F", "UCU" = "S", "UAU" = "Y", "UGU" = "C",
                              "UUC" = "F", "UCC" = "S", "UAC" = "Y", "UGC" = "C",
                              "UUA" = "L", "UCA" = "S", "UAA" = "_", "UGA" = "_",
                              "UUG" = "L", "UCG" = "S", "UAG" = "_", "UGG" = "W",
                              "CUU" = "L", "CCU" = "P", "CAU" = "H", "CGU" = "R",
                              "CUC" = "L", "CCC" = "P", "CAC" = "H", "CGC" = "R",
                              "CUA" = "L", "CCA" = "P", "CAA" = "Q", "CGA" = "R",
                              "CUG" = "L", "CCG" = "P", "CAG" = "Q", "CGG" = "R",
                              "AUU" = "I", "ACU" = "T", "AAU" = "N", "AGU" = "S",
                              "AUC" = "I", "ACC" = "T", "AAC" = "N", "AGC" = "S",
                              "AUA" = "I", "ACA" = "T", "AAA" = "K", "AGA" = "R",
                              "AUG" = "M", "ACG" = "T", "AAG" = "K", "AGG" = "R",
                              "GUU" = "V", "GCU" = "A", "GAU" = "D", "GGU" = "G",
                              "GUC" = "V", "GCC" = "A", "GAC" = "D", "GGC" = "G",
                              "GUA" = "V", "GCA" = "A", "GAA" = "E", "GGA" = "G",
                              "GUG" = "V", "GCG" = "A", "GAG" = "E", "GGG" = "G")
