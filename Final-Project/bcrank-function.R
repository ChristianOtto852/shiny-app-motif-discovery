if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("BCRANK")
library(BCRANK)
library(Biostrings)


Bcrank_motif_finder <- function(DNA, k) {
  dna <- DNAStringSet(unlist(DNA))
  fasta <- tempfile(fileext = ".fa")
  writeXStringSet(dna, fasta)
  
  bcrank_results <- bcrank(fafile = fasta,
                   length = k,
                   restarts = 20)
  top_motif <- toptable(bcrank_results, 1)
  
  pwm_mat <- pwm(top_motif, normalize = TRUE)
  
  pseudo_prob <- .01
  
  pwm_mat_with_pseudo <- (pwm_mat + pseudo_prob) / (1 + pseudo_prob * 4)
  
  pwm_width <- ncol(pwm_mat_with_pseudo)
  
  best_motifs <- c()
  
  for(seq in unlist(DNA)){
    best_prob <- -Inf
    for (i in 1:(nchar(seq) - pwm_width + 1)) {
      kmer <- substr(seq, i, i + pwm_width - 1)
      prob <- 0
      for (j in 1:(nchar(kmer))) {
        letter <- substr(kmer, j, j)
        prob <- prob + log(pwm_mat_with_pseudo[letter, j])
      }
      if (prob > best_prob) {
        best_prob <- prob
        best_motif <- kmer
      }
    }
    best_motifs <- c(best_motifs, best_motif)
  }
  return(best_motifs)
}