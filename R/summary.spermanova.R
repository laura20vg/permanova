summary.spermanova <- function(perm, latex=FALSE){
  cat(" ###### PERMANOVA : MANOVA basado en distancias #######\n")
  cat("________________________________________________\n")
  print(perm$call)
  cat("Number of permutations:", perm$nperm, "\n\n")
  SC<- c(round(perm$Inicial$BSS,3), round(perm$Inicial$WSS,2), round(perm$Inicial$TSS,3))
  grados <- c(perm$Inicial$glb, perm$Inicial$glw, perm$Inicial$glt)
  Fexp <- c(round(perm$Inicial$Fexp,3), "", "")
  pvalor <- c(round(perm$pval,10), "", "")
  res <- data.frame(SC, grados, Fexp, pvalor, row.names = c("Entre grupos", "Dentro de los grupos", "TOTAL"))
  print(res)

  if (latex){
    print(xtable(res, caption="Resultados del PERMANOVA"))
  }
}



