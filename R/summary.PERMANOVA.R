summary.PERMANOVA <- function(perm, latex=FALSE){
  cat(" ###### PERMANOVA : MANOVA basado en distancias #######\n")
  cat("________________________________________________\n")
  print(perm$call)
  cat("\n\nNumber of permutations:", perm$nperm, "\n\n")
  cat("Matriz de Contrastes:\n\n")
  print(perm$C)
  cat("\n\nContraste Global:\n")
  print(perm$Inicial$Global)
  cat("\n\nContrastes Individuales:\n")
  print(perm$Inicial$Contrastes)
  if (!is.null(perm$Inicial$Efectos))
  print(perm$Inicial$Efectos)

  if (latex){
    print(xtable(perm$C, caption="Matriz de Contrastes"))
    print(xtable(perm$Inicial$Global,  digits=c(0, 3, 3, 0, 0, 3, 6), caption="Contraste global"))
    print(xtable(perm$Inicial$Contrastes,  digits=c(0, 3, 3, 0, 0, 3, 6), caption="Contrastes individuales"))
    if (!is.null(perm$Inicial$Efectos))
      print(xtable(perm$Inicial$Efectos,  digits=c(0, 3, 3, 0, 0, 3, -3), caption="Contrastes para los efectos"))

    print(xtable(perm$Inercias, caption="Varianza explicada por las coordenadas principales de las medias"))
  }
}
