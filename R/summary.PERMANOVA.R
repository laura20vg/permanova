summary.PERMANOVA <- function(x, Latex=FALSE){
  cat(" ###### PERMANOVA Analysis #######\n\n")
  cat("Call\n")
  print(x$call)
  cat("________________________________________________\n\n")
  cat("Contrast Matrix\n")
  print(x$C)
  cat("________________________________________________\n\n")
  if (!is.null(x$Effects)){
    cat("Effects vector\n")
    print(x$Effects)
    cat("________________________________________________\n\n")
  }
  cat(paste("PerMANOVA\n"))
  print(x$Initial$Global)
  cat("________________________________________________\n\n")
  cat(paste("Contrasts\n"))
  print(rbind(x$Initial$Contrastes, x$Initial$Global))
  cat("________________________________________________\n\n")
  if (!is.null(x$Effects)){
    cat("Effects \n")
    print(rbind(x$Initial$Efectos, x$Initial$Global))
    cat("________________________________________________\n\n")
  }
  if (Latex){
    xtable(round(x$C, digits=0), caption="Contrast Matrix")
    cat("________________________________________________\n\n")
    xtable(x$Initial$Global, caption="PerMANOVA")
    cat("________________________________________________\n\n")
    xtable(rbind(x$Initial$Contrastes, x$Initial$Global), caption="perMANOVA with contrasts")
    
  }
}

