summary.BootDisMANOVA <- function(x, Latex=TRUE){
  cat(" ###### Bootstrap Distance Based MANOVA Analysis #######\n\n")
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
  cat(paste("MANOVA\n"))
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
    xtable(x$Initial$Global)
    xtable(rbind(x$Initial$Contrastes, x$Initial$Global))
    xtable(rbind(x$Initial$Efectos, x$Initial$Global))
  }
}