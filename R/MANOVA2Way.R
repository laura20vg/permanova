MANOVA2Way <- function(Y, f1, f2, Interacyion=TRUE, M=NULL, InitialTransform = 5, AddOnes=FALSE){
  p1=length(levels(f1))
  c1=t(contr.helmert(p1))
  p2=length(levels(f2))
  c1=t(contr.helmert(p2))
  ng=p1*p2
  C=matrix(0,ng-1, ng)
  k=0
  for (i in 1:p1)
    for (j in 1:p2){
      k=k+1
    }

}



