MANOVA.Estimacion <- function(Y, X, C, M) {

  Man=list()
  L = dim(X)[2] #Número de grupos
  I = dim(Y)[1] #Número de de individuos
  J = dim(Y)[2] #Número de variables
  r = min(c(L - 1, J)) # Rango de la solucion
   C=matrix(C, ncol=L)

  B=ginv(t(X) %*% X) %*% t(X) %*% Y

  O = C %*% B %*% M

  R=C %*% ginv(t(X) %*% X) %*% t(C)

  E = t(M) %*% t(Y) %*% (diag(I)- (X %*% ginv(t(X) %*% X) %*% t(X))) %*% Y %*% M
  H = t(O) %*% ginv(R) %*% O

  E2=matrixsqrtinv(E)

  YB2= E2 %*% H %*% E2

  dvs=svd(YB2, 0)
  VP=dvs$d[1:r]
  Inercia=(VP/sum(VP))*100
  V=dvs$v
  A=E2 %*% V

  Man$B=B
  Man$E=E
  Man$H=H
  Man$VP=VP
  Man$Inercia=Inercia
  Man$V=V
  Man$A=A

  Wilks = det(diag(1/(1 + VP)))
  glh = rankMatrix(H)[1]
  J=rankMatrix(H+E)[1]

  gle = I - L
  t = ((glh^2 * J^2 - 4)/(J^2 + glh^2 - 5))^0.5
  w = gle + glh - 0.5 * (J + glh + 1)
  df1w = J * glh
  df2w = w * t - 0.5 * (J * glh - 2)
  s=min(c(J, glh))
  Wilksf = ((1 - Wilks^(1/t))/(Wilks^(1/t))) * (df2w/df1w)
  Wilksp = 1 - pf(Wilksf, df1w, df2w)
  Wilsr=c(Wilks, Wilksf, df1w, df2w,  Wilksp)
  # Esto es para el caso en el que no haya matriz de contrastes. Hay que ver la adaptación

  Pillai=sum(VP/(1 + VP))
  b=max(c(J,glh))
  m=(abs(J-glh)-1)/2
  n=(gle-J-1)/2
  Pillaif= ((2*n+s+1)*Pillai)/((2*m+s+1)*(s-Pillai))
  df1p=(2*m+s+1)*s
  df2p=(2*n+s+1)*s
  Pillaip = 1 - pf(Pillaif, df1p, df2p)

  Pillair=c(Pillai, Pillaif, df1p, df2p,  Pillaip)

  Hotteling=sum(VP)
  a=J*glh
  a=s*(2*m+s+1)
  Hottelingf=(Hotteling)*((2*(s*n+1))/(s^2*(2*m+s+1)))
  df1h=s*(2*m+s+1)
  df2h=2*(s*n+1)
  Hottelingp = 1 - pf(Hottelingf, df1h, df2h)
  Hottelingr=c(Hotteling, Hottelingf, df1h, df2h,  Hottelingp)


  Roy=VP[1]
  v1= (abs(J-glh)-1)/2
  v2= (gle-J-1)/2
  Royf=Roy*(2*v2+2)/(2*v1+2)
  df1r=(2*v1+2)
  df2r=(2*v2+2)
  Royp = 1 - pf(Royf, df1p, df2p)
  Royr=c(Roy, Royf, df1r, df2r,  Royp)

  ManovaR=rbind(Wilsr, Pillair, Hottelingr, Royr)
  rownames(ManovaR) =c("Wilk's lambda", "Pillai Trace", "Lawley-Hotteling", "Roy's greatest root")
  colnames(ManovaR)=c("Estadístico", "F aprox", "GL. Num", "GL. Den", "p-valor")
  Man$MANOVA=ManovaR
  return(Man)

}
