PERMANOVA <- function(Distance, grupo, C=NULL, Efectos=NULL, nperm = 1000, seed=NULL, CoordPrinc=FALSE, dimens=2, PCoA="Standard", ProjectInd=TRUE, tol=1e-4, DatosIni=TRUE) {

  D=Distance$D
  Coefficient=Distance$Coefficient
  cl <- match.call()
  
  PCoAs= c("Standard", "Weighted")
  if (is.numeric(PCoA)) PcoA=PCoAs(PCoA)
  
  if (!is.factor(grupo)) stop("The grouping variable must be a factor")
  
  if (!is.null(seed)) set.seed(seed)
  Perm = list() #Container for the solution
  Perm$call=cl
  # Setting the properties of data
  if (is.null(rownames(D)))
    rownames(D) <- rownames(D, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(D)
  
  Perm$Title = "MANOVA BASED ON  based DISTANCES"
  Perm$Type = "PERMANOVA"
  Perm$Distances = D
  
  if (is.factor(grupo)) {
    GroupNames = levels(grupo)
  }
  
  L = length(levels(grupo)) #Número de grupos
  I = dim(D)[1] #Número de de individuos
  X = Factor2Binary(grupo)
  colnames(X)=levels(grupo)
  rownames(X)=rownames(D)
  
  if (is.null(C)){
    C=diag(L)
    rownames(C)=paste("C",GroupNames)
    colnames(C)=GroupNames
    Cwasnull=TRUE
  }
  
  nc = dim(C)[1]
  
  Perm$C=C
  
  if (!is.null(Efectos)){
    Perm$Effects=Efectos
  }
  
  print("PerMANOVA Inicial")
  # PERMANOVA INICIAL
  Perm$Initial=PERMANOVA.Estimacion(D, X, C, Efectos)
  
  # Permutaciones
  Ftotal=matrix(0, 1,  nperm)
  FContrastes= matrix(0, nc, nperm)
  if (!is.null(Efectos)){
    ne=length(levels(Efectos))
    FEfectos=matrix(0, ne, nperm)
  }
 
  if (is.null(seed))
    set.seed(1)
  else
    set.seed(seed)
  
  for (i in 1:nperm){
    if (i%%100==0)
    print(i)
    muestra=sample.int(I)
    Man=PERMANOVA.Estimacion(D[muestra, muestra], X, C, Efectos)
    Ftotal[i]=Man$Global[5]
    FContrastes[,i]=Man$Contrastes[,5]
    if (!is.null(Efectos)){
      FEfectos[,i]=Man$Efectos[,5]
    }
  }
  
  Perm$DistMuestral=Ftotal
  Perm$pvalue=((sum(Ftotal >= Perm$Initial$Global[1,5]) +1 )/ (nperm+1))
  Perm$Initial$Global= cbind(Perm$Initial$Global, Perm$pvalue)
  
  colnames(Perm$Initial$Global)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")
  rownames(Perm$Initial$Global)="Total"
  
  Perm$DistMuestralC=FContrastes
  pvalC=(apply(FContrastes > matrix(Perm$Initial$Contrastes[,5], ncol=1) %*% matrix(1, 1, nperm), 1, sum)+1)/ (nperm+1)
  Perm$Initial$Contrastes=cbind(Perm$Initial$Contrastes, pvalC)
  colnames(Perm$Initial$Contrastes)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")
  
  if (!is.null(Efectos)){
    Perm$DistMuestralEffects=FEfectos
    pvalE=(apply(FEfectos > matrix(Perm$Initial$Efectos[,5], ncol=1) %*% matrix(1, 1, nperm), 1, sum)+1)/ (nperm+1)
    Perm$Initial$Efectos=cbind(Perm$Initial$Efectos, pvalE)
    colnames(Perm$Initial$Efectos)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")
  }

  
  Perm$CoordPrinc=CoordPrinc
  
  if (CoordPrinc){
   
    N=t(X) %*% X
    #D=0.5*D^2
    FS=solve(N)%*% (t(X) %*% (0.5*D^2) %*% X) %*% solve(N)
    f=matrix(diag(FS), L, 1)
    
    DB=2*FS - f %*% matrix(1, 1, L) - t(f %*% matrix(1, 1, L))
    DB=0.5*DB
    
    switch(PCoA, Standard = {
      H=(diag(L) - matrix(1, L, L)/L)
      B =  H %*% DB %*% H
    },Weighted = {
      H=(diag(L) - matrix(1, L, 1) %*% matrix(diag(N), 1, L)/I)
      B = H %*% DB  %*% H
    })
    
    solut <- svd(B)
    b=diag(B)
    vp=solut$d[1:dimens]
    Inertia=(solut$d[1:dimens]/sum(solut$d)) * 100
    Perm$ExplainedVariance = Inertia
    Perm$Inertias= cbind(vp, Inertia, cumsum(Inertia))
    rownames(Perm$Inertias)=paste("PCo", 1:dimens)
    colnames(Perm$Inertias)=c("Own value", "Explained Variance", "Cumulative")
    Y = solut$u %*% diag(sqrt(solut$d))
    d0=apply(Y^2,1, sum)
    
    rownames(Y)=levels(grupo)
    
    st <- apply(Y^2, 1, sum)
    Perm$MeanCoordinates=Y[,1:dimens]
    colnames(Perm$MeanCoordinates)=paste("Dim",1:dimens)
    colnames(Perm$MeanCoordinates)=paste("Dim",1:dimens)
    qlr <- diag(1/st) %*% (Perm$MeanCoordinates^2)
    Perm$Qualities=round(qlr[, 1:dimens]*100, digits=2)
    rownames(Perm$Qualities)=levels(grupo)
    colnames(Perm$Qualities)=paste("Dim",1:dim(Perm$Qualities)[2])
    
    Perm$CummulativeQualities=t(apply(Perm$Qualities,1,cumsum))
    
    if (ProjectInd){
      Y = Y[,1:dimens]
      if (DatosIni){
        YY=as.matrix(Distance$Data)
        Means= solve(N) %*% t(X) %*% YY
        Di=DistContinuous(Means, y=YY, coef = Distance$Coefficient)$D^2
        Yi=-0.5 * ginv(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,I,1) %*% matrix(d0,1,L))
        Yi=t(Yi)
      }
      else{
        G <- (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I) %*% (-0.5 * D^2) %*% (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I)
        soluc <- svd(G)
        dimefec=sum(soluc$d > tol)
        YY=soluc$u[,1:dimefec] %*% diag(sqrt(soluc$d[1:dimefec]))
        Means= solve(N) %*% t(X) %*% YY
        Di=DistContinuous(Means, y=YY, coef = 1)$D
        Yi=-0.5 * solve(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,I,1) %*% matrix(d0,1,L))
        Yi=t(Yi)
      }
      
      rownames(Yi)=rownames(D)
      Perm$RowCoordinates=Yi
      
    }
    
    Perm=AddCluster2Biplot(Perm, ClusterType="us", Groups=grupo )
  }
  class(Perm)=c("PERMANOVA")
  return(Perm)
}





