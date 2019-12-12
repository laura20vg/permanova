\name{PERMANOVA}
\alias{PERMANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PERMANOVA: MANOVA based on distances
}
\description{
The correct application of MANOVA needs normal and homocedastic data and the number of variables be much smaller than the number of individuals, but for many applications the conditions do not hold. To extend the application to this data Anderson develops PERMANOVA. This non-parametric test based on distances uses permutation to approximate the sampling distribution of the test statistic.

Contrasts and Effects can be added to the calculations.
}
\usage{
PERMANOVA(Distance, grupo, C = NULL, Efectos = NULL, nperm = 1000, seed = NULL, CoordPrinc=FALSE, dimens = 2, PCoA = "Standard", ProjectInd = TRUE, tol = 1e-04, DatosIni = TRUE)
}

\arguments{
  \item{Distance}{
  A list of three elements containing the data, the distances between individuals and type of distance used.
}
  \item{grupo}{
  A factor containing the groups to compare.
}
  \item{C}{
  Contrast matrix. By default it is null and the identity is used.
}
  \item{Efectos}{
  A factor with the effects in the rows of C. By default it is null and each row is considered as an effect.
}
  \item{nperm}{
  Number of permutations to perform. By default is 1000.
}
  \item{seed}{
  Seed to start permutations. By default is null.
}

  \item{CoordPrinc}{
  Should the principal coordinates be calculated?. By default is FALSE.
}
  \item{dimens}{
  Number of dimensions to choose in the PCoA. By default is 2.
}
  \item{PCoA}{
  The type of Principal Coordinates Analysis. 
  
  There are two possibilities:
  
   * Standard
   
   * Weighted
   
  By default is Standard
}
  \item{ProjectInd}{
  Show the row coordinates. By default is TRUE.
}
  \item{tol}{
  Tolerance
}
  \item{DatosIni}{
  It contains the initial data. By default is TRUE.
}
}
\details{
The function performs a PERMANOVA Analysis.
}
\value{
The PERMANOVA function create a list that return:
\item{call}{Function}
\item{Title}{MANOVA BASADO EN PERMUTACIONES}
\item{Type}{PERMANOVA}
\item{Distances}{A matrix containing the distances between individuals.}
\item{C}{Contrasts Matrix.}
\item{Initial}{ Containing two matrices:

    * Global -> Global contrast.
  
    * Contrastes ->Contrar for groups.
}
\item{DistMuestral}{Sample distribution of F-exp from permutations.}
\item{pvalue}{Estimate p-valor for PERMANOVA.}
\item{ExplainedVariance}{Explained variance by Principal Coordinates selected.}
\item{Inertias}{Own value, Explained variance, Cumulative explained variance.}
\item{MeanCoordinates}{Mean Coordinates by groups for the dimensions obtained in the Principal Coordinates Analysis.}
\item{Qualities}{Qualities representation by groups for the dimensions of PCoA.}
\item{CummulativeQualities}{Cummulative qualities representation.}
\item{ClusterType}{Cluster type selected.}
\item{Clusters}{Clusters created.}
\item{ClusterNames}{Names of clusters}
\item{ClusterColors}{Colors of clusters, color name and HTML code.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46.

Anderson, M. J. (2005). Permanova: a fortran computer program for permutational multivariate analysis of variance. Department of Statistics, University of Auckland, New Zealand, 24.
% hay que colocar a
}
\author{
Laura Vicente González and José Luis Vicente Villardón

}


\examples{
data(wine)
X = wine[,4:21]
D = DistanciasContinuas (X)
perwine=PERMANOVA(D, wine$Group)
perwine

% Separation of main and interaction effects. We have 4 groups (RD86, RD87, TO86, TO87)
% We have three independent contrasts to separate the main and interaction effects
% The contrast to compare the two origins RD and TO is (1, 1, -1, 1)
% The contrast to compare the two years 86 and 87 is (1, -1, 1, 1)
% The interaction is obtained from the product of both (1, -1, -1, 1)
% The three contrasts are orthogonal and  togeteher form the matrix C
C = matrix(c(1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1), nrow=3, byrow=TRUE)
rownames(C)=c("C1", "C2", "C3")
colnames(C)=levels(wine$Group)
% The sets of contrasts allows for the separation of the explaines sum of squares
% into three ortogonal parts corresponding to the main effects of origin and year
% and the interaction. The vector for effects is
effects=factor(c(1,2,3))
levels(effects)=c("Origin", "Year", "Interaction")
perwine2=PERMANOVA(D, wine$Group, C=C, Efectos=effects)
summary(perwine2)
}

