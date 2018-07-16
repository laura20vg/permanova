\name{MANOVA}
\alias{MANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Analysis of Variance (MANOVA)
}
\description{
Performs a Multivariate Analysis of Variance (MANOVA).
}
\usage{
MANOVA(Y, grupo, C = NULL, M = NULL, Efectos = NULL, InitialTransform = 5, AddOnes = FALSE, Contrastes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A matrix containing the response variables
}
  \item{Group}{
A factor containing the groups to compare (or the trantments)
}
  \item{C}{
A matrix of contrasts
}
  \item{M}{
A matrix with the linear combinations of the variables
}
  \item{Effects}{
A factor with the effects in the rows of C
}
  \item{InitialTransform}{
Initial transformation of the response variables
}
  \item{AddOnes}{
Add a column of ones to the design matrix. By default is false.
}
  \item{Contrasts}{
Shoud each contrast be tested separately?. By default is TRUE.
}
}
\details{
Performs a general MANOVA to compare several groups or treatments. Aditional contrasts can be tested using the contrasts matrix \eqn{\bf{C}}. This can also be separated into several effects.
}
\value{
An object of class "MANOVA"  and "Canonical.Biplot"
}
\references{
Amaro, I. R., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Gabriel KR (1971) The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3):453-467.

Gabriel, K. R. (1995). MANOVA biplots for two-way contingency tables. WJ Krzanowski (Ed.), Recent advances in descriptive multivariate analysis, Oxford University Press, Toronto. 227-268.

Gower y Hand (1996): Biplots. Chapman & Hall.

Varas, M. J., Vicente-Tavera, S., Molina, E., & Vicente-Villardon, J. L. (2005). Role of canonical biplot method in the study of building stones: an example from Spanish monumental heritage. Environmetrics, 16(4), 405-419.

Santana, M. A., Romay, G., Matehus, J., Vcente-Villardon, J. L., & Demey, J. R. (2009). simple and low-cost strategy for micropropagation of cassava (Manihot esculenta Crantz). African Journal of Biotechnology, 8(16).

}
\author{
Laura Vicente González
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(wine)
Y=as.matrix(wine[,4:21])
group=wine$Group
manvin=MANOVA(Y, Group=group)
plot(manvin, PlotClus=TRUE, PlotCircle=FALSE, PlotVars=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
