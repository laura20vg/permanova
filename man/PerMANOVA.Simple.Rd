\name{PerMANOVA.Simple}
\alias{PerMANOVA.Simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PERMANOVA from distancies matrix
%%  ~~function to do ... ~~
}
\description{This function makes a PERMANOVA from distancies matrix and factor for groups.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PerMANOVA.Simple(D, grupo, nperm = 999, seed = NULL, C = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
  A matrix containing the distances between individuals.
%%     ~~Describe \code{D} here~~
}
  \item{grupo}{
  A factor containing the groups to compare.
%%     ~~Describe \code{grupo} here~~
}
  \item{nperm}{
  Number of permutation that want to perform. By default it is 999.
%%     ~~Describe \code{nperm} here~~
}
  \item{seed}{
  Seed to start permutations. By default it is null.
%%     ~~Describe \code{seed} here~~
}
  \item{C}{
  Contrast matrix. By default it is null.
%%     ~~Describe \code{C} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The PERMANOVA.Simple function create a list that return:
\item{call}{Function}
\item{nperm}{Number of permutation.}
\item{Inicial}{Containing:

* TSS -> Total sum of squares.

* BSS -> Between groups sum of squares.

* WSS -> Within groups sum of squares.

* glt -> Total degrees of freedom.

* glb -> Between groups degrees of freedom.

* glw -> Within groups degrees of freedom.

* Fexp -> Experimental F.
}
\item{Fvals}{F values of the permutation.}
\item{pval}{Estimate p-valor for PERMANOVA.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46.

Anderson, M. J. (2005). Permanova: a fortran computer program for permutational multivariate analysis of variance. Department of Statistics, University of Auckland, New Zealand, 24.
%% ~put references to the literature/web site here ~
}
\author{
Laura Vicente González

José Luis Vicente Villardón
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(wine)
X = wine[,4:21]
Dist = DistanciasContinuas (X)
PerMANOVA.Simple(Dist$D, wine$Group)
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
