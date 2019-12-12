\name{PostHocComp}
\alias{PostHocComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post Hoc pairwise comparisons}
\description{Creates a contrast C matrix of post hoc comparisons among groups.}
\usage{
PostHocComp(grupo)
}

\arguments{
  \item{grupo}{A factor containing the groups or treatments.
}
}
\details{
Creates a contrast C matrix of post hoc comparisons among groups. The rows of the contrast matrix are not orthogonal.
}
\value{The PostHocComp function return a matrix to compare the levels of groups introduced.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Laura Vicente González

José Luis Vicente Villardón
}

\examples{
data(wine)
PostHocComp(wine$Group)
}

