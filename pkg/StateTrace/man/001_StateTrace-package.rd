\name{StateTrace-package}
\alias{StateTrace-package}
\alias{StateTrace}
\docType{package}
\title{
Bayesian Ordinal Analysis for State-Trace Data
}
\description{
State-trace analysis (Bamber, 1979) is a graphical method for determining whether one or more than one latent (i.e., not observable) variable mediates an apparent dissociation between the effect of two experimental manipulations. State-trace analysis makes only ordinal assumptions and so is not confounded by range effects that plague alternative methods, especially when performance is measured on a bounded response scale (e.g., accuracy). \code{StateTrace} automates many aspects of a state-trace analysis of accuracy and other binary response data, including implementing Bayesian methods quantifying evidence about the outcomes of a state-trace experiment and the creation of customisable graphs.
}
\details{
\tabular{ll}{
Package: \tab StateTrace\cr
Type: \tab Package\cr
Version: \tab 1.0-4\cr
Date: \tab 2012-05-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package offers users the option of using either a non-GUI (i.e., command line input) or GUI version for the majority of the functions. Although examples for the non-GUI version are provided in this help documentation for each function, users should see Prince, Hawkins, Love and Heathcote (submitted) for an abridged description of the GUI; alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace", package="StateTrace")}.
}
\seealso{
\code{\link{guista}}, for the main user interface for the \code{StateTrace} package.\cr
To generate an sta object see \code{\link{staMake}} or alternatively \code{\link{stFirst}}, which will also run the first pass of sampling. See \code{\link{stSample}} to either begin the sampling process or further refine the sampling and posterior estimates. Results can then be examined using \code{\link{stSummary}} and \code{\link{stProbplot}}. See \code{\link{stBootav}} to obtain bootstrap participant averages, which are utilised in generating the state-trace plots in \code{\link{stPlot}}. See \code{\link{staManage}} for joining sta objects and managing the stored samples.\cr
\code{\link{DFIE}}, provides example data suitable for \code{StateTrace}.
}
\author{
Melissa Prince, Guy Hawkins, Jonathon Love and Andrew Heathcote\cr
Maintainer: David Elliott <David.Elliott@newcastle.edu.au>
}
\references{
Bamber, D. (1979). State-trace analysis: A method of testing simple theories of causation. \emph{Journal of Mathematical Psychology, 19,} 137-181.\cr
Prince, M., Brown, S., & Heathcote, A. (2011). The design and analysis of state-trace experiments. \emph{Psychological Methods}.\cr
Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
\cr
\cr
}
