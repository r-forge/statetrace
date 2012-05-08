\name{DFIE}
\alias{DFIE}
\docType{data}
\title{
Disproportionate Face Inversion Effect choice accuracy data
}
\description{
This data concerns recognition memory judgements from an experiment that aimed to examine evidence for the Disproportionate Face Inversion Effect (DFIE); that is, the finding that memory performance is disproportionately affected by inversion for faces than non-face stimuli. The DFIE is one of the primary pieces of evidence to suggest that faces are encoded on an extra latent dimension (commonly called configural encoding), which is not available to non-face stimuli. Note this is also the main example 2AFC data used by Prince, Hawkins, Love and Heathcote (2011) and the vignette provided in this package. 
}
\usage{
data(DFIE)
}
\details{
A data frame containing one row per design cell per participant
\describe{
  \item{\code{P    }}{Participant identifier; 1:18.}
  \item{\code{S    }}{State factor levels; 'F' = Face stimuli, 'H' = House stimuli.}
  \item{\code{D    }}{Dimension factor levels; 'U' = Upright study presentation, 'I' = Inverted study presentation.}
  \item{\code{T    }}{Trace factor levels; three levels of study presentation time, for upright (66, 100, 200ms) and for inverted (200, 600, 1800ms).}
  \item{\code{n    }}{number of successes per design cell; i.e., number of correct identifications.}
  \item{\code{N    }}{total number of observations per design cell.}
}
In this experiment, images of unfamiliar faces and houses were studied either upright or inverted for varying presentation times and then tested in the same "matched" orientation using a two alternative forced choice recognition task (i.e., participants were shown two test images simultaneously, one previously studied and one new, and were asked to indicate which of the two images had been presented at study). Participants attended three 1 hour sessions over which they completed 52 blocks of 18 study images and 18 test pairs, yielding 78 observations per design cell.\cr
\cr
Note this data set has been presented in a "summary" format, which has a row for each design cell containing the summed frequencies. However, \code{StateTrace} is very flexible in it's data input capabilities, and so the data may alternatively be presented in a "trial" format, which has one row of data for each trial. In contrast to the summary format, trial files have a single column after the trace factor levels that contains a binary response indicator (e.g., correct and error responses which are coded as 1 and 0 respectively). For a full description of the data formats suitable for \code{StateTrace} see Prince et al., (2011) or the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}.. 
}
\references{
Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr 
}
\examples{
## load DFIE data set
data(DFIE)
## Average proportion correct
tapply(DFIE$n,list(DFIE$T,DFIE$D,DFIE$S),mean)
}
      