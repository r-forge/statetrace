\name{stBootav}
\alias{stBootav}
\alias{guistBootav}
\title{
Generate bootstrap averages for state-trace plots
}
\description{
  \code{stBootav} is used to obtain bootstrap participant averages, which are utilised in generating the state-trace plots. Two methods of calculating the bootstrap average are offered, where participants can either be treated as a fixed effect or as a sample from the population.
}
\usage{
guistBootav()
stBootav(staname = "", exclude = NULL, 
         eav = TRUE, tav = TRUE, mav = TRUE, 
         nbsamp = 10000, acc = TRUE, bootparticipants = FALSE)
}
\arguments{
  \item{staname}{
     character string indicating the name of the sta object.
}
  \item{exclude}{
     boolean vector indicating participants to be excluded.
}
  \item{eav}{
     logical argument specifying if the bootstrap average should be calculated for the encompassing model.
}
  \item{tav}{
     logical argument specifying if the bootstrap average should be calculated for the trace model.
}
  \item{mav}{
     logical argument specifying if the bootstrap average should be calculated for the monotonic model.
}
  \item{nbsamp}{
     number of bootstrap samples to draw from each model.
}
  \item{acc}{
     logical value specifying whether a probability measure should be used when calculating the bootstrap averages. When \code{TRUE}, the proportion correct is used for B0 designs and the Hit minus False Alarm rate is used for B2 designs. When \code{FALSE}, the corresponding inverse cumulative normal (\emph{z}) transformation is used; i.e., \emph{z}(proportion correct) for B0 designs and the signal detection \emph{d'} measure for B2 designs.
}
  \item{bootparticipants}{
     logical value indicating the method for calculating the bootstrap average. If \code{TRUE}, will treat the participants as a sample from a population and so re-samples them as well as re-sampling the posterior estimates. However, the user also has the option of calculating the bootstrap average by treating the participants as a fixed effect and so only does the latter re-sampling of the posterior estimates.
}
}
\details{
  Note that both a GUI (i.e., \code{guistBootav}) and non-GUI (i.e., \code{stBootav}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}. 
}
\value{
  \code{stBootav} will return an updated version of the sta object, with the bootstrap samples stored separately depending on the accuracy measure selected (i.e., \code{$all$p} or \code{$all$z}) and the model for which they belong to (\code{$e} for the encompassing model, \code{$t} for the trace model, \code{$m} for the monotonic model).
}
\references{
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
}
\seealso{
 \code{\link{stPlot}}, to generate the state-trace plot using the bootstrap averages.
}
\examples{
\dontrun{
  ## stBootav must be given an existing sta object which has completed at  
  ## least one pass of sampling.
  ## see stFirst for an example to produce an sta object that can 
  ## then be used to run the following:

  ## To calculate the bootstrap average treating the participants as a sample 
  ## from a population, use:
    
  stBootav(staname = "DFIE.sta", bootparticipants = FALSE)
}
}