\name{stSummary}
\alias{stSummary}
\alias{guistSummary}
\title{
Summary model selection results
}
\description{
  Summary results can be extracted from an existing sta object using \code{stSummary}, which controls output of tabular model selection results to the R console. \code{stSummary} also provides information about the status of sampling for an sta object, as well as a range of additional results including prior probabilities and raw sample counts, which can be used to instantiate different model selection strategies without obtaining new samples.
}
\usage{
guistSummary()
stSummary(bosname = "", BF = FALSE, traceTrue = FALSE, eachs = TRUE,  
          nsig = 3, exclude = NULL, sorts = "", extras = "", 
          guiarg = NULL)
}
\arguments{
  \item{bosname}{
     character string indicating the name of the sta object from which to extract the summary results.
}
  \item{BF}{
       logical value indicating whether the model selection results should be reported as Bayes factors (\code{TRUE}) or posterior model probabilities (\code{FALSE}).
}
  \item{traceTrue}{
        logical value (conditional on \code{BF = FALSE}) specifying if probabilities should be calculated based on all four models (i.e., an exhaustive selection strategy) or by excluding the non-trace model from the set (i.e., the trace-true strategy).
}
  \item{eachs}{
     logical value indicating if the results should be reported for individual participants.
}
  \item{nsig}{
     numeric value specifying the numerical rounding of the group (and individual participant) results.
}
  \item{exclude}{
     boolean vector indicating participants to be excluded.
}
  \item{sorts}{
     character string specifying the sort order for individual participant results (i.e., \code{"Non-Trace model"}, \code{"No-Overlap model"}, \code{"Uni-dimensional model"} or \code{"Multi-dimensional model"}). If \code{sorts=""}, individual results will be sorted by the participant identifier.
}
  \item{extras}{
     character vector specifying additional summary results stored in the sta object that should be printed. These include the \code{"Data sources"}, prior probabilities (i.e., \code{"Trace model prior probability"}, \code{"No-Overlap model prior probability"}, \code{"Uni-dimensional model prior probability"} and \code{"Multi-dimensional model prior probability"}) and raw sample counts (i.e., \code{"Non-Trace model samples from encompassing model"}, \code{"Total samples from encompassing model"}, \code{"Non-overlapping monotonic samples from trace model"}, \code{"Overlapping monotonic samples from trace model"}, \code{"Non-monotonic samples from trace model"}, \code{"Monotonic model counts"} and \code{"Total samples from trace model"}).
}
  \item{guiarg}{
     hidden argument relating to the multi-option list display available for the \code{sorts} and \code{extras} arguments in the gui version of \code{stSummary}.
}
}
\details{
 Note that both a GUI (i.e., \code{guistSummary}) and non-GUI (i.e., \code{stSummary}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}.
}
\value{
 \code{stSummary} will print the selected summary results in the R console. It will first report whether the sampling is complete and the accuracy criteria (credible interval and precision) used for estimating the posterior proportions. If the sampling is incomplete, it will print the accuracy criteria and the estimated computation time remaining to satisfy these parameters. \code{stSummary} will then report the number of participants included in the summary results, followed by the group aggregate model selection results (as either Bayes factors or posterior model probabilities). If any additional summary values were selected (e.g., individual participant model selection results, prior probabilities or raw sample counts) these are printed after the group level results.
}
\references{
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.
}
\seealso{
 \code{\link{stProbplot}}, for a graphical method of displaying individual participant posterior model probabilities.
}
\examples{
\dontrun{
  ## stSummary must be given an existing sta object which has completed at least 
  ## one pass of sampling
  ## see stFirst for an example to produce an sta object that can 
  ## then be used to run the following.

  ## To extract the posterior model probabilities for the "DFIE.sta" object and 
  ## also report the posterior model probabilities for the individual participants 
  ## (sorted in order of the "No-Overlap model" results) as well as the prior 
  ## probabilities for the Trace model, No-Overlap model and Uni-dimensional model use:
    
  stSummary(bosname="DFIE.sta", sorts = "No-Overlap model", 
            extras = c("Trace model prior probability",
            "No-Overlap model prior probability", 
            "Uni-dimensional model prior probability"))
}
}