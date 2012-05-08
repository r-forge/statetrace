\name{guista}
\alias{guista}
\title{
Start GUI for StateTrace. 
}
\description{
\code{guista} displays the main user interface for the \code{StateTrace} package. This interface contains a button linking to each of the seven functions (with an available GUI display), which are used to manage the entire analysis process: from entering the raw data, to monitoring and refining the sampling estimates, as well as for obtaining summary statistics and plots. Note there is no non-GUI equivalent for this function.
}
\usage{
guista()
}
\value{
This function returns nothing. The raw data can be loaded, sampling initiated and refined, and summary results extracted by selecting one of the buttons on the displayed GUI. 
}
\seealso{
To generate an sta object and run the first pass of sampling, use \code{\link{stFirst}}. Once a first pass is complete, use \code{\link{stSample}} to obtain accurate posterior model probability estimates. Then examine results with \code{\link{stSummary}} and \code{\link{stProbplot}}. Run \code{\link{stBootav}} to obtain participant results for average trace and monotonic plots generated with \code{\link{stPlot}}. Use \code{\link{staManage}} to join sta objects and manage stored samples.
}
\examples{
## Open main GUI
guista()
}