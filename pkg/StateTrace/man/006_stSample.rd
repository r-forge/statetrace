\name{stSample}
\alias{stSample}
\alias{guistSample}
\title{
Refine Sampling Estimates
}
\description{
The sampling estimates can be further refined using \code{stSample}, which will repeatedly draw samples from the encompassing and trace models until the required precision (which can be set independently for each model) is satisfied or the maximum computation time allowed has elapsed.
}
\usage{
guistSample()
stSample(bosname = "", refresh = FALSE, maxt = 0, ci = 95, 
         BFe = TRUE, sampe = 1e+05, ed = 5e-04, nkeepe = 10000, 
         BFt = TRUE, sampt = 5000, burn = 100, td = 0.005,  
         nkeept = 10000, nkeepm = Inf,
         verbose = 1)
}
\arguments{
  \item{bosname}{
     character string indicating the name of the sta object.
}
  \item{refresh}{
     logical value corresponding to two modes available for running \code{stSample}: \code{Refresh = TRUE} is used to update or "refresh" the control parameters given \code{d} and \code{ci}. This mode is fast (i.e., no sampling) unless no sampling has been done yet (e.g., when using \code{staMake} to create the sta object), in which case it will complete one pass to get the necessary timing information. \code{Refresh = FALSE} will continue sampling and updating the sta object for at most time = maxt (hours). If however, \code{maxt = 0} just one pass will be completed.
}
  \item{maxt}{
     maximum run time (hours).
}
  \item{ci}{
     desired size of the credible interval for the estimated posterior proportions (0 - 100). Note this parameter may also be specified as a scalar or vector with one entry per participant.
}
  \item{BFe}{
     logical value indicating if \code{stSample} should continue sampling from the encompassing model.
}
  \item{sampe}{
     number of samples to draw from the encompassing model per run.
}
  \item{ed}{
     desired precision for the posterior estimates from the encompassing samples (0 - 1). Note this parameter may also be specified as a scalar or vector with one entry per participant.
}
  \item{nkeepe}{
     number of encompassing model samples to keep stored in the sta object, which will be used in plotting (e.g., \code{\link{stBootav}} and \code{\link{stPlot}}).
}
  \item{BFt}{
     logical value indicating if \code{stSample} should continue sampling from the trace model.
}
  \item{sampt}{
     number of samples to draw from the trace model per run.
}
  \item{burn}{
     number of samples from the trace model per run to discard as "burn-in".
}
  \item{td}{
     desired precision for the posterior estimates from the trace samples (0 - 1). Note this parameter may also be specified as a scalar or vector with one entry per participant.
}
  \item{nkeept}{
     number of trace model samples to keep stored in the sta object, which will be used in plotting (e.g., \code{\link{stBootav}} and \code{\link{stPlot}}).
}
  \item{nkeepm}{
     number of monotonic model samples to keep stored in the sta object, which will be used in plotting (e.g., \code{\link{stBootav}} and \code{\link{stPlot}}).
}
  \item{verbose}{
     numeric value (\code{0, 1, 2}) controlling the information printed to the R console during sampling: \code{0} is silent, \code{1} prints the estimated total time remaining after each run for all participants and \code{2} adds estimated timings per participant.
}
}
\details{
\code{stSample} may be required to complete further sampling if: (a) the analysis did not reach completion during the previous pass; (b) additional subjects were added to the sta object or (c) the user wishes to alter any of the parameter values (e.g., credible interval and corresponding precision, or the number of samples drawn per run) that were used in the previous pass. This process of running additional passes can be repeated as many times as is necessary for all argument values to be satisfied. \cr
\cr
Note that both a GUI (i.e., \code{guistSample}) and non-GUI (i.e., \code{stSample}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}.
}
\value{
\code{stSample} will return an updated version of the sta object and depending on the \code{verbose} value specified, will also print updates in the R console for the user to monitor the progress of the sampling.
}
\references{
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.
}
\seealso{
 \code{\link{staMake}}, for a detailed break-down of the sta object.
 }
\examples{
\dontrun{
## stSample must be given as existing sta object: 
      ## load DFIE data set 
      data(DFIE)
      ## save data as a text file in the current working directory
      write.table(DFIE,file="DFIE.txt",sep="\t",row.names=FALSE)
      ## generate an empty sta object called "DFIE.sta"

      staMake(staname="DFIE.sta",fnams="DFIE.txt",multiparticipant=TRUE)

## To complete a single pass of the data (maxt=0) and receive maximal timing updates
    ## Note this takes approximately 15 minutes to run (Netbook with 1GB RAM)
 stSample(bosname = "DFIE.sta", verbose = 2)
}
}