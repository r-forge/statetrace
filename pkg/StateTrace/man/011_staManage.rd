\name{staManage}
\alias{staManage}
\alias{guistaManage}
\title{
Manage sta objects
}
\description{
  \code{staManage} allows the user to bind multiple sta objects, check MCMC convergence and save the mcmc object for further analysis by coda, as well as manage the number of stored samples in the sta object (and hence its size).
}
\usage{
guistaManage()
staManage(stanames = "", staname = stanames[1], 
          checkConverge = FALSE, nmcmc = 5000, mcmcSave = "", 
          nkeepe = Inf, nkeept = Inf, nkeepm = Inf, 
          keepbestm = FALSE)
}
\arguments{
  \item{stanames}{
     either a single character name, in which case stored samples are managed or alternatively a vector of character names, in which case the corresponding objects are joined and saved to a new object.
}
  \item{staname}{
     optional character name for the new, combined sta object. By default, \code{staManage} will save the combined object back into the first element of \code{stanames}. 
}
  \item{checkConverge}{
     logical argument specifying if the convergence of the MCMC trace chains should be checked using the \code{coda} package.
}
  \item{nmcmc}{
     length of the number of trace samples drawn per run. This value needs to be some multiple of the total number of samples drawn and should produce at least two chains.
}
  \item{mcmcSave}{
     character name for the mcmc object if the user wants to save the MCMC samples. This object (which is a list of \code{mcmc.list} objects) can then be further analysed using \code{coda}.
}
  \item{nkeepe}{
     number of encompassing model samples to keep stored in the sta object.
}
  \item{nkeept}{
     number of trace model samples to keep stored in the sta object.
}
  \item{nkeepm}{
     number of monotonic model samples to keep stored in the sta object.
}
  \item{keepbestm}{
     logical argument which can specify if only the best (i.e., most frequently occurring, and hence most probable) monotonic model samples should remain stored in the sta object. If \code{TRUE}, removes all but the best monotonic samples, as recorded in \code{$m$control$nm} rather than in the actual stored samples.
}
}
\details{
  Note that both a GUI (i.e., \code{guistaManage}) and non-GUI (i.e., \code{staManage}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}.
}
\value{
  \code{staManage} will return an updated version of the sta object if asked to combine multiple objects or if asked to modify the number of samples stored in an existing sta object. If asked to check the convergence of the MCMC trace chains, Gelman's multivariate "R-hat" statistic for MCMC convergence is reported in the R console. Finally if asked to save the MCMC samples as a separate object, \code{staManage} will create a new object in the R workspace which is a list of \code{mcmc.list} objects (i.e., one list per participant).
}
\references{
 Plummer, M., Best, N., Cowles, K., & Vines, K. (2006). CODA: Convergence Diagnosis and Output Analysis for MCMC. \emph{R News, 6}, 7-11.\cr
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
}
\examples{
\dontrun{
#### Example 1####
  ## staManage must be given an existing sta object which has completed at least 
  ## one pass of sampling from the trace model in order to check the convergence 
  ## of these chains.
  ## see stFirst for an example to produce an sta object that can 
  ## then be used to run the following

  ## To check the convergence of the trace chains and save the mcmc object:

  staManage(stanames = "DFIE.sta", checkConverge = TRUE, nmcmc = 5000, 
            mcmcSave = "DFIEmcmc")


          
#### Example 2 ####
  ## If the user had a very large sample to analyse (such as when running a 
  ## simulation study), the data could be separated into more manageable  
  ## portions and run simultaneously on separate machines.
  ## the separate sta objects (e.g., DFIE-1, DFIE-2 etc.) can then be 
  ## combined using:

  staManage(stanames = c("DFIE-1","DFIE-2"), staname = "DFIE")
         
}
}