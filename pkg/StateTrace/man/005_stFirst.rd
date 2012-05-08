\name{stFirst}
\alias{stFirst}
\alias{guistFirst}
\title{
Generate sta object and complete the first pass
}
\description{
The \code{stFirst} function provides a wrapper for several other functions that have arguments set at values minimally sufficient for initiating analysis.
}
\usage{
guistFirst()
stFirst(staname = "", fnams = NULL, folder = "", extension = "txt",  
        sep = "tab",multiparticipant = FALSE, header = TRUE,  
        usecols = NULL, na.strings = "NA", acc = TRUE)
}
\arguments{
  \item{staname}{
     character string assigning a name to the sta object
}
  \item{fnams}{
     character vector with either a single element specifying the name of a file in the working directory or the path + name of a file if included in another directory. This argument value may also be multiple elements containing such specifications. Alternatively, use the following \code{folder} and \code{extension} arguments to specify the files to load.
}
  \item{folder}{
     character string of the directory containing the data file/s
}
  \item{extension}{
     character string of the file extension of the data file/s. The \code{folder} and \code{extension} arguments can be used to load all of the data files contained in the specified folder by indicating the directory to search in (i.e., \code{folder}) and the type of file to load (i.e., \code{extension}).
}
  \item{sep}{
     field separator character used in the data files.
}
  \item{multiparticipant}{
     logical value specifying if the data files contain data for multiple participants (\code{TRUE}) or only a single participant (\code{FALSE}).
}
  \item{header}{
     logical argument to specify if the data file contains a header row.
}
  \item{usecols}{
     vector of integers or character strings specifying the columns to use from the data file.
}
  \item{na.strings}{
     string used to specify empty cells in the data files. Where this string is located in the fifth column of the file (for trial files, column "C" and for summary files column "n"), these rows will be stripped from further analysis.
}
  \item{acc}{
     logical value specifying whether a probability measure should be used when calculating the bootstrap averages and generating the state-trace plot. When \code{TRUE}, the proportion correct is used for B0 designs and the Hit minus False Alarm rate is used for B2 designs. When \code{FALSE}, the corresponding inverse cumulative normal (\emph{z}) transformation is used; i.e., \emph{z}(proportion correct) for B0 designs and the signal detection \emph{d'} measure for B2 designs. 
}
}
\details{
   The \code{stFirst} function will first call \code{staMake} to read-in the raw data and create an empty sta object. \code{stSample} is then called to complete two runs sampling from the encompassing and trace models for each participant and the convergence of the MCMC trace samples is assessed using \code{staManage}. Finally, \code{stBootav} is called to obtain the bootstrap average results from the encompassing posterior and a state-trace plot is created for each participant and the group aggregate using \code{stPlot}.\cr
   \cr
   Note that both a GUI (i.e., \code{guistFirst}) and non-GUI (i.e., \code{stFirst}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}. 
}
\value{
\code{stFirst} will first create an sta object, which will hold information extracted from the raw data and will be updated as the sampling process is initiated. It will also print a number of updates in the R console to allow the user to monitor the progress of the analysis. First, the 'datasources' of each file being read-in are reported, followed by a summary of the design extracted from the raw data files. Next updates relating to the sampling from the encompassing posterior are printed for each participant and overall: (a) the current precision of the credible interval of the trace posterior estimate for each independent chain as well as (b) the estimated time remaining (in minutes). Equivalent updates are then reported for the sampling from the trace model, including (a) the current precision of the credible intervals for the no-overlap, uni-dimensional and multi-dimensional posterior estimates, and (b) the estimated time remaining (in minutes). Finally, once the first pass of sampling is complete \code{stFirst} reports Gelman's multivariate "R-hat" statistic for MCMC convergence, draws 10,000 bootstrap average samples from the encompassing model, and creates a preliminary state-trace plot for each participant and the group aggregate, which are all displayed in individual graphics devices.
}
\references{
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
}
\seealso{
 \code{\link{staMake}}, to read in raw data and create an empty sta object but not complete any of the extras done by \code{stFirst}.\cr
 For a detailed description of the default parameter values and the value returned by each function called by \code{stFirst} see, \code{\link{staMake}}, \code{\link{stSample}}, \code{\link{staManage}}, \code{\link{stBootav}}, \code{\link{stPlot}}, and see \code{\link{DFIE}}, for example data suitable for \code{StateTrace}.
}
\examples{
\dontrun{
## load DFIE data set
data(DFIE)
## save data as a text file in the current working directory
write.table(DFIE,file="DFIE.txt",sep="\t",row.names=FALSE)

## generate an sta object called "DFIE.sta" and run first pass
stFirst(staname="DFIE.sta",fnams="DFIE.txt",multiparticipant=TRUE)
}
}