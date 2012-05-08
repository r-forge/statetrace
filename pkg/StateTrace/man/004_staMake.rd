\name{staMake}
\alias{staMake}
\title{
Generate sta object
}
\description{
\code{staMake} allows the user to read-in the appropriate data file/s, which will then be used to create the sta object that holds the generated output. 
}
\usage{
staMake(staname = "", fnams = NULL, folder = "", extension = "txt", 
        usecols = NULL, header = TRUE, sep = "\t", na.strings = "NA",
        multiparticipant = FALSE, a = 1)
}
\arguments{
  \item{staname}{
     character string assigning a name to the sta object.
}
  \item{fnams}{
     character vector with either a single element specifying the name of a file in the working directory or the path + name of a file if included in another directory. This argument value may also be multiple elements containing such specifications. Alternatively, the following \code{folder} and \code{extension} arguments can be used to specify the files to load.
}
  \item{folder}{
     character string of the directory containing the data file/s.
}
  \item{extension}{
     character string of the file extension of the data file/s. The \code{folder} and \code{extension} arguments can be used to load all of the data files contained in the specified folder by indicating the directory to search in (i.e., \code{folder}) and the type of file to load (i.e., \code{extension}).
}
  \item{usecols}{
     vector of integers or character strings specifying the columns to use from the data file/s.
}
  \item{header}{
     logical argument to specify if the data file/s contains a header row.
}
  \item{sep}{
     field separator character used in the data file/s.
}
  \item{na.strings}{
     string used to specify empty cells in the data file/s. Where this string is located in the fifth column of the file (for trial files column "C", and for summary files column "n"), these rows will be stripped from further analysis.
}
  \item{multiparticipant}{
     logical value specifying if the data file contains data for multiple participants (\code{TRUE}) or only a single participant (\code{FALSE}) 
}
  \item{a}{
     specifies the first beta shape parameter for the prior distribution. Here we assume an encompassing prior of uniform distributions, with the uniform being a special case of the Beta distribution, Beta(a, b), where a = b = 1 for the uniform.
}
}
\details{
   Note there is no GUI equivalent for this function. However, if a GUI is desired by the user \code{\link{guistFirst}} can be used as an alternative for reading in raw data file/s and creating the sta object. However, unlike \code{staMake} which will only read-in the data and create an empty sta object, \code{guistFirst} will also sample from the encompassing and trace models, check the convergence of the MCMC trace chains, calculate bootstrap average results and generate preliminary state-trace plots at both the individual participant and group aggregate level.
}
\value{
  Creates the sta object, which is a list of lists. At the highest level:
\item{all }{results summarising over participants}
\item{ss }{list with one entry per participant, with attributes \code{"slevs"}, \code{"dlevs"} and \code{"nt"} giving the names of the state and dimension factor levels and the number of trace levels, common to all participants. Each entry of \code{ss} is named '1', '2', ... and has an attribute \code{"datasouce"} giving the data file name from which the data was obtained (or the "P" column level corresponding to that participant when data was all obtained from one 'multiparticipant' file).}
  Each entry in \code{ss} is a list with components:
\item{d }{augmented representation of the participant's raw data. Note this is the only component filled in on initial creation}
\item{m }{statistics for each model: \code{t} = trace model, \code{m} = monotonic model}
\item{p }{plotting information for each model: \code{e} = encompassing model, \code{t} = trace model, \code{m} = monotonic model}
 Contained within the data representation, \code{d}, list is:
\item{sname }{extracts the participant's identifier; e.g., 'P1', 'P2'...}
\item{m }{matrix containing the number of successes for each design cell obtained from the raw data file (rows corresponds to each state x dimension combination, while the columns represent the trace factor levels)}
\item{n }{matrix containing the total number of observations for each design cell, obtained from the raw data file}
\item{s1 }{first Beta shape parameter matrix for the encompassing posterior, \code{m+a}, where a=1 for the uniform prior}
\item{s2 }{second Beta shape parameter matrix for the encompassing posterior, \code{n-m+b}, where b=1 for the uniform prior}
\item{odr }{list containing the order expected for trace factor levels, with one entry for each independent chain named 'c1', 'c2'...}
\item{priorp }{prior probability for observing the order specified by the trace factor for each chain. This prior is calculated analytically from Beta (for details see Prince, Brown & Heathcote, 2011)}
\item{dsn }{design of the sta object; either 'B0' (i.e., no baseline) or 'B2' (i.e., state baseline)}
\item{chain.names }{names of the chains that can be run independently, where 'D1' = Dimension level 1, 'D2' = Dimension level 2, 'S1' = State level 1, 'S2' = State level 2}
  Each element in the model statistics, \code{m}, list (for both the trace, \code{$t}, and monotonic, \code{$m}, models) has \code{control} (i.e., ancillary) information used in obtaining the statistics. This information includes:
\item{models }{the models under examination; for \code{$m$t}, 'M1'=trace and 'M2'=encompassing and for \code{$m$m}, 'M1'=1D and 'M2'=MultiD}
\item{d }{precision of the credible interval specified by the user; this value is set independently for \code{$m$t} and for \code{$m$m}}
\item{doned }{logical vector indicating whether the sampling is completed to the specified precision; for \code{$m$t} this is reported for each chain}
\item{ci }{credible interval specified by the user}
\item{tpers }{average time taken per sample, in hours}
\item{moret }{optimal estimate for how much more computation time will be required, in hours}
\item{moren }{optimal estimate for how many more samples are required}
\item{cis }{current observed credible interval; for \code{$m$t} this is reported for each independent chain and for \code{$m$m} is reported for the no-overlap, uni-dimensional and multi-dimensional models}
\item{dcis }{current observed precision of the credible interval; as above, this is reported for each independent chain for \code{$m$t} and for the no-overlap, uni-dimensional and multi-dimensional models for \code{$m$m}}
\item{ntraces }{for \code{$m$t} only; reports the number of samples drawn from the encompassing model that respect the trace order per chain}
\item{Ntraces }{for \code{$m$t} only; reports the total number of samples drawn from the encompassing model per chain}
\item{pp1d }{for \code{$m$m} only; prior probability for the 'NoLap' (i.e., non-overlapping portion of the monotonic model) and 'OLap' (i.e., overlapping portion of the monotonic model)}
\item{nm }{for \code{$m$m} only; vector of the monotonic models observed from the trace model samples, and the number of times each monotonic model was observed}
\item{n1d }{for \code{$m$m} only; count of the observed monotonic models where data traces were non-overlapping ('NoLap') and overlapping ('OLap')}
\item{Nm }{for \code{$m$m} only; total number of samples drawn from the trace model}
 Contained within the plotting, \code{p}, list is:
\item{e }{list of the samples drawn from the encompassing model}
\item{t }{list of the samples respecting the trace model}
\item{m }{list of the observed monotonic model samples}
}
\references{
Prince, M., Brown, S., & Heathcote, A. (2011). The design and analysis of state-trace experiments. \emph{Psychological Methods.}
}
\seealso{
 \code{\link{guistFirst}}, for a GUI alternative for reading in raw data file/s and creating the sta object.\cr
 \code{\link{DFIE}}, for example data suitable for \code{StateTrace}.
}
\examples{
\dontrun{
## load DFIE data set
data(DFIE)
## save data as a text file in the current working directory
write.table(DFIE,file="DFIE.txt",sep="\t",row.names=FALSE)

## generate an empty sta object called "DFIE.sta"
staMake(staname="DFIE.sta",fnams="DFIE.txt",multiparticipant=TRUE)
}
}