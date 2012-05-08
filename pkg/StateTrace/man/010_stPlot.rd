\name{stPlot}
\alias{stPlot}
\alias{guistPlot}
\title{
State-Trace plots for the group average and individual participants
}
\description{
  \code{stPlot} offers a method for generating the state-trace plot, which can be used to visualise the estimates obtained by sampling from the encompassing posterior.
}
\usage{
guistPlot()
stPlot(bosname = "", exclude = NULL, main = NA, xlab = NA, ylab = NA, 
       legnamd1 = NA, legnamd2 = NA, acc = TRUE, alls = TRUE, 
       symbd1 = 1, symbd2 = 2, stat = "mode", line = "e", lpts = TRUE, 
       preg = TRUE, p = 0.68, smoothfac = 5, linew = 4, 
       xmin = NA, xmax = NA, ymin = NA, ymax = NA, guiarg = NULL)
}
\arguments{
  \item{bosname}{
     character string indicating the name of the sta object
}
  \item{exclude}{
     boolean vector indicating participants to be excluded.
}
  \item{main}{
     title to display above the state-trace plot.
}
  \item{xlab}{
     x-axis label corresponding to the first level of the state factor. If left as the default \code{NA} value, this argument will be filled with the code provided in the raw data file to denote the first level of the state factor. (NB; this is also true for the following \code{ylab}, \code{legnamd1} and \code{legnamd2} arguments). Note also that the accuracy measure specified for the plot will also be included in parentheses after the axis label. 
}
  \item{ylab}{
     y-axis label corresponding to the second level of the state factor. 
}
  \item{legnamd1}{
     first level of the dimension factor to appear in the legend of the plot.  
}
  \item{legnamd2}{
     second level of the dimension factor to appear in the legend of the plot. 
}
  \item{acc}{
     logical value specifying whether a probability measure should be used when calculating the bootstrap averages and generating the state-trace plot. When \code{TRUE}, the proportion correct is used for B0 designs and the Hit minus False Alarm rate is used for B2 designs. When \code{FALSE}, the corresponding inverse cumulative normal (\emph{z}) transformation is used; i.e., \emph{z}(proportion correct) for B0 designs and the signal detection \emph{d'} measure for B2 designs.  Note this value must match the \code{acc} value used for \code{stBootav}.
}
  \item{alls}{
     logical argument indicating if the state-trace plot should be generated for the group average (\code{TRUE}) or for each individual participant (\code{FALSE}).
}
  \item{symbd1}{
    plotting character (i.e., symbol) to use for the first level of the dimension factor. This value may either be assigned using an integer code (see \code{?points} for the available symbols in R) or by using a character string to define the symbol (i.e., \code{"Unfilled circles"}, \code{"Unfilled upright triangles"}, \code{"Unfilled inverted triangles"}, \code{"Unfilled squares"} or \code{"Unfilled diamonds"}). Note that by default only unfilled symbols are used as \code{stPlot} will label each shape with a number corresponding to the level of the trace factor.
}
  \item{symbd2}{
    plotting character (i.e., symbol) to use for the second level of the dimension factor. As above, this value may either be assigned using an integer code or by using a character string to define the symbol.
}
  \item{stat}{
     statistic to use for model plotting. This value can be set as either the \code{"mean"}, \code{"median"}, or \code{"mode"} of the encompassing posterior.
}
  \item{line}{
     type of line to display on the state-trace plot. This value can either be set to display no lines (\code{"n"}), the data traces (\code{"e"}), or the best (i.e., most frequent) trace (\code{"t"}) or monotonic (\code{"m"}) model.
}
  \item{lpts}{
     logical argument specifying if half-size versions of the data plot symbols should be included on the selected line (except for \code{line="e"}, data traces, as by definition these lines join the data points).
}
  \item{preg}{
     logical argument indicating if the credible p regions for each posterior estimate should be displayed on the plot. These credible regions are obtained using a linear binned two-dimensional kernal smoother (see Wand & Jones, 1995).
}
  \item{p}{
     integer value specifying the width of the credible p region.
}
  \item{smoothfac}{
     integer value indicating the kernal smoothing factor used to find the center of the posterior distribution.
}
  \item{linew}{
     width of the plot line/s (1-10).  
}
  \item{xmin}{
    minimum x-axis value (0 - 1). Note that \code{stPlot} will automatically scale the axes of the plot to fit the data when left at the default "NA" value.
}
  \item{xmax}{
     maximum x-axis value (0 - 1).
}
  \item{ymin}{
     minimum y-axis value (0 - 1).
}
  \item{ymax}{
     maximum y-axis value (0 - 1).
}
  \item{guiarg}{
     hidden argument relating to the multi-option list available for the \code{stat}, \code{line}, \code{symbd1} and \code{symbd2} arguments in the GUI version of \code{stProbplot}.
}
}
\details{
  Note that both a GUI (i.e., \code{guistPlot}) and non-GUI (i.e., \code{stPlot}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}. \cr
  If using \code{guistPlot} note that the GUI will remain open after the function has been executed allowing the user to progressively customise the plot without re-calling the function each time a parameter value is altered. Note that clicking the 'Cancel' button will dismiss the GUI but does not cancel the most recent parameter values assigned or the output reproduced.
}
\value{
  \code{stPlot} will return a single plot if \code{alls = TRUE} or one plot for each individual participant if \code{alls = FALSE}. In both cases the plots will appear in separate graphics devices.
}
\section{known issues}{
  Mac users please note that there is a known bug in the GUI version of this function. This issue seems to be related to the progressive customisation available to both of the plotting functions (\code{stProbplot} and \code{stPlot}) when using the GUI. There are no issues with either plotting function when all input is provided through the command line and so we recommend this method to Mac users (a command line example is provided below). However, if a GUI is desired, Mac users should note that when directly calling \code{guistPlot()} it seems that clicking the OK button twice will produce the required output. Alternatively, when first calling \code{guista()} and then clicking the button for \code{stPlot}, the only way to view the generated plot is to click the "Cancel" button. We are working on addressing this issue, but any feedback or suggestions would be welcome.
}
\references{
 Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
 Wand, M.P., & Jones, M.C. (1995). \emph{Kernal Smoothing.} Chapman and Hall, London. 
}
\seealso{
 \code{\link{stBootav}}, to obtain the bootstrap averages used in generating the state-trace plots.
}
\examples{
\dontrun{
  ## stPlot must be given an existing sta object which has completed at least 
  ## one pass of sampling and has obtained the bootstrap samples from at least 
  ## the encompassing model.
  ## see stFirst for an example to produce an sta object that can 
  ## then be used to run the following

  ## To create a state-trace plot for the group aggregate, use:
    
  stPlot(bosname = "DFIE.sta", xlab = "Face Accuracy", ylab = "House Accuracy", 
         legnamd1 = "Upight", legnamd2 = "Inverted", 
         symbd1 = "Unfilled upright triangles", 
         symbd2 = "Unfilled inverted triangles")
}
}