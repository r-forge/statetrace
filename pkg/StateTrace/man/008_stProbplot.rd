\name{stProbplot}
\alias{stProbplot}
\alias{guistProbplot}
\title{
Plot individual participant model selection probabilities
}
\description{
  The individual participant posterior model probabilities can be assessed graphically using \code{stProbplot}, which will generate a four panelled plot (one for each model) that displays the posterior model probabilities for each participant as well as optional information about the corresponding group aggregate results.
}
\usage{
guistProbplot()
stProbplot(bosname = "", exclude = NULL,  maint = "", maino = "", 
          mainu = "", mainm = "", ylabt = "p(Non-Trace)", 
          ylabo = "p(No-Overlap)", ylabu = "p(Uni-dimensional)", 
          ylabm = "p(Multi-dimensional)", xlab = "Sorted Participants",
          symb = 1, weakl = TRUE, strongl = FALSE, plines = FALSE,
          pnames = TRUE, ymin = 0, ymax = 1, guiarg = NULL)
}
\arguments{
  \item{bosname}{
     character string indicating the name of the sta object from which results should be extracted.
}
  \item{exclude}{
     boolean vector indicating participants to be excluded.
}
  \item{maint}{
     Title for the Non-Trace panel.
}
  \item{maino}{
     Title for the No-Overlap panel.
}
  \item{mainu}{
     Title for the Uni-dimensional panel.
}
  \item{mainm}{
     Title for the Multi-dimensional panel.
}
  \item{ylabt}{
     y-axis label for the Non-Trace panel.
}
  \item{ylabo}{
     y-axis label for the No-Overlap panel.
}
  \item{ylabu}{
     y-axis label for the Uni-dimensional panel.
}
  \item{ylabm}{
     y-axis label for the Multi-dimensional panel.
}
  \item{xlab}{
     x-axis label, which is common to all four panels.
}
  \item{symb}{
     plotting character (i.e., symbol) to use. This value may either be assigned using an integer code (see \code{?points} for the available symbols in R) or by using a character string to define the symbol (i.e., \code{"Unfilled circles"}, \code{"Filled circles"}, \code{"Unfilled triangles"}, \code{"Filled triangles"}, \code{"Unfilled squares"}, \code{"Filled squares"}, \code{"Lower case letters"}, \code{"Upper case LETTERS"} or \code{"Participant numbers"}).
} 
  \item{weakl}{
     logical value indicating whether a dashed line should be included at Raftery's (1995) criteria for equivocal, positive and strong evidence; namely p(0.05, 0.25, 0.5, 0.75, 0.95).
}
  \item{strongl}{
     logical value indicating whether a dashed line should show Raftery's (1995) criterion for very strong evidence; namely p(0.01, 0.99).
}
   \item{plines}{
     logical value specifying if a heavy dashed line should be included to illustrate the group aggregate result for each respective model.
}
  \item{pnames}{
     logical value specifying if the group posterior model probability should be included in the panel title.
}
  \item{ymin}{
    minimum y-axis value (0 - 1). 
}
  \item{ymax}{
     maximum y-axis value (0 - 1).
}
  \item{guiarg}{
     hidden argument relating to the multi-option list available for the \code{symb} argument in the GUI version of \code{stProbplot}.
}
}
\details{
  Note that both a GUI (i.e., \code{guistProbplot}) and non-GUI (i.e., \code{stProbplot}) version is available for this function. An abridged description of the GUI is available in Prince, Hawkins, Love and Heathcote (2011). Alternatively, a detailed example is provided in the \code{StateTrace} vignette, which can be accessed using \code{vignette(topic="StateTrace",package="StateTrace")}.\cr 
  If using \code{guistProbplot} note that the GUI will remain open after the function has been executed allowing the user to progressively customise the plot without re-calling the function each time a parameter value is altered. Note that clicking the 'Cancel' button will dismiss the GUI but does not cancel the most recent parameter values assigned or the output reproduced.
}
\value{
  \code{stProbplot} will return a four panelled plot in a separate graphics device. Each panel corresponds to one of the four models (Non-Trace, No-Overlap, Uni-dimensional, Multi-dimensional) and plots the posterior model probabilities for each individual participant.
}
\section{known issues}{
  Mac users please note that there is a known bug in the GUI version of this function. This issue seems to be related to the progressive customisation available to both of the plotting functions (\code{stProbplot} and \code{stPlot}) when using the GUI. There are no issues with either plotting function when all input is provided through the command line and so we recommend this method to Mac users (a command line example is provided below). However, if a GUI is desired, Mac users should note that when directly calling \code{guistProbplot()} it seems that clicking the OK button twice will produce the required output. Alternatively, when first calling \code{guista()} and then clicking the button for \code{stProbplot}, the only way to view the generated plot is to click the "Cancel" button. We are working on addressing this issue, but any feedback or suggestions would be welcome.
}
\references{
Prince, M., Hawkins, G., Love, J., & Heathcote, A. (2011). An R package for state-trace analysis. Manuscript submitted for publication.\cr
Raftery, A.E. (1995). Bayesian model selection in social research. \emph{Sociological Methodology, 25,} 111-163.
}
\seealso{
 \code{\link{stSummary}}, for a non-graphical method of assessing individual participant posterior model probabilities.
}
\examples{
\dontrun{
  ## stProbplot must be given an existing sta object which has completed at least 
  ## one pass of sampling
  ## see stFirst for an example to produce an sta object that can 
  ## then be used to run the following:

  ## To create a plot where each subjects' data is represented by their participant 
  ## number.
    
  stProbplot(bosname = "DFIE.sta",  maint = "Non-Trace Model", 
             maino = "No-Overlap Model", mainu = "Uni-dimensional Model", 
             mainm = "Multi-dimensional Model", symb = "Participant numbers")
}
}
