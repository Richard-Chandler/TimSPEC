\name{SSP585data}
\alias{SSP585data}
\alias{GlobalTemps585}
\alias{ERF585}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data for use in \pkg{TimSPEC} package examples}
\description{The datasets provided here contain annual estimates of global mean surface temperature from 1950 to 2022, together with projected temperatures from five global climate models under shared socioeconomic pathway 5.8-5 and the corresponding global net radiative forcings.}
%- maybe also 'usage' for other objects documented here.
\details{The \code{SSP585data} dataset contains two data frames, as follows:
\describe{
\item{\code{GlobalTemps585}:}{This contains the temperature data. Column 1 is \code{Year}, from 1950 to 2100; Column 2 is \code{Obs} and contains estimated global temperatures (in degrees Celsius) to 2022, obtained from the HadCRUT5 dataset (\cite{Morice et al. 2020}). The remaining columns are produced using five of the global climate models (GCMs) contributing to Phase 6 of the Coupled Model Intercomparison Project (CMIP6; \cite{Eyring et al. 2016}): the column names are of the form \code{$$$$$$$$.###} where \code{$$$$$$$$} is a code representing the GCM that was used and \code{###} is the run number, starting from zero (see the Examples section below). 

Tha HadCRUT data represent estimates derived from land surface stations, ships and buoys: data from 1850 onwards are publicly available from \url{https://crudata.uea.ac.uk/cru/data/temperature/}, but the early part of the record is more uncertain due to a relative sparsity of measurements. The underlying measurements have at least 90\% global coverage continuously from 1950 onwards (the period considered here).

The series produced by the climate models were produced by running each model conditioned on a variety of inputs, including atmospheric greenhouse gas concentrations and other forcings. These atmospheric forcings are based on historical estimates until 2014, and subsequently on the 'SSP5-8.5' socioeconomic scenario of fossil-fueled global development and high greenhouse gas emissions (\cite{IPCC 2023}, Box SPM.1). Where multiple runs are available from the same model, they differ in the state of the climate that was used to initialise them - which usually corresponds to a plausible configuration at some point in the mid-19th century. 
}
\item{\code{ERF585}:}{This contains two columns: \code{Year} and \code{NetERF}. The latter is the net 'effective radiative forcing' (ERF) series, in watts per square metre, corresponding to the complex set of atmospheric forcings underpinning the climate model outputs in \code{GlobalTemps585}. The data are from \url{https://github.com/IPCC-WG1/Chapter-7}.}
}
}

\note{Usage of these data must be accompanied by appropriate acknowledgements as follows:

\describe{
\item{HadCRUT5:}{\emph{The HadCRUT5 global temperature data are provided by the Climatic Research Unit (University of East Anglia) and Met Office, and were calculated from information provided at \url{https://crudata.uea.ac.uk/cru/data/temperature/}.}}
\item{CMIP6 climate model outputs:}{\emph{We acknowledge the World Climate Research Programme, which, through its Working Group on Coupled Modelling, coordinated and promoted CMIP6. We thank the climate modeling groups for producing and making available their model output, the Earth System Grid Federation (ESGF) for archiving the data and providing access, and the multiple funding agencies who support CMIP6 and ESGF.}}
}
}

\author{Richard E. Chandler <r.chandler@ucl.ac.uk>}
\references{
Eyring, V., Bony, S., Meehl, G. A., Senior, C. A., Stevens, B., Stouffer, R. J., and Taylor, K. E. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. \emph{Geoscientific Model Development}, \bold{9(5)}:1937--1958. \doi{10.5194/gmd-9-1937-2016}.

IPCC (2023). Summary for policymakers. In Lee, H. and Romero, J., editors, \emph{Climate Change 2023:
Synthesis Report. A Report of the Intergovernmental Panel on Climate Change. Contribution of Working
Groups I, II and III to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change.}
Intergovernmental Panel on Climate Change, Geneva, Switzerland. Available from \url{https://www.ipcc.
ch/report/ar6/syr/}.

Morice, C. P., Kennedy, J. J., Rayner, N. A., Winn, J. P., Hogan, E., Killick, R. E., Dunn, R. J. H., Osborn, T. J., Jones, P. D., and Simpson, I. R. (2020). An updated assessment of near-surface temperature change from 1850: the HadCRUT5 dataset. \emph{J. Geophys. Res.}, \bold{126}:e2019JD032361. \doi{10.1029/2019JD032361}.
}
\examples{
## Load the data 
data(SSP585data)
## Plot the forcings
plot(ERF585$Year, ERF585$NetERF, type="l", 
     xlab="Year", ylab=expression("W ms"^-1), 
     main="Net radiative forcing, SSP5.8-5, 1950-2100")

## Tabulate numbers of runs from each climate model
print(names(GlobalTemps585)) # Year, Obs, Model1.000, Model1,001, ...
GCMIDs <- 
  unlist(lapply(strsplit(names(GlobalTemps585)[-(1:2)], "\\\\."), 
                FUN=function(x) x[1]))
print(table(GCMIDs))
}
