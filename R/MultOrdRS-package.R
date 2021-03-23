
#' Model Multivariate Ordinal Responses Including Response Styles
#' 
#' A model for multivariate ordinal responses. The response is modelled 
#' using a mixed model approach that is also capable of the inclusion 
#' of response style effects of the respondents.
#' 
#' 
#' @name MultOrdRS-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://orcid.org/0000-0002-0392-1580}
#' @seealso \code{\link{multordRS}}  \code{\link{ctrl.multordRS}}  \code{\link{plot.MultOrdRS}}
#' @keywords multivariate ordinal response style adjacent categories cumulative 
#' @references Schauberger, Gunther and Tutz, Gerhard (2021): Multivariate Ordinal Random Effects Models Including Subject and Group Specific Response Style Effects, 
#' \emph{Statistical Modelling}, \url{https://journals.sagepub.com/doi/10.1177/1471082X20978034}
#' @examples
#' data(tenseness)
#' 
#' ## create a small subset of the data to speed up calculations
#' set.seed(1860)
#' tenseness <- tenseness[sample(1:nrow(tenseness), 300),]
#' 
#' ## scale all metric variables to get comparable parameter estimates
#' tenseness$Age <- scale(tenseness$Age)
#' tenseness$Income <- scale(tenseness$Income)
#' 
#' ## two formulas, one without and one with explanatory variables (gender and age)
#' f.tense0 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ 1"))
#' f.tense1 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ Gender + Age"))
#' 
#' 
#' 
#' ####
#' ## Adjacent Categories Models
#' ####
#' 
#' ## Multivariate adjacent categories model, without response style, without explanatory variables
#' m.tense0 <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE, cores = 2))
#' m.tense0
#' 
#' \donttest{
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1 <- multordRS(f.tense0, data = tenseness)
#' m.tense1
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2 <- multordRS(f.tense1, data = tenseness, control = ctrl.multordRS(XforRS = FALSE))
#' m.tense2
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3 <- multordRS(f.tense1, data = tenseness)
#' m.tense3
#' 
#' plot(m.tense3)
#' 
#' 
#' 
#' ####
#' ## Cumulative Models
#' ####
#' 
#' ## Multivariate cumulative model, without response style, without explanatory variables
#' m.tense0.cumul <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE), 
#'   model = "cumulative")
#' m.tense0.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1.cumul <- multordRS(f.tense0, data = tenseness, model = "cumulative")
#' m.tense1.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2.cumul <- multordRS(f.tense1, data = tenseness, control = ctrl.multordRS(XforRS = FALSE), 
#'   model = "cumulative")
#' m.tense2.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3.cumul <- multordRS(f.tense1, data = tenseness, model = "cumulative")
#' m.tense3.cumul
#' 
#' plot(m.tense3.cumul)
#' 
#' ################################################################
#' ## Examples from Schauberger and Tutz (2020) 
#' ## Data from the German Longitudinal Election Study (GLES) 2017
#' ################################################################
#' 
#' ####
#' ## Source: German Longitudinal Election Study 2017 
#' ## Rossteutscher et al. 2017, https://doi.org/10.4232/1.12927
#' ####
#' 
#' ## load GLES data
#' data(GLES17)
#' 
#' ## scale data
#' GLES17[,7:11] <- scale(GLES17[,7:11])
#' 
#' ## define formula
#' f.GLES <- as.formula(cbind(RefugeeCrisis, ClimateChange, Terrorism, 
#'                        Globalization, Turkey, NuclearEnergy) ~ 
#'                        Age + Gender + Unemployment + EastWest + Abitur)
#' 
#' ## fit adjacent categories model without and with response style parameters
#' m.GLES0 <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(RS = FALSE, cores = 6))
#' m.GLES <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(cores = 6))
#' 
#' m.GLES0
#' m.GLES
#' 
#' plot(m.GLES, main = "Adjacent categories model")
#' 
#' 
#' ## fit cumulative model without and with response style parameters (takes pretty long!!!)
#' m.GLES20 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6, RS = FALSE))
#' 
#' m.GLES2 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6))
#' 
#' m.GLES20
#' m.GLES2
#' 
#' plot(m.GLES2, main = "Cumulative model")
#' 
#'}
NULL


#' Tenseness data from the Freiburg Complaint Checklist (tenseness)
#' 
#' Data from the Freiburg Complaint Checklist. The data contain all 8 items corresponding to the scale 
#' \emph{Tenseness} for 1847 participants of the standardization sample of the Freiburg Complaint Checklist. 
#' Additionally, several person characteristics are available. 
#' 
#' @name tenseness
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist with 1847 observations. 
#' All items refer to the scale \emph{Tenseness} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Clammy_hands}{Do you have clammy hands?}
#' \item{Sweat_attacks}{Do you have sudden attacks of sweating?}
#' \item{Clumsiness}{Do you notice that you behave clumsy?}
#' \item{Wavering_hands}{Are your hands wavering frequently, e.g. when lightning a cigarette or when holding a cup?}
#' \item{Restless_hands}{Do you notice that your hands are restless?}
#' \item{Restless_feet}{Do you notice that your feet are restless?}
#' \item{Twitching_eyes}{Do you notice unvoluntary twitching of your eyes?}
#' \item{Twitching_mouth}{Do you notice unvoluntary twitching of your mouth?}
#' \item{Gender}{Gender of the participant}
#' \item{Household}{Does participant live alone in a houshold or together with others?}
#' \item{WestEast}{is the participant from East Germany (former GDR) or West Germany?}
#' \item{Age}{Age in 15 categories, treated as continuous variable}
#' \item{Abitur}{Does the participant have Abitur (a-levels)?}
#' \item{Income}{Income in 11 categories, treated as continuous variable}
#'  }
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(tenseness)
#' 
#' ## create a small subset of the data to speed up calculations
#' set.seed(1860)
#' tenseness <- tenseness[sample(1:nrow(tenseness), 300),]
#' 
#' ## scale all metric variables to get comparable parameter estimates
#' tenseness$Age <- scale(tenseness$Age)
#' tenseness$Income <- scale(tenseness$Income)
#' 
#' ## two formulas, one without and one with explanatory variables (gender and age)
#' f.tense0 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ 1"))
#' f.tense1 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ Gender + Age"))
#' 
#' 
#' 
#' ####
#' ## Adjacent Categories Models
#' ####
#' 
#' ## Multivariate adjacent categories model, without response style, without explanatory variables
#' m.tense0 <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE))
#' m.tense0
#'
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1 <- multordRS(f.tense0, data = tenseness)
#' m.tense1
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2 <- multordRS(f.tense1, data = tenseness, control = ctrl.multordRS(XforRS = FALSE))
#' m.tense2
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3 <- multordRS(f.tense1, data = tenseness)
#' m.tense3
#' 
#' plot(m.tense3)
#' 
#' 
#' 
#' ####
#' ## Cumulative Models
#' ####
#' 
#' ## Multivariate cumulative model, without response style, without explanatory variables
#' m.tense0.cumul <- multordRS(f.tense0, data = tenseness, control = 
#'   ctrl.multordRS(RS = FALSE), model = "cumulative")
#' 
#' m.tense0.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1.cumul <- multordRS(f.tense0, data = tenseness, model = "cumulative")
#' m.tense1.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2.cumul <- multordRS(f.tense1, data = tenseness, 
#'   control = ctrl.multordRS(XforRS = FALSE), model = "cumulative")
#'   
#' m.tense2.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## with explanatory variables 
#' ## for location AND response style
#' m.tense3.cumul <- multordRS(f.tense1, data = tenseness, model = "cumulative")
#' m.tense3.cumul
#' 
#' plot(m.tense3.cumul)
#'}
NULL



#' German Longitudinal Election Study 2017 (GLES17)
#' 
#' Data from the German Longitudinal Election Study (GLES) from 2017 (Rossteutscher  et
#' al., 2017, https://doi.org/10.4232/1.12927). The GLES is a long-term study of the German electoral process.
#' It collects pre- and post-election data for several federal elections, the
#' data used here originate from the pre-election study for 2017.
#' 
#' @name GLES17
#' @docType data
#' @format A data frame containing data from the German Longitudinal Election Study with 2036 observations. 
#' The data contain socio-demographic information about the participants as well as their responses to items about specific political fears.
#' \describe{ 
#' \item{RefugeeCrisis}{How afraid are you due to the refugee crisis? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{ClimateChange}{How afraid are you due to the global climate change? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{Terrorism}{How afraid are you due to the international terrorism? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{Globalization}{How afraid are you due to the globalization? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{Turkey}{How afraid are you due to the political developments in Turkey? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{NuclearEnergy}{How afraid are you due to the use of nuclear energy? (Likert scale from 1 (not afraid at all) to 7 (very afraid))}
#' \item{Age}{Age in years} 
#' \item{Gender}{0: male, 1: female}
#' \item{EastWest}{0: West Germany, 1: East Germany}
#' \item{Abitur}{High School Diploma, 1: Abitur/A levels, 0: else} 
#' \item{Unemployment}{1: currently unemployed, 0: else} 
#' }
#' @references Rossteutscher, S., Schmitt-Beck, R., Schoen, H., Wessels, B., Wolf, C., Bieber, I., Stovsand, L.-C., Dietz, M., and Scherer, P. (2017). 
#' Pre-election cross section (GLES 2017). \emph{GESIS Data Archive, Cologne, ZA6800 Data file Version 2.0.0.}, \doi{10.4232/1.12927}.
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2021): Multivariate Ordinal Random Effects Models Including Subject and Group Specific Response Style Effects, 
#' \emph{Statistical Modelling}, \url{https://journals.sagepub.com/doi/10.1177/1471082X20978034}
#' @source
#' \url{https://gles-en.eu/} and 
#' \doi{10.4232/1.12927}
#' @keywords datasets
#' @examples
#' \donttest{
#' ###############################################################
#' ## Examples from Schauberger and Tutz (2020) 
#' ## Data from the German Longitudinal Election Study (GLES) 2017
#' ###############################################################
#' 
#' ####
#' ## Source: German Longitudinal Election Study 2017 
#' ## Rossteutscher et al. 2017, https://doi.org/10.4232/1.12927
#' ####
#' 
#' ## load GLES data
#' data(GLES17)
#' 
#' ## scale data
#' GLES17[,7:11] <- scale(GLES17[,7:11])
#' 
#' ## define formula
#' f.GLES <- as.formula(cbind(RefugeeCrisis, ClimateChange, Terrorism, 
#'                        Globalization, Turkey, NuclearEnergy) ~ 
#'                        Age + Gender + Unemployment + EastWest + Abitur)
#' 
#' ## fit adjacent categories model without and with response style parameters
#' m.GLES0 <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(RS = FALSE, cores = 6))
#' m.GLES <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(cores = 6))
#' 
#' m.GLES0
#' m.GLES
#' 
#' plot(m.GLES, main = "Adjacent categories model")
#' 
#' 
#' ## fit cumulative model without and with response style parameters (takes pretty long!!!)
#' m.GLES20 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6, RS = FALSE))
#' 
#' m.GLES2 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6))
#' 
#' m.GLES20
#' m.GLES2
#' 
#' plot(m.GLES2, main = "Cumulative model")
#' 
#'}
NULL

