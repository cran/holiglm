#' Bike Sharing Dataset
#'
#' This data set contains the daily count of rented bikes from the the Capital Bikeshare system
#' in Washington D.C., USA, for the years 2011 and 2012. The dataset is already prepared
#' (correct types + factor encodings) for model building.
#'
#' 
#' @docType data
#' @keywords datasets
#' @format A data.frame of dimension 731 x 12 containing daily data related
#'	to related bikes.
#' 
#' \describe{
#'  \item{dteday}{a date vector giving the date of the rental.}
#'  \item{season}{a factor with levels `spring`, `summer`, `fall` and `winter`.}
#'  \item{year}{a factor with levels `2011` and `2012`.}
#'  \item{mnth}{a factor with levels `Jan`, `Feb`, `Mar`, `Apr`, `May`, `Jun`, `Jul`, `Aug`, `Sep`, `Oct`, `Nov` and `Dec`.}
#'  \item{holiday}{a boolean vector indicating if the day is a holiday.}
#'  \item{weathersit}{a factor with levels `good`, `neutral`, `bad` and `very bad` giving the weather situation.}
#'  \item{temp}{a numeric vector containing max-normalized temperature in Celsius with 41 as maximum.}
#'  \item{atemp}{a numeric vector containing max-normalized feeling temperature in Celsius with 50 as maximum.}
#'  \item{hum}{a numeric vector containing max-normalized humidity with 100 as maximum.}
#'  \item{windspeed}{a numeric vector containing max-normalized windspeed with 67 as maximum.}
#'  \item{cnt}{an integer vector containing counts of rented bikes.}
#' }
#' 
#' @references Fanaee-T, Hadi. (2013). Bike Sharing Dataset. UCI Machine Learning Repository.
#'
#' @source \url{https://www.ics.uci.edu/~mlearn/MLRepository.html}
#' @name bike
#' @examples
#' data("bike")
#' hglm(formula = cnt ~ ., data=bike, family="poisson")
#' 
NULL
