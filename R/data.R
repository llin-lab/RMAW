#' GMM fitting component data coming from mouse atlas 
#'
#'
#' @format A list containing the GMM info of 4 subjects onto 2 dimensions:
#' \describe{
#' \item{ww}{Concatenated weights of each component across all subjects}
#' \item{supp}{Mean and Variance info of each component across all subjects}
#' \item{stride}{Number of components for each GMM corresponding to each subject}
#' \item{m}{Targeted number of components for barycenter computation}
#' }
"mouse_2"


#' GMM fitting component data coming from COVID study
#'
#'
#' @format A list containing the GMM info of 3 subjects onto 2 dimensions:
#' \describe{
#' \item{ww}{Concatenated weights of each component between both subjects}
#' \item{supp}{Mean and Variance info of each component between the two subjects}
#' \item{stride}{Number of components for each GMM corresponding to each subject}
#' \item{m}{Targeted number of components for barycenter computation}
#' \item{d}{Number of dimensions}
#' }
"mouse_2"