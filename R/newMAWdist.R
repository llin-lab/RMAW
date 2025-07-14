#' Compute the signed Minimized Aggregated Wasserstein (MAW) Distance
#'
#' Calculates the signed MAW distance between two Gaussian Mixture Models (GMMs),
#' where the components' Gaussian parameters are provided as supports.
#' The distance is computed using Gaussian Wasserstein pairwise distances,
#' optimal transport, and cosine similarity-based adjustment.
#'
#' @param d Integer. The dimension of the Gaussian components.
#' @param supp1 Matrix. A \code{d x n1} matrix where each column is a Gaussian component (mean vector) of the first GMM.
#' @param supp2 Matrix. A \code{d x n2} matrix where each column is a Gaussian component (mean vector) of the second GMM (typically the barycenter).
#' @param w1 Numeric vector. A vector of length \code{n1} containing the mixture weights for \code{supp1}. Must sum to 1.
#' @param w2 Numeric vector. A vector of length \code{n2} containing the mixture weights for \code{supp2}. Must sum to 1.
#' @param method Character. The method to compute pairwise Gaussian Wasserstein distances. Default is \code{"Schur"}.
#'
#' @return A list with:
#' \describe{
#'   \item{dist}{Numeric. The minimized aggregated Wasserstein distance with cosine similarity adjustment.}
#'   \item{gammaij}{Matrix. The optimal transport plan (coupling matrix) between the two GMMs.}
#' }
#'
#' @details
#' The function first computes pairwise Gaussian Wasserstein distances between components of \code{supp1} and \code{supp2},
#' solves an optimal transport problem based on the weights \code{w1} and \code{w2}, and applies cosine similarity
#' between shifted vectors to modulate the final distance.
#'
#' @examples
#' # Example usage:
#' d <- 2
#' supp1 <- matrix(rnorm(2 * 3), nrow = 2)
#' supp2 <- matrix(rnorm(2 * 4), nrow = 2)
#' w1 <- rep(1/3, 3)
#' w2 <- rep(1/4, 4)
#' result <- newMawdist(d, supp1, supp2, w1, w2)
#' result$dist
#'
#' @export
newMawdist <- function(d, supp1, supp2, w1, w2, method = "Schur") {
  pairdist <- GaussWasserstein(d, supp1, supp2, method)
  result <- ot(pairdist, w1, w2)
  gammaij <- result$gammaij

  cos_similariy <- matrix(0, nrow = length(w1), ncol = length(w2))
  for (i in seq_along(w1)) {
    for (j in seq_along(w2)) {
      if (sum(supp1[1:d, i] - supp2[1:d, j]) != 0) {
        cos_similariy[i, j] <- lsa::cosine(supp1[1:d, i] - supp2[1:d, j], rep(1, d))
      }
    }
  }

  dist <- sum(pairdist * gammaij * cos_similariy)
  return(list(dist = dist, gammaij = gammaij))
}
