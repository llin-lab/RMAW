#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Solve complex optimization problems,
#' especailly linear programs (LPs) constrained optimization.
#'
#'
#' @param c: Objective function coefficients (vector)
#' @param A: Constraints matrix
#' @param blc: Lower bounds on constraints
#' @param buc: Upper bounds on constraints
#' @param blx: Lower bounds on variables
#' @param bux: Upper bounds on variables
#'
#' @return Returns a list that contains MAW distance between two GMMs as well as
#' the optimal transport between them.
#'
#'
#' @import lpSolve
#' @export
#' @concept MAW ot
#'
#' @examples
#' \dontrun{
#' library(RMAW)
#' system.file("extdata", "mouse_2.mat", package = "myfirstpackage")
#'
#' # mouse_2 is a an example data file containing 8 GMMs mean and variance
#' stride <- mouse_2$stride
#' supp1 <- mouse_2$supp
#' result <- Mawdist(d = 2,
#'                   supp1 = mouse$supp[,1:stride[1]],
#'                   supp2 = mouse$supp[,(stride[1]+1):(stride[1]+stride[2])],
#'                   w1 = mouse$ww[1:stride[1]],
#'                   w2 = mouse$ww[(stride[1]+1):(stride[1]+stride[2])])
#'
#'
#' }
#'


opt <- function(c, A, blc, buc, blx, bux) {

  #
  #
  #

  # Set up the direction of the optimization (minimize)
  direction <- "min"

  # Define the constraint directions
  const_dir <- rep("=", length(blc))  # For equality constraints

  # If there are inequalities in buc, adjust constraint directions accordingly
  for (i in 1:length(buc)) {
    if (buc[i] < Inf) {
      const_dir[i] <- "<="
    }
    if (blc[i] > -Inf) {
      const_dir[i] <- ">="
    }
  }

  # Solve the linear program
  result <- lpSolve::lp(direction = direction,
               objective.in = c,
               const.mat = A,
               const.dir = const_dir,
               const.rhs = pmax(blc, buc), # Use the tighter bound as rhs
               lower = blx,
               upper = bux,
               all.int = FALSE) # If integer variables are not required

  # Extract results
  if (result$status == 0) {
    x <- result$solution
    obj_val <- result$objval
  } else {
    stop("Optimization did not converge.")
  }

  # Return solution and objective value
  return(list(xx = x, obj_val = obj_val))
}





