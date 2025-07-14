#' @include MAWdist.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find barycenter across GMMs
#'
#' Find the barycenter across a set of GMMs.
#' The barycenter can later be used to integrate the GMMs.
#'
#' The main steps of this procedure are outlined below. For a more detailed
#' description of the methodology, please see Lin L, Shi W, Ye J, Li J. (2023):
#' \doi{10.1111/biom.13630}
#'
#'
#'
#' @param stride An integer vector of the form (J_1,...,J_N) that indicates
#' the number of components of each GMMs.
#' @param instanceW A vector of 1 with length of stride
#' @param supp A matrix that stores the mean and covariance for all the GMMs.
#' More specifically, the dimension is (d+d*d)x(J_1+J_2+...+J_N),
#' where d is the dimension of the original data.
#' For each column of supp, the first d rows store the mean vector,
#' and the last d*d rows store the vectorized covariance matrix for each component.
#' @param w A joined vector of weights across GMMs
#' @param c0 GMM barycenter initialization. Usually, it can be set the same as the first GMM
#' or SCT
#' @param options list controlling pre-set parameters. \itemize{'tau'} default is 10.
#' \itemize{'badmm_tol'} default is 1e-4. \itemize{'badmm_max_iters'} default is 2000.
#' \itemize{'badmm_rho'} default is 10*mean(median(C))
#' where C is the total Wasserstain distance across all GMMs.
#'
#' @return Returns a list that contains MAW barycneter as well as the optimal transport
#' between barycenter and each GMM.
#'
#' @references Lin L, Shi W, Ye J, Li J. (2023)
#' Multisource single-cell data integration by MAW barycenter for Gaussian mixture models.
#' Biometrics, 79, 866–877. https://doi.org/10.1111/biom.13630
#'
#' @importFrom pracma meshgrid repmat pdist2 eps
#' @importFrom Matrix sparseMatrix Schur
#'
#' @export
#' @concept MAW
#'
#' @examples
#' \dontrun{
#' library(RMAW)
#' data("mouse_2")
#'
#' # mouse_2 is a an example data file containing 8 GMMs mean and variance
#' m <- 15 # user should specify the number of components for MAW barycenter
#' stride <- t(mouse_2$stride)
#' mouse_2$ww <- t(mouse_2$ww)
#' N <- length(stride)
#' instanceW <- matrix(1, nrow = 1, ncol = N)
#' c0 <- list()
#' c0$supp <- mouse_2$supp[, 1:m]
#' c0$w <- mouse_2$ww[1:m] / sum(mouse_2$ww[1:m])
# set the number of iterations and badmm_rho (no need to change)
#' options <- list()
#' options$badmm_max_iters <- 1000
#' options$badmm_rho <- 10
#'
#' result <- centroid_sphBregman_GMM(stride, instanceW, mouse_2$supp, mouse_2$ww, c0, options)
#'
#' }
#'


centroid_sphBregman_GMM <- function(stride, instanceW, supp, w, c0, options){
  # The algorithmic prototype of Wasserstein Barycenter using Bregman ADMM
  #
  # The Matlab code has been created by Jianbo Ye (jxy198 [AT] ist.psu.edu).
  # and R code is tranlated from Matlab code.

  d <- floor(sqrt(dim(supp)[1]))
  n <- length(stride)
  m <- length(w)
  posvec <- c(1,cumsum(stride)+1)

  if (is.null(c0)) {
    stop('Please give a GMM barycenter initialization.')
    c <- centroid_init(stride, supp, w, options)
  }
  else {
    c <- c0
  }
  support_size <- length(c$w)
  #load(['cstart' num2str(n) '-' num2str(support_size) '.mat']);

  X <- matrix(0, nrow = support_size, ncol = m)
  Y <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
  Z <- X
  spIDX_rows <- matrix(0, nrow = support_size * m, ncol = 1)
  spIDX_cols <- matrix(0, nrow = support_size * m, ncol = 1)
  for (i in 1:n) {
    temp <- pracma::meshgrid((i-1)*support_size + (1:support_size), posvec[i]:(posvec[i+1]-1))
    xx <- temp$X
    yy <- temp$Y
    ii <- support_size*(posvec[i]-1) + (1:(support_size*stride[i]))
    spIDX_rows[ii,1] <- t(xx)
    spIDX_cols[ii,1] = t(yy)
  }
  spIDX <- pracma::repmat(diag(x = 1, nrow = support_size, ncol = support_size), 1, n)

  # initialization
  suppW <- matrix(0, nrow = 1, ncol = m)
  for (i in 1:n) {
    suppW[posvec[i]:(posvec[i+1]-1)] <- instanceW[i]
    Z[,posvec[i]:(posvec[i+1]-1)] <- 1 / (support_size * stride[i])
  }
  if (d > 1) {
    C <- pracma::pdist2(t(c$supp[1:d,]), t(supp[1:d,]))^2
  }
  else {
    C <- pracma::pdist2(t(t(c$supp[1:d,])), t(t(supp[1:d,])))^2
  }
  C <- C + gaussian_wd(c$supp[(d+1):nrow(c$supp),], supp[(d+1):nrow(supp),])
  nIter <- 2000
  if (!is.null(options$badmm_max_iters)) {
    nIter <- options$badmm_max_iters
  }

  if (!is.null(options$badmm_rho)) {
    rho <- options$badmm_rho*mean(apply(C,2,median))
  }
  else {
    rho <- 10 * mean(apply(C,2,median))
  }

  C <- C * as.vector(suppW)
  rho <- rho * median(instanceW)

  if (!is.null(options$tau)) {
    tau <- options$tau
  }

  else {
    tau <- 10
  }


  if (!is.null(options$badmm_tol)) {
    badmm_tol <- options$badmm_tol
  }
  else {
    badmm_tol <- 1e-4
  }


  for (iter in 1:nIter) {
    # update X
    X <- Z * exp((C+Y)/(-rho)) + pracma::eps(x = 1)
    X <-  t(t(X) * (as.vector(w) / colSums(X)))

    # update Z
    Z0 <- Z
    Z <- X * exp(Y/rho) + pracma::eps(x = 1)
    spZ <- Matrix::sparseMatrix(i = as.vector(spIDX_rows), j = as.vector(spIDX_cols), x = as.vector(Z), dims = c(support_size * n, m))
    tmp <- rowSums(as.matrix(spZ))
    tmp <- matrix(tmp, nrow = support_size, ncol = n)
    dg <-  1/tmp * as.vector(c$w)
    dg <- Matrix::sparseMatrix(i = 1:(support_size*n), j = 1:(support_size*n), x = as.vector(dg))
    Z <- as.matrix(spIDX %*% dg %*% spZ)

    # update Y
    Y = Y + rho * (X - Z);

    # update c.w
    tmp <- t(t(tmp) * 1/colSums(tmp))
    sumW <- rowSums(sqrt(tmp))^2 # (R2)
    #sumW = sum(tmp,2)' # (R1)
    c$w <- sumW / sum(sumW)
    #c.w = Fisher_Rao_center(tmp')

    # update c.supp and compute C (lazy)
    if (iter %% tau == 0 & is.null(options$support_points)) {
      tmpX <- t(t(X) * as.vector(suppW))
      c$supp[1:d,] <- (supp[1:d,] %*% t(tmpX)) * 1 / pracma::repmat(rowSums(tmpX), d, 1)
      c$supp[(d+1):nrow(c$supp),] <- gaussian_mean(supp[(d+1):nrow(supp),], tmpX, c$supp[(d+1):nrow(c$supp),])
      if (d > 1) {
        C <- pracma::pdist2(t(c$supp[1:d,]), t(supp[1:d,]))^2
      }
      else {
        C <- pracma::pdist2(t(t(c$supp[1:d,])), t(t(supp[1:d,])))^2
      }
      C <- C + gaussian_wd(c$supp[(d+1):nrow(c$supp),], supp[(d+1):nrow(supp),])
      C = t(t(C) * as.vector(suppW))
    }



    # The constraint X=Z are not necessarily strongly enforced
    # during the update of w, which makes it suitable to reset
    # lagrangian multipler after a few iterations
    # if (mod(iter, 100) == 0)
      #          Y(:,:) = 0;
      #           if primres > 10*dualres
      #             rho = 2 * rho;
      #             fprintf(' *2');
      #           elseif 10*primres < dualres
      #             rho = rho / 2;
      #             fprintf(' /2');
      #           end
      #end

    # output
    if (iter %% 100 == 0 | iter == 10) {
      primres <- norm(X - Z, type = "F") / norm(Z, type = "F")
      dualres <- norm(Z - Z0, type = "F") / norm(Z, type = "F")
      cat(sprintf('\t %d %f %f %f \n', iter, sum(t(t(C * X) * as.vector(suppW))) / sum(instanceW),
                  primres, dualres))
      if (sqrt(dualres * primres) < badmm_tol) {
        break
      }
    }

  }
  return(list(c = c, X = X))
}


#' Update for the GMM variance part of the MAW barycenter
#
# @param object.1 First Seurat object to merge
# @param object.2 Second Seurat object to merge
# @param reduction Name of DimReduc to use. Must be an SVD-based DimReduc (eg, PCA or LSI)
# so that the loadings can be used to project new embeddings. Must be present
# in both input objects, with a substantial overlap in the features use to construct
# the SVDs.
# @param dims dimensions used for rpca
# @param projected.name Name to store projected SVDs under (eg, "projectedpca")
# @param features Features to use. Will subset the SVD loadings to use these features
# before performing projection. Typically uses the anchor.features for integration.
# @param do.center Center projected values (subtract mean)
# @param do.scale Scale projected values (divide by SD)
# @param slot Name of slot to pull data from. Should be scale.data for PCA and data for LSI
# @param verbose Display messages
# @return Returns a merged Seurat object with two projected SVDs (object.1 -> object.2, object.2 -> object.1)
# and a merged SVD (needed for within-dataset neighbors)


gaussian_mean <- function(V, w, Sigma0) {
  # size(V) = [d*d, n]
  # size(w) = [m, n]
  # size(Sigma) = [d*d, m]
  if (is.vector(V)) V <-  matrix(V, nrow = 1)
  d <- sqrt(dim(V)[1])
  stopifnot(dim(V)[1] == as.integer(d*d))
  n <- dim(V)[2]
  stopifnot(n == dim(w)[2])
  m <- dim(w)[1]
  w <- w * (1 / rowSums(w))

  Sigma <- Sigma0
  if (d > 1) {
    Sigma <- array(Sigma, dim = c(d, d, m))
    old_Sigma <- array(0, c(d, d, m))
    V <- array(V, dim = c(d, d, n))
    while (max(abs(old_Sigma - Sigma)) > 1e-5 * max(abs(Sigma))) {
      old_Sigma <- Sigma
      Sigma <- array(0, dim = c(d, d, m))

      for (j in 1:m) {
        mem <- sqrtm_old(old_Sigma[,,j])
        for (i in 1:n) {
          Sigma[,,j] <- Sigma[,,j] + w[j,i] * sqrtm_old(mem %*% V[,,i] %*% mem)
        }
        # Sigma(:,:,j) = sqrtm_batch_it(V, w(j,:), mem);
      }
    }
  }
  else if (d == 1) {
    V <- as.vector(V)
    Sigma = (w %*% sqrt(V))^2
  }

  Sigma = matrix(Sigma, nrow = d*d, ncol = m)

  return(Sigma)
}


#' Compute the Wasserstain distance between two GMMs' components on the variance part
#'
#' @param V1 First GMM variance part
#' @param V2 Second GMM variance part
#' @return Returns matrix that contains the Wasserstain distance between two GMMs' components
#' (V1 -> V2 with rows correspond to components in the first GMM
#' and columns correspond to components in the second GMM)

gaussian_wd <- function(V1, V2) {
  # size(V1) = [d*d, n1]
  # size(V2) = [d*d, n2]
  # size(D) = [n1, n2]
  if (is.vector(V1)) V1 <- matrix(V1, nrow = 1)
  if (is.vector(V2)) V2 <- matrix(V2, nrow = 1)
  d <- sqrt(dim(V1)[1])
  stopifnot(dim(V1)[1] == as.integer(d*d))
  stopifnot(dim(V2)[1] == as.integer(d*d))
  n1 <- dim(V1)[2]
  n2 <- dim(V2)[2]

  D <- matrix(0, nrow = n1, ncol = n2)
  V1 <- array(V1, dim = c(d, d, n1))
  V2 <- array(V2, dim = c(d, d, n2))

  for (i in 1:n1) {
    if (d > 1) {
      d1 <- sqrtm_old(V1[,,i])
      t1 <- sum(diag(V1[,,i]))
    }
    else {
      d1 <- sqrt(V1[,,i])
      t1 <- V1[,,i]
    }

    for (j in 1:n2) {
      if (d > 1) {
        v2 <- V2[,,j]
        d2 <- sqrtm_old(d1 %*% v2 %*% d1)
        D[i,j] = t1 + sum(diag(v2)) - 2 * sum(diag(d2))
      }
      else {
        v2 = V2[j]
        d2 = sqrt(d1 * v2 * d1)
        D[i,j] = t1 + v2 - 2 * d2
      }
    }
    # D(i,:) = t1*ones(1,n2) + sqrtm_batch_ud(V2,d1);
  }
  return(D)
}

#' Compute the square root of a matrix
#'
#' @param A the matrix that is about to be computed for square root.
#' Usually, it should be a positive definite symmetric matrix, e.g. covariance matrix
#' @return Returns the square root of the input matrix

sqrtm_old <- function(A) {
  n <- dim(A)[1]
  result <- Matrix::Schur(A)
  Q <- result$Q
  Tr <- result$`T`
  R <- matrix(0, nrow = n, ncol = n)

  if (sum(diag(Tr) == 0) <= 1) {
    if (sum(diag(Tr) < 0) == 0) {
      # Compute the square root of an upper triangular matrix with at most 1 zero element on diagonal
      for (j in 1:n) {
        R[j,j] <- sqrt(Tr[j,j])
        if (j > 1) {
          for (i in (j-1):1) {
            s <- 0
            if ((i+1) <= (j-1)) {
              for (k in (i+1):(j-1)) {
                s <- s + R[i,k] * R[k,j]
              }
            }
            if (R[i,i] + R[j,j] == 0) R[i,j] <- (Tr[i,j] - s) / (R[i,i] + R[j,j] + pracma::eps(x = 1))
            else R[i,j] <- (Tr[i,j] - s) / (R[i,i] + R[j,j])
          }
        }

      }
      X <- Q %*% R %*% t(Q)
    }
    else {
      # Compute the square root of an upper triangular matrix with at most 1 zero element but some small negative values (close to 0) on diagonal.
      #if (Tr[j,j] < 0 & Tr[j,j] > -1e-4) {
        #warning("Due to the computational limit, the covariance matrix is positive semi-definite and the eigenvalues are negative in format.")
      #}
      for (j in 1:n) {
        R[j,j] <- sqrt(as.complex(Tr[j,j]))
        if (j > 1) {
          for (i in (j-1):1) {
            s <- 0
            if ((i+1) <= (j-1)) {
              for (k in (i+1):(j-1)) {
                s <- s + R[i,k] * R[k,j]
              }
            }
            R[i,j] <- (Tr[i,j] - s) / (R[i,i] + R[j,j])
          }
        }

      }
      X <- Re(Q %*% R %*% t(Q))
    }
  }
  else stop(paste("In sqrm_old function, the symmetric matrix may be negative definite,",
  "or the upper triangular matrix in Schur decomposition has at least 2 zero elements on diagonal,",
  "which may indicate that the square root of the matrix may not exist."))

  return(X)
}


















