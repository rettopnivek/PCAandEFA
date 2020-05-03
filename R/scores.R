# The 'scores' function
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2020-05-02

# Table of contents
# 1) scores
# 2) Dispatch methods
#   2.1) scores.matrix
#   2.2) scores.PCA

###
### 1) Generic function 'score'
###

#' Computes Principal Component Scores
#'
#' A function to compute principal component
#' scores, either via a N x K matrix of N
#' standardized observations for K measures
#' and a K x K weights matrix, or by extracting
#' the scores from an object of class \code{PCA}.
#'
#' @param x Either a N x K matrix of N standardized
#'   observations for K measures, or an object of
#'   class \code{PCA}.
#' @param w An optional K x K weights matrix (the loadings
#'   returned after conducting PCA).
#' @param k An optional number indicating
#'   the number of component scores to retain.
#'
#' @return Description
#'
#' @examples
#' # Example of function
#'
#'
#' @export
scores = function( x, ... ) {
  UseMethod( "scores", x )
}

###
### 2) Dispatch methods
###

# 2.1)
#' @rdname scores
#' @export

scores.matrix = function( x, w, k = NULL ) {

  if ( is.null( k ) ) {
    k = ncol( x )
  }

  out = x %*% w[,1:k]

  return( out )
}

# 2.2)
#' @rdname scores
#' @export

scores.PCA = function( x, k = NULL ) {

  if ( is.null( k ) ) {
    k = length( names( x ) )
  }

  out = results$princomp[[1]]$scores[,1:k]

  return( out )
}

