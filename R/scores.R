# The 'scores' function
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2020-02-22

# Table of contents
# 1) scores
# 2) Dispatch methods
#   2.1) scores.matrix

###
### 1) Generic function 'score'
###

#' Title Case
#'
#' Description.
#'
#' @param variable Description.
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
#' #' @rdname scores
#' @export

scores.matrix = function( x, w, k ) {

  out = x %*% w[,1:k]

  return( out )
}
