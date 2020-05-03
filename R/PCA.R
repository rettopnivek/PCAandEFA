# The 'PCA' function
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you
# have any questions or comments
# Last updated 2020-02-22

# Table of contents
# 1) PCA
# 2) S3 methods
#   2.1) is.PCA
#   2.2) names.PCA
#   2.3) print.PCA
#   2.4) loadings.PCA
#   2.5) plot.PCA

###
### 1) PCA
###

#' Principal Component Analysis
#'
#' A convenience function that conducts a
#' prinicipal components analysis (PCA) on
#' the specified rows and columns of a
#' data frame using R's base
#' function \code{\link[stats]{princomp}}.
#'
#' Once a PCA has been conducted,
#' prinicipal component scores are
#' computed for all observations.
#'
#' @param df A data frame.
#' @param rows An optional logical vector matching in
#'   length to the number of rows in 'df' indicating
#'   which rows to use for the PCA.
#' @param columns An optional vector of column names
#'   in 'df' to use for the PCA.
#' @param patterns An optional vector of character
#'   strings, incomplete patterns to use to identify
#'   either columns to include in or columns to exclude
#'   from the PCA.
#' @param exclude Logical; if TRUE, the 'patterns' vector
#'   are strings to identify columns to exclude; otherwise
#'   the 'patterns' vector are strings to identify columns
#'   to include.
#' @param retain Description.
#'
#' @return An object of class 'PCA'.
#'
#' @examples
#' # Example of function
#' results = PCA( mtcars, patterns = c( 'p', 'r' ) )
#'
#' @export

PCA = function( df,
                rows = NULL,
                columns = NULL,
                patterns = NULL,
                exclude = T,
                retain = NULL ) {

  ### Prep data ###

  # If specific columns to include were
  # not specified
  if ( is.null( columns ) ) {

    # If no options specified, include
    # all columns in data frame
    if ( is.null( patterns ) ) {
      columns = colnames( df )
    } else {
      # If pattern completion is specified

      # Matrix of matches
      matches = sapply( patterns, function( string ) {
        out = grepl( string, colnames( df ), fixed = T )
        return( out )
      } )
      if ( !is.null( dim( matches ) ) ) {
        selected_columns = rowSums( matches ) > 0
      } else {
        selected_columns = matches
      }
      names( selected_columns ) = colnames( df )

      if ( exclude ) {
        # When patterns refer to columns to exclude
        selected_columns = !selected_columns
      }

      columns = colnames(df)[ selected_columns ]


    }

  }

  # By default use all rows in data frame
  if ( is.null( rows ) ) {
    rows = rep( TRUE, nrow( df ) )
  }

  # Extract data to be fitted and
  # convert to matrix
  X = as.matrix( df[ rows, columns ] )

  ### Prep variables for PCA ###

  # Number of observations
  N = nrow( X )
  # Number of measures
  K = ncol( X )

  # Means and standard deviations for
  # rescaling data
  M = colMeans( X )
  SD = apply( X, 2, sd )

  # Standardize input for PCA
  Z = scale( X, center = M, scale = SD )

  ### PCA results ###

  # Z = Pw

  # Singular value decomposition
  udv = svd( Z )
  d = udv$d # The singular values
  U = udv$u # The left singular vectors
  V = udv$v # The right singular vectors
  D = diag( d )
  # udv$d = The singlular values of Z
  # udv$u = Left singular vectors of Z
  # udv$v = Right singular vectors of Z

  # Z = U %*% D %*% t( V )

  # Use R's built-in function for PCA
  pca = princomp( Z, scores = T )

  # Extract eigenvalues and eigenvectors
  evv = eigen( var( Z ) )

  # Proportion of variance each component accounts for
  prop_of_var = pca$sdev^2 / sum( pca$sdev^2 )

  # Convert loadings to simple matrix
  L = matrix( pca$loadings, K, K )

  # Compute correlations between component scores and
  # loadings
  R = matrix( NA, K, K )
  rownames( R ) = columns
  colnames( R ) = paste0( 'Comp.', 1:K )
  for ( rw in 1:K ) {
    for ( cl in 1:K ) {
      R[rw,cl] = cor( Z[,rw], pca$scores[,cl] )
    }
  }

  if ( all( rows ) ) {
    pcs = pca$scores
    colnames( pcs ) = paste0( 'Comp.', 1:K )
  } else {
    # Take all columns
    NX = as.matrix( df[,columns] )
    pcs = scores( NX, L, K )
  }

  ### Output ###

  # List of outputs
  if ( FALSE ) {
    out = list(
      data = list(
        input = df,
        raw = X,
        standardized = Z,
        M = M, SD = SD
      ),
      princomp = pca,
      SD = pca$sdev,
      loadings = L,
      correlations = R,
      proportion = prop_of_var,
      rows = rows,
      columns = columns,
      K = K,
      N = N,
      scores = pcs
    )
  }
  out = list(
    data = list(
      input = df,
      raw = X,
      standardized = Z,
      M = M, SD = SD
    ),
    princomp = list( pca ),
    SD = pca$sdev,
    loadings = L,
    correlations = R,
    proportion = prop_of_var,
    K = K,
    N = N,
    measures = colnames( X ),
    eigenvalues = evv$values,
    eigenvectors = evv$vectors,
    SVD = list(
      d = d,
      U = U,
      D = D,
      V = V
    )
  )

  class( out ) = 'PCA'

  return( out )
}

###
### 2) S3 methods
###

# 2.1)
#' #' @rdname PCA
#' @export

is.PCA = function( x ) {
  inherits( x, 'PCA' )
}

# 2.2)
#' #' @rdname PCA
#' @export

names.PCA = function( x ) {
  x$measures
}

# 2.3)
#' #' @rdname PCA
#' @export

print.PCA = function( x, digits = 2 ) {

  tbl = data.frame(
    Component = paste0( 'Comp.', 1:( x$K ) ),
    Proportion.Variance = round( x$proportion, digits ),
    Cumulative.Proportion = round( cumsum( x$proportion ), digits ),
    SD = round( x$SD, digits ),
    Cut.off = ''
  )
  row.names( tbl ) = 1:nrow( tbl )

  print( tbl )
}

# 2.4)
#' #' @rdname PCA
#' @export

loadings.PCA = function( x ) {
  x$loadings
}

# 2.5)
#' #' @rdname PCA
#' @export

plot.PCA = function( x, type = 'scree',
                     lnSz = 2,
                     ptSz = 1.25,
                     axSz = 1.25,
                     axPos = -1.25,
                     lbSz = 1.25,
                     lbPos = 1.75,
                     pos_col = colors()[499],
                     neg_col = colors()[122],
                     R_interval = .1
                     ) {

  types = list(
    scree = c( 'Scree', 'scree',
               'Eigenvalues', 'eigenvalues',
               'Eigen', 'eigen' ),
    prop = c( 'Proportions', 'proportions',
              'Proportion', 'proportion',
              'Prop', 'prop' ),
    cumul = c( 'Cumulative', 'cumulative',
               'Cumul', 'cumul',
               'Cum', 'cum' ),
    heatmap = c( 'Heatmap', 'heatmap', 'Heat', 'heat' ),
    R = c( 'Correlations', 'correlations',
           'Correlation', 'correlation',
           'Corr', 'corr',
           'Cor', 'cor',
           'R', 'r' ),
    R2 = c( 'R-squared', 'r-squared',
            'Squared', 'squared',
            'R squared', 'r squared',
            'R2', 'r2' )
  )

  ### Scree plot ###

  if ( type %in% types$scree ) {

    K = x$K
    SD = x$SD

    xl = c( 0, K+1 )
    yl = c( 0, ceiling( SD[1] ) )

    # Blank plot
    plot( xl, yl,
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = ' ',
          xlab = ' ',
          bty = 'n' )

    # Grid lines
    segments( xl[1], 1,
              xl[2], 1,
              col = 'grey80',
              lwd = lnSz,
              lty = 2 )
    pos = unique( round( seq( 1, x$K, length.out = 5 ) ) )
    segments( pos, yl[1],
              pos, yl[2],
              lwd = lnSz, col = 'grey80' )

    # Plot eigenvalues
    lines( 1:K, SD, lwd = lnSz )
    points( 1:K, SD, pch = 19, cex = ptSz )

    # Axes
    axis( 2, round( seq( 0, yl[2], length.out = 5 ), 1 ),
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Standard deviations',
           side = 2, line = lbPos, cex = lbSz )
    axis( 1, pos,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Component',
           side = 1, line = lbPos, cex = lbSz )

    # Border
    segments( xl[ c( 1, 1, 1, 2 ) ],
              yl[ c( 1, 2, 1, 1 ) ],
              xl[ c( 2, 2, 1, 2 ) ],
              yl[ c( 1, 2, 2, 2 ) ],
              lwd = lnSz, col = 'black' )

  }

  ### Proportion of variance accounted for by each component ###

  if ( type %in% types$prop ) {

    K = x$K
    p = x$proportion

    xl = c( 0, K+1 )
    mx = 10*max( p )
    mx = ceiling( mx )
    yl = c( 0, mx/10 )

    # Blank plot
    plot( xl, yl,
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = ' ',
          xlab = ' ',
          bty = 'n' )

    # Grid lines
    pos = unique( round( seq( 1, x$K, length.out = 5 ) ) )
    segments( pos, yl[1],
              pos, yl[2],
              lwd = lnSz, col = 'grey80' )

    # Plot proportion of variance
    lines( 1:K, p, lwd = lnSz )
    points( 1:K, p, pch = 19, cex = ptSz )

    # Axes
    axis( 2, round( seq( 0, yl[2], length.out = 5 ), 2 ),
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Proportion of variance',
           side = 2, line = lbPos, cex = lbSz )
    axis( 1, pos,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Component',
           side = 1, line = lbPos, cex = lbSz )

    # Border
    segments( xl[ c( 1, 1, 1, 2 ) ],
              yl[ c( 1, 2, 1, 1 ) ],
              xl[ c( 2, 2, 1, 2 ) ],
              yl[ c( 1, 2, 2, 2 ) ],
              lwd = lnSz, col = 'black' )

  }

  ### Cumulative proportion of variance accounted for by each component ###

  if ( type %in% types$cumul ) {

    K = x$K
    p = cumsum( x$proportion )

    xl = c( 0, K+1 )
    mx = 10*max( p )
    mx = ceiling( mx )
    yl = c( 0, 1 )

    # Blank plot
    plot( xl, yl,
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = ' ',
          xlab = ' ',
          bty = 'n' )

    # Grid lines
    pos = unique( round( seq( 1, x$K, length.out = 5 ) ) )
    segments( pos, yl[1],
              pos, yl[2],
              lwd = lnSz, col = 'grey80' )

    segments( xl[1], seq( 0, 1, .1 ),
              xl[2], seq( 0, 1, .1 ),
              lwd = lnSz, col = 'grey80' )

    # Plot proportion of variance
    lines( 1:K, p, lwd = lnSz )
    points( 1:K, p, pch = 19, cex = ptSz )

    # Axes
    axis( 2, seq( 0, 1, .2 ),
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Cumulative proportion of variance',
           side = 2, line = lbPos, cex = lbSz )
    axis( 1, pos,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Component',
           side = 1, line = lbPos, cex = lbSz )

    # Border
    segments( xl[ c( 1, 1, 1, 2 ) ],
              yl[ c( 1, 2, 1, 1 ) ],
              xl[ c( 2, 2, 1, 2 ) ],
              yl[ c( 1, 2, 2, 2 ) ],
              lwd = lnSz, col = 'black' )

  }

  if ( type %in% types$heatmap ) {

    # Specify range of colors for the
    # positive correlation values

    # Specify the two colors
    color_0 = colors()[1] # White for 0
    color_1 = pos_col

    # Define a function that will create the gradient
    color_gradient_function = colorRampPalette( c( color_0, color_1 ) )

    pos_col_assign = seq( 0 + R_interval/2, 1 - R_interval/2, R_interval )

    # Generate a vector with a specified number of steps between the colors
    pos_colors <- color_gradient_function( length(pos_col_assign) + 1 )

    # Specify range of colors for the
    # negative correlation values

    # Specify the two colors
    color_0 = colors()[1] # White for 0
    color_1 = neg_col

    # Define a function that will create the gradient
    color_gradient_function <- colorRampPalette( c( color_1, color_0 ) )

    neg_col_assign = -pos_col_assign

    # Generate a vector with a specified number of steps between the colors
    neg_colors = color_gradient_function( length(pos_col_assign) + 1 )

    all_colors = c( neg_colors, pos_colors[-1] )
    all_assign = c( -1, rev( neg_col_assign ), pos_col_assign, 1 )

    match_color = function( R, all_assign, all_colors ) {

      index = rep( 0, 2 )

      index = max( which( all_assign <= R ) )
      out = all_colors[ index[1] ]

      return( out )
    }

    # Number of components
    K = x$K

    # Grid
    xy = expand.grid( 1:K, 1:K )

    # Create  blank plot
    xl = c( 0, K+2 )
    yl = c( 0, K+2 )

    plot( xl, yl,
          type = 'n',
          xaxt = 'n',
          yaxt = 'n',
          ylab = ' ',
          xlab = ' ',
          bty = 'n' )

    # Loop over cells
    xa = rep( 0, 4 )
    ya = rep( 0, 4 )

    for ( j in 1:nrow( xy ) ) {

      # xy[1,1:2] = <1,1>
      # xy[2,1:2] = <2,1>

      # Matrix starts at <1,1> top, left

      # Plot <1,1> is matrix <30,1>
      # Plot <2,1> is matrix <29,1>

      xa[1:2] = xy[j,2] - 1
      xa[3:4] = xy[j,2]
      ya[1] = (K+1) - (xy[j,1] - 1)
      ya[2:3] = (K+1) - (xy[j,1])
      ya[4] = (K+1) - (xy[j,1] - 1)

      clr = match_color( x$correlations[ xy[j,1], xy[j,2] ],
                         all_assign,
                         all_colors )
      polygon( xa, ya, border = NA, col = clr )

    }

    # Add color key
    col_key_y = seq( 1, K+1, length.out = length( all_colors ) + 1 )


    for ( j in 1:( length( col_key_y ) - 1 ) ) {

      xa[1:2] = K + 1
      xa[3:4] = K + 2
      ya[1] = col_key_y[j]
      ya[2:3] = col_key_y[j+1]
      ya[4] = col_key_y[j]

      clr = all_colors[j]

      polygon( xa, ya, col = clr, border = NA )

    }

    segments( 0:K, 1, 0:K, K+1, col = 'grey80' )

    # Plotting border
    segments( 0, 1, 0, K+1, lwd = lnSz )
    segments( 0, 1, K, 1, lwd = lnSz )
    segments( K, 1, K, K+1, lwd = lnSz )
    segments( K, K+1, 0, K+1, lwd = lnSz )


    axis( 4,
          seq( 1, K+1, length.out = 4*2 + 1 ),
          c( seq( -1, 0, .25 ),
             seq( .25, 1, .25 ) ),
          tick = F, line = axPos, cex.axis = axSz*.75 )

    pos = unique( round( seq( 1, x$K, length.out = 5 ) ) )
    axis( 1, pos - .5, pos,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( 'Component',
           side = 1, line = lbPos, cex = lbSz )

    mtext( 'Standardized variable',
           side = 2, line = lbPos, cex = lbSz )

    mtext( 'Correlations between observed and component scores',
           side = 3, line = lbPos, cex = lbSz )

  }

  if ( type %in% types$R ) {

  }

  if ( type %in% types$R2 ) {

  }

}

