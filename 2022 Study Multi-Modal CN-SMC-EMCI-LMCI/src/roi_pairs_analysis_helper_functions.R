#
# roi_pairs_analysis_helper_functions.R
#

#
# extract_anchor_roi
#
extract_anchor_roi <- function( roi_name ) {
  
  tmp <- paste( unlist( str_split( roi_name, pattern = "\\." ) )[1:5], collapse = "." )
  
} # extract_anchor_roi <- function( roi_name ) {

#
# Read in the name list files
#
get_roi_time_series <- function( rid_list, roi_list, exp_dir, ts_len = 187 ) {
  # roi_list <- roi_list_anchors_fname;  k <- 1; irid <- 1; iroi <- 1; ts_len = 187; exp_dir <- paste0( exp_dir, adni3_protocol )
  

  roi_list <- read.table( roi_list, header = FALSE, sep = " ", stringsAsFactors = FALSE )
  names( roi_list ) <- c( "roi" )
  
  df_res <- list()
  rid_missing <- NULL
  
  tryCatch (
    {
      k <- 1
      for ( irid in 1:nrow( rid_list ) ) {
        
        for ( iroi in 1:nrow( roi_list ) ) {
          
          ts_name <- paste( rid_list[ irid, ] %>% select( - dxbl ) %>% unlist(), collapse = ".", sep = "" )
          ts_name <- gsub( " ", "", paste( exp_dir, ts_name, "/timeCourse.", roi_list$roi[iroi], ".txt", sep = "" ) )
          if ( file.exists( ts_name ) ) {
            tmp <- read.table( ts_name, header = FALSE, sep = " ", colClasses = c( "numeric" ), stringsAsFactors = FALSE ) %>% t()
            rownames( tmp ) <- NULL
            if ( length(tmp) == ts_len ) {
              df_res[[ k ]] <-
                data.frame( rid = as.character( rid_list$rid[irid] ),
                            dx = as.character( rid_list$dx[irid] ),
                            dxbl = as.character( rid_list$dxbl[irid] ),
                            roi = as.character( roi_list$roi[iroi] ),
                            tmp, stringsAsFactors = FALSE )
              k <- k + 1
            } # if ( nrow(tmp) == ts_len ) {
          } else {
            rid_missing <- c( rid_missing, k )
          } # if ( file.exists( ) ) {
          
        } # for ( iroi in roi_list ) {
        
      } # for ( irid in nrow( rid_list ) ) {
      
      df_res <- rbindlist( df_res )
      cat( " get_roi_time_series: Expected RID and ROI Count = ", nrow(rid_list), nrow(roi_list),
           " ; Observed = ", uniqueN( df_res$rid ), uniqueN( df_res$roi ),
           " ; Time Series Length Expected & Observed = ", ts_len, ncol( df_res ) - 2, "\n" )
      
    },
    error = function( x ) {
      message( " ** get_roi_time_series: ERROR ISSUED: ", x, "\n" );
    } # error = function( x ) {
  ) # tryCatch
  
  return( df_res )
  
} # get_roi_time_series <- function( roi_anchor_list ) {

#
# permute_labels
#
permute_labels <- function ( x ) {
  
  x[ sample( length(x), length(x), replace=FALSE ) ] <- x
  return( x )
  
} # permute_labels <- function ( x ) {

#
# t_test_with_cohens_d
#
t_test_with_cohens_d <- function( x, labels ) {
  # UnitTest: val = df[ 1, ]; val = df[ itest, ];
  
  res = rep( NA, 4 )
    
    if ( ( sum( is.na( x ) ) < length( x ) )  & ( var( x, na.rm = TRUE ) > 0 ) ) {
      
      t_res <- t.test( x ~ labels )
      d_res <- effsize::cohen.d( x ~ labels )
      res = c( t_res$statistic, t_res$p.value, t_res$parameter, d_res$estimate );
      
    } # if ( ) {
    
  names( res ) = c( "t_stat", "t_p", "t_df", "d" );
  
  return( res );
  
} # t_test_with_cohens_d = function( val, labels ) {
