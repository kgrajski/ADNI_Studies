#
# ADNI3.HelperFunctions.R
#

#
# Helper Functions
#

#
# ClassCheck
#
ClassCheck = function ( df ) {
  
  numFields = dim( df )[2];
  
  cOut = NULL;
  
  for ( iField in 1:numFields ) {
    
    cOut = c( cOut, class( df[ 1, iField ] ) );
    #cat( " ClassCheck: Field =", names( df )[ iField ], class( df[ , names( df )[ iField ] ] ), "\n" );
    
  } # for ( iField in 1:numFields ) {
  
  return( cOut );
  
} # ClassCheck = function ( df ) {

ConvertDateType0 = function( x ) {
  #  Convert: YYYY-MM-DD to MM/DD/YYYY
  # x = adnimerge_tmp$EXAMDATE;
  
  y = rep( "", length( x ) );
  for ( i in 1:length( x ) ) {
    
    itmp = unlist( str_split( x[ i ], pattern = "-" ) );
    y[ i ] = paste( itmp[2], itmp[3], itmp[1], sep = "/" );
    
  } # for ( i in 1:nrow( x ) ) {
  
  return( y );
  
} # ConvertDateType0 = function( ) {

ConvertDateType1 = function( x ) {
  #  Convert: YYYYMMDD:XXXXXX to MM/DD/YYYY
  # x = df_study$SERIES_DATE.mprage;
  
  y = rep( "", length( x ) );
  for ( i in 1:length( x ) ) {
    
    itmp = unlist( str_split( x[ i ], pattern = ":" ) )[1];
    y[ i ] = paste( substr( itmp, 7, 8 ), substr( itmp, 5, 6 ), substr( itmp, 1, 4 ), sep = "/" );
    
  } # for ( i in 1:nrow( x ) ) {
  
  return( y );
  
} # ConvertDateType0 = function( ) {

#
#	DoANOVA
#
if ( 0 ) {
  
  DoANOVA = function ( var.CN, var.MCI, var.Dementia ) {
    # Testing: var.CN = f.CN$AgeOnSeries; var.MCI = f.MCI$AgeOnSeries; var.Dementia=f.Dementia$AgeOnSeries
    
    z.Var = c ( var.CN, var.MCI, var.Dementia );
    z.Lab = c ( rep("CN", length(var.CN)), rep("MCI", length(var.MCI)), rep("ZD", length(var.Dementia)) );
    zSet = data.frame ( z.Var, z.Lab );
    return ( aov ( z.Var ~ z.Lab, data = zSet ) );
    
  } # DoANOVA = function ( var.CN, var.MCI, var.Dementia ) {
  
} # if ( 0 )

#
#	DoANOVA
#
DoANOVA = function ( vlist ) {
  
  num_levels = length( vlist );
  z.Var = NULL;
  z.Lab = NULL;
  z.Mean = rep( NA, num_levels ); names( z.Mean ) = names( vlist );
  z.SD = rep( NA, num_levels ); names( z.Mean ) = names( vlist );
  z.Min = rep( NA, num_levels );  names( z.Mean ) = names( vlist );
  z.Max = rep( NA, num_levels );  names( z.Mean ) = names( vlist );
  z.Median = rep( NA, num_levels );  names( z.Mean ) = names( vlist );
  for ( i in 1:num_levels ) {
    
    x = unlist( vlist[[ names( vlist )[ i ] ]] );
    z.Mean[ i ] = mean( x );
    z.SD[ i ] = sqrt( var( x ) );
    z.Min[ i ] = min( x );
    z.Max[ i ] = max( x );
    z.Median[ i ] = median( x );
    
    z.Var = c( z.Var, x );
    z.Lab = c( z.Lab, rep( names( vlist )[i], nrow( vlist[[ names( vlist )[ i ] ]] ) ) );
    
  } # for ( i in 1:num_levels ) {
  
  zSet = data.frame ( z.Var, z.Lab );
  res_anova = aov ( z.Var ~ z.Lab, data = zSet );
  
  res = list( res_anova = res_anova, mean = z.Mean, min = z.Min, max = z.Max, median = z.Median, sd = z.SD );
  
  return ( res );
  
} # DoANOVA = function ( var.CN, var.MCI, var.Dementia ) {

#
#	Do a Paired T Test
#
DoPairedTTest = function ( zScores, zLabels ) {
  return ( pairwise.t.test ( zScores, zLabels ) );
} # DoPairedTTest = function ( ) {

#
#	DoPairWiseW
#
DoPairWiseW = function ( var.CN, var.MCI, var.Dementia ) {
  # Testing: var.CN = f.CN$AgeOnSeries; var.MCI = f.MCI$AgeOnSeries; var.Dementia=f.Dementia$AgeOnSeries
  alternative = "two.sided";
  paired = FALSE;
  t.cn.v.mci = as.numeric(unlist(wilcox.test ( var.MCI, var.CN, alternative=alternative, paired=paired ))["p.value"]) 
  t.cn.v.zd = as.numeric(unlist(wilcox.test ( var.Dementia, var.CN, alternative=alternative, paired=paired  ))["p.value"]) 
  t.mci.v.zd = as.numeric(unlist(wilcox.test ( var.Dementia, var.MCI, alternative=alternative, paired=paired ))["p.value"]) 
  return ( c(t.cn.v.mci, t.cn.v.zd, t.mci.v.zd ) );
  
} # DoPairWiseT = function ( var.CN, var.MCI, var.Dementia ) {

#
#	DoPairWiseT
#
DoPairWiseT = function ( var.CN, var.MCI, var.Dementia ) {
  # Testing: var.CN = f.CN$AgeOnSeries; var.MCI = f.MCI$AgeOnSeries; var.Dementia=f.Dementia$AgeOnSeries
  alternative = "two.sided";
  paired = FALSE;
  var.equal = FALSE;
  t.cn.v.mci = as.numeric(unlist(t.test ( var.MCI, var.CN, alternative=alternative, paired=paired, var.equal=var.equal ))["p.value"]) 
  t.cn.v.zd = as.numeric(unlist(t.test ( var.Dementia, var.CN, alternative=alternative, paired=paired, var.equal=var.equal ))["p.value"]) 
  t.mci.v.zd = as.numeric(unlist(t.test ( var.Dementia, var.MCI, alternative=alternative, paired=paired, var.equal=var.equal ))["p.value"]) 
  return ( c(t.cn.v.mci, t.cn.v.zd, t.mci.v.zd ) );
  
} # DoPairWiseT = function ( var.CN, var.MCI, var.Dementia ) {

#
#	FDRCorrectedPValue
#
FDRCorrectedPValue = function ( pV, q ) {
  #	Test: pV = pValues[, iStat, iPair ];
  
  pStar <- NA;
  pV <- pV[ !is.na(pV) ];
  spvals <- sort ( as.vector(pV), decreasing=FALSE );
  numTests <- length( spvals );
  #plot( spvals , type="p", pch="*" );
  tmp <- spvals <= ( seq ( 1, numTests  ) * q / numTests );
  imax <- which ( tmp );
  if ( length(imax) > 0 ) {
    imax <- max(imax);
    pStar <- spvals[imax];
  } # if ( imax > 0 )
  return ( pStar );
  
} # FDRCorrectedPValue = function ( pValues, q ) {

#
#	GetVecVar
#
GetVecVar = function ( varPlusKey, keyList ) {
  # Testing varPlusKey = pRes[, c("IMAGE_ID", statNames[1]) ]; keyList = full.CN$LONI_IMAGE; iRec = 1;
  numRec = length ( keyList );
  res = rep ( NA, numRec );
  for ( iRec in 1:numRec ) {
    tmp = which ( varPlusKey[,1] == as.character(keyList[iRec]) );
    if ( length(tmp) ) {
      res[iRec] = varPlusKey[tmp,2];
    } else {
      stop( "Not Found iRec=", iRec );
    } # if ( length(tmp)>0 ) {
  } # for ( iRec in 1:numRec ) {
  return ( res );
} # GetVecVar = function ( varPlusKey, key ) {

#
# NiceHistogramPlot
#
NiceHistogramPlot = function( xvec, title_text, xlab_text, ylab_text ) {
  # xvec = table( cdf_matched[ ! cdf_matched[ , "Vehicle.Make" ] == "",  "Vehicle.Make" ] );
  # title_text = "Title"; xlab_text = "Attribute Name"; ylab_text = "Count";
  
  tmp = data.frame( xvec, stringsAsFactors = FALSE ); names( tmp ) = "value";
  p1 = ggplot( tmp, aes( x = value ) ) + 
    geom_histogram( colour = "black", fill = "steelblue", binwidth = 1 ) + theme_minimal() +
    ggtitle( title_text ) + xlab( xlab_text ) + ylab( ylab_text );
  
  plot( p1 );
  
} # NiceHistogramPlot = function( ) {

#
# OneHotConvert
#   The df should contain only factor and numeric classes
#
OneHotConvert = function ( df ) {
  
  return( data.frame( predict( dummyVars( " ~ . ", data = df, fullrank = TRUE ), newdata = df ) ) );
  
} # OneHotConvert = function ( df ) {

#
#	PairwiseCohenD
#	Expect data to contain numerical values and labels to contain 3 label values.
#	Perform the three pair-wise computations in some sequence as t-test.
#	e.g., CN vs MCI, CN v ZDementia, MCI v ZDementia
#
PairwiseCohenDTests = function ( data, labels, resampleData ) {
  # Testing: data = c( stats.CN[,1], stats.MCI[,1], stats.Dementia[,1]); labels = DX; resampleData = FALSE;
  # Testing:  data = c ( f.CN$Entorhinal, f.MCI$Entorhinal, f.Dementia$Entorhinal ); labels = DX.label; resampleData = FALSE;
  dVals = rep ( NA, 3 ) ;
  pooled = TRUE;
  paired = FALSE;
  na.rm = TRUE;
  iloc.CN = which ( labels == "CN" );
  iloc.MCI = which ( labels == "MCI" );
  iloc.ZDementia = which ( labels == "Dementia" );
  
  dVals[1] = PairwiseDVal( data[iloc.MCI], data[iloc.CN], "MCI", "CN", pooled, paired, na.rm, resampleData );
  dVals[2] = PairwiseDVal( data[iloc.ZDementia], data[iloc.CN], "Dementia", "CN", pooled, paired, na.rm, resampleData );
  dVals[3] = PairwiseDVal( data[iloc.ZDementia], data[iloc.MCI], "Dementia", "MCI", pooled, paired, na.rm, resampleData );
  
  return ( abs ( as.numeric ( dVals ) ) );
  
} # PairwiseCohenDTests = function ( z.CompScores, z.Labels ) {

#
#	PairwiseDVal
#
PairwiseDVal = function ( x, y, xLabel, yLabel, pooled, paired, na.rm, nullHypCondition ) {
  #	Test: x = data[iloc.MCI]; y = data[iloc.CN];
  naCut = 5;
  dVal = NA;
  nx = length(x); ny = length(y);
  if ( nullHypCondition ) {
    x = x[sample(seq(1,nx), nx, replace=TRUE)];
    y = y[sample(seq(1,ny), ny, replace=TRUE)];
  } # if ( nullHypCondition ) {
  if ( sum(!is.na(x))>=naCut & sum(!is.na(y))>=naCut ) {
    dVal = unlist( cohen.d( c(x,y), as.factor(c(rep(xLabel,nx), rep(yLabel,ny))), pooled=pooled, paired=paired, na.rm=na.rm ) )[3];
  } # if ( sum(!is.na(x))>=naCut & sum(!is.na(y))>=naCut ) {
  return ( dVal );
} # PairwiseTVal = function ( ) {

#
# QuickPairWiseVisual
#
#   This routine will produce two plots
#   1) Histogram comparing the two distributions, with x% trimming from both ends
#   2) Scatter plot
#
QuickPairWiseVisual = function( x, y, z, trim, xText, yText, zText, scatter_plots_flag = FALSE, main_title = NULL,
                                xaxis_text = "xlab", yaxis_text = "ylab" ) {
  # x = vlist[[1]]; y = vlist[[2]]; z = vlist[[3]]; trim = 0;
  # xText = "CN"; yText = "MCI"; zText = "DEMENTIA"; scatter_plots_flag = FALSE; main_title = paste( ifield, " vs DX", sep = "" ); xaxis_text = ifield; yaxis_text = "Density";
  
  
  xLimits = quantile( x, probs = c( trim, 1-trim ), na.rm = TRUE );
  ix = !is.na( x ) & ( x >= xLimits[1] ) & ( x <= xLimits[2] );
  
  yLimits = quantile( y, probs = c( trim, 1-trim ), na.rm = TRUE );
  iy = !is.na( y ) & ( y >= yLimits[1] ) & ( y <= yLimits[2] );
  
  if ( length(z) ) {
    
    zLimits = quantile( z, probs = c( trim, 1-trim ), na.rm = TRUE );
    iz = !is.na( z ) & ( z >= zLimits[1] ) & ( z <= zLimits[2] );
    tz = z[ iz ];
    
  } # if ( length(z) ) {
  
  tx = x[ ix ];
  ty = y[ iy ];
  
  #
  # Do density and boxplots
  #
  if ( length( z ) ) {
    
    data = data.frame( value=c( tx, ty, tz ), source=factor( c( rep( xText,length(tx) ),
                                                                rep( yText,length(ty) ),
                                                                rep( zText,length(tz) ) ), levels = c( xText, yText, zText ), ordered = TRUE )
    );
    
  } else {
    
    data = data.frame( value=c( tx, ty ), source=factor( c( rep( xText,length(tx) ), rep( yText,length(ty) ) ), levels = c( xText, yText, zText ), ordered = TRUE ) );
    
  } #  if z != NULL ) {
  
  p1 = ggplot( data, aes( x = value, fill = source ) ) + geom_density( alpha=0.25 ) + theme_minimal() +theme( legend.position="bottom" ) +
    ggtitle( main_title ) + xlab( xaxis_text ) + ylab( yaxis_text );
  
  p2 = ggplot( data, aes( x = source, y = value ) ) + geom_boxplot() + coord_flip() + theme_minimal() + theme( legend.position="bottomQ" ) +
    ggtitle( main_title ) + xlab( xaxis_text ) + ylab( yaxis_text );
  
  #
  # Do scatter plot
  #
  if ( scatter_plots_flag ) {
    if ( length( z ) ) {
      
      data = data.frame( x = x[ ix & iy ], y = y[ ix & iy ] );
      p3 = ggplot( data, aes( x = x, y = y ) ) + geom_point( shape=1 ) + geom_smooth() + labs( x = xText ) + labs( y = yText );
      cat( " QuickPairWiseVisual: %_Fill(x) =", 100.0 * round( sum(!is.na(x)) / length(x), 2 ),
           ": %_Fill(y) =", 100.0 * round( sum(!is.na(y)) / length(y), 2 ),
           ": %_Full(x,y)=", 100.0 * round( sum( ix & iy) / length(x), 2 ),
           ": rho =", cor ( x[ ix & iy ],y[ ix & iy ] ),
           ": med =", median( x[ ix & iy ] ), median( y[ ix & iy ] ),
           "\n" );
      
      data = data.frame( x = x[ ix & iz ], y = z[ ix & iz ] );
      p4 = ggplot( data, aes( x = x, y = y ) ) + geom_point( shape=1 ) + geom_smooth() + labs( x = xText ) + labs( y = zText );
      cat( " QuickPairWiseVisual: %_Fill(x) =", 100.0 * round( sum(!is.na(x)) / length(x), 2 ),
           ": %_Fill(z) =", 100.0 * round( sum(!is.na(z)) / length(z), 2 ),
           ": %_Full(x,z)=", 100.0 * round( sum( ix & iz) / length(x), 2 ),
           ": rho =", cor ( x[ ix & iz ],z[ ix & iz ] ),
           ": med =", median( x[ ix & iz ] ), median( z[ ix & iz ] ),
           "\n" );
      
      data = data.frame( x = y[ iy & iz ], y = z[ iy & iz ] );
      p5 = ggplot( data, aes( x = x, y = y ) ) + geom_point( shape=1 ) + geom_smooth() + labs( x = yText ) + labs( y = zText );
      cat( " QuickPairWiseVisual: %_Fill(y) =", 100.0 * round( sum(!is.na(y)) / length(y), 2 ),
           ": %_Fill(z) =", 100.0 * round( sum(!is.na(z)) / length(z), 2 ),
           ": %_Full(y,z)=", 100.0 * round( sum( iy & iz) / length(x), 2 ),
           ": rho =", cor ( y[ iy & iz ],z[ iy & iz ] ),
           ": med =", median( y[ iy & iz ] ), median( z[ iy & iz ] ),
           "\n" );
      
      multiplot( p1, p2, p3, p4, p5, cols = 2 );
      
    } else {
      
      data = data.frame( x = x[ ix & iy ], y = y[ ix & iy ] );
      p3 = ggplot( data, aes( x = x, y = y ) ) + geom_point( shape=1 ) + geom_smooth() + labs( x = xText ) + labs( y = yText )
      
      multiplot( p1, p2, p3, cols = 2 );
      
      cat( " QuickPairWiseVisual: %_Fill(x) =", 100.0 * round( sum(!is.na(x)) / length(x), 2 ),
           ": %_Fill(y) =", 100.0 * round( sum(!is.na(y)) / length(y), 2 ),
           ": %_Full(x,y)=", 100.0 * round( sum( ix & iy) / length(x), 2 ),
           ": rho =", cor ( x[ ix & iy ],y[ ix & iy ] ),
           "\n" );
      
    } #   if ( length( z ) ) {
  } else {
    
    plot( p1 );
    plot( p2 );
    
  } #   if ( scatter_plots_flag ) {
  
} # QuickPairWiseVisual = function( ) {

#
# QuickPairWiseVisual_Categorical
#
QuickPairWiseVisual_Categorical = function ( x, y, xlabel, ylabel, title_text ) {
  # Unit_Test: x = unlist( df_all ); y = unlist( df_sample ); xlabel = "All Retail Accounts"; ylabel = "12-Month Cohort"; title_text = "PTYPE";
  
  code = factor( c( x, y ) );
  label = c( rep(xlabel, length(x) ), rep( ylabel, length(y) ) );
  
  x_table = table( code[ label == xlabel ] );
  x_table = x_table / sum( x_table );
  
  y_table = table( code[ label == ylabel ] );
  y_table = y_table / sum( y_table );
  
  df = data.frame( x_table, y_table ); df = df[ , c(1,2,4 ) ]; names( df ) = c( "Types", "ALL", "12-Month" );
  df.m = melt( df, id.vars = "Types" );
  
  p = ggplot( df.m, aes( Types, value ) ) +   
    geom_bar( aes( fill = variable ), position = "dodge", stat="identity" ) + theme_minimal() + theme( legend.position="bottom" ) +
    ggtitle( title_text ) + xlab( title_text ) + ylab( "Density" );
  
  plot( p );
  
} # QuickPairWiseVisual = function ( ) {

#
# QuickANOVAScan
#
QuickANOVAScan <- function( df ) {
  
  #
  # Return will be a list of lists!
  #
  res = list()
  
  #
  # Do ANOVA on the population
  #
  field_list_anova <- unique( tolower( c( "AGE_ON_EXAMDATE", "PTEDUCAT", "CDRSB", "MMSE", "MOCA", "weight_mprage",
                                         names( df )[ grep( "adas", names( df ) ) ],
                                         names( df )[ grep( "rav", names( df ) ) ],
                                         names( df )[ grep( "ecog", names( df ) ) ],
                                         names( df )[ grep( "_aseg_norm", names( df ) ) ],
                                         names( df )[ grep( "_subhipp_norm", names( df ) ) ],
                                         names( df )[ grep( "_amyg_norm", names( df ) ) ],
                                         names( df )[ grep( "_thickness", names( df ) ) ] ) ) )
  field_list_anova <- c( field_list_anova, c( "avg_censored_motion", "max_censored_displ", "df_left", "censor_fraction",
                                              "avg_tsnr", "global_corr", "dice_coeff", "dfused" ) )
  field_list_anova = field_list_anova[ ! field_list_anova %in% c ( "rid_mprage_rsfmri_key" ) ];
  res_anova = list();
  for ( ifield in field_list_anova ) {
    
    vlist = list( CN = df %>% dplyr::filter( dx_bl == "CN" ) %>% dplyr::select( one_of( ifield ) ),
                  SMC = df %>% dplyr::filter( dx_bl == "SMC" ) %>% dplyr::select( one_of( ifield ) ),
                  EMCI = df %>% dplyr::filter( dx_bl == "EMCI" ) %>% dplyr::select( one_of( ifield ) ),
                  LMCI = df %>% dplyr::filter( dx_bl == "LMCI" ) %>% dplyr::select( one_of( ifield ) ) );
    
    df_tmp <- DoANOVA( vlist );
    res_anova[[ ifield ]] <- unlist( summary(df_tmp$res_anova ) )
    
    cat( " ** ANOVA: ", ifield, ":", res_anova[[ ifield ]][ "Pr(>F)1" ], "\n" );
    
  } # for ( ifield in field_list_anova ) {
  
  #
  # Chi-Sq on Clinical Populations
  #
  
  #
  # Gender
  #
  M = as.table( rbind( 	c( nrow( df %>% dplyr::filter( dx_bl == "CN", ptgender == "Female" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "SMC", ptgender == "Female" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "EMCI", ptgender == "Female" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "LMCI", ptgender == "Female" ) ) ),
                        
                        c( nrow( df %>% dplyr::filter( dx_bl == "CN", ptgender == "Male" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "SMC", ptgender == "Male" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "EMCI", ptgender == "Male" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "LMCI", ptgender == "Male" ) ) )
  ) );
  dimnames( M ) = list( PTGENDER = c("F", "M"), DX = c( "CN","SMC", "EMCI", "LMCI") );
  ggplot( data = data.frame( M ), aes( x = DX, y = Freq, fill = PTGENDER ) ) +
    geom_bar( stat = "identity", position = position_dodge() ) +
    ggtitle( "GENDER vs DX" ) + xlab( "DX" ) + ylab( "Count" );
  gender_chisq <- chisq.test( M );
  
  #
  # APOE4
  #
  M = as.table( rbind( 	c( nrow( df %>% dplyr::filter( dx_bl == "CN", apoe4 == "0" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "SMC", apoe4 == "0" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "EMCI", apoe4 == "0" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "LMCI", apoe4 == "0" ) ) ),
                        
                        c( nrow( df %>% dplyr::filter( dx_bl == "CN", apoe4 == "1" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "SMC", apoe4 == "1" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "EMCI", apoe4 == "1" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "LMCI", apoe4 == "1" ) ) ),
                        
                        c( nrow( df %>% dplyr::filter( dx_bl == "CN", apoe4 == "2" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "SMC", apoe4 == "2" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "EMCI", apoe4 == "2" ) ),
                           nrow( df %>% dplyr::filter( dx_bl == "LMCI", apoe4 == "2" ) ) )
  ) );
  dimnames( M ) = list( APOE4 = c( "0", "1", "2" ), DX = c( "CN","SMC", "EMCI", "LMCI") );
  ggplot( data = data.frame( M ), aes( x = DX, y = Freq, fill = APOE4 ) ) +
    geom_bar( stat = "identity", position = position_dodge() ) +
    ggtitle( "APOE4 vs DX" ) + xlab( "DX" ) + ylab( "Count" );
  apoe4_chisq <- chisq.test( M );
  
  #
  # Do some FDR correction analysis for the different groupings of ANOVA tests:
  #   aseg; sub_hip; aparc; 
  #
  q = 0.05;
  
  #
  # Thickness
  #
  ifield_thick <- field_list_anova[ grep( "_thickness", field_list_anova ) ]
  pV <- NULL;
  for ( ifield in ifield_thick ) {
    
    pV <- c( pV, res_anova[[ ifield ]][ "Pr(>F)1" ] );
    
  } # for ( ifield in ifield_thick ) {
  df_thickness <- data.frame( field = ifield_thick, p = pV );
  df_thickness$p_star_thickness = FDRCorrectedPValue( pV, q );
  df_thickness <- df_thickness %>% arrange( p );
  cat( "** Thickness\n" )
  
  #
  # Aseg Norm
  #
  ifield_aseg_norm <- field_list_anova[ grep( "_aseg_norm", field_list_anova ) ]
  ifield_aseg_norm <- ifield_aseg_norm[ ! grepl("cerebellum", ifield_aseg_norm) ]
  ifield_aseg_norm <- ifield_aseg_norm[ ! grepl("_vent_", ifield_aseg_norm) ]
  ifield_aseg_norm <- ifield_aseg_norm[ ! grepl("_ventricle_", ifield_aseg_norm) ]
  pV <- NULL;
  for ( ifield in ifield_aseg_norm ) {
    
    pV <- c( pV, res_anova[[ ifield ]][ "Pr(>F)1" ] );
    
  } # for ( ifield in ifield_aseg_norm ) {
  df_aseg_norm <- data.frame( field = ifield_aseg_norm, p = pV );
  df_aseg_norm$p_star_aseg_norm <- FDRCorrectedPValue( pV, q );
  df_aseg_norm <- df_aseg_norm %>% arrange( p );
  cat( "** Aseg_norm\n" )
  
  #
  # Sub Hipp
  #
  ifield_sub_hipp_norm = field_list_anova[ grep( "_subhipp_norm", field_list_anova ) ]
  pV = NULL;
  for ( ifield in ifield_sub_hipp_norm ) {
    
    pV = c( pV, res_anova[[ ifield ]][ "Pr(>F)1" ] );
    
  } # for ( ifield in ifield_sub_hipp_norm ) {
  df_sub_hipp_norm = data.frame( field = ifield_sub_hipp_norm, p = pV );
  df_sub_hipp_norm$p_star_sub_hipp_norm = FDRCorrectedPValue( pV, q );
  df_sub_hipp_norm <- df_sub_hipp_norm %>% arrange( p );
  cat( "** SubHipp\n" )
  
  #
  # Amyg
  #
  ifield_amyg_norm = field_list_anova[ grep( "_amyg_norm", field_list_anova ) ]
  pV = NULL;
  for ( ifield in ifield_amyg_norm ) {
    
    pV = c( pV, res_anova[[ ifield ]][ "Pr(>F)1" ] );
    
  } # for ( ifield in ifield_amyg_norm ) {
  df_amyg_norm = data.frame( field = ifield_amyg_norm, p = pV );
  df_amyg_norm$p_star_amyg_norm = FDRCorrectedPValue( pV, q );
  df_amyg_norm <- df_amyg_norm %>% arrange( p );
  cat( "** Amyg\n" )
  
  return ( list( res_anova = res_anova, gender_chisq = gender_chisq, apoe4_chisq = apoe4_chisq,
                 df_thickness = df_thickness, df_aseg_norm = df_aseg_norm, 
                 df_sub_hipp_norm = df_sub_hipp_norm, df_amyg_norm = df_amyg_norm ) )
  
} # QuickANOVAScan = function( df )

#
# ReadCSVFile
#
ReadCSVFile = function( rawSourceDirName, fName ) {
  # Test: iZipcode = zipCodeList[ iZip ]; df.out = df.fused$df.fused;
  # fName = paste( rawSourceDirName, "fused.", iZipcode, ".csv", sep="" );
  
  cat( " ReadCSVFile: Reading File=", fName, "\n" );
  #write.table( df.out, file = fName, col.names = TRUE, eol = '\n', na = "NA", quote=TRUE, sep = "," );
  #tmp = read.table( file = fName, header = TRUE, quote="\"'", row.names = NULL, sep = ",", stringsAsFactors=FALSE );
  df.in = fread( fName, header=TRUE, stringsAsFactors=FALSE );
  
  #
  # Debug support - on occasion may read in fewer lines than written out.  That means there are some
  #   funky characters lurking in the data.
  #
  if ( 0 ) {
    
    tmp = count.fields( file = fName, quote="\"", sep = "," );
    length( tmp ); table( tmp );
    tmp = readLines( fName );
    tmp.strsplit = strsplit( tmp, split="\",\"" );
    tmp.strsplit.length = sapply( tmp.strsplit, length );
    bad.apple.866 = which( tmp.strsplit.length == 866 );
    bad.apple.873 = which( tmp.strsplit.length == 873 )
    
    QuickCompPeek( df.fused$df.fused, df.fused$df.fused, "CATRINABAILEY", "FIRST_LAST_NAME_KEY", varList, varList );
    
  } #   if ( 0 ) {
  
  # 
  # read.csv automatically adds a "row.number" column; get rid of it.
  #
  
  return( df.in );
  
} # ReadCSVFile = function( df.fused, rawSourceDirName, iZipcode ) {

#
# WriteCSVFile
#
WriteCSVFile = function( df.out, rawSourceDirName, fName ) {
  # Test: iZipcode = zipCodeList[ iZip ]; df.out = df.fused$df.fused;
  # fName = paste( rawSourceDirName, "fused.", iZipcode, ".csv", sep="" );
  # df.out = res_study[[ "base_indv_pc_bev" ]]; fName = report_file_name; df.out = dtmp; fName = file_name; rawSourceDirName = NULL;
  
  cat( " WriteCSVFile: Writing File=", fName, "\n" );
  write.table( df.out, file = fName, col.names = TRUE, na = "NA", quote=TRUE, row.names=FALSE, sep = "," );
  
  # tmp = read.csv( file = fName, header = TRUE, quote="\"", row.names = NULL, sep = ',', stringsAsFactors=FALSE );
  # df.in = tmp[ , 2:dim( tmp )[2] ];
  
} # WriteCSVFile = function( df.fused, rawSourceDirName, iZipcode ) {

#
# WriteExcellAllSheets
#
WriteExcelReport = function( x, fName ) {
  # Test: x = list( Training = gain_lift_train, Validation = gain_lift_validation ); fName = "test.xlsx";
  
  WriteXLS( x, ExcelFileName = fName, SheetNames = NULL, col.names = TRUE,
            AdjWidth = TRUE, AutoFilter = FALSE, BoldHeaderRow = TRUE );
  
} # WriteExcellAllSheets = function( fName ) {
#