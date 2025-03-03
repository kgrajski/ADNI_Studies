#
#	ADNI3.GetToKnow.R
#

#
# 2Oct2018 - Load and explore ADNIMERGE.
#

#
# Set the stage.
#
rm ( list = ls() );
assign( "last.warning", NULL, envir=baseenv() );

#
# Remove (previous) adnimerge.
#
rm ( list = ls() );
assign( "last.warning", NULL, envir=baseenv() );
remove.packages("adnimerge");

#
# Reboot (so to speak)
#
library( plyr );
library( dplyr );

library( caret );
library( corrplot );
library( data.table );
library( doParallel );
library( effsize );
library( Hmisc );
library( psych );
library( stringr );
library( WriteXLS );

#
# Read Helper Functions
#
source( "ADNI3.HelperFunctions.R" );

#
# Install ADNIMERGE.  Get this from the ADNI Download Study Data:
#  https://ida.loni.usc.edu/pages/access/studyData.jsp?categoryId=16&subCategoryId=43
#
if ( 0 ) {
  install.packages( "//Users/kag/Documents/kag/ADNI3/data/ADNIMERGE/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source" );
  help( package = "ADNIMERGE" );
  library( ADNIMERGE );
} # if ( 0 )

#
# Get adni3 from adnimerge.
#
adni3 = adnimerge %>% dplyr::filter( COLPROT == "ADNI3" );
summary ( adni3 );

#
# Get ADNI3 image database search results.
#
loni_list = read.csv( file =  "//Users/kag/Documents/kag/ADNI3/data/IDA_SEARCH/idaSearch_11_21_2021.csv", header = TRUE, stringsAsFactors = FALSE );
loni_list$Description = str_squish( toupper( loni_list$Description ) );
descriptions = unique( loni_list$Description );

#
# MPRAGE
#   May need to do further subsetting.  Let's see where we end up first.
#
contains_mprage = descriptions[ grep( "MPRAGE", descriptions ) ];
df_mprage = loni_list %>% dplyr::filter( Description %in% contains_mprage );
nrow( df_mprage );
sort( table( df_mprage$Description ), decreasing = TRUE );

#
# fMRI
#
df_fmri = loni_list %>% dplyr::filter( Modality == "fMRI" );
sort( table( df_fmri$Description ), decreasing = TRUE );
contains_rsfmri = c( "AXIAL RSFMRI (EYES OPEN)", "AXIAL FCMRI (EYES OPEN)", "AXIAL MB RSFMRI (EYES OPEN)",
                     "AXIAL RESTING FCMRI (EYES OPEN)", "AXIAL FCMRI (EYES OPEN)_REPEAT" );
df_rsfmri = df_fmri %>% dplyr::filter( Description %in% contains_rsfmri );
nrow( df_rsfmri );
sort( table( df_rsfmri$Imaging.Protocol ), decreasing = TRUE );
sort( table( df_rsfmri$Description ), decreasing = TRUE );

#
# Winnow the fMRI list by QC results.
#

#
# Read in the QC data, but we'll probably just use what comes with ADNIMERGE.
#
qc = read.csv( file =  "//Users/kag/Documents/kag/ADNI3/data/MR_IMAGE_QUALITY/MAYOADIRL_MRI_QUALITY_ADNI3.csv", header = TRUE, stringsAsFactors = FALSE );
dim( qc );

#
# Condition the mayoadirl_mri_imageqc data frame
#
data( mayoadirl_mri_imageqc );
df_qc = mayoadirl_mri_imageqc;
df_qc$LONI_IMAGE = as.integer( df_qc$LONI_IMAGE );
df_qc = df_qc %>% dplyr::filter( LONI_IMAGE %in% c( df_mprage$Image.ID, df_rsfmri$Image.ID ) );

df_qc$STUDY_OVERALLPASS = as.integer( df_qc$STUDY_OVERALLPASS );
df_qc$STUDY_COMMENTS = ifelse( ( is.na( df_qc$STUDY_COMMENTS ) ), ( "" ), ( df_qc$STUDY_COMMENTS ) );
df_qc$STUDY_PROTOCOL_STATUS = as.integer( df_qc$STUDY_PROTOCOL_STATUS );
df_qc$STUDY_PROTOCOL_COMMENT = ifelse( ( is.na( df_qc$STUDY_PROTOCOL_COMMENT ) ), ( "" ), ( df_qc$STUDY_PROTOCOL_COMMENT ) );
df_qc$SERIES_QUALITY = as.integer( df_qc$SERIES_QUALITY );
df_qc$SERIES_COMMENTS = ifelse( ( is.na( df_qc$SERIES_COMMENTS ) ), ( "" ), ( df_qc$SERIES_COMMENTS ) );
df_qc$PROTOCOL_STATUS = as.integer( df_qc$PROTOCOL_STATUS );
df_qc$PROTOCOL_COMMENTS = ifelse( ( is.na( df_qc$PROTOCOL_COMMENTS ) ), ( "" ), ( df_qc$PROTOCOL_COMMENTS ) );
df_qc$STUDY_MEDICAL_ABNORMALITIES = ifelse( ( is.na( df_qc$STUDY_MEDICAL_ABNORMALITIES ) ), ( "" ), ( df_qc$STUDY_MEDICAL_ABNORMALITIES ) );
df_qc$MEDICAL_EXCLUSION = as.integer( df_qc$MEDICAL_EXCLUSION );
df_qc$QUALIFICATION = as.integer( df_qc$QUALIFICATION );

#
# Pre-join checks for intersection of image IDs
#
length( intersect( df_rsfmri$Image.ID, df_qc$LONI_IMAGE ) );
length( intersect( df_mprage$Image.ID, df_qc$LONI_IMAGE ) );

#
# Join on df_rsfmri$Image.ID and df_qc$LONI_IMAGE
#
nrow( df_rsfmri );
df_rsfmri = df_rsfmri %>% left_join( df_qc, by = c( "Image.ID" = "LONI_IMAGE" ) );
nrow( df_rsfmri );

df_rsfmri = df_rsfmri %>% dplyr::filter( STUDY_COMMENTS == "", STUDY_PROTOCOL_STATUS == 1, STUDY_PROTOCOL_COMMENT == "QC'd",
                                         SERIES_QUALITY < 4, SERIES_COMMENTS == "", PROTOCOL_STATUS == 1,
                                         PROTOCOL_COMMENTS == "", SERIES_SELECTED == TRUE, STUDY_MEDICAL_ABNORMALITIES == 0,
                                         MEDICAL_EXCLUSION == 0, QUALIFICATION == 0 ) %>% group_by( Subject.ID, Study.Date );

#
# Join on df_mprage$Image.ID and df_qc$LONI_IMAGE
#
nrow( df_mprage );
df_mprage = df_mprage %>% left_join( df_qc, by = c( "Image.ID" = "LONI_IMAGE" ) );
nrow( df_mprage );

df_mprage = df_mprage %>% dplyr::filter( STUDY_COMMENTS == "", STUDY_PROTOCOL_STATUS == 1, STUDY_PROTOCOL_COMMENT == "QC'd",
                                         SERIES_QUALITY < 4, SERIES_COMMENTS == "", PROTOCOL_STATUS == 1,
                                         PROTOCOL_COMMENTS == "", SERIES_SELECTED == TRUE, STUDY_MEDICAL_ABNORMALITIES == 0,
                                         MEDICAL_EXCLUSION == 0, QUALIFICATION == 0 ) %>% group_by( Subject.ID, Study.Date );

#
# Join mprage and rsfmri into final candidate list
#
df_study = inner_join( df_mprage, df_rsfmri, by = c( "RID", "Visit", "Study.Date" ), suffix = c( ".mprage", ".rsfmri" ) );

#
# Attempt a join with adnimerge (brings in all the adnimerge goodies).
#
adnimerge_tmp = adnimerge %>% dplyr::filter( RID %in% df_study$RID, COLPROT == "ADNI3" );
adnimerge_tmp$EXAMDATE = as.character( adnimerge_tmp$EXAMDATE );
adnimerge_tmp$EXAMDATE.bl = as.character( adnimerge_tmp$EXAMDATE.bl );
adnimerge_tmp$RID = as.integer( adnimerge_tmp$RID );

#
# Quite a few RID have > entries in adnimerge.  Need to look at these.
#
head( sort( table( adnimerge_tmp$RID ), decreasing = TRUE ), 200 );
table( adnimerge_tmp$DX );
table( adnimerge_tmp$MMSE ); sum( is.na( adnimerge_tmp$MMSE ) );
table( adnimerge_tmp$CDRSB ); sum( is.na( adnimerge_tmp$CDRSB ) );
head( sort( table( adnimerge_tmp$RID ), decreasing = TRUE ), 100 );

#
# Try with EXAMDATE
#
tmp1 = df_study %>% left_join( adnimerge_tmp, by = c( "RID" = "RID", "Study.Date" = "EXAMDATE" ) );
head( sort( table( tmp1$RID ), decreasing = TRUE ), 20 );
table( ( tmp1$DX ) );
tmp2 = df_study %>% left_join( adnimerge_tmp, by = c( "RID" = "RID", "Study.Date" = "EXAMDATE.bl" ) );

#
# Evaluate for multiple scans per RID.
#
dim( df_study );
uniqueN( df_study$RID );
sort( table( df_study$RID ), decreasing = TRUE );
df_study %>% ungroup() %>% select( Visit ) %>% table()

#
# Generate the initial candidate clinical study population.
# The goal is one record per RID.  Some downside here, because the images may fail
# in later stages of the pipeline.
#
tmp = list();
tmp[[1]] = df_study %>% dplyr::filter( Visit == "ADNI3 Year 2 Visit" );
tmp[[2]] = df_study %>% dplyr::filter( ! ( RID %in% unique( tmp[[1]]$RID ) ) & Visit == "ADNI3 Year 1 Visit" );
tmp[[3]] = df_study %>% dplyr::filter( ! RID %in% unique( c( tmp[[1]]$RID, tmp[[2]]$RID ) ) );
df_study = rbindlist( tmp, use.names = TRUE );
nrow( df_study );
unique( df_study$RID );
uniqueN( df_study$RID );
sort( table( df_study$RID ), decreasing = TRUE );

#
# Need to join with adnimerge
#     Might be able to do it on image id? - NOPE - Weird.
#     By Exam date?
#
if ( 1 ) {
  
  #intersect( as.integer( df_study$Image.ID.mprage ), as.integer( adnimerge$IMAGEUID[ ! is.na( adnimerge$IMAGEUID ) ] ) );
  #intersect( as.integer( df_study$Image.ID.rfmri ), as.integer( adnimerge$IMAGEUID[ ! is.na( adnimerge$IMAGEUID ) ] ) );
  
  names( adnimerge )[ grep( "DATE", toupper( names( adnimerge ) ) ) ];
  names( df_study )[ grep( "DATE", toupper( names( df_study ) ) ) ];
  head( df_study$Study.Date ); head( df_study$SERIES_DATE.mprage ); head( df_study$SERIES_DATE.rsfmri );
  class( df_study$Study.Date ); class( df_study$SERIES_DATE.mprage ); class( df_study$SERIES_DATE.rsfmri );
  head( adnimerge_tmp$EXAMDATE );
  class( adnimerge_tmp$EXAMDATE );
  
  adnimerge_tmp$EXAMDATE_CONV = ConvertDateType0( adnimerge_tmp$EXAMDATE );
  class( df_study$Study.Date ); class( adnimerge_tmp$EXAMDATE_CONV );
  uniqueN( df_study$Study.Date );
  uniqueN( adnimerge_tmp$EXAMDATE_CONV );
  intersect( df_study$Study.Date, adnimerge_tmp$EXAMDATE_CONV );
  
  df_study$SERIES_DATE_CONV = ConvertDateType1( df_study$SERIES_DATE.mprage );
  head(  df_study$SERIES_DATE_CONV ); head( adnimerge_tmp$EXAMDATE_CONV );
  class( df_study$SERIES_DATE_CONV ); class( adnimerge_tmp$EXAMDATE_CONV );
  uniqueN( df_study$SERIES_DATE_CONV ); uniqueN( adnimerge_tmp$EXAMDATE_CONV );
  intersect( df_study$SERIES_DATE_CONV, adnimerge_tmp$EXAMDATE_CONV );
  
  head( df_study$RID );
  adnimerge_tmp %>% dplyr::filter( RID == 1261 );
  df_study %>% dplyr::filter( RID == 1261 );
  
  #
  # convert df_study study date and adnimerge exam date to a common date format and find closest.
  #
  save_df_study = df_study;
  df_study = save_df_study;
  
  df_study$STUDY_DATE_AS_DATE = as.Date(df_study$Study.Date, format = "%m/%d/%Y" );
  adnimerge_tmp$EXAMDATE_AS_DATE = as.Date( adnimerge_tmp$EXAMDATE );
  adnimerge_tmp$recloc = seq( 1, nrow( adnimerge_tmp ) );
  recloc = data.frame( recloc = rep( NA, nrow( df_study ) ), dist = rep( NA, nrow( df_study ) ) );
  df_study$adnimerge_rec = rep( NA, nrow( df_study ) );
  df_study$adnimerge_rec_chk = rep( NA, nrow( df_study ) );
  df_study$adnimerge_rec_dist = rep( NA, nrow( df_study ) );
  
  #
  # Final check to set expectation that we will find df_study RID in adnimerge_tmp
  #
  rid_list_common = intersect( adnimerge_tmp$RID, df_study$RID );
  df_study = df_study %>% dplyr::filter( RID %in% rid_list_common );
  
  # as.Date("2019-06-01") - as.Date("2019-07-01")
  
  for ( irec in 1:nrow( df_study ) ) {
    
    df_tmp = adnimerge_tmp %>% dplyr::filter( RID == df_study$RID[ irec ] ) %>% dplyr::select( EXAMDATE_AS_DATE, recloc );
    tmp_dist = as.Date( df_tmp$EXAMDATE_AS_DATE ) - df_study$STUDY_DATE_AS_DATE[ irec ];
    iloc = which.min( abs( tmp_dist ) );
    df_study$adnimerge_rec[ irec ] = df_tmp$recloc[ which.min( abs( tmp_dist ) ) ];
    df_study$adnimerge_rec_chk[ irec ] = sum( tmp_dist < 0 ) == length( tmp_dist );
    df_study$adnimerge_rec_dist[ irec ] = tmp_dist[ iloc ];
    
  } # for ( iRID in uniqueN( df_study$RID ) ) {
  
  #
  # Check results before committing.
  #
  table( df_study$adnimerge_rec_dist );
  df_study_vis_adj = df_study %>% dplyr::filter( abs( adnimerge_rec_dist ) <= 365 );
  adnimerge_tmp = adnimerge_tmp[ df_study_vis_adj$adnimerge_rec, ] %>% dplyr::select( - RID );
  df = cbind( df_study_vis_adj, adnimerge_tmp, stringsAsFactors = FALSE );
  names( df ) = toupper( names( df ) );
  
  #
  # Need minimum field fill rate to proceed.
  #   Have to have AGE, GENDER, PTEDUCAT, MMSE, CDRSB, DX
  #
  nrow( df );
  field_list_req = c( "AGE", "PTGENDER", "PTEDUCAT", "MMSE", "CDRSB", "DX", "YEARS.BL" );
  df = df[ df %>% select( all_of ( field_list_req ) ) %>% complete.cases(), ];
  df$AGE_ON_EXAMDATE = df$AGE + df$YEARS.BL;
  
  #
  # Subset on AGE_ON_EXAMDATE?
  #
  if ( 0 ) {
    
    age_min = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( AGE_ON_EXAMDATE ) %>% unlist() %>% min();
    age_max = df %>% dplyr::filter( DX == "Dementia" ) %>% dplyr::select( AGE_ON_EXAMDATE ) %>% unlist() %>% max();
    #df = df %>% dplyr::filter( AGE_ON_EXAMDATE >= age_mean - age_sd, AGE_ON_EXAMDATE <= age_mean + age_sd );
    
  } # if ( 0 )
  
  #
  # Write the data frame.
  #
  WriteExcelReport( df, "//Users/kag/Documents/kag/ADNI3/data//df_study_adnimerge_21Nov2021.xlsx" );
  
  #
  # Write the text file used to control follow-on batch runs.
  #
  df_tmp = df %>% select( SUBJECT.ID.RSFMRI, IMAGE.ID.MPRAGE, IMAGE.ID.RSFMRI, DX );
  df_tmp$IMAGE.ID.MPRAGE = paste( "I", df_tmp$IMAGE.ID.MPRAGE, sep="" );
  df_tmp$IMAGE.ID.RSFMRI = paste( "I", df_tmp$IMAGE.ID.RSFMRI, sep="" );
  write.table( df_tmp, file = "/Users/kag/Documents/kag/ADNI3/data///df_study_adnimerge_21Nov2021.csv", row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE );
  
} # if ( 1 ) {

#
# Generate a string of image IDs that can be used to download from LONI.
#
if ( 0 ) {
  
  tmp1 = paste( "I", df$IMAGE.ID.MPRAGE, sep = "", collapse = "," );
  tmp2 = paste( "I", df$IMAGE.ID.RSFMRI, sep = "", collapse = "," );
  tmp3 = paste( tmp1, tmp2, sep = "", collapse = "," );
  
} # if ( 0 )

#
# Generate Study Participant Demographic and Metadata Analysis
#

if ( 1 ) {
  
  if ( 0 ) {
    
    field_list_quant = c( "FDG", "PIB", "AV45", "ABETA", "TAU", "PTAU" );
    for ( ifield in field_list_quant ) {
      
      cat( " ** ", ifield, "\n\n" );
      print( summary( df[, ifield ] ) );
      
    } # for ( ifield in field_list_cog ) {
    cor( df %>% select( one_of( field_list_quant ) ), use = "complete.obs" );
    
  } # if ( 0 )
  
  #
  # ANOVAs on Clinical Populations
  #
  field_list_anova = c( "AGE_ON_EXAMDATE", "PTEDUCAT", "CDRSB", "MMSE", "MOCA",
                        names( df )[ grep( "ADAS", names( df ) ) ], names(df)[ grep( "RAV", names( df ) ) ], names(df)[ grep( "ECOG", names( df ) ) ] );
  
  field_list_anova = c( "AGE_ON_EXAMDATE", "PTEDUCAT", "CDRSB", "MMSE", "MOCA" );
  
  anova_res = list();
  for ( ifield in field_list_anova ) {
    
    vlist = list( CN = df %>% dplyr::filter( DX == "CN" ) %>% dplyr::select( one_of( ifield ) ),
                  MCI = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( one_of( ifield ) ),
                  Dementia = df %>% dplyr::filter( DX == "Dementia" ) %>% dplyr::select( one_of( ifield ) ) );
    #QuickPairWiseVisual( vlist[[1]], vlist[[2]], vlist[[3]], 0, "CN", "MCI", "DEMENTIA", FALSE, paste( ifield, " vs DX", sep = "" ), ifield, "Density" );
    anova_res[[ ifield ]] = DoANOVA( vlist );
    cat( " ** ANOVA: ", ifield, ":", unlist( summary(anova_res[[ ifield ]]$res_anova ) )[ "Pr(>F)1" ], "\n" );
    
  } # for ( ifield in field_list_anova ) {
  
  #
  # Corrplot
  #
  if ( 1 ) {
    
    corrplot( cor( cbind( OneHotConvert( data.frame( df$DX ) ), df[ , field_list_anova ] ), use = "complete" ), method = "circle", tl.cex = 0.45, tl.col = 1 );
    
  } # if ( 0 )
  
  #
  # Chi-Sq on Clinical Populations
  #
  #
  # Gender
  #
  M = as.table( rbind( 	c( nrow( df %>% dplyr::filter( DX == "CN", PTGENDER == "Female" ) ),
                           nrow( df %>% dplyr::filter( DX == "MCI", PTGENDER == "Female" ) ),
                           nrow( df %>% dplyr::filter( DX == "Dementia", PTGENDER == "Female" ) ) ),
                        c( nrow( df %>% dplyr::filter( DX == "CN", PTGENDER == "Male" ) ),
                           nrow( df %>% dplyr::filter( DX == "MCI", PTGENDER == "Male" ) ),
                           nrow( df %>% dplyr::filter( DX == "Dementia", PTGENDER == "Male" ) ) )
  ) );
  dimnames( M ) = list( PTGENDER = c("F", "M"), DX = c( "CN","MCI", "Dementia") );
  ggplot( data = data.frame( M ), aes( x = DX, y = Freq, fill = PTGENDER ) ) +
    geom_bar( stat = "identity", position = position_dodge() ) +
    ggtitle( "GENDER vs DX" ) + xlab( "DX" ) + ylab( "Count" );
  chisq.test( M );
  
  #
  # APOE4
  #
  df_tmp = df[ ! is.na( df$APOE4 ), ];
  M = as.table( rbind( 	c( nrow( df_tmp %>% dplyr::filter( DX == "CN", APOE4 == "0" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "MCI", APOE4 == "0" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "Dementia", APOE4 == "0" ) ) ),
                        c( nrow( df_tmp %>% dplyr::filter( DX == "CN", APOE4 == "1" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "MCI", APOE4 == "1" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "Dementia", APOE4 == "1" ) ) ),
                        c( nrow( df_tmp %>% dplyr::filter( DX == "CN", APOE4 == "2" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "MCI", APOE4 == "2" ) ),
                           nrow( df_tmp %>% dplyr::filter( DX == "Dementia", APOE4 == "2" ) ) )
  ) );
  dimnames( M ) = list( APOE4 = c( "0", "1", "2" ), DX = c( "CN","MCI", "Dementia") );
  ggplot( data = data.frame( M ), aes( x = DX, y = Freq, fill = APOE4 ) ) +
    geom_bar( stat = "identity", position = position_dodge() ) +
    ggtitle( "APOE4 vs DX" ) + xlab( "DX" ) + ylab( "Count" );;
  chisq.test( M );
  
} # if ( 1 )
