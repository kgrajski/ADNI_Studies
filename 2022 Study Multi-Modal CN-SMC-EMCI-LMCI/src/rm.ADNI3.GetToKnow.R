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

library( plyr );
library( dplyr );

library( corrplot );
library( data.table );
library( Hmisc );
library( psych );
library( stringr );
library( WriteXLS );

#
# Install ADNIMERGE.
#
#remove.packages("ADNIMERGE");
#install.packages( "/home/rstudio/ADNI3/dev-adni3/05July2019_Downloads/05July2019_ADNIMERGE_0.0.1.tar.gz ", repos = NULL, type = "source" );
#help( package = "ADNIMERGE" );
library( ADNIMERGE );

#
# Read Helper Functions
#
source( "ADNI3.HelperFunctions.R" );

#
# Get adni3 from adnimerge.
#
adni3 = adnimerge %>% dplyr::filter( COLPROT == "ADNI3" );
summary ( adni3 );

#
# Get ADNI3 image database search results.
#
loni_list = read.csv( file =  "/home/rstudio/ADNI3/dev-adni3/05July2019_Downloads/05July2019_idaSearch.csv", header = TRUE, stringsAsFactors = FALSE );
loni_list$Description = str_squish( toupper( loni_list$Description ) );
descriptions = unique( loni_list$Description );
    
    #
    # MPRAGE
    #
contains_mprage = descriptions[ grep( "MPRAGE", descriptions ) ];
contains_mprage = contains_mprage[ grep( "SAGITTAL", contains_mprage ) ];
df_mprage = loni_list %>% dplyr::filter( Description %in% contains_mprage );
nrow( df_mprage );
head( sort( table( df_mprage$Imaging.Protocol ), decreasing = TRUE ) );

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
  # Read in the QC data
  #
qc = read.csv( file =  "/home/rstudio/ADNI3/dev-adni3/05July2019_Downloads/05July2019_MAYOADIRL_MRI_QUALITY_ADNI3.csv", header = TRUE, stringsAsFactors = FALSE );
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
nrow( df_rsfmri );

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

nrow( df_mprage );


#
# Join mprage and rsfmri into final candidate list
#
df_study = inner_join( df_mprage, df_rsfmri, by = c( "RID", "Visit", "Study.Date" ), suffix = c( ".mprage", ".rsfmri" ) );

#
# Attempt a join with adnimerge
#
  #
  # Prep a subset of adnimerge for further work
  #
adnimerge_tmp = adnimerge %>% dplyr::filter( RID %in% df_study$RID, COLPROT == "ADNI3" );
adnimerge_tmp$EXAMDATE = as.character( adnimerge_tmp$EXAMDATE );
adnimerge_tmp$EXAMDATE.bl = as.character( adnimerge_tmp$EXAMDATE.bl );
adnimerge_tmp$RID = as.integer( adnimerge_tmp$RID );

    #
    # Some RID have 2 entries in adnimerge.  Need to look at these.
    #
head( sort( table( adnimerge_tmp$RID ), decreasing = TRUE ), 50 );
table( adnimerge_tmp$DX );
table( adnimerge_tmp$MMSE ); sum( is.na( adnimerge_tmp$MMSE ) );
table( adnimerge_tmp$CDRSB ); sum( is.na( adnimerge_tmp$CDRSB ) );

    #
    # Explore: adnimerge to image data merge using EXAMDATE
    #
tmp1 = df_study %>% left_join( adnimerge_tmp, by = c( "RID" = "RID", "Study.Date" = "EXAMDATE" ) );
head( sort( table( tmp1$RID ), decreasing = TRUE ), 20 );
table( ( tmp1$DX ) );
tmp2 = df_study %>% left_join( adnimerge_tmp, by = c( "RID" = "RID", "Study.Date" = "EXAMDATE.bl" ) );

    #
    # Explore: look for multiple scans per RID.  May want to keep latest, but as ADNI3 proceeds, maybe not.
    #
dim( df_study );
uniqueN( df_study$RID );
head( sort( table( df_study$RID ), decreasing = TRUE ), 25 );

    #
    # From previous iterations of the download/investigate cycle, at this stage we know that there are
    # problems with certain of the rsfmri images.  Remove these from df_study.
    # On 5-July-2019, the problematic rsfmri consisted of those sessions that produced a single dcm file 9456 dcm files - odd, but true!
    #
bad_rsfmri_images = read.table( file = "/home/rstudio/ADNI3/loading_dock/bad_rsfmri_images.txt", header = FALSE, stringsAsFactors = FALSE );
names( bad_rsfmri_images ) = c( "RID", "IMAGE.ID.RSFMRI" );
#df_study_save = df_study;
df_study = df_study_save;
df_study = df_study %>% dplyr::filter( ! Image.ID.rsfmri %in% bad_rsfmri_images$IMAGE.ID.RSFMRI );

#
# Generate the initial candidate clinical study population
#
table( df_study$Visit );
tmp = list();
tmp[[1]] = df_study %>% dplyr::filter( Visit == "ADNI3 Year 1 Visit" );
tmp[[2]] = df_study %>% dplyr::filter( ! RID %in% unique( tmp[[1]]$RID ) );
df_study = rbindlist( tmp, use.names = TRUE );
nrow( df_study );

    #
    # Need to join with adnimerge
    #     Might be able to do it on image id? - NOPE - Weird.
    #     By Exam date?
    #
if ( 1 ) {
  
  #intersect( as.integer( df_study$Image.ID.mprage ), as.integer( adnimerge$IMAGEUID[ ! is.na( adnimerge$IMAGEUID ) ] ) );
  #intersect( as.integer( df_study$Image.ID.rfmri ), as.integer( adnimerge$IMAGEUID[ ! is.na( adnimerge$IMAGEUID ) ] ) );
  
    #
    # Explore only.
    #
  names( adnimerge )[ grep( "DATE", toupper( names( adnimerge ) ) ) ];
  names( df_study )[ grep( "DATE", toupper( names( df_study ) ) ) ];
  head( df_study$Study.Date ); head( df_study$SERIES_DATE.mprage ); head( df_study$SERIES_DATE.rsfmri );
  class( df_study$Study.Date ); class( df_study$SERIES_DATE.mprage ); class( df_study$SERIES_DATE.rsfmri );
  head( adnimerge_tmp$EXAMDATE );
  class( adnimerge_tmp$EXAMDATE );
  
    #
    # Owing to Explore section above take action.
    #
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
  
    #
    # Explore only.
    #
  head( df_study$RID );
  adnimerge_tmp %>% dplyr::filter( RID == 1261 );
  df_study %>% dplyr::filter( RID == 1261 );
  
    #
    # Conert convert df_study study date and adnimerge exam date to a common date format and find closest.
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
  df_study_vis_adj = df_study %>% dplyr::filter( abs(adnimerge_rec_dist) < 9*30 );
  adnimerge_tmp = adnimerge_tmp[ df_study_vis_adj$adnimerge_rec, ] %>% dplyr::select( - RID );
  
      #
      # Switch from df_study to df as main dataframe
      #
  df = cbind( df_study_vis_adj, adnimerge_tmp, stringsAsFactors = FALSE );
  names( df ) = toupper( names( df ) );
  
      #
      # Need minimum field fill rate to proceed.
      #   Have to have AGE, GENDER, PTEDUCAT, MMSE, CDRSB, DX
      #
  nrow( df );
  field_list_req = c( "AGE", "PTGENDER", "PTEDUCAT", "MMSE", "CDRSB", "DX", "YEARS.BL" );
  df = df[ complete.cases( df[ , field_list_req ] ), ];
  df$AGE_ON_EXAMDATE = df$AGE + df$YEARS.BL;
  
      #
      # Subset on AGE_ON_EXAMDATE - this may (or may not) happen when seeking to mtach classes.
      #
  if ( 0 ) {
    
    age_min_mci = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( AGE_ON_EXAMDATE ) %>% min();
    age_max_cn = df %>% dplyr::filter( DX == "CN" ) %>% dplyr::select( AGE_ON_EXAMDATE ) %>% max();
    df = df %>% dplyr::filter( AGE_ON_EXAMDATE >= age_min_mci, AGE_ON_EXAMDATE <= age_max_cn );
    
  } # if ( 0 )

      #
      # Write the data frame.
      #
  WriteExcelReport( df, "df_study_adnimerge_05July2019.xlsx" );
  
} # if ( 1 ) {

#
# This section should be done in interactive mode.
# Generate a string of image IDs that can be used to download from LONI.
#
if ( 0 ) {
  
    #
    # From earlier iterations, we know which are the "bad" RID, so just focus on those.
    #
  df_bad = df %>% dplyr::filter( SUBJECT.ID.RSFMRI %in% bad_rsfmri_images$RID );
  
  tmp1 = paste( "I", df_bad$IMAGE.ID.MPRAGE, sep = "" );
  write.table( tmp1, file = "/home/rstudio/ADNI3/re_sched.mprage.txt", row.names = FALSE, col.names = FALSE, quote = FALSE );
  
  tmp2 = paste( "I", df_bad$IMAGE.ID.RSFMRI, sep = "" );
  write.table( tmp2, file = "/home/rstudio/ADNI3/re_sched.rsfmri.txt", row.names = FALSE, col.names = FALSE, quote = FALSE );
  
  tmp3 = data.frame( RID = df_bad$SUBJECT.ID.MPRAGE, MPRAGE = tmp1, RSFMRI = tmp2, DX = df_bad$DX, stringsAsFactors = FALSE ) %>% arrange( RID );
  write.table( tmp3, file = "/home/rstudio/ADNI3/re_sched.mprage.rsfmri.dx.txt", row.names = FALSE, col.names = FALSE, quote = FALSE );
  
  paste( tmp1, collapse = "," );
  paste( tmp2, collapse = "," );
  
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
  anova_res = list();
  for ( ifield in field_list_anova ) {
    
    vlist = list( CN = df %>% dplyr::filter( DX == "CN" ) %>% dplyr::select( one_of( ifield ) ),
                  MCI = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( one_of( ifield ) ),
                  Dementia = df %>% dplyr::filter( DX == "Dementia" ) %>% dplyr::select( one_of( ifield ) ) );
    QuickPairWiseVisual( vlist[[1]], vlist[[2]], vlist[[3]], 0, "CN", "MCI", "DEMENTIA", FALSE, paste( ifield, " vs DX", sep = "" ), ifield, "Density" );
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




#####################
#####################
#
#
# Once the rsfmri files are converted to nii and json ingest the json files and inspect
#
#
#####################
#####################

library( rjson );

  #
  # Get the list of jsons to read.
  #
image_list = read.table( file = "/home/rstudio/ADNI3/sched.mprage.rsfmri.dx.txt", header = FALSE, stringsAsFactors = FALSE );
names( image_list ) = c( "SUBJ_ID", "MPRAGE_ID", "RSFMRI_ID", "DX" );

  #
  # Read the jsons
  #
json = list();
name_list = NULL;
for ( irec in 1:nrow( image_list ) ) {
  
  im = image_list$RSFMRI_ID[ irec ];
  json[[ im ]] = fromJSON( file = paste( "/home/rstudio/ADNI3/rsfmri/", im, "/", im, ".json", sep = "" ) );
  json[[ im ]][ "SUBJ_ID" ] = image_list$SUBJ_ID[ irec ];
  json[[ im ]][ "MPRAGE_ID" ] = as.character( gsub( "I", "", image_list$MPRAGE_ID[ irec ] ) );
  json[[ im ]][ "RSFMRI_ID" ] = as.character( gsub( "I", "", image_list$RSFMRI_ID[ irec ] ) );
  json[[ im ]][ "DXX" ] = image_list$DX[ irec ];
  name_list = c( name_list, names( json[[ im ]] ) );
  
} # for ( im in image_list$RSFMRI_ID ) {
tmp_name_stats = sort( table( name_list ), decreasing = TRUE );

  #
  # Assemble a data table.  Need to do this way since not every json file contains the same fields.
  #
df = data.frame( matrix( "", nrow = length( json ), ncol = length( tmp ) ), stringsAsFactors = FALSE ); names( df ) = names( tmp );
for ( irec in 1:length( json ) ) {
  
  for ( ifield in names( json[[ irec ]] ) ) {
    
    df[ irec, ifield ] = as.character( json[[ irec ]][ ifield ] );
    
  } # for ( ifield in names( json[[ irec ]] ) ) {
  
} # for ( iname in unique( name_list ) ) {

  #
  # Merge the JSON data with the schedule of images and other data.
  #
df_adni3 = read_excel( "../../ADNI3/df_study_adnimerge_20March2019.xlsx" );
df_adni3$IMAGE.ID.MPRAGE = as.character( df_adni3$IMAGE.ID.MPRAGE );
df_adni3$IMAGE.ID.RSFMRI = as.character( df_adni3$IMAGE.ID.RSFMRI );

  class( df_adni3$SUBJECT.ID.MPRAGE ); class( df_adni3$IMAGE.ID.MPRAGE ); class( df_adni3$IMAGE.ID.RSFMRI ); 
  class( df$SUBJ_ID ); class( df$MPRAGE_ID ); class( df$RSFMRI_ID );

df_tmp = df %>% select( SUBJ_ID, MPRAGE_ID, RSFMRI_ID, DXX );
df_adni3_tmp = df_adni3 %>% select( SUBJECT.ID.MPRAGE, IMAGE.ID.MPRAGE, IMAGE.ID.RSFMRI );

df_x = df_adni3 %>% left_join( df, by = c( "SUBJECT.ID.MPRAGE" = "SUBJ_ID", "IMAGE.ID.MPRAGE" = "MPRAGE_ID", "IMAGE.ID.RSFMRI" = "RSFMRI_ID" ) );

  #
  # Write the data frame.
  #

  #
  # After doing some viewing we know we don't want to proceed with: 019_S_4293
  #
df_x = df_x %>% dplyr::filter( SUBJECT.ID.MPRAGE != "019_S_4293" );
WriteExcelReport( df_x, "/home/rstudio/ADNI3/df_study_adnimerge_json_02April2019.xlsx" );
WriteCSVFile( df_x, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_02April2019.csv" );

  #
  # Explore subsetting the data, for example, on manufacturer.
  #
df_siemens = df_x %>% dplyr::filter( Manufacturer == "Siemens" );
table( df_siemens$DX ); table( df_siemens$DXX )

df_philips = df_x %>% dplyr::filter( Manufacturer == "Philips" );
table( df_philips$DX ); table( df_philips$DXX )

df_mb = df_x %>% dplyr::filter( ProtocolName == "Axial_MB_rsfMRI_(Eyes_Open)" );
table( df_mb$DX ); table( df_mb$DXX ); table( df_mb$Manufacturer );

df_x$SliceThickness = as.character( round( as.numeric( df_x$SliceThickness ), 1 ) );
table( df_x$ProtocolName, df_x$SliceThickness );

#
# Again generate Study Participant Demographic and Metadata Analysis
#
if ( 1 ) {
  
  df = data.frame( df_x );
  
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
                        names( df )[ grep( "ADAS", names( df ) ) ],
                        names(df)[ grep( "RAV", names( df ) ) ],
                        names(df)[ grep( "ECOG", names( df ) ) ] );
  anova_res = list();
  for ( ifield in field_list_anova ) {
    
    vlist = list( CN = df %>% dplyr::filter( DX == "CN" ) %>% dplyr::select( one_of( ifield ) ),
                  MCI = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( one_of( ifield ) ),
                  Dementia = df %>% dplyr::filter( DX == "Dementia" ) %>% dplyr::select( one_of( ifield ) ) );
    QuickPairWiseVisual( vlist[[1]], vlist[[2]], vlist[[3]], 0, "CN", "MCI", "DEMENTIA", FALSE, paste( ifield, " vs DX", sep = "" ), ifield, "Density" );
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

WriteExcelReport( df_x, "/home/rstudio/ADNI3/df_study_adnimerge_json_08April2019.xlsx" );
WriteCSVFile( df_x, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_08April2019.csv" );

#
# Process the stats from AFNI
#
afni_stats = ReadCSVFile( NULL, "//home/rstudio/ADNI3/Exp.2019.04.06/ssReview.Summary.txt" );
summary( afni_stats );

df = data.frame( df_x %>% left_join( afni_stats, by = c( "SUBJECT.ID.MPRAGE" = "RID" ) ), stringsAsFactors = FALSE );
df$DX = df$DX.y;

WriteExcelReport( df, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_all_08April2019.xlsx" );
WriteCSVFile( df, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_all_08April2019.csv" );

df = df %>% dplyr::filter( SUBJECT.ID.MPRAGE != "018_S_2133" );

WriteExcelReport( df, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_filt_08April2019.xlsx" );
WriteCSVFile( df, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_filt_08April2019.csv" );

save_df = df = ReadCSVFile( NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_filt_08April2019.csv" );

  #
  # Can we get CDRGlobal?
  #
df_cdr = cdr %>% dplyr::filter( COLPROT == "ADNI3", RID %in% df$RID ) %>% group_by( RID );
df_cdr$EXAMDATE_YEAR = as.numeric( substr( df_cdr$USERDATE, 1, 4 ) );
df_cdr$EXAMDATE_MONTH = as.numeric( substr( df_cdr$USERDATE, 6, 7 ) );
df_cdr$EXAMDATE_DAY = as.numeric( substr( df_cdr$USERDATE, 9, 10 ) );
df$EXAMDATE_YEAR = as.numeric( substr( df$EXAMDATE_AS_DATE, 1, 4 ) );
df$EXAMDATE_MONTH = as.numeric( substr( df$EXAMDATE_AS_DATE, 6, 7 ) );
df$EXAMDATE_DAY = as.numeric( substr( df$EXAMDATE_AS_DATE, 9, 10 ) );

df = df %>% inner_join( df_cdr, by = c( "RID", "EXAMDATE_YEAR" ) );

WriteExcelReport( df, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_filt_cdr_09April2019.xlsx" );
WriteCSVFile( df, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_filt_cdr_09April2019.csv" );

df = df %>% dplyr::filter( DICE_COEFF > 0.77, CENSOR_FRACTION < 0.25, MAX_CENSORED_DISPL < 3.0, AVG_TSNR > 100, DF_LEFT > 80, GLOBAL_CORR < 0.1,
                                AVG_CENSORED_MOTION < 1.0 );

df = df %>% dplyr::filter( ! ( DX == "CN" & CDGLOBAL != 0 ) );
df = df %>% dplyr::filter( ! ( DX == "CN" & CDRSB.x != 0 ) );
df = df %>% dplyr::filter( ! ( DX == "MCI" & CDGLOBAL != 0.5 ) );
df = df %>% dplyr::filter( ! ( DX == "MCI" & CDGLOBAL < 0.5 ) );

table( df$CDGLOBAL, df$CDRSB.x );
table( df$CDRSB.x, df$DX );
table( df$CDGLOBAL, df$DX );
table( df$MMSE, df$DX );

#
# Study Participant Demographic and Metadata Analysis
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
  field_list_anova = c( "AGE_ON_EXAMDATE", "PTEDUCAT", "CDRSB.x", "CDGLOBAL", "MMSE", "MOCA",
                        "AVG_CENSORED_MOTION", "MAX_CENSORED_DISPL", "DF_LEFT", "CENSOR_FRACTION", "AVG_TSNR", "GLOBAL_CORR", "DICE_COEFF",
                        names( df )[ grep( "ADAS", names( df ) ) ],
                        names(df)[ grep( "RAV", names( df ) ) ],
                        names(df)[ grep( "ECOG", names( df ) ) ] );
  anova_res = list();
  for ( ifield in field_list_anova ) {
    
    vlist = list( CN = df %>% dplyr::filter( DX == "CN" ) %>% dplyr::select( one_of( ifield ) ),
                  MCI = df %>% dplyr::filter( DX == "MCI" ) %>% dplyr::select( one_of( ifield ) ),
                  Dementia = df %>% dplyr::filter( DX == "Dementia" ) %>% dplyr::select( one_of( ifield ) ) );
    QuickPairWiseVisual( vlist[[1]], vlist[[2]], vlist[[3]], 0, "CN", "MCI", "DEMENTIA", FALSE, paste( ifield, " vs DX", sep = "" ), ifield, "Density" );
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

  #
  # Final clean up.
  #
df$ImageType = rep( "", nrow( df ) );
WriteExcelReport( df, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_09April2019.xlsx" );
WriteCSVFile( df, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_09April2019.csv" );

#df_test = ReadCSVFile( NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_09April2019.csv" );

tmp = df %>% dplyr::filter( DX == "CN" ) %>% select( SUBJECT.ID.MPRAGE, IMAGE_ID_MPRAGE, IMAGE_ID_RSFMRI, DX );
tmp = paste( tmp$SUBJECT.ID.MPRAGE, tmp$IMAGE_ID_MPRAGE, tmp$IMAGE_ID_RSFMRI, tmp$DX, sep = "." );
write.table( tmp, file = "/home/rstudio/ADNI3/subj_mprage_rsfmri_dx_09Apr2019.CN", row.names = FALSE, col.names = FALSE, quote = FALSE );

tmp = df %>% dplyr::filter( DX == "MCI" ) %>% select( SUBJECT.ID.MPRAGE, IMAGE_ID_MPRAGE, IMAGE_ID_MPRAGE, IMAGE_ID_RSFMRI, DX );
tmp = paste( tmp$SUBJECT.ID.MPRAGE, tmp$IMAGE_ID_MPRAGE, tmp$IMAGE_ID_RSFMRI, tmp$DX, sep = "." );
write.table( tmp, file = "/home/rstudio/ADNI3/subj_mprage_rsfmri_dx_09Apr2019.MCI", row.names = FALSE, col.names = FALSE, quote = FALSE );

tmp = df %>% dplyr::filter( DX == "Dementia" ) %>% select( SUBJECT.ID.MPRAGE, IMAGE_ID_MPRAGE, IMAGE_ID_MPRAGE, IMAGE_ID_RSFMRI, DX );
tmp = paste( tmp$SUBJECT.ID.MPRAGE, tmp$IMAGE_ID_MPRAGE, tmp$IMAGE_ID_RSFMRI, tmp$DX, sep = "." );
write.table( tmp, file = "/home/rstudio/ADNI3/subj_mprage_rsfmri_dx_09Apr2019.Dementia", row.names = FALSE, col.names = FALSE, quote = FALSE );

tmp = paste( df$SUBJECT.ID.MPRAGE, df$IMAGE_ID_MPRAGE, df$IMAGE_ID_RSFMRI, df$DX, sep = "." );
write.table( tmp, file = "/home/rstudio/ADNI3/subj_mprage_rsfmri_dx_09Apr2019.ALL", row.names = FALSE, col.names = FALSE, quote = FALSE );


#
# Bring in the FreeSurfer data.
#














