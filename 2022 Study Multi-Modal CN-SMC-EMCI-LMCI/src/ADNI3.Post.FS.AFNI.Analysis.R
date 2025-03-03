#
# ADNI3.Post.FS.AFNI.Analysis.R
#

#
# Modified 11-Dec-2021 for updated ADNI3
#

rm ( list = ls() );
assign( "last.warning", NULL, envir=baseenv() );

library( data.table );
library( ggplot2 );

library( plyr );
library( dplyr );

library( readxl );

source( "ADNI3.HelperFunctions.R" );

  #
  # Read in the starting study data schedule.
  #
df_study <- read_excel( "/adni3/exp0/df_study_adnimerge_21Nov2021.xlsx" )
head( df_study[ , c( "SUBJECT.ID.MPRAGE", "IMAGE.ID.MPRAGE", "IMAGE.ID.RSFMRI", "DX" ) ] )
df_study$rid_mprage_rsfmri_key <- paste( df_study$SUBJECT.ID.MPRAGE,
                                        paste( "I", df_study$IMAGE.ID.MPRAGE, sep = "" ),
                                        paste( "I", df_study$IMAGE.ID.RSFMRI, sep = "" ),
                                        df_study$DX, sep = "." )
names( df_study ) <- gsub( "\\.", "_", names( df_study ) )

  #
  # Read in the Freesurfer Segmentation table.
  #
df_aseg = read.table( "/adni3/exp0/aseg_stats_table.csv", header = TRUE, row.names = 1 );
names( df_aseg ) = gsub( "\\.", "_", names( df_aseg ) );
names( df_aseg ) = tolower( names( df_aseg ) );
df_aseg$rid_mprage_rsfmri_key = rownames( df_aseg );

      #
      # Do the aseg normalizations (left)
      #
ifield_list = names( df_aseg )[ grep( "left_", names( df_aseg ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "left_wm_hypointensities", "left_non_wm_hypointensities" ) ];
df_aseg_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_aseg_norm = cbind( df_aseg_norm, df_aseg[ , ifield_list[ ifield ] ] / df_aseg[ , "brainsegvolnotvent" ] );
  
} # for ( ifield in ifield_list ) {
df_aseg_norm = data.frame( df_aseg_norm );
names( df_aseg_norm ) = paste( ifield_list, "_aseg_norm", sep = "" );
head( df_aseg_norm );
df_aseg = cbind( df_aseg, df_aseg_norm );

    #
    # Do the aseg normalizations (right)
    #
ifield_list = names( df_aseg )[ grep( "right_", names( df_aseg ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "right_wm_hypointensities", "right_non_wm_hypointensities" ) ];
df_aseg_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_aseg_norm = cbind( df_aseg_norm, df_aseg[ , ifield_list[ ifield ] ] / df_aseg[ , "brainsegvolnotvent" ] );
  
} # for ( ifield in ifield_list ) {
df_aseg_norm = data.frame( df_aseg_norm );
names( df_aseg_norm ) = paste( ifield_list, "_aseg_norm", sep = "" );
head( df_aseg_norm );
df_aseg = cbind( df_aseg, df_aseg_norm );

  #
  # Read in the Freesurfer Parcellation tables.
  #
df_aparc_lh = read.table( "/adni3/exp0/aparc_stats_table_lh.txt", header = TRUE, row.names = 1 );
names( df_aparc_lh ) = gsub( "\\.", "_", names( df_aparc_lh ) );
names( df_aparc_lh ) = tolower( names( df_aparc_lh ) );
df_aparc_lh$rid_mprage_rsfmri_key = rownames( df_aparc_lh );

df_aparc_rh = read.table( "/adni3/exp0/aparc_stats_table_rh.txt", header = TRUE, row.names = 1 );
names( df_aparc_rh ) = gsub( "\\.", "_", names( df_aparc_rh ) );
names( df_aparc_rh ) = tolower( names( df_aparc_rh ) );
df_aparc_rh$rid_mprage_rsfmri_key = rownames( df_aparc_rh );
  
    #
    # Do a simple plot
    #
df = data.frame( left = df_aparc_lh$etiv, right = df_aparc_rh$etiv );
ggplot( df, aes( x=left, y=right) ) + geom_point()
identical( df$left, df$right )

df = data.frame( left = df_aparc_lh$brainsegvolnotvent, right = df_aparc_rh$brainsegvolnotvent );
ggplot( df, aes( x=left, y=right) ) + geom_point()
identical( df$left, df$right )

  #
  # Join df_aparc_lh, df_aparc_rh.
  #
df_aparc = full_join( df_aparc_lh, df_aparc_rh, by = c( "rid_mprage_rsfmri_key", "brainsegvolnotvent", "etiv" ) );

  #
  # Join df_aseg, df_aparc
  #   (First time - do some quick checks on names.)
  #
ncol( df_aseg ); ncol( df_aparc ); intersect( names( df_aseg ), names( df_aparc ) );
df_fs = full_join( df_aseg, df_aparc, by = c( "rid_mprage_rsfmri_key" ) );

  #
  # Read in the hippocampal sub-field results.
  #
df_subhipp = read.table( "/adni3/exp0/hipposf.csv", header = TRUE, row.names = 1 );
names( df_subhipp ) = gsub( "\\.", "_", names( df_subhipp ) );
names( df_subhipp ) = tolower( names( df_subhipp ) );
df_subhipp$rid_mprage_rsfmri_key = rownames( df_subhipp );

      #
      # Do the hippocampal sub-field normalizations (left)
      #
ifield_list = names( df_subhipp )[ grep( "left_", names( df_subhipp ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "left_whole_hippocampus", "left_hippocampal_fissure" ) ];
df_subhipp_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_subhipp_norm = cbind( df_subhipp_norm, df_subhipp[ , ifield_list[ ifield ] ] / df_subhipp[ , "left_whole_hippocampus" ] );
  
} # for ( ifield in ifield_list ) {
df_subhipp_norm = data.frame( df_subhipp_norm );
names( df_subhipp_norm ) = paste( ifield_list, "_subhipp_norm", sep = "" );
head( df_subhipp_norm );
df_subhipp = cbind( df_subhipp, df_subhipp_norm );

      #
      # Do the hippocampal sub-field normalizations (right)
      #
ifield_list = names( df_subhipp )[ grep( "right_", names( df_subhipp ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "right_whole_hippocampus", "right_hippocampal_fissure" ) ];
df_subhipp_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_subhipp_norm = cbind( df_subhipp_norm, df_subhipp[ , ifield_list[ ifield ] ] / df_subhipp[ , "right_whole_hippocampus" ] );
  
} # for ( ifield in ifield_list ) {
df_subhipp_norm = data.frame( df_subhipp_norm );
names( df_subhipp_norm ) = paste( ifield_list, "_subhipp_norm", sep = "" );
head( df_subhipp_norm );
df_subhipp = cbind( df_subhipp, df_subhipp_norm );

  #
  # Join df_fs (aparc + aseg) with df_subhipp
  #   Check keys first, just to be sure.
  #
df_fs_save = df_fs;
df_subhipp_save = df_subhipp
ncol( df_fs ); ncol( df_subhipp ); intersect( names( df_fs ), names( df_subhipp ) );
head(df_fs$rid_mprage_rsfmri_key); head(df_subhipp$rid_mprage_rsfmri_key);
df_fs$rid_only = substr( df_fs$rid_mprage_rsfmri_key, 1, 10 )
df_subhipp$rid_only = substr( df_subhipp$rid_mprage_rsfmri_key, 1, 10 )
icommon = intersect( df_fs$rid_only, df_subhipp$rid_only )
df_subhipp <- df_subhipp %>% filter( rid_only %in% icommon )
df_fs = full_join( df_fs, df_subhipp, by = c( "rid_only" ) );

  #
  # Read in the amygdala sub-field results.
  #
df_amyg = read.table( "/adni3/exp0/amygnucu.csv", header = TRUE, row.names = 1 );
names( df_amyg ) = gsub( "\\.", "_", names( df_amyg ) );
names( df_amyg ) = tolower( names( df_amyg ) );
df_amyg$rid_mprage_rsfmri_key = rownames( df_amyg );
df_amyg$rid_only = substr( df_amyg$rid_mprage_rsfmri_key, 1, 10 )

  #
  # Do the amygdala sub-field normalizations (left)
  #
ifield_list = names( df_amyg )[ grep( "left_", names( df_amyg ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "left_whole_amygdala" ) ];
df_amyg_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_amyg_norm = cbind( df_amyg_norm, df_amyg[ , ifield_list[ ifield ] ] / df_amyg[ , "left_whole_amygdala" ] );
  
} # for ( ifield in ifield_list ) {
df_amyg_norm = data.frame( df_amyg_norm );
names( df_amyg_norm ) = paste( ifield_list, "_amyg_norm", sep = "" );
head( df_amyg_norm );
df_amyg = cbind( df_amyg, df_amyg_norm );

#
# Do the amygdala sub-field normalizations (right)
#
ifield_list = names( df_amyg )[ grep( "right_", names( df_amyg ) ) ];
ifield_list = ifield_list[ ! ifield_list %in% c( "right_whole_amygdala" ) ];
df_amyg_norm = NULL;
for ( ifield in 1:length(ifield_list) ) {
  
  df_amyg_norm = cbind( df_amyg_norm, df_amyg[ , ifield_list[ ifield ] ] / df_amyg[ , "right_whole_amygdala" ] );
  
} # for ( ifield in ifield_list ) {
df_amyg_norm = data.frame( df_amyg_norm );
names( df_amyg_norm ) = paste( ifield_list, "_amyg_norm", sep = "" );
head( df_amyg_norm );
df_amyg = cbind( df_amyg, df_amyg_norm );

#
# Join df_fs (aparc + aseg) with df_amyg
#   Check keys first, just to be sure.
#
df_amyg_save = df_amyg;
ncol( df_fs ); ncol( df_amyg ); intersect( names( df_fs ), names( df_amyg ) );
head(df_fs$rid_mprage_rsfmri_key); head(df_subhipp$rid_mprage_rsfmri_key);
df_amyg$rid_only = substr( df_amyg$rid_mprage_rsfmri_key, 1, 10 )
icommon = intersect( df_fs$rid_only, df_amyg$rid_only )
df_amyg <- df_amyg %>% filter( rid_only %in% icommon )
df_fs = full_join( df_fs, df_amyg, by = c( "rid_only" ) );

  #
  # Read in AFNI results.  Adjust field names as prep for subsequent joins.
  #
df_afni = fread( "/adni3/exp0/ssreview_summary.txt", header = TRUE, quote = "" );
names( df_afni ) = gsub( "\\.", "_", names( df_afni ) );
names( df_afni ) = tolower( names( df_afni ) );
df_afni$rid_mprage_rsfmri_key = paste( df_afni$rid, df_afni$image_id_mprage, df_afni$image_id_rsfmri, df_afni$dx, sep = "." );
df_afni$rid_only = substr( df_afni$rid_mprage_rsfmri_key, 1, 10 )
df_afni = df_afni %>% select( -rid, -image_id_mprage, -image_id_rsfmri, -dx );

  #
  # Join df_fs and df_afni
  #   Do a final check to make sure no common field names.
  #
intersect( names( df_fs ), names( df_afni ) );
icommon = intersect( df_fs$rid_only, df_afni$rid_only )
df = full_join( df_fs, df_afni, by = c( "rid_only" ) );
dim( df );
df_save <- df;

  #
  # Join df_study
  #  Here we fix some sloppiness in the handling of rid_mprage_rsfmri_key
  #
df$rid_mprage_rsfmri_key = df$rid_mprage_rsfmri_key.x
df = left_join( df, df_study, by = c( "rid_mprage_rsfmri_key" ) );
dim( df );

  #
  # Write out the join results.
  #
names( df ) = tolower( names( df ) );
write.table( df, "/adni3/exp0//df_11Dec2021.csv", col.names = TRUE, quote = TRUE, row.names = FALSE );

  #
  #
  # BEGIN ANALYSIS
  #
  #
ilocs = is.na( df$df_left ); sum( ilocs );
if ( sum(ilocs) ) {
  df_missing = df[ ilocs, ]; nrow( df_missing );
  df = df[ ! ilocs, ]; nrow( df );
  write.table( df, "/adni3/exp0//df_11Dec2021_filtered.csv", col.names = TRUE, quote = TRUE, row.names = FALSE );
} # ilocs = is.na( df$df_left ); sum( ilocs );

    #
    # Do a simple plot
    #
df %>% select( contains( "presubiculum" ) ) %>%
  ggplot( aes( x=left_presubiculum_head_subhipp_norm, y=right_presubiculum_head_subhipp_norm ) ) + geom_point()

    #
    # Determine which to drop based on AFNI stats
    #
ifield_afni = c( "AVG_CENSORED_MOTION", "MAX_CENSORED_DISPL", "DF_LEFT", "CENSOR_FRACTION",
                 "AVG_TSNR", "GLOBAL_CORR", "DICE_COEFF", "DFUSED" );
ifield_afni = tolower( ifield_afni );

      #
      # Do a summary & plot
      #
df %>% select( all_of( ifield_afni ) )  %>% summary()
df %>% select( all_of( ifield_afni ) ) %>%
  reshape2::melt() %>%
  ggplot(aes(x=value)) + geom_histogram() + facet_wrap(~ variable, ncol=2, scales='free')
  
      #
      # Set some thresholds
      #
max_censored_displ_thresh = 3.0;
censor_fraction_thresh = 0.25;
global_corr_thresh = 0.05;
df_go = df %>% dplyr::filter( max_censored_displ < max_censored_displ_thresh &
                                censor_fraction < censor_fraction_thresh &
                                global_corr < global_corr_thresh &
                                dx != "Dementia");

      #
      # Now let's deal with comparison of dx.bl and dx.
      # Could have done this at earlier points, but OK.
      # Look at the table then set some rules for dx.bl vs dx
      #
df_go %>% select( dx_bl, dx ) %>% table()
df_go %>% select( dx_bl, dx ) %>% View()

      #
      # The rule:
      # Keep only: CN, CN; EMCI, MCI; LMCI, MCI; SMC, CN
      #
df_go <- df_go %>% dplyr::filter( dx_bl == "CN" & dx == "CN" |
                                    dx_bl == "EMCI" & dx == "MCI" |
                                    dx_bl == "LMCI" & dx == "MCI" |
                                    dx_bl == "SMC" & dx == "CN" )

    #
    # Do a summary & plot
    #
df_go %>% select( all_of( ifield_afni ) )  %>% summary()
df_go %>% select( all_of( ifield_afni ) ) %>%
  reshape2::melt() %>%
  ggplot(aes(x=value)) + geom_histogram() + facet_wrap(~ variable, ncol=2, scales='free')

df_go %>% select( dx_bl, dx ) %>% table()

      #
      # Write out the join results.
      #
write.table( df_go, "/adni3/exp0/df_go_11Dec2021.csv", col.names = TRUE, quote = TRUE, row.names = FALSE );


      #
      # NOTICE
      # RETURNING Analysis Entry Point
      # NOTICE
      #
df <- fread( "/adni3/exp0/df_go_11Dec2021.csv", header = TRUE, stringsAsFactors = FALSE );

      #
      # Check dx_bl, dx, etc.
      #
head( df$dx_bl ); table( df$dx_bl ); sum( table( df$dx_bl ) ); sum( is.na( df$dx_bl ) );
ilocs  = which( is.na( df$dx_bl ) );
labels_dx_bl = levels( factor( df$dx_bl ) );

head( df$dx ); table( df$dx ); sum( table( df$dx ) ); sum( is.na( df$dx ) );
ilocs  = which( is.na( df$dx ) );
labels_dx = levels( factor( df$dx ) );

table( df$dx_bl, df$dx );

#
# Clean up and then summarize
#  Get rid of columns that are all NA
#
df <- df %>% select( - all_of( names( which( apply( is.na( df ), 2, sum ) == nrow( df ) ) ) ) )
#df_summary <- df %>% mutate_all( min )
#write.table( df_summary %>% t(), "/adni3/exp0/df_go_summary_11Dec2021.csv", col.names = FALSE, quote = TRUE, row.names = TRUE );

#
# Do a quick ANOVA scan on everything.
#
all_qanova <- QuickANOVAScan( df )

#
# In this section, we iterate, because we may trim some subjects.
#
df_save <- df
df <- df %>% filter( ! rid_only == "130_S_6019" )

#
# We know that in ADNI3 there is a BASIC and an ADVANCED.
# We will want to see breakdowns that way, for reference.
#
df_time <- fread( "/adni3/exp0/time_stats.csv", header = FALSE, stringsAsFactors = FALSE )
names(df_time) <- c( "rid", "image_id_mprage", "image_id_rsfmri", "dx", "dx_bl", "num_time_steps", "time_step" )
df_time %>% select( num_time_steps, time_step ) %>% table()

  # For those with 3.0 second time step...how does QuickANOVAScan look? (N=134)
rid_list <- df_time %>% select( rid, time_step ) %>% filter( time_step == 3.0 ) %>% select( rid ) %>% unlist()
df_basic <- df %>% filter( rid_only %in% rid_list )
anova_basic <- df_basic %>% QuickANOVAScan()
df_basic %>% select( dx_bl, dx ) %>% table()
df_tmp = df_basic %>% select( subject_id_mprage, image_id_mprage, image_id_rsfmri, dx, dx_bl ) %>%
  mutate( image_id_mprage = paste( "I", image_id_mprage, sep = "" ),
          image_id_rsfmri = paste( "I", image_id_rsfmri, sep = "" ) );
write.table( df_tmp, file = "/adni3/exp0/df_go_basic_22Dec2021_id_list.csv", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE );

  # For those with 3.0 second time step...how does QuickANOVAScan look? (N=49)
rid_list <- df_time %>% select( rid, time_step ) %>% filter( time_step == 0.607 ) %>% select( rid ) %>% unlist()
df_adv <- df %>% filter( rid_only %in% rid_list )
anova_adv <- df_adv %>% QuickANOVAScan()
df_adv %>% select( dx_bl, dx ) %>% table()
df_tmp = df_adv %>% select( subject_id_mprage, image_id_mprage, image_id_rsfmri, dx, dx_bl ) %>%
  mutate( image_id_mprage = paste( "I", image_id_mprage, sep = "" ),
          image_id_rsfmri = paste( "I", image_id_rsfmri, sep = "" ) );
write.table( df_tmp, file = "/adni3/exp0/df_go_adv_22Dec2021_id_list.csv", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE );

#
# Confirm the final "schedule" of RID, MPRAGE, RSFMRI, DX, DX_BL
# Comes in handy when doing AFNI ROI
#
df_tmp = df %>% select( subject_id_mprage, image_id_mprage, image_id_rsfmri, dx, dx_bl ) %>%
            mutate( image_id_mprage = paste( "I", image_id_mprage, sep = "" ),
                    image_id_rsfmri = paste( "I", image_id_rsfmri, sep = "" ) );
write.table( df_tmp, file = "/adni3/exp0/df_go_22Dec2021_id_list.csv", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE );

##################################
##################################
#
# FINAL PUB Time Analyses
#
##################################
##################################

# Read in the base database of results
df <- fread( "/adni3/exp0/df_go_11Dec2021.csv", header = TRUE, stringsAsFactors = FALSE )
df_save <- df

#   Process basic
# Read in the final basic list - for which we found usable results
df_go_basic <- fread( "/adni3/exp0/df_go_basic_22Dec2021_id_list.csv", header = FALSE, stringsAsFactors = FALSE )
names( df_go_basic ) <- c( "subject_id_mprage", "image_id_mprage", "image_id_rsfmri", "dx", "dx_bl" )
anova_basic <- df %>% dplyr::filter( subject_id_mprage %in% df_go_basic$subject_id_mprage ) %>% QuickANOVAScan()
anova_basic[["df"]] <- df %>% dplyr::filter( subject_id_mprage %in% df_go_basic$subject_id_mprage )
save( anova_basic, file = "/adni3/exp0/df_go_basic_22Dec2021_id_list_anova.Rdat")

#   Process adv
# Read in the final adv list - for which we found usable results
df_go_adv <- fread( "/adni3/exp0/df_go_adv_22Dec2021_id_list.csv", header = FALSE, stringsAsFactors = FALSE )
names( df_go_adv ) <- c( "subject_id_mprage", "image_id_mprage", "image_id_rsfmri", "dx", "dx_bl" )
anova_adv <- df %>% dplyr::filter( subject_id_mprage %in% df_go_adv$subject_id_mprage ) %>% QuickANOVAScan()
anova_adv[["df"]] <- df %>% dplyr::filter( subject_id_mprage %in% df_go_adv$subject_id_mprage )
save( anova_adv, file = "/adni3/exp0/df_go_adv_22Dec2021_id_list_anova.Rdat")

