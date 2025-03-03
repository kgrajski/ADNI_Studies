#
# ADNI3.MPRAGEFileManagement.R
#

#
# More or less one-time up-front data staging script.
#

#
# 22March2019
# 29Nov2021 - modified
#

#
# 1.  Setting up the mprage directory structure
#
# In the shell run: find | grep nii
#   This locates all the nii files.
# Generate a list of image IDs accordingly.
#   Visually inspect the list to make sure there aren't any strays, such as "center/sample-001.nii.gz"
#

sched_mprage = read.csv( file = "/home/rstudio/ADNI3/sched.mprage.txt", header = FALSE, stringsAsFactors = FALSE );
df_tmp = read.csv( file = "/home/rstudio/ADNI3/rm.txt", header = FALSE, stringsAsFactors = FALSE );
rid_list_mprage = rep( "", nrow( df_tmp ) );
image_list_mprage = rep( "", nrow( df_tmp ) );
cmd_line_args = rep( "", nrow( df_tmp ) );
icount = 1;
for ( irec in 1:nrow( df_tmp ) ) {
  
      #
      # Get the image id
      #
  tmp = str_split( as.character( df_tmp[ irec, 1 ] ), pattern = "_" );
  tmp = gsub( ".nii", "", tmp[[ 1 ]][ length( tmp[[ 1 ]]) ] );
  if ( tmp %in% unlist( sched_mprage ) ) {
    
    image_list_mprage[ icount ] = paste( tmp[[ 1 ]][ length( tmp[[ 1 ]]) ], ".nii", sep = "" );
    
        #
        # Get the RID
        #
    tmp = str_split( as.character( df_tmp[ irec, 1 ] ), pattern = "/" );
    rid_list_mprage[ icount ] = tmp[[ 1 ]][ 3 ];
    
        #
        # Form the command line input arguments
        #
    cmd_line_args[ icount ] = paste( df_tmp[ irec, 1 ], rid_list_mprage[ icount ], image_list_mprage[ icount ], sep = " " );
    
    icount = icount + 1;
    
  } # if ( tmp %in% sched_mprage ) {

} # for ( irec in 1:nrow(tmp) ) {

rid_list_mprage = rid_list_mprage[ 1:( icount - 1) ];
image_list_mprage = image_list_mprage[ 1:( icount - 1) ];
cmd_line_args = cmd_line_args[ 1:( icount - 1) ];

write.table( cmd_line_args, file = "/home/rstudio/ADNI3/mprage.setupdir.txt", row.names = FALSE, col.names = FALSE, quote = FALSE );

#
# 2.  Shell script to rename and consolidate all mprage files into a single diretory (that FreeSurfer will read from).
# tcsh script here
#   Copy & Paste these commands in the directory from which the find above was run.
#
foreach line ( "`cat /home/rstudio/ADNI3/mprage.setupdir.txt`" )

  set argv = ( $line )
  set RAW = $1
  set RID = $2
  set IMAGE_ID = $3
  echo $RAW $RID $IMAGE_ID
  
  rm -rf mprage/$RID
  cp $RAW mprage/$RID.$IMAGE_ID

end
ls mprage | grep nii > file_list_recon_all.txt

#
# 3.  Run recon-all
#   Copy & Paste these commands in the directory from which the find above was run.
#
rm -f batch.txt
foreach line ( "`cat /adni3/exp0/df_study_adnimerge_21Nov2021.csv`" )

    set argv = ( $line )
    set ID_SUBJ = $1
    set ID_MPRAGE = $2
    set ID_RSFMRI = $3
    set DX = $4
  
    echo "recon-all -all -i /adni3/data/mprage/$ID_SUBJ.$ID_MPRAGE.nii -subject $ID_SUBJ.$ID_MPRAGE -sd /adni3/data/fsres &" >> batch.txt

end
cat batch.txt | sed -e 's/\.nii/_shft.nii/' > rm.txt
split -n 4 rm.txt batch.

#
# 4. Gather the results for freesurfer
#
quantifyHAsubregions.sh hippoSf T1 /adni3/exp0/hipposf.csv /adni3/data/fsres
quantifyHAsubregions.sh amygNuc T1 /adni3/exp0/amygnucu.csv /adni3/data/fsres

#
# For the rsfMRI data do similar operation to consolidate images into a simpler diretory.
#







