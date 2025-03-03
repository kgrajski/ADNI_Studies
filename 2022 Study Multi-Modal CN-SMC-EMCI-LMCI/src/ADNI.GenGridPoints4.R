#
#	ADNI.GenGridPoints4.R
#

#	13July2017 KAG
#	17July2017 KAG - File driven.
# 01May2018.  Revisiting for reproducibility...
#

#	Best to do setwd into an experiment directory since a lot of files could be generated here.

rm ( list = ls() );
assign( "last.warning", NULL, envir=baseenv() );

#


#
#
#
VecDiffEuc = function ( x, y ) {
  #	Testing: x = c(0, 0, 0); y = c(1, 1, 1);
  return( sqrt( sum( (x - y)^2 ) ) );
} # VecDiffEuc = function ( x, y ) {

#
#
#
GenGrid0 = function ( v0, gmax, gstep ) {
  
  #	Test: v0 = c(0, 0, 0); gmax = 36; gstep = 6;
  
  x0 = v0[1]; y0 = v0[2]; z0 = v0[3];
  xyz = NULL;
  
  for ( z in (z0 + seq( -gmax, gmax, gstep ) ) ) {
    
    for ( y in (y0 + seq( -gmax, gmax, gstep ) ) ) {
      
      for ( x in (x0 + seq( -gmax, gmax, gstep ) ) ) {
        
        xyz = rbind( xyz, c( x, y, z ) );
        
      } # for ( x in (x0 + grange) ) {
      
    } # for ( y in (y0 + grange) ) {
    
  } # for ( z in (z0 + grange) ) {
  
  return( xyz );
  
} # GenGrid0 = function ( v0, gmax, rmin ) {

GenGrid = function ( xmin, xmax, ymin, ymax, zmin, zmax, gstep ) {
  
  #	Test: v0 = c(0, 0, 0); gmax = 36; gstep = 6;
  
  xyz = NULL;
  
  for ( z in seq( zmin, zmax, gstep ) ) {
    
    for ( y in seq( ymin, ymax, gstep ) ) {
      
      for ( x in seq( xmin, xmax, gstep ) ) {
        
        xyz = rbind( xyz, c( x, y, z ) );
        
      } # for ( x in seq( xmin, xmax, gstep ) ) ) {
      
    } # for ( y in seq( ymin, ymax, gstep ) ) ) {
    
  } # for ( z in seq( zmin, zmax, gstep ) ) ) {
  
  return( xyz );
  
} # GenGrid = function ( v0, gmax, rmin ) {


#
#	Begin
#

#
# Switch into the working directory.
#

setwd("/home/rstudio/ADNI3/Exp_30June2019/Masks/ROI.anchor.candidates")

#	Constants
gstep = 3;


#	Read in the table of cluster results.  These are NOT in MNI.
#cRpt = read.table( file="v55.1.Grid.4")

file_list = c( "bbox.mtl.list.txt", "bbox.amyg.list.txt" )
for ( ifile in file_list ) {
  
  bboxes = read.table( file=paste( "../", ifile, sep = "" ) )
  colnames( bboxes ) = c("roi", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax" );
  
  #	Generate candidate ROI seeds
  
  numPoints = 0;
  xyz_all = NULL;
  for ( iROI in 1:nrow(bboxes ) ) {
    
    xyz = round( GenGrid( bboxes$xmin[iROI], bboxes$xmax[iROI], bboxes$ymin[iROI], bboxes$ymax[iROI], bboxes$zmin[iROI], bboxes$zmax[iROI], gstep ), 0 );
    numPoints = numPoints + nrow( xyz );
    cat( " ** ROI = ", levels(bboxes[ ,1 ])[ iROI ], " N ROI candidate anchors =", nrow( xyz ), "\n" );
    xyz_all = rbind( xyz_all, xyz );
    
    if ( 1 ) {
      for ( iPt in 1:nrow(xyz) ) {
        
        tmp = data.frame( bboxes$roi[ iROI ], xyz[ iPt, 1 ], xyz[ iPt, 2 ], xyz[ iPt, 3 ] );
        write.table( tmp, file=paste( levels(bboxes[ ,1 ])[ iROI ], xyz[ iPt, 1 ], xyz[ iPt, 2 ], xyz[ iPt, 3 ], sep="." ), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t" );
        
      } # for ( iPt = 1:numPoints )
    } # if ( 0 )
    
  } # for ( iROI in roiNameList ) {
  
  
} # for ( ifile in file_list ) {



