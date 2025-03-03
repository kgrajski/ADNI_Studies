#
#	File: ADNI.ProcFreeSurferASeg.R
#
rm ( list = ls() );
assign( "last.warning", NULL, envir=baseenv() );

library ( plyr );
library ( dplyr );
library ( effsize );
library ( Hmisc );
library ( psych );

source("ADNI.ROI.PairAnalysis.HelperFunctions10.R");
source( "ADNI3.HelperFunctions.R" );

	######################################
	#
	#	Helper Functions
	#
	######################################
	
	#
	#	ANOVABatch
	#	Assumes that col #1 is RID and removes it.
	#

ANOVABatch = function ( xDat.CN, xDat.MCI, xDat.Dementia, q ) {
	#   Testing: ANOVABatch ( aseg.CN, aseg.MCI, aseg.Dementia );
	#	Testing: xDat.CN = aseg.CN; xDat.MCI = aseg.MCI; xDat.Dementia = aseg.Dementia
	#	Testing: xDat.CN = aparc.CN; xDat.MCI = aparc.MCI; xDat.Dementia = aparc.Dementia
  # Testing: xDat.CN = hpcsf.CN; xDat.MCI = hpcsf.MCI; xDat.Dementia = hpcsf.Dementia;
	
	num.CN = nrow( xDat.CN ); num.MCI = nrow( xDat.MCI ); num.Dementia = nrow( xDat.Dementia );
	numRec = num.CN + num.MCI + num.Dementia;

	xDat.CN = xDat.CN %>% select( - RID ); xDat.MCI = xDat.MCI %>% select( - RID ); xDat.Dementia = xDat.Dementia %>% select( - RID );
	xDat.ALL = rbind( xDat.CN, xDat.MCI, xDat.Dementia );
	
	numTests = ncol( xDat.ALL );
	numDXPairs = 3;
	xDat.anova = list();
	xDat.anova.Fstat = rep ( NA, numTests );
	xDat.anova.pVal = rep ( NA, numTests );
	xDat.m = array( NA, c ( numTests, numDXPairs ) );
	xDat.se = array( NA, c ( numTests, numDXPairs ) );
	xDat.t = array( NA, c ( numTests, numDXPairs ) );
	xDat.p = array( NA, c ( numTests, numDXPairs ) );
	xDat.d = array( NA, c ( numTests, numDXPairs ) );
	z.Labels = c ( rep ( AddressDXLabel(1), num.CN ) );
	z.Labels = c ( z.Labels, rep ( AddressDXLabel(2), num.MCI ) );
	z.Labels = c ( z.Labels, rep ( AddressDXLabel(3), num.Dementia ) );

	iCount = 1;
	for ( iTest in 1:numTests ) {
	
		z.xDat = xDat.CN[ , iTest ];
		z.xDat = c( z.xDat, xDat.MCI[ , iTest ] );
		z.xDat = c( z.xDat, xDat.Dementia[ , iTest ] );
		
		xDat.m[ iTest, 1 ] = mean( xDat.CN[ , iTest ] );
		xDat.m[ iTest, 2 ] = mean( xDat.MCI[ , iTest ] );
		xDat.m[ iTest, 3 ] = mean( xDat.Dementia[ , iTest ] );
	
		xDat.se[ iTest, 1 ] = sqrt( var( xDat.CN[ , iTest ] ) ) / sqrt( num.CN );
		xDat.se[ iTest, 2 ] = sqrt( var( xDat.MCI[ , iTest ] ) ) / sqrt( num.MCI );
		xDat.se[ iTest, 3 ] = sqrt( var( xDat.Dementia[ , iTest ] ) ) / sqrt( num.Dementia );
	
		tmp = PairwiseTTests ( z.xDat, z.Labels, FALSE );
		xDat.t[ iTest, ] = tmp$tStats;
		xDat.p[ iTest, ] = c( unlist((tmp$t1)["p.value"] ), unlist((tmp$t2)["p.value"] ), unlist((tmp$t3)["p.value"] ) );
		xDat.d[ iTest, ] = PairwiseCohenDTests ( z.xDat, z.Labels, FALSE );
	
		xDat.anova[[iCount]] = aov ( z.xDat ~ z.Labels, data = data.frame( z.xDat, z.Labels ) );
		xDat.anova.Fstat[iCount] = unlist(summary(xDat.anova[[iCount]]))["F value1"];
		xDat.anova.pVal[iCount] = unlist(summary(xDat.anova[[iCount]]))["Pr(>F)1"];
	
		iCount = iCount + 1;
			
	} # for ( iTest in 1:numTests ) {
	anovaRes = round( data.frame ( xDat.anova.Fstat, xDat.anova.pVal, xDat.m, xDat.t, xDat.p, xDat.d, xDat.se ), 5 );
	rownames(anovaRes) = names(xDat.ALL);
	colnames(anovaRes) = c("Fstat", "pVal", "m.CN", "m.MCI", "m.D",
	                       "t.MCI.CN", "t.D.CN", "t.D.MCI",
	                       "p.MCI.CN", "p.D.CN", "p.D.MCI",
	                       "d.MCI.CN", "d.D.CN", "d.D.MCI",
	                       "se.CN", "se.MCI", "se.D" );
	tmp = sort.int( anovaRes[,"Fstat"], decreasing=TRUE, index.return=TRUE );
	anovaRes = anovaRes[ tmp$ix, ];

	#	Apply correction if any.
	
	fdr.pVal = FDRCorrectedPValue ( anovaRes$pVal, q );
	cat(" FDR pVal=", fdr.pVal, "\n" );
	sigRes = which( anovaRes$pVal <= fdr.pVal );
	anovaRes[ sigRes, ]
	
	return ( list( anovaRes=anovaRes, sigRes=sigRes, fdr.pVal=fdr.pVal ) );

} # ANOVABatch = function ( x.CN, x.MCI, x.Dementia ) {

	
GroupAndScatterComboPlot = function ( sigVarName, pRes, d.CN, d.MCI, d.Dementia ) {
	
	yLabText = "Volume (mm^3)"
	numDX = 3;
	mainText = paste( sigVarName, sep=" " );	
	tmp = round( pRes[9:11], 4 );
	pValText = paste( "P: ", tmp[1], "; ", tmp[2], "; ", tmp[3], sep="" );
	tmp = round( pRes[12:14], 3 );
	dValText = paste( "d: ", tmp[1], "; ", tmp[2], "; ", tmp[3], sep="" );
	mainText = paste ( mainText, pValText, dValText, sep="\n" );

		#	Figure 1.0.  Group Averages and SEs

	ga.d = data.frame( x  = c(1:numDX), y = c(pRes[3:5]), se = c(pRes[15:17]) );
	ga.xlim = c ( 1, 3.5 );
	ga.ylim = c ( 	min( pRes[3:5] ) - max( pRes[15:17] ), max( pRes[3:5] ) + max( pRes[15:17] ) )
	with ( data = ga.d, expr = plot ( x, y, xlim=ga.xlim, ylim=ga.ylim, xaxp = c(1, numDX, 2 ), type="l", lty=1, xaxt="n", col=1,
				main=mainText, xlab="Diagnosis", ylab=yLabText ) );
	axis( 1, at=c(1:3), labels=c("CN", "MCI", "D") );
	with ( data = ga.d, expr = ( errbar( x, y, y+se, y-se, add=T, cap=.015 ) ) );

		#	Figure 2.0.  Group scatter plots.
		
	tmp = c ( d.CN, d.MCI, d.Dementia);
	rid.scat.xlim = c ( 0.5, 3.5 );
	rid.scat.ylim =  c ( min(tmp, na.rm=TRUE), max(tmp, na.rm=TRUE) );
	rid.d = list(); rid.d[[1]] = d.CN; rid.d[[2]] = d.MCI; rid.d[[3]] = d.Dementia;

	plot ( c(rid.scat.ylim[1], rid.scat.ylim[2]), type="n", xlim=rid.scat.xlim, ylim=rid.scat.ylim, xaxt="n", main=mainText, xlab="Diagnosis", ylab=yLabText );

	for ( iDX in 1:numDX ) {
		tmp.rid.d = data.frame( x = rep( iDX, length(rid.d[[iDX]])), y  = rid.d[[iDX]]);
		with ( data = tmp.rid.d, expr = points ( jitter(x), y, ylim = rid.scat.ylim, xaxp = c(1, numDX, 2 ), type="p", pch=iDX, xaxt="n" ) );
		axis( 1, at=c(1:3), labels=c("CN", "MCI", "D") );
	} # 	for ( iDX in 1:iDX ) {
	
} # GroupAndScatterComboPlot = function ( sigVarName, pRes, d.CN, d.MCI, d.Dementia ) {
	
	
	
SpecialGroupAndScatterComboPlot = function ( sigVarName, pRes, d.CN, d.MCI, d.Dementia, sigVarName.Plus, pRes.Plus, d.CN.Plus, d.MCI.Plus, d.Dementia.Plus,
												flexTitleText1=NULL, flexTitleText2=NULL, flexYLabText=NULL  ) {
		# Test: sigVarName=sigVarName1; pRes=unlist(anovaRes[ sigVarLoc1, ]); d.CN=unlist(aseg.CN[,sigVarName1]); d.MCI=unlist(aseg.MCI[,sigVarName1]); d.Dementia=unlist(aseg.Dementia[,sigVarName1]);
		# Test: sigVarName.Plus=sigVarName2; pRes.Plus=unlist(anovaRes[ sigVarLoc2, ]); d.CN.Plus=unlist(aseg.CN[,sigVarName2]); d.MCI.Plus=unlist(aseg.MCI[,sigVarName2]); d.Dementia.Plus=unlist(aseg.Dementia[,sigVarName2]);

	yLabText = "Volume (mm^3)"
	if ( length( flexYLabText) ) { yLabText = flexYLabText; }
	numDX = 3;
	tmpText = "Group Avg Vol Est vs CDR\n";
	if ( length( flexTitleText1 ) ) { tmpText = flexTitleText1; }
	mainText = paste( tmpText, sigVarName, " (Solid), ", sigVarName.Plus, " (Dash)", sep="" );

		#	Figure 1.0.  Group Averages and SEs

	ga.d = data.frame( x  = c(1:numDX), y = c(pRes[3:5]), se = c(pRes[15:17]) );
	ga.d.Plus = data.frame( x  = 0.25 + c(1:numDX), y = c(pRes.Plus[3:5]), se = c(pRes.Plus[15:17]) );
	ga.xlim = c ( 1, 3.25 );
	ga.ylim = c ( min( pRes[3:5], pRes.Plus[3:5] ) - max( pRes[15:17], pRes.Plus[15:17] ), max( pRes[3:5], pRes.Plus[3:5] ) + max( pRes[15:17], pRes.Plus[15:17] ) );
	
	with ( data = ga.d, expr = plot ( x, y, xlim=ga.xlim, ylim=ga.ylim, xaxp = c(1, numDX, 2 ), type="l", lty=1, xaxt="n", col=1,
				main=mainText, xlab="", ylab=yLabText ) );
	with ( data = ga.d, expr = ( errbar( x, y, y+se, y-se, add=T, cap=.015 ) ) );
	
	with ( data = ga.d.Plus, expr = lines ( x, y, xlim=ga.xlim, ylim=ga.ylim, xaxp = c(1, numDX, 2 ), type="l", lty=2, xaxt="n", col=1,
				main=mainText, xlab="", ylab=yLabText ) );
	with ( data = ga.d.Plus, expr = ( errbar( x, y, y+se, y-se, add=T, cap=.015 ) ) );
	axis( 1, at=c(1:3), labels=c("CN", "MCI", "D") );

		#	Figure 2.0.  Group scatter plots.
		
	tmp = c ( d.CN, d.MCI, d.Dementia);
	tmp.Plus = c ( d.CN.Plus, d.MCI.Plus, d.Dementia.Plus );
	
	rid.scat.xlim = c ( 0.5, 3.5 );
	rid.scat.ylim =  c ( min(tmp, tmp.Plus, na.rm=TRUE), max(tmp, tmp.Plus, na.rm=TRUE) );
	
	rid.d = list(); rid.d[[1]] = d.CN; rid.d[[2]] = d.MCI; rid.d[[3]] = d.Dementia;
	rid.d.Plus = list(); rid.d.Plus[[1]] = d.CN.Plus; rid.d.Plus[[2]] = d.MCI.Plus; rid.d.Plus[[3]] = d.Dementia.Plus;

	tmpText = "Subject Vol Est vs CDR\n";
	if ( length( flexTitleText2) ) { tmpText = flexTitleText2; }
	mainText = paste( tmpText, sigVarName, " (Black), ", sigVarName.Plus, " (Blue)", sep="" );
	plot ( c(rid.scat.ylim[1], rid.scat.ylim[2]), type="n", xlim=rid.scat.xlim, ylim=rid.scat.ylim, xaxt="n", main=mainText, xlab="", ylab="" );

	for ( iDX in 1:numDX ) {
		tmp.rid.d = data.frame( x = rep( iDX, length(rid.d[[iDX]])), y  = rid.d[[iDX]]);
		with ( data = tmp.rid.d, expr = points ( jitter(x), y, ylim = rid.scat.ylim, xaxp = c(1, numDX, 2 ), type="p", pch=iDX, xaxt="n" ) );
		axis( 1, at=c(1:3), labels=c("CN", "MCI", "D") );
	} # 	for ( iDX in 1:iDX ) {
		
	for ( iDX in 1:numDX ) {
		tmp.rid.d = data.frame( x = rep( iDX, length(rid.d.Plus[[iDX]])), y  = rid.d.Plus[[iDX]]);
		with ( data = tmp.rid.d, expr = points ( 0.25 + jitter(x), y, ylim = rid.scat.ylim, xaxp = c(1, numDX, 2 ), type="p", pch=iDX, xaxt="n", col=4 ) );
	} # 	for ( iDX in 1:iDX ) {
	
} # SpecialGroupAndScatterComboPlot = function ( sigVarName, pRes, d.CN, d.MCI, d.Dementia ) {

	######################################
	#
	#	MAIN
	#
	######################################

q.fdr = 0.10;
df_study = ReadCSVFile( NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_09April2019.csv" );
df_study$ID = paste( df_study$SUBJECT.ID.MPRAGE, df_study$IMAGE_ID_MPRAGE, df_study$IMAGE_ID_RSFMRI, df_study$DX, sep = "." );
save_df_study = df_study;

	####################
	#	
	#	HPC Subfields
	#
	####################

df_hpcsf = ReadCSVFile( NULL, "/home/rstudio/ADNI3/Exp.2019.04.06/hippocampal_subfield_table.csv" );

df_study = df_study %>% left_join( df_hpcsf, by = c( "ID" = "Subject" ) );
WriteExcelReport( df_study, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_09April2019.xlsx" );
WriteCSVFile( df_study, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_09April2019.csv" );

hpcsf.CN = df_study %>% dplyr::filter( DX == "CN" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.CN = dim(hpcsf.CN)[1];
numVar = dim(hpcsf.CN)[2];
for ( iRec in 1:nrow( hpcsf.CN ) ) {
	hpcsf.CN[ iRec, 2:13 ] = 100.0 * ( hpcsf.CN[ iRec,2:13 ] / ( hpcsf.CN$left_Whole_hippocampus[ iRec ] ) );
	hpcsf.CN[ iRec, 15:26 ] = 100.0 * ( hpcsf.CN[ iRec,15:26 ] / ( hpcsf.CN$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

hpcsf.MCI = df_study %>% dplyr::filter( DX == "MCI" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.MCI = dim(hpcsf.MCI)[1];
numVar = dim(hpcsf.MCI)[2];
for ( iRec in 1:nrow( hpcsf.MCI ) ) {
  hpcsf.MCI[ iRec, 2:13 ] = 100.0 * ( hpcsf.MCI[ iRec,2:13 ] / ( hpcsf.MCI$left_Whole_hippocampus[ iRec ] ) );
  hpcsf.MCI[ iRec, 15:26 ] = 100.0 * ( hpcsf.MCI[ iRec,15:26 ] / ( hpcsf.MCI$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

hpcsf.Dementia = df_study %>% dplyr::filter( DX == "Dementia" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.Dementia = dim(hpcsf.Dementia)[1];
numVar = dim(hpcsf.Dementia)[2];
for ( iRec in 1:nrow( hpcsf.Dementia ) ) {
  hpcsf.Dementia[ iRec, 2:13 ] = 100.0 * ( hpcsf.Dementia[ iRec,2:13 ] / ( hpcsf.Dementia$left_Whole_hippocampus[ iRec ] ) );
  hpcsf.Dementia[ iRec, 15:26 ] = 100.0 * ( hpcsf.Dementia[ iRec,15:26 ] / ( hpcsf.Dementia$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

hpcsf.CN.save = hpcsf.CN;
hpcsf.MCI.save = hpcsf.MCI;
hpcsf.Dementia.save = hpcsf.Dementia;
hpcsf.ALL.save = rbind( hpcsf.CN, hpcsf.MCI, hpcsf.Dementia );

field_list_drop = c( "left_Whole_hippocampus", "right_Whole_hippocampus", "left_hippocampal-fissure", "right_hippocampal-fissure" );
hpcsf.CN = hpcsf.CN %>% select( - one_of( field_list_drop ) );
hpcsf.MCI = hpcsf.MCI %>% select( - one_of( field_list_drop ) );
hpcsf.Dementia = hpcsf.Dementia %>% select( - one_of( field_list_drop ) );
hpcsf.ALL = rbind( hpcsf.CN, hpcsf.MCI, hpcsf.Dementia );

tmp = ANOVABatch ( hpcsf.CN, hpcsf.MCI, hpcsf.Dementia, q.fdr );
tmp$anovaRes[ tmp$sigRes, ];

tmp = ANOVABatch ( hpcsf.CN %>% select( RID, contains( "left_" ) ),
                   hpcsf.MCI  %>% select( RID, contains( "left_" ) ),
                   hpcsf.Dementia  %>% select( RID, contains( "left_" ) ), q.fdr );
tmp$anovaRes[ tmp$sigRes, ];

tmp = ANOVABatch ( hpcsf.CN %>% select( RID, contains( "right_" ) ),
                   hpcsf.MCI  %>% select( RID, contains( "right_" ) ),
                   hpcsf.Dementia  %>% select( RID, contains( "right_" ) ), q.fdr );
tmp$anovaRes[ tmp$sigRes, ];


	####################
	#	
	#	ASEG
	#
	####################

df_aseg = rbind( read.table( NULL, "/home/rstudio/ADNI3/Exp.2019.04.06/aseg_table.CN" ),
                 ReadCSVFile( NULL, "/home/rstudio/ADNI3/Exp.2019.04.06/aseg_table.MCI" ),
                 ReadCSVFile( NULL, "/home/rstudio/ADNI3/Exp.2019.04.06/aseg_table.Dementia" ) );

df_study = df_study %>% left_join( df_aseg, by = c( "ID" = "Subject" ) );
WriteExcelReport( df_study, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_aesg_10April2019.xlsx" );
WriteCSVFile( df_study, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_aseg_10April2019.csv" );

aseg.CN = df_study %>% dplyr::filter( DX == "CN" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.CN = dim(aseg.CN)[1];
numVar = dim(aseg.CN)[2];
for ( iRec in 1:nrow( aseg.CN ) ) {
  aseg.CN[ iRec, 2:13 ] = 100.0 * ( aseg.CN[ iRec,2:13 ] / ( aseg.CN$left_Whole_hippocampus[ iRec ] ) );
  aseg.CN[ iRec, 15:26 ] = 100.0 * ( aseg.CN[ iRec,15:26 ] / ( aseg.CN$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

aseg.MCI = df_study %>% dplyr::filter( DX == "MCI" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.MCI = dim(aseg.MCI)[1];
numVar = dim(aseg.MCI)[2];
for ( iRec in 1:nrow( aseg.MCI ) ) {
  aseg.MCI[ iRec, 2:13 ] = 100.0 * ( aseg.MCI[ iRec,2:13 ] / ( aseg.MCI$left_Whole_hippocampus[ iRec ] ) );
  aseg.MCI[ iRec, 15:26 ] = 100.0 * ( aseg.MCI[ iRec,15:26 ] / ( aseg.MCI$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

aseg.Dementia = df_study %>% dplyr::filter( DX == "Dementia" ) %>% select( RID, starts_with( "left_" ), starts_with( "right_" ) );
num.Dementia = dim(aseg.Dementia)[1];
numVar = dim(aseg.Dementia)[2];
for ( iRec in 1:nrow( aseg.Dementia ) ) {
  aseg.Dementia[ iRec, 2:13 ] = 100.0 * ( aseg.Dementia[ iRec,2:13 ] / ( aseg.Dementia$left_Whole_hippocampus[ iRec ] ) );
  aseg.Dementia[ iRec, 15:26 ] = 100.0 * ( aseg.Dementia[ iRec,15:26 ] / ( aseg.Dementia$right_Whole_hippocampus[ iRec ] ) );
} # for ( iRec in 1:numRec ) {

aseg.CN.save = aseg.CN;
aseg.MCI.save = aseg.MCI;
aseg.Dementia.save = aseg.Dementia;
aseg.ALL.save = rbind( aseg.CN, aseg.MCI, aseg.Dementia );

field_list_drop = c( "left_Whole_hippocampus", "right_Whole_hippocampus", "left_hippocampal-fissure", "right_hippocampal-fissure" );
aseg.CN = aseg.CN %>% select( - one_of( field_list_drop ) );
aseg.MCI = aseg.MCI %>% select( - one_of( field_list_drop ) );
aseg.Dementia = aseg.Dementia %>% select( - one_of( field_list_drop ) );
aseg.ALL = rbind( aseg.CN, aseg.MCI, aseg.Dementia );

tmp = ANOVABatch ( aseg.CN, aseg.MCI, aseg.Dementia, q.fdr );
tmp$anovaRes[ tmp$sigRes, ];

tmp = ANOVABatch ( aseg.CN %>% select( RID, contains( "left_" ) ),
                   aseg.MCI  %>% select( RID, contains( "left_" ) ),
                   aseg.Dementia  %>% select( RID, contains( "left_" ) ), q.fdr );
tmp$anovaRes[ tmp$sigRes, ];

tmp = ANOVABatch ( aseg.CN %>% select( RID, contains( "right_" ) ),
                   aseg.MCI  %>% select( RID, contains( "right_" ) ),
                   aseg.Dementia  %>% select( RID, contains( "right_" ) ), q.fdr );
tmp$anovaRes[ tmp$sigRes, ];


	####################
	#	
	#	APARC
	#
	####################
	
normalizeFlag = TRUE;
tiffFlag = FALSE;
iMeasure = "area"
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.CN", sep="." );
	aparc.CN.lh = read.table( file=fName.lh, header=TRUE );
	aparc.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.CN.lh); tmp[1] = "RID"; names(aparc.CN.lh) = tmp;
	tmp = names(aparc.CN.rh); tmp[1] = "RID"; names(aparc.CN.rh) = tmp;
	dim(aparc.CN.lh); dim(aparc.CN.rh); 
	num.CN = dim(aparc.CN.lh)[1];
	numVar = dim(aparc.CN.lh)[2];
	aparc.CN = merge ( aparc.CN.lh, aparc.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.MCI", sep="." );
	aparc.MCI.lh = read.table( file=fName.lh, header=TRUE );
	aparc.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.MCI.lh); tmp[1] = "RID"; names(aparc.MCI.lh) = tmp;
	tmp = names(aparc.MCI.rh); tmp[1] = "RID"; names(aparc.MCI.rh) = tmp;
	dim(aparc.MCI.lh); dim(aparc.MCI.rh); 
	num.MCI = dim(aparc.MCI.lh)[1];
	numVar = dim(aparc.MCI.lh)[2];
	aparc.MCI = merge ( aparc.MCI.lh, aparc.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.Dementia", sep="." );
	aparc.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	aparc.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.Dementia.lh); tmp[1] = "RID"; names(aparc.Dementia.lh) = tmp;
	tmp = names(aparc.Dementia.rh); tmp[1] = "RID"; names(aparc.Dementia.rh) = tmp;
	dim(aparc.Dementia.lh); dim(aparc.Dementia.rh); 
	num.Dementia = dim(aparc.Dementia.lh)[1];
	numVar = dim(aparc.Dementia.lh)[2];
	aparc.Dementia = merge ( aparc.Dementia.lh, aparc.Dementia.rh, by="RID" );
	
	useTheseCols = c( 1:35, 37:70 );	#	Do not use the rh_WhiteSurfArea_area and lh_WhiteSurfArea_area values.
	aparc.CN = aparc.CN[ , useTheseCols ]; aparc.MCI = aparc.MCI[ , useTheseCols ]; aparc.Dementia = aparc.Dementia[ , useTheseCols ];
	set.lh = 2:35; set.rh = 36:69;
	if ( 0 ) {
		for ( iRec in 1:num.CN ) {
			aparc.CN[ iRec, set.lh ]  = 100.0 * ( aparc.CN[ iRec, set.lh ] / sum( aparc.CN[ iRec, set.lh ] ) );
			aparc.CN[ iRec, set.rh ] = 100.0 * ( aparc.CN[ iRec, set.rh ] / sum( aparc.CN[ iRec, set.rh ] ) );
		}
		for ( iRec in 1:num.MCI ) {
			aparc.MCI[ iRec, set.lh ]  = 100.0 * ( aparc.MCI[ iRec, set.lh ] / sum( aparc.MCI[ iRec, set.lh ] ) );
			aparc.MCI[ iRec, set.rh ] = 100.0 * ( aparc.MCI[ iRec, set.rh ] / sum( aparc.MCI[ iRec, set.rh ] ) );
		}
		for ( iRec in 1:num.Dementia ) {
			aparc.Dementia[ iRec, set.lh ]  = 100.0 * ( aparc.Dementia[ iRec, set.lh ] / sum( aparc.Dementia[ iRec, set.lh ] ) );
			aparc.Dementia[ iRec, set.rh ] = 100.0 * ( aparc.Dementia[ iRec, set.rh ] / sum( aparc.Dementia[ iRec, set.rh ] ) );
		}
	} # if ( normalizeFlag ) {
		
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			aparc.CN[ iRec, 2:69 ]  = 100.0 * ( aparc.CN[ iRec, 2:69  ] / sum( aparc.CN[ iRec, 2:69  ] ) );
		}
		for ( iRec in 1:num.MCI ) {
			aparc.MCI[ iRec, 2:69  ]  = 100.0 * ( aparc.MCI[ iRec, 2:69  ] / sum( aparc.MCI[ iRec, 2:69  ] ) );
		}
		for ( iRec in 1:num.Dementia ) {
			aparc.Dementia[ iRec, 2:69  ]  = 100.0 * ( aparc.Dementia[ iRec, 2:69  ] / sum( aparc.Dementia[ iRec, 2:69  ] ) );
		}
	} # if ( normalizeFlag ) {
	
	tmp = ANOVABatch ( aparc.CN, aparc.MCI, aparc.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
iMeasure = "volume"
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.CN", sep="." );
	aparc.CN.lh = read.table( file=fName.lh, header=TRUE );
	aparc.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.CN.lh); tmp[1] = "RID"; names(aparc.CN.lh) = tmp;
	tmp = names(aparc.CN.rh); tmp[1] = "RID"; names(aparc.CN.rh) = tmp;
	dim(aparc.CN.lh); dim(aparc.CN.rh); 
	num.CN = dim(aparc.CN.lh)[1];
	numVar = dim(aparc.CN.lh)[2];
	aparc.CN = merge ( aparc.CN.lh, aparc.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.MCI", sep="." );
	aparc.MCI.lh = read.table( file=fName.lh, header=TRUE );
	aparc.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.MCI.lh); tmp[1] = "RID"; names(aparc.MCI.lh) = tmp;
	tmp = names(aparc.MCI.rh); tmp[1] = "RID"; names(aparc.MCI.rh) = tmp;
	dim(aparc.MCI.lh); dim(aparc.MCI.rh); 
	num.MCI = dim(aparc.MCI.lh)[1];
	numVar = dim(aparc.MCI.lh)[2];
	aparc.MCI = merge ( aparc.MCI.lh, aparc.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.Dementia", sep="." );
	aparc.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	aparc.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.Dementia.lh); tmp[1] = "RID"; names(aparc.Dementia.lh) = tmp;
	tmp = names(aparc.Dementia.rh); tmp[1] = "RID"; names(aparc.Dementia.rh) = tmp;
	dim(aparc.Dementia.lh); dim(aparc.Dementia.rh); 
	num.Dementia = dim(aparc.Dementia.lh)[1];
	numVar = dim(aparc.Dementia.lh)[2];
	aparc.Dementia = merge ( aparc.Dementia.lh, aparc.Dementia.rh, by="RID" );
	
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			aparc.CN[ iRec, 2:69 ]  = 100.0 * ( aparc.CN[ iRec, 2:69  ] / sum( aparc.CN[ iRec, 2:69  ] ) );
		}
		for ( iRec in 1:num.MCI ) {
			aparc.MCI[ iRec, 2:69  ]  = 100.0 * ( aparc.MCI[ iRec, 2:69  ] / sum( aparc.MCI[ iRec, 2:69  ] ) );
		}
		for ( iRec in 1:num.Dementia ) {
			aparc.Dementia[ iRec, 2:69  ]  = 100.0 * ( aparc.Dementia[ iRec, 2:69  ] / sum( aparc.Dementia[ iRec, 2:69  ] ) );
		}
	} # if ( normalizeFlag ) {
	
	tmp = ANOVABatch ( aparc.CN, aparc.MCI, aparc.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
iMeasure = "thickness"
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.CN", sep="." );
	aparc.CN.lh = read.table( file=fName.lh, header=TRUE );
	aparc.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.CN.lh); tmp[1] = "RID"; names(aparc.CN.lh) = tmp;
	tmp = names(aparc.CN.rh); tmp[1] = "RID"; names(aparc.CN.rh) = tmp;
	dim(aparc.CN.lh); dim(aparc.CN.rh); 
	num.CN = dim(aparc.CN.lh)[1];
	numVar = dim(aparc.CN.lh)[2];
	aparc.CN = merge ( aparc.CN.lh, aparc.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.MCI", sep="." );
	aparc.MCI.lh = read.table( file=fName.lh, header=TRUE );
	aparc.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.MCI.lh); tmp[1] = "RID"; names(aparc.MCI.lh) = tmp;
	tmp = names(aparc.MCI.rh); tmp[1] = "RID"; names(aparc.MCI.rh) = tmp;
	dim(aparc.MCI.lh); dim(aparc.MCI.rh); 
	num.MCI = dim(aparc.MCI.lh)[1];
	numVar = dim(aparc.MCI.lh)[2];
	aparc.MCI = merge ( aparc.MCI.lh, aparc.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.Dementia", sep="." );
	aparc.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	aparc.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.Dementia.lh); tmp[1] = "RID"; names(aparc.Dementia.lh) = tmp;
	tmp = names(aparc.Dementia.rh); tmp[1] = "RID"; names(aparc.Dementia.rh) = tmp;
	dim(aparc.Dementia.lh); dim(aparc.Dementia.rh); 
	num.Dementia = dim(aparc.Dementia.lh)[1];
	numVar = dim(aparc.Dementia.lh)[2];
	aparc.Dementia = merge ( aparc.Dementia.lh, aparc.Dementia.rh, by="RID" );
	
	#colsOfInterest = c (1:dim());
	#aparc.CN = aparc.CN[ , colsOfInterest ];
	#aparc.MCI = aparc.MCI[ , colsOfInterest ];
	#aparc.Dementia = aparc.Dementia[ , colsOfInterest ];
	
	tmp = ANOVABatch ( aparc.CN, aparc.MCI, aparc.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
	anovaRes = tmp$anovaRes;
	sigVarName1 = "lh_entorhinal_thickness";
	sigVarLoc1 = which( rownames( anovaRes ) == sigVarName1 );
	sigVarName2 = "rh_entorhinal_thickness";
	sigVarLoc2 = which( rownames( anovaRes ) == sigVarName2 );
	if ( tiffFlag ) {
		tiff ( filename=paste("v44.4.FreeSurfer.aparc.ECSummary", "tiff", sep="."), compression="lzw" );
	} else {
		x11(); 
	} # if ( tiffFlag )
	par(mfcol=c(1,2), pty="s");
	SpecialGroupAndScatterComboPlot ( "L-EC", unlist(anovaRes[ sigVarLoc1, ]), unlist(aparc.CN[,sigVarName1]),
										unlist(aparc.MCI[,sigVarName1]), unlist(aparc.Dementia[,sigVarName1]),
									"R-EC", unlist(anovaRes[ sigVarLoc2, ]), unlist(aparc.CN[,sigVarName2]),
										unlist(aparc.MCI[,sigVarName2]), unlist(aparc.Dementia[,sigVarName2]),
										"Group Avg EC Thickness vs CDR\n", "Subject EC Thickness vs CDR\n", "Thickness (mm)" );
	if ( tiffFlag ) { dev.off(); } # if ( tiffFlag )

    # Fstat    pVal    m.CN   m.MCI     m.D p.MCI.CN  p.D.CN p.D.MCI d.MCI.CN  d.D.CN d.D.MCI   se.CN  sd.MCI    sd.D
	print(tmp$anovaRes[ tmp$sigRes, c("Fstat", "pVal", "m.CN", "p.MCI.CN", "d.MCI.CN", "m.MCI", "p.D.CN", "d.D.CN", "m.D", "p.D.MCI", "d.D.MCI" )] )
	
	itmp = which(tmp$anovaRes[ tmp$sigRes, "m.MCI" ] > tmp$anovaRes[ tmp$sigRes, "m.CN" ])
	print(tmp$anovaRes[ itmp, c("Fstat", "pVal", "m.CN", "p.MCI.CN", "d.MCI.CN", "m.MCI", "p.D.CN", "d.D.CN", "m.D", "p.D.MCI", "d.D.MCI" )] )

	itmp = which(tmp$anovaRes[ tmp$sigRes, "m.D" ] > tmp$anovaRes[ tmp$sigRes, "m.CN" ])
	print(tmp$anovaRes[ itmp, c("Fstat", "pVal", "m.CN", "p.MCI.CN", "d.MCI.CN", "m.MCI", "p.D.CN", "d.D.CN", "m.D", "p.D.MCI", "d.D.MCI" )] )
	

	####################
	#	
	#	aparc.DKTatlas40
	#
	####################
measureList = c ( "volume", "area", "thickness" );
for ( iMeasure in measureList ) {
	
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	normalizeFlag = FALSE;
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.DKTatlas40.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.DKTatlas40.CN", sep="." );
	aparc.DKTatlas40.CN.lh = read.table( file=fName.lh, header=TRUE );
	aparc.DKTatlas40.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.DKTatlas40.CN.lh); tmp[1] = "RID"; names(aparc.DKTatlas40.CN.lh) = tmp;
	tmp = names(aparc.DKTatlas40.CN.rh); tmp[1] = "RID"; names(aparc.DKTatlas40.CN.rh) = tmp;
	dim(aparc.DKTatlas40.CN.lh); dim(aparc.DKTatlas40.CN.rh); 
	num.CN = dim(aparc.DKTatlas40.CN.lh)[1];
	numVar = dim(aparc.DKTatlas40.CN.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			tmp.lh = sum( aparc.DKTatlas40.CN.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.DKTatlas40.CN.rh[2:(numVar-1)]);
			aparc.DKTatlas40.CN.lh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.CN.lh[iRec,2:numVar] / tmp.lh );
			aparc.DKTatlas40.CN.rh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.CN.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.DKTatlas40.CN = merge ( aparc.DKTatlas40.CN.lh, aparc.DKTatlas40.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.DKTatlas40.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.DKTatlas40.MCI", sep="." );
	aparc.DKTatlas40.MCI.lh = read.table( file=fName.lh, header=TRUE );
	aparc.DKTatlas40.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.DKTatlas40.MCI.lh); tmp[1] = "RID"; names(aparc.DKTatlas40.MCI.lh) = tmp;
	tmp = names(aparc.DKTatlas40.MCI.rh); tmp[1] = "RID"; names(aparc.DKTatlas40.MCI.rh) = tmp;
	dim(aparc.DKTatlas40.MCI.lh); dim(aparc.DKTatlas40.MCI.rh); 
	num.MCI = dim(aparc.DKTatlas40.MCI.lh)[1];
	numVar = dim(aparc.DKTatlas40.MCI.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.MCI ) {
			tmp.lh = sum( aparc.DKTatlas40.MCI.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.DKTatlas40.MCI.rh[2:(numVar-1)]);
			aparc.DKTatlas40.MCI.lh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.MCI.lh[iRec,2:numVar] / tmp.lh );
			aparc.DKTatlas40.MCI.rh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.MCI.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.DKTatlas40.MCI = merge ( aparc.DKTatlas40.MCI.lh, aparc.DKTatlas40.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.DKTatlas40.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.DKTatlas40.Dementia", sep="." );
	aparc.DKTatlas40.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	aparc.DKTatlas40.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.DKTatlas40.Dementia.lh); tmp[1] = "RID"; names(aparc.DKTatlas40.Dementia.lh) = tmp;
	tmp = names(aparc.DKTatlas40.Dementia.rh); tmp[1] = "RID"; names(aparc.DKTatlas40.Dementia.rh) = tmp;
	dim(aparc.DKTatlas40.Dementia.lh); dim(aparc.DKTatlas40.Dementia.rh); 
	num.Dementia = dim(aparc.DKTatlas40.Dementia.lh)[1];
	numVar = dim(aparc.DKTatlas40.Dementia.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.Dementia ) {
			tmp.lh = sum( aparc.DKTatlas40.Dementia.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.DKTatlas40.Dementia.rh[2:(numVar-1)]);
			aparc.DKTatlas40.Dementia.lh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.Dementia.lh[iRec,2:numVar] / tmp.lh );
			aparc.DKTatlas40.Dementia.rh[iRec,2:numVar] = 100.0 * ( aparc.DKTatlas40.Dementia.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.DKTatlas40.Dementia = merge ( aparc.DKTatlas40.Dementia.lh, aparc.DKTatlas40.Dementia.rh, by="RID" );
	tmp = ANOVABatch ( aparc.DKTatlas40.CN, aparc.DKTatlas40.MCI, aparc.DKTatlas40.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );

} # for ( iMeasure in measureList )

	####################
	#	
	#	aparc.a2009s
	#
	####################
	
measureList = c ( "volume", "thickness", "area" );
for ( iMeasure in measureList ) {
	
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	normalizeFlag = FALSE;
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.a2009s.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.a2009s.CN", sep="." );
	aparc.a2009s.CN.lh = read.table( file=fName.lh, header=TRUE );
	aparc.a2009s.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.a2009s.CN.lh); tmp[1] = "RID"; names(aparc.a2009s.CN.lh) = tmp;
	tmp = names(aparc.a2009s.CN.rh); tmp[1] = "RID"; names(aparc.a2009s.CN.rh) = tmp;
	dim(aparc.a2009s.CN.lh); dim(aparc.a2009s.CN.rh); 
	num.CN = dim(aparc.a2009s.CN.lh)[1];
	numVar = dim(aparc.a2009s.CN.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			tmp.lh = sum( aparc.a2009s.CN.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.a2009s.CN.rh[2:(numVar-1)]);
			aparc.a2009s.CN.lh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.CN.lh[iRec,2:numVar] / tmp.lh );
			aparc.a2009s.CN.rh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.CN.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.a2009s.CN = merge ( aparc.a2009s.CN.lh, aparc.a2009s.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.a2009s.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.a2009s.MCI", sep="." );
	aparc.a2009s.MCI.lh = read.table( file=fName.lh, header=TRUE );
	aparc.a2009s.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.a2009s.MCI.lh); tmp[1] = "RID"; names(aparc.a2009s.MCI.lh) = tmp;
	tmp = names(aparc.a2009s.MCI.rh); tmp[1] = "RID"; names(aparc.a2009s.MCI.rh) = tmp;
	dim(aparc.a2009s.MCI.lh); dim(aparc.a2009s.MCI.rh); 
	num.MCI = dim(aparc.a2009s.MCI.lh)[1];
	numVar = dim(aparc.a2009s.MCI.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.MCI ) {
			tmp.lh = sum( aparc.a2009s.MCI.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.a2009s.MCI.rh[2:(numVar-1)]);
			aparc.a2009s.MCI.lh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.MCI.lh[iRec,2:numVar] / tmp.lh );
			aparc.a2009s.MCI.rh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.MCI.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.a2009s.MCI = merge ( aparc.a2009s.MCI.lh, aparc.a2009s.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "aparc.a2009s.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "aparc.a2009s.Dementia", sep="." );
	aparc.a2009s.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	aparc.a2009s.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(aparc.a2009s.Dementia.lh); tmp[1] = "RID"; names(aparc.a2009s.Dementia.lh) = tmp;
	tmp = names(aparc.a2009s.Dementia.rh); tmp[1] = "RID"; names(aparc.a2009s.Dementia.rh) = tmp;
	dim(aparc.a2009s.Dementia.lh); dim(aparc.a2009s.Dementia.rh); 
	num.Dementia = dim(aparc.a2009s.Dementia.lh)[1];
	numVar = dim(aparc.a2009s.Dementia.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.Dementia ) {
			tmp.lh = sum( aparc.a2009s.Dementia.lh[2:(numVar-1)]);
			tmp.rh = sum( aparc.a2009s.Dementia.rh[2:(numVar-1)]);
			aparc.a2009s.Dementia.lh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.Dementia.lh[iRec,2:numVar] / tmp.lh );
			aparc.a2009s.Dementia.rh[iRec,2:numVar] = 100.0 * ( aparc.a2009s.Dementia.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	aparc.a2009s.Dementia = merge ( aparc.a2009s.Dementia.lh, aparc.a2009s.Dementia.rh, by="RID" );
	tmp = ANOVABatch ( aparc.a2009s.CN, aparc.a2009s.MCI, aparc.a2009s.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
} # for ( iMeasure in measureList )

	####################
	#	
	#	entorhinal_exvivo
	#
	####################

measureList = c ( "volume", "area", "thickness" );
for ( iMeasure in measureList ) {
	
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );
	normalizeFlag = FALSE;
	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "entorhinal_exvivo.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "entorhinal_exvivo.CN", sep="." );
	entorhinal_exvivo.CN.lh = read.table( file=fName.lh, header=TRUE );
	entorhinal_exvivo.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(entorhinal_exvivo.CN.lh); tmp[1] = "RID"; names(entorhinal_exvivo.CN.lh) = tmp;
	tmp = names(entorhinal_exvivo.CN.rh); tmp[1] = "RID"; names(entorhinal_exvivo.CN.rh) = tmp;
	dim(entorhinal_exvivo.CN.lh); dim(entorhinal_exvivo.CN.rh); 
	num.CN = dim(entorhinal_exvivo.CN.lh)[1];
	numVar = dim(entorhinal_exvivo.CN.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			tmp.lh = sum( entorhinal_exvivo.CN.lh[2:(numVar-1)]);
			tmp.rh = sum( entorhinal_exvivo.CN.rh[2:(numVar-1)]);
			entorhinal_exvivo.CN.lh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.CN.lh[iRec,2:numVar] / tmp.lh );
			entorhinal_exvivo.CN.rh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.CN.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	entorhinal_exvivo.CN = merge ( entorhinal_exvivo.CN.lh, entorhinal_exvivo.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "entorhinal_exvivo.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "entorhinal_exvivo.MCI", sep="." );
	entorhinal_exvivo.MCI.lh = read.table( file=fName.lh, header=TRUE );
	entorhinal_exvivo.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(entorhinal_exvivo.MCI.lh); tmp[1] = "RID"; names(entorhinal_exvivo.MCI.lh) = tmp;
	tmp = names(entorhinal_exvivo.MCI.rh); tmp[1] = "RID"; names(entorhinal_exvivo.MCI.rh) = tmp;
	dim(entorhinal_exvivo.MCI.lh); dim(entorhinal_exvivo.MCI.rh); 
	num.MCI = dim(entorhinal_exvivo.MCI.lh)[1];
	numVar = dim(entorhinal_exvivo.MCI.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.MCI ) {
			tmp.lh = sum( entorhinal_exvivo.MCI.lh[2:(numVar-1)]);
			tmp.rh = sum( entorhinal_exvivo.MCI.rh[2:(numVar-1)]);
			entorhinal_exvivo.MCI.lh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.MCI.lh[iRec,2:numVar] / tmp.lh );
			entorhinal_exvivo.MCI.rh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.MCI.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	entorhinal_exvivo.MCI = merge ( entorhinal_exvivo.MCI.lh, entorhinal_exvivo.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "entorhinal_exvivo.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "entorhinal_exvivo.Dementia", sep="." );
	entorhinal_exvivo.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	entorhinal_exvivo.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(entorhinal_exvivo.Dementia.lh); tmp[1] = "RID"; names(entorhinal_exvivo.Dementia.lh) = tmp;
	tmp = names(entorhinal_exvivo.Dementia.rh); tmp[1] = "RID"; names(entorhinal_exvivo.Dementia.rh) = tmp;
	dim(entorhinal_exvivo.Dementia.lh); dim(entorhinal_exvivo.Dementia.rh); 
	num.Dementia = dim(entorhinal_exvivo.Dementia.lh)[1];
	numVar = dim(entorhinal_exvivo.Dementia.lh)[2];
	if ( normalizeFlag ) {
		for ( iRec in 1:num.Dementia ) {
			tmp.lh = sum( entorhinal_exvivo.Dementia.lh[2:(numVar-1)]);
			tmp.rh = sum( entorhinal_exvivo.Dementia.rh[2:(numVar-1)]);
			entorhinal_exvivo.Dementia.lh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.Dementia.lh[iRec,2:numVar] / tmp.lh );
			entorhinal_exvivo.Dementia.rh[iRec,2:numVar] = 100.0 * ( entorhinal_exvivo.Dementia.rh[iRec,2:numVar] / tmp.rh );
		} # for ( iRec in 1:numRec ) {
	} # if ( normalizeFlag )
	entorhinal_exvivo.Dementia = merge ( entorhinal_exvivo.Dementia.lh, entorhinal_exvivo.Dementia.rh, by="RID" );
	tmp = ANOVABatch ( entorhinal_exvivo.CN, entorhinal_exvivo.MCI, entorhinal_exvivo.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
} # for ( iMeasure in measureList )

anovaRes = tmp$anovaRes;
sigVarName1 = "lh_lh.entorhinal_exvivo.label_thickness";
sigVarLoc1 = which( rownames( anovaRes ) == sigVarName1 );
sigVarName2 = "rh_rh.entorhinal_exvivo.label_thickness";
sigVarLoc2 = which( rownames( anovaRes ) == sigVarName2 );
if ( tiffFlag ) {
		tiff ( filename=paste("v44.4.FreeSurfer.ECSummary", "tiff", sep="."), compression="lzw" );
} else {
		x11(); 
} # if ( tiffFlag )
par(mfcol=c(1,2), pty="s");
SpecialGroupAndScatterComboPlot ( "L-EC", unlist(anovaRes[ sigVarLoc1, ]), unlist(entorhinal_exvivo.CN[,sigVarName1]),
										unlist(entorhinal_exvivo.MCI[,sigVarName1]), unlist(entorhinal_exvivo.Dementia[,sigVarName1]),
									"R-EC", unlist(anovaRes[ sigVarLoc2, ]), unlist(entorhinal_exvivo.CN[,sigVarName2]),
										unlist(entorhinal_exvivo.MCI[,sigVarName2]), unlist(entorhinal_exvivo.Dementia[,sigVarName2]),
										"Group Avg EC Thickness vs CDR\n", "Subject EC Thickness vs CDR\n", "Thickness (mm)" );
if ( tiffFlag ) { dev.off(); } # if ( tiffFlag )



		#	Junk yard
		
if ( 0 ) {
	
	
	####################
	#	
	#	BA
	#
	####################

iMeasure = c ( "volume" );
normalizeFlag = TRUE;
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.CN", sep="." );
	BA.CN.lh = read.table( file=fName.lh, header=TRUE );
	BA.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.CN.lh); tmp[1] = "RID"; names(BA.CN.lh) = tmp;
	tmp = names(BA.CN.rh); tmp[1] = "RID"; names(BA.CN.rh) = tmp;
	dim(BA.CN.lh); dim(BA.CN.rh); 
	num.CN = dim(BA.CN.lh)[1];
	BA.CN = merge ( BA.CN.lh, BA.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.MCI", sep="." );
	BA.MCI.lh = read.table( file=fName.lh, header=TRUE );
	BA.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.MCI.lh); tmp[1] = "RID"; names(BA.MCI.lh) = tmp;
	tmp = names(BA.MCI.rh); tmp[1] = "RID"; names(BA.MCI.rh) = tmp;
	dim(BA.MCI.lh); dim(BA.MCI.rh); 
	num.MCI = dim(BA.MCI.lh)[1];
	BA.MCI = merge ( BA.MCI.lh, BA.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.Dementia", sep="." );
	BA.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	BA.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.Dementia.lh); tmp[1] = "RID"; names(BA.Dementia.lh) = tmp;
	tmp = names(BA.Dementia.rh); tmp[1] = "RID"; names(BA.Dementia.rh) = tmp;
	dim(BA.Dementia.lh); dim(BA.Dementia.rh); 
	num.Dementia = dim(BA.Dementia.lh)[1];
	BA.Dementia = merge ( BA.Dementia.lh, BA.Dementia.rh, by="RID" );

	
	tmp = ANOVABatch ( BA.CN, BA.MCI, BA.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
iMeasure = c ( "area" );
normalizeFlag = TRUE;
	cat( "\n\n\n*************", iMeasure, "*************", "\n" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.CN", sep="." );
	BA.CN.lh = read.table( file=fName.lh, header=TRUE );
	BA.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.CN.lh); tmp[1] = "RID"; names(BA.CN.lh) = tmp;
	tmp = names(BA.CN.rh); tmp[1] = "RID"; names(BA.CN.rh) = tmp;
	dim(BA.CN.lh); dim(BA.CN.rh); 
	num.CN = dim(BA.CN.lh)[1];
	BA.CN = merge ( BA.CN.lh, BA.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.MCI", sep="." );
	BA.MCI.lh = read.table( file=fName.lh, header=TRUE );
	BA.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.MCI.lh); tmp[1] = "RID"; names(BA.MCI.lh) = tmp;
	tmp = names(BA.MCI.rh); tmp[1] = "RID"; names(BA.MCI.rh) = tmp;
	dim(BA.MCI.lh); dim(BA.MCI.rh); 
	num.MCI = dim(BA.MCI.lh)[1];
	BA.MCI = merge ( BA.MCI.lh, BA.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.Dementia", sep="." );
	BA.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	BA.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.Dementia.lh); tmp[1] = "RID"; names(BA.Dementia.lh) = tmp;
	tmp = names(BA.Dementia.rh); tmp[1] = "RID"; names(BA.Dementia.rh) = tmp;
	dim(BA.Dementia.lh); dim(BA.Dementia.rh); 
	num.Dementia = dim(BA.Dementia.lh)[1];
	BA.Dementia = merge ( BA.Dementia.lh, BA.Dementia.rh, by="RID" );
	
	useTheseCols = c ( 1:14, 16:28 );
	BA.CN = BA.CN[ , useTheseCols ];
	BA.MCI = BA.MCI[ , useTheseCols ];
	BA.Dementia = BA.Dementia[ , useTheseCols ];
	numVar = 27;
	if ( normalizeFlag ) {
		for ( iRec in 1:num.CN ) {
			BA.CN[ iRec, 2:numVar ] = 100.0 * ( BA.CN[ iRec, 2:numVar ] / sum( BA.CN[ iRec, 2:numVar ] ) );
		}
		for ( iRec in 1:num.MCI ) {
			BA.MCI[ iRec, 2:numVar ] = 100.0 * ( BA.MCI[ iRec, 2:numVar ] / sum( BA.MCI[ iRec, 2:numVar ] ) );
		}
		for ( iRec in 1:num.Dementia ) {
			BA.Dementia[ iRec, 2:numVar ] = 100.0 * ( BA.Dementia[ iRec, 2:numVar ] / sum( BA.Dementia[ iRec, 2:numVar ] ) );
		}
	} # if ( normalizeFlag ) {
	
	tmp = ANOVABatch ( BA.CN, BA.MCI, BA.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
iMeasure = c ( "thickness" );

	cat( "\n\n\n*************", iMeasure, "*************", "\n" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.CN", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.CN", sep="." );
	BA.CN.lh = read.table( file=fName.lh, header=TRUE );
	BA.CN.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.CN.lh); tmp[1] = "RID"; names(BA.CN.lh) = tmp;
	tmp = names(BA.CN.rh); tmp[1] = "RID"; names(BA.CN.rh) = tmp;
	dim(BA.CN.lh); dim(BA.CN.rh); 
	num.CN = dim(BA.CN.lh)[1];
	BA.CN = merge ( BA.CN.lh, BA.CN.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.MCI", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.MCI", sep="." );
	BA.MCI.lh = read.table( file=fName.lh, header=TRUE );
	BA.MCI.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.MCI.lh); tmp[1] = "RID"; names(BA.MCI.lh) = tmp;
	tmp = names(BA.MCI.rh); tmp[1] = "RID"; names(BA.MCI.rh) = tmp;
	dim(BA.MCI.lh); dim(BA.MCI.rh); 
	num.MCI = dim(BA.MCI.lh)[1];
	BA.MCI = merge ( BA.MCI.lh, BA.MCI.rh, by="RID" );

	fName.lh = paste( "v44.4.matched.rid.image.lh", iMeasure, "BA.Dementia", sep="." );
	fName.rh = paste( "v44.4.matched.rid.image.rh", iMeasure, "BA.Dementia", sep="." );
	BA.Dementia.lh = read.table( file=fName.lh, header=TRUE );
	BA.Dementia.rh = read.table( file=fName.rh, header=TRUE );
	tmp = names(BA.Dementia.lh); tmp[1] = "RID"; names(BA.Dementia.lh) = tmp;
	tmp = names(BA.Dementia.rh); tmp[1] = "RID"; names(BA.Dementia.rh) = tmp;
	dim(BA.Dementia.lh); dim(BA.Dementia.rh); 
	num.Dementia = dim(BA.Dementia.lh)[1];
	BA.Dementia = merge ( BA.Dementia.lh, BA.Dementia.rh, by="RID" );
	
	tmp = ANOVABatch ( BA.CN, BA.MCI, BA.Dementia, q.fdr );
	print( tmp$anovaRes );
	print(tmp$anovaRes[ tmp$sigRes, ] );
	
	anovaRes = tmp$anovaRes;
	sigVarName1 = "lh_perirhinal_thickness";
	sigVarLoc1 = which( rownames( anovaRes ) == sigVarName1 );
	sigVarName2 = "rh_perirhinal_thickness";
	sigVarLoc2 = which( rownames( anovaRes ) == sigVarName2 );
	if ( tiffFlag ) {
		tiff ( filename=paste("v44.4.FreeSurfer.PerRhSummary", "tiff", sep="."), compression="lzw" );
	} else {
		x11(); 
	} # if ( tiffFlag )
	par(mfcol=c(1,2), pty="s");
	SpecialGroupAndScatterComboPlot ( "L-PRhC", unlist(anovaRes[ sigVarLoc1, ]), unlist(BA.CN[,sigVarName1]),
										unlist(BA.MCI[,sigVarName1]), unlist(BA.Dementia[,sigVarName1]),
									"R-PRhC", unlist(anovaRes[ sigVarLoc2, ]), unlist(BA.CN[,sigVarName2]),
										unlist(BA.MCI[,sigVarName2]), unlist(BA.Dementia[,sigVarName2]),
										"Group Avg PRhC Thickness vs CDR\n", "Subject PRhC Thickness vs CDR\n", "Thickness (mm)" );
	if ( tiffFlag ) { dev.off(); } # if ( tiffFlag )
	
	sigVarName1 = "lh_BA3a_thickness";
	sigVarLoc1 = which( rownames( anovaRes ) == sigVarName1 );
	sigVarName2 = "rh_BA3a_thickness";
	sigVarLoc2 = which( rownames( anovaRes ) == sigVarName2 );
	if ( tiffFlag ) {
		tiff ( filename=paste("v44.4.FreeSurfer.A3aSummary", "tiff", sep="."), compression="lzw" );
	} else {
		x11(); 
	} # if ( tiffFlag )
	par(mfcol=c(1,2), pty="s");
	SpecialGroupAndScatterComboPlot ( "L-3a", unlist(anovaRes[ sigVarLoc1, ]), unlist(BA.CN[,sigVarName1]),
										unlist(BA.MCI[,sigVarName1]), unlist(BA.Dementia[,sigVarName1]),
									"R-3a", unlist(anovaRes[ sigVarLoc2, ]), unlist(BA.CN[,sigVarName2]),
										unlist(BA.MCI[,sigVarName2]), unlist(BA.Dementia[,sigVarName2]),
										"Group Avg PRhC Thickness vs CDR\n", "Subject PRhC Thickness vs CDR\n", "Thickness (mm)" );
	if ( tiffFlag ) { dev.off(); } # if ( tiffFlag )

}


