df_aseg = ReadCSVFile( NULL, "/home/rstudio/ADNI3/Exp.2019.04.06/hippocampal_subfield_table.csv" );

df_study = df_study %>% left_join( df_aseg, by = c( "ID" = "Subject" ) );
WriteExcelReport( df_study, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_09April2019.xlsx" );
WriteCSVFile( df_study, NULL, "/home/rstudio/ADNI3/df_study_adnimerge_json_afni_final_fsdat_09April2019.csv" );

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
