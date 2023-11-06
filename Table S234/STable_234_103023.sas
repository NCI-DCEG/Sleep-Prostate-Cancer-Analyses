/***************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: S_Table_234_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s): Supplementary Tables 2, 3, & 4.

 PROGRAM FUNCTION: Generate summary statistics of characteristics presented in Supplementary Tables 2, 3, & 4 of published manuscript.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Data have no missingness and are thus based on pre-imputed data.
				   No datasets are generated in this step.
***************************************************************************************************************************************/

******************************************************************************;
*  SYSTEM OPTIONS                                                            *;
******************************************************************************;
options nofmterr;
******************************************************************************;
*  LIBNAMES                                                                  *;
******************************************************************************;
libname i '';/*ADD IN LIBRARY AS APPROPRIATE FOR DATASOURCE*/

****************************************************************************;
*  PROGRAM                                                                 *;
****************************************************************************;
/*READ IN ANALYTIC DATASET; n=34,260*/
data a77; set i.X; run; /*Read in pre-imputed dataset*/

/*GENERATE OUTPUT FOR SUPPLEMENTARY TABLES 2, 3, & 4*/

/*SUPPLEMENTARY TABLE 2 RANGES FOR SLEEP CHARACTERISTICS*/

/*Sleep Duration by categorical levels*/
proc sort data=a77; by sleepduration_acc_65; run;
proc means data=a77 n nmiss min max maxdec=2;
by sleepduration_acc_65;
var SLEEPDURATION_H_AVG; 
run;

/*Sleep Onset by Quartiles*/
proc sort data=a77; by rank_onset; run;
proc means data=a77 n nmiss min max maxdec=2;
by rank_onset;
var SLEEPONSET_H_AVG; 
run;

/*Sleep Midpoint by categorical levels*/
proc sort data=a77; by midsleep_avg_h_cat; run;
proc means data=a77 n nmiss min max maxdec=2;
by midsleep_avg_h_cat;
var SLEEPMIDPOINT_T_IAVG; 
run;

/*Wakeup Time by Quartiles*/
proc sort data=a77; by rank_wakeup; run;
proc means data=a77 n nmiss min max maxdec=2;
by rank_wakeup;
var WAKEUP_H_AVG; 
run;

/*Social Jetlag by Levels*/
proc sort data=a77; by SOCIALJETLAG_H_CAT; run;
proc means data=a77 n nmiss min max maxdec=2;
by SOCIALJETLAG_H_CAT;
var SOCIALJETLAG_H; 
run;

/*Sleep Efficiency by Quartiles*/
proc sort data=a77; by rank_se; run;
proc means data=a77 n nmiss min max maxdec=2;
by rank_se;
var SE_AVG; 
run;

/*WASO by Levels*/
proc sort data=a77; by WASO_M_CAT3; run;
proc means data=a77 n nmiss min max maxdec=2;
by WASO_M_CAT3;
var WASO_H_AVG; 
run;

/*WASO by Quartiles*/
proc sort data=a77; by rank_WASO; run;
proc means data=a77 n nmiss min max maxdec=2;
by rank_WASO;
var WASO_H_AVG;
run;

/*Wake_Active by Quartiles*/
proc sort data=a77; by rank_active; run;
proc means data=a77 n nmiss min max maxdec=2;
by rank_active;
var WAKE_ACTIVE_H_AVG;
run;

/*Frequency of WASO by Levels*/
proc sort data=a77; by wasocountv7; run;
proc means data=a77 n nmiss min max maxdec=7;
by wasocountv7;
var WASOC30M; 
run;


/*SUPPLEMENTARY TABLE 3 STATISTICS*/
proc sort data=a77; by n_eid; run;
%macro means (v1);
proc means data=a77 mean std q1 median q3 maxdec=2;
var &v1;
run;
%mend;
%means(v1=SLEEPDURATION_H_AVG); /*Sleep Duration (continuous)*/
%means(v1=SLEEPONSET_H_AVG); /*Sleep Onset (continuous)*/
%means(v1=SLEEPMIDPOINT_T_IAVG); /*Sleep Midpoint (continuous)*/
%means(v1=WAKEUP_H_AVG); /*Wakeup Time (continuous)*/
%means(v1=SOCIALJETLAG_H); /*Social Jetlag (continuous)*/
%means(v1=SE_AVG); /*Sleep Efficiency (continuous)*/
%means(v1=WASO_H_AVG); /*WASO (continuous)*/
%means(v1=WASOC30M); /*Frequency of WASO (continuous)*/

/*SUPPLEMENTARY TABLE 4 FREQUENCIES*/
%macro t1fq (v1);
proc freq data=a77;
table &v1;
run;
%mend;
%t1fq(v1=sleepduration_acc_65); 
%t1fq(v1=rank_onset);
%t1fq(v1=midsleep_avg_h_cat);
%t1fq(v1=rank_wakeup);
%t1fq(v1=SOCIALJETLAG_H_CAT);
%t1fq(v1=rank_se);
%t1fq(v1=waso_m_cat3);
%t1fq(v1=rank_WASO);
%t1fq(v1=rank_active);
%t1fq(v1=wasocountv7);
