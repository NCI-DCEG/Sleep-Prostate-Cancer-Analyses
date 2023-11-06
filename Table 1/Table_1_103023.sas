/*******************************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: Table_1_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s): Table 1.; Results Paragraph 1

 PROGRAM FUNCTION: Generate summary statistics of characteristics presented in Table 1 of published manuscript and in first paragraph of results.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Data are pre-imputed data and based on complete case in Table 1.
				   No datasets are generated in this step.

*********************************************************************************************************************************************************/
******************************************************************************;
*  OPTIONS                                                                   *;
******************************************************************************;
options nofmterr;
******************************************************************************;
*  INPUT FILES                                                               *;
******************************************************************************;
libname i''; /*Read in library for dataset*/

****************************************************************************;
*  PROGRAM                                                                 *;
****************************************************************************;

/*READ IN ANALYTIC DATASET; n=34,260*/
data a77; set i.X/*Read in dataset*/; run;

/*TABLE 1 SLEEP DURATION-FREQUENCIES*/
proc freq data=a77;
table sleepduration_acc_65/missing;
run;

/*TABLE 1 MEANS*/
proc sort data=a77; by sleepduration_acc_65; run;
%macro means (v1);
proc means data=a77 mean std maxdec=8 nolabels;
by sleepduration_acc_65;
var &v1;
run;

%mend;
%means(v1=accel_entry_age); 
%means(v1=MVPA_HRS_0M_AVG);
%means(v1=SLEEPONSET_H_AVG);
%means(v1=SLEEPMIDPOINT_T_IAVG);
%means(v1=WAKEUP_H_AVG);
%means(v1=SE_AVG);
%means(v1=SOCIALJETLAG_H);
%means(v1=WASO_H_AVG);
%means(v1=WASOC30M);
%means(v1=accel_followup_days);

/*TABLE 1 FREQUENCIES*/
%macro t1fq (v1);
proc freq data=a77;
table &v1*sleepduration_acc_65; 
run;
%mend;
%t1fq(v1=midsleep_avg_h_cat);
%t1fq(v1=bmicat2);
%t1fq(v1=HEALTH);
%t1fq(v1=alc3);
%t1fq(v1=ever_smoke);
%t1fq(v1=educ);
%t1fq(v1=INCOME);
%t1fq(v1=nswork_v7);
%t1fq(v1=coffee);
%t1fq(v1=tea);
%t1fq(v1=rank_townsend);
%t1fq(v1=race_v1);
%t1fq(v1=psatest);
%t1fq(v1=prostate_fam);
%t1fq(v1=diabetes);

proc freq data=a77;
table 
p_disorders*sleepduration_acc_65
p_hyperplasia*sleepduration_acc_65
p_inflammation*sleepduration_acc_65/missing;
run;

/*Townsend Deprivation Index Levels by Quartiles presented in Table 1*/
proc sort data=a77; by rank_townsend; run;
proc means data=a77 n nmiss mean std min q1 median q3 max maxdec=2;
by rank_townsend;
var n_189_0_0; 
run;

/*DAYS CONTRIBUTED--RESULTS 1st Paragraph*/
proc means data=a77 n nmiss mean std;
var day_count_men;
run;
