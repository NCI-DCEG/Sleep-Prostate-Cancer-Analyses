/***********************************************************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: Table_2_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s): Table 2.

 PROGRAM FUNCTION: Generate HR and 95% confidence intervals for associations in Table 2 of published manuscript.
				   Generates total counts, person-years, and case numbers by counts as well.
				   Contains log-rank statistical tests for sleep duration and sleep midpoint.
				   Contains statistical tests for linearity for other exposure-outcome associations.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Unadjusted models, log rank tests, case counts, person-years, and case numbers based on pre-imputed data given the associated variables have no missing.
				   Adjusted model output based on imputed data given covariates have missingness.
*************************************************************************************************************************************************************************************/
******************************************************************************;
*  SYSTEM OPTIONS                                                            *;
******************************************************************************;
options nofmterr;

******************************************************************************;
*  INPUT FILES                                                               *;
******************************************************************************;
libname i ''; /*Read in input library*/
******************************************************************************;
*  OUTPUT FILES                                                              *;
******************************************************************************;
libname o '';/*Set output location*/
****************************************************************************;
*  PROGRAM                                                                 *;
****************************************************************************;

/*READ IN ANALYTIC DATASET; n=34,260*/
/*Recoding social jetlag post imputation so that analyses will be subset to those with missing data*/
/*Results are similar whether using complete case social jetlag or missing indicator social jetlag in analyses*/
data a5; 
set i.X; /*Set pre-imputed dataset*/
if socialjetlag_h_cat=3 then socialjetlag_h_cat=.; 
run;

/*TABLE 2 PERSON-TIME CALCULATION*/
%macro ptime (v1);
proc sort data=a5; by &v1;run;
proc means data=a5 n nmiss mean min max sum maxdec=2 nolabels;
by &v1;
var accel_followup_yrs;
run;
%mend;
%ptime(v1=sleepduration_acc_65); 
%ptime(v1=rank_onset);
%ptime(v1=midsleep_avg_h_cat);
%ptime(v1=rank_wakeup);
%ptime(v1=SOCIALJETLAG_H_CAT);


/*TABLE 2 EXPOSURE-OUTCOME FREQUENCIES*/ 
%macro fq (v1);
proc freq data=a5; 
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq(v1=sleepduration_acc_65); 
%fq(v1=rank_onset);
%fq(v1=midsleep_avg_h_cat);
%fq(v1=rank_wakeup);
%fq(v1=SOCIALJETLAG_H_CAT);

/*Unadjusted analyses-categorical exposures-no sleep variables, age, or outcome variables are missing and thus valid inference is made prior to imputation*/
/*Results are equivalent when associations are evaluated in a multiply imputed model, though SAS will return a warning that no data is missing*/
%macro uns(v0, v1);
proc phreg data=a5; 
class &v1(ref=&v0);
model accel_followup_days*pca(0)=&v1/rl;
run;
%mend;
%uns(v0='3', v1=sleepduration_acc_65);
%uns(v0='0', v1=rank_onset);
%uns(v0='1', v1=midsleep_avg_h_cat);
%uns(v0='0', v1=rank_wakeup);
%uns(v0='0', v1=SOCIALJETLAG_H_CAT);

/*Unadjusted analyses-linear exposures-no sleep variables, age, or outcome variables are missing and thus valid inference is made prior to imputation*/
/*Results are equivalent when associations are evaluated in a multiply imputed model, though SAS will return a warning that no data is missing*/
%macro unsl(v1);
proc phreg data=a5; 
model accel_followup_days*pca(0)=&v1/rl;
run;
%mend;
%unsl(v1=SLEEPONSET_H_AVG);
%unsl(v1=WAKEUP_H_AVG);
%unsl(v1=SOCIALJETLAG_H);

/*Read in imputed dataset for covariate-adjusted analyses*/
data dat_impute2; set i.X; run; 

/*Rename variables due to issues with variable length or characters in subsequent analyses*/
data a6mi; 
rename sleepduration_acc_65=DURATION nswork_v7=nswork; 
set dat_impute2;
run;

/*Sort data prior to conducting imputed analyses*/
proc sort data=a6mi; by n_eid _imputation_; run;
proc sort data=a6mi; by _imputation_; run;

/*Duration*/
proc phreg data=a6mi; 
  by _imputation_;
  class DURATION(ref='3') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=DURATION 
nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_DURATION CovB=o.cov_E_DURATION; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_DURATION CovB=o.cov_E_DURATION;
class DURATION nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects DURATION
nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_DURATION; 
run;

data o.pa_DURATION; set o.va_DURATION; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Rank_Onset*/
proc phreg data=a6mi; 
  by _imputation_;
  class rank_onset(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=rank_onset nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_onset CovB=o.cov_E_rankonset; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_onset CovB=o.cov_E_rankonset;
class rank_onset nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects rank_onset nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_rank_onset; 
run;

data o.pa_rank_onset; set o.va_rank_onset; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Sleep Midpoint*/
proc phreg data=a6mi; 
  by _imputation_;
  class midsleep_avg_h_cat(ref='1') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=midsleep_avg_h_cat nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_midsleep_avg_h_cat CovB=o.cov_E_midpoint; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_midsleep_avg_h_cat CovB=o.cov_E_midpoint;
class midsleep_avg_h_cat nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects midsleep_avg_h_cat nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_midsleep_avg_h_cat; 
run;

data o.pa_midsleep_avg_h_cat; set o.va_midsleep_avg_h_cat; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;


/*Rank_Wakeup*/
proc phreg data=a6mi; 
  by _imputation_;
  class rank_wakeup(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=rank_wakeup nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_wakeup CovB=o.cov_E_rankwakeup; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_wakeup CovB=o.cov_E_rankwakeup;
class rank_wakeup nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects rank_wakeup nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_rank_wakeup; 
run;

data o.pa_rank_wakeup; set o.va_rank_wakeup; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Social Jetlag*/
proc phreg data=a6mi; 
where socialjetlag_h>.;
  by _imputation_;
  class SOCIALJETLAG_H_CAT(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=SOCIALJETLAG_H_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_SOCIALJETLAG_H_CAT CovB=cov_E_sjl; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_SOCIALJETLAG_H_CAT CovB=cov_E_sjl;
class SOCIALJETLAG_H_CAT nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects SOCIALJETLAG_H_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_SOCIALJETLAG_H_CAT; 
run;

data o.pa_SOCIALJETLAG_H_CAT; set o.va_SOCIALJETLAG_H_CAT; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Sleep Onset (continuous)*/
proc phreg data=a6mi; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=SLEEPONSET_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_SLEEPONSET_H_AVG CovB=o.cov_SLEEPONSET_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_SLEEPONSET_H_AVG CovB=o.cov_SLEEPONSET_H_AVG;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects SLEEPONSET_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_SLEEPONSET_H_AVG; 
run;

data o.pa_SLEEPONSET_H_AVG; set o.va_SLEEPONSET_H_AVG; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Wakeup Time Continuous*/
proc phreg data=a6mi; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WAKEUP_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WAKEUP_H_AVG CovB=o.cov_WAKEUP_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WAKEUP_H_AVG covB=o.cov_WAKEUP_H_AVG;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WAKEUP_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WAKEUP_H_AVG; 
run;

data o.pa_WAKEUP_H_AVG; set o.va_WAKEUP_H_AVG; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Social jetlag (continuous)*/
proc phreg data=a6mi; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=SOCIALJETLAG_H nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_SOCIALJETLAG_H CovB=o.cov_SOCIALJETLAG_H; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_SOCIALJETLAG_H CovB=o.cov_SOCIALJETLAG_H;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects SOCIALJETLAG_H nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_SOCIALJETLAG_H; 
run;

data o.pa_SOCIALJETLAG_H; set o.va_SOCIALJETLAG_H; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*CREATE INDICATOR VARIABLES FOR TEST FOR LINEARITY ANALYSES*/
/*PROC MIANALYZE cannot have a class statement and test statement at the same time, so all class variables must be recoded as dummy variables*/
data a7; set a6mi;
NSWORK_0=0;
if nswork=0 then NSWORK_0=1;
NSWORK_1=0;
if nswork=1 then NSWORK_1=1;
NSWORK_2=0;
if nswork=2 then NSWORK_2=1;
NSWORK_3=0;
if nswork=3 then NSWORK_3=1;
NSWORK_4=0;
if nswork=4 then NSWORK_4=1;

BMICAT2_1=0;
if bmicat2=1 then BMICAT2_1=1;
BMICAT2_2=0;
if bmicat2=2 then BMICAT2_2=1;
BMICAT2_3=0;
if bmicat2=3 then BMICAT2_3=1;
BMICAT2_4=0;
if bmicat2=4 then BMICAT2_4=1;
BMICAT2_5=0;
if bmicat2=5 then BMICAT2_5=1;

EVER_SMOKE_0=0;
if ever_smoke=0 then EVER_SMOKE_0=1;
EVER_SMOKE_1=0;
if ever_smoke=1 then EVER_SMOKE_1=1;
EVER_SMOKE_2=0;
if ever_smoke=2 then EVER_SMOKE_2=1;

HEALTH_1=0;
if HEALTH=1 then HEALTH_1=1;
HEALTH_2=0;
if HEALTH=2 then HEALTH_2=1;
HEALTH_3=0;
if HEALTH=3 then HEALTH_3=1;
HEALTH_4=0;
if HEALTH=4 then HEALTH_4=1;

ALC3_0=0;
if alc3=0 then ALC3_0=1;
ALC3_1=0;
if alc3=1 then ALC3_1=1;
ALC3_2=0;
if alc3=2 then ALC3_2=1;
ALC3_3=0;
if alc3=3 then ALC3_3=1;
ALC3_4=0;
if alc3=4 then ALC3_4=1;

educ_1=0;
if educ=1 then educ_1=1;
educ_2=0;
if educ=2 then educ_2=1;
educ_3=0;
if educ=3 then educ_3=1;

INCOME_1=0;
if income=1 then INCOME_1=1;
INCOME_2=0;
if income=2 then INCOME_2=1;
INCOME_3=0;
if income=3 then INCOME_3=1;
INCOME_4=0;
if income=4 then INCOME_4=1;
INCOME_5=0;
if income=5 then INCOME_5=1;

RANK_ACC_PA_0=0;
if rank_acc_pa=0 then RANK_ACC_PA_0=1;
RANK_ACC_PA_1=0;
if rank_acc_pa=1 then RANK_ACC_PA_1=1;
RANK_ACC_PA_2=0;
if rank_acc_pa=2 then RANK_ACC_PA_2=1;
RANK_ACC_PA_3=0;
if rank_acc_pa=3 then RANK_ACC_PA_3=1;

RACE_BI_0=0;
if race_bi=0 then RACE_BI_0=1;
RACE_BI_1=0;
if race_bi=1 then RACE_BI_1=1;

RANK_TOWNSEND_0=0;
if rank_townsend=0 then RANK_TOWNSEND_0=1;
RANK_TOWNSEND_1=0;
if rank_townsend=1 then RANK_TOWNSEND_1=1;
RANK_TOWNSEND_2=0;
if rank_townsend=2 then RANK_TOWNSEND_2=1;
RANK_TOWNSEND_3=0;
if rank_townsend=3 then RANK_TOWNSEND_3=1;

PSATEST_0=0;
if psatest=0 then PSATEST_0=1;
PSATEST_1=0;
if psatest=1 then PSATEST_1=1;

PROSTATE_FAM_0=0;
if prostate_fam=0 then PROSTATE_FAM_0=1;
PROSTATE_FAM_1=0;
if prostate_fam=1 then PROSTATE_FAM_1=1;

COFFEE_0=0;
if coffee=0 then COFFEE_0=1;
COFFEE_1=0;
if coffee=1 then COFFEE_1=1;
COFFEE_2=0;
if coffee=2 then COFFEE_2=1;
COFFEE_3=0;
if coffee=3 then COFFEE_3=1;
COFFEE_4=0;
if coffee=4 then COFFEE_4=1;

TEA_0=0;
if tea=0 then TEA_0=1;
TEA_1=0;
if tea=1 then TEA_1=1;
TEA_2=0;
if tea=2 then TEA_2=1;
TEA_3=0;
if tea=3 then TEA_3=1;
TEA_4=0;
if tea=4 then TEA_4=1;

DIABETES_0=0;
if diabetes=0 then DIABETES_0=1;
DIABETES_1=0;
if diabetes=1 then DIABETES_1=1;
run;

/*TEST FOR LINEARITY*/

/*Sleep Onset (continuous)*/
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=SLEEPONSET_H_AVG 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_SLEEPONSET_H_AVG CovB=o.cov_LT_SLEEPONSET_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_SLEEPONSET_H_AVG CovB=o.cov_LT_SLEEPONSET_H_AVG;
modeleffects SLEEPONSET_H_AVG
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test SLEEPONSET_H_AVG/mult;
ods output TestMultStat=o.test_1LAD_SLEEPONSET_H_AVG TestParameterEstimates=o.test_1LAD_p_SLEEPONSET_H_AVG;
run;

/*Wakeup Time (continuous)*/
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WAKEUP_H_AVG 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_WAKEUP_H_AVG CovB=o.cov_LT_WAKEUP_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_WAKEUP_H_AVG CovB=o.cov_LT_WAKEUP_H_AVG;
modeleffects WAKEUP_H_AVG
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WAKEUP_H_AVG/mult;
ods output TestMultStat=o.test_1LAD_WAKEUP_H_AVG TestParameterEstimates=o.test_1LAD_p_WAKEUP_H_AVG;
run;

/*Social jetlag (continuous)*/
proc phreg data=a7; 
where SOCIALJETLAG_H>.;
by _imputation_;
model accel_followup_days*pca(0)=SOCIALJETLAG_H 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_SOCIALJETLAG_H CovB=o.cov_LT_SOCIALJETLAG_H; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_SOCIALJETLAG_H CovB=o.cov_LT_SOCIALJETLAG_H;
modeleffects SOCIALJETLAG_H
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test SOCIALJETLAG_H/mult;
ods output TestMultStat=o.test_1LAD_SOCIALJETLAG_H TestParameterEstimates=o.test_1LAD_p_SOCIALJETLAG_H;
run;


/*PROC LIFETEST FOR TABLE 2*/
proc sort data=a5; by n_eid;run;
%macro lrt (v1);
proc lifetest data=a5;
strata &v1;
time accel_followup_days*pca(0);
test &v1;
run;
%mend;
%lrt(v1=sleepduration_acc_65);
%lrt(v1=midsleep_avg_h_cat);
