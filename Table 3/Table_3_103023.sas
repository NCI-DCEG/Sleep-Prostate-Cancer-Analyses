/***********************************************************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: Table_3_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s): Table 3.

 PROGRAM FUNCTION: Generate HR and 95% confidence intervals for associations in Table 3 of published manuscript.
				   Generates total counts, person-years, and case numbers by counts as well.
				   Contains statistical tests for linearity for sleep quality exposure-outcome associations.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Unadjusted models, log rank tests, case counts, person-years, and case numbers based on pre-imputed data given the associated variables have no missing.
				   Adjusted model output based on imputed data given covariates have missingness.
				   Tests for linearity completed using adjusted models.
				   No datasets are generated in this script.
*************************************************************************************************************************************************************************************/
******************************************************************************;
*     SYSTEM OPTIONS                                                         *;
******************************************************************************;
options nofmterr;
******************************************************************************;
* INPUT FILES                                                                *;
******************************************************************************;
libname i '';/*Read in input library*/
******************************************************************************;
* OUTPUT FILES                                                               *;
******************************************************************************;
libname o '';/*Set output library*/
****************************************************************************;
*  PROGRAM                                                                  *;
****************************************************************************;

/*READ IN ANALYTIC DATASET; n=34,260*/
data a5; 
set i.X;
run;

/*TABLE 2-3 PERSON-TIME CALCULATION*/
%macro ptime (v1);
proc sort data=a5; by &v1;run;
proc means data=a5 n nmiss mean min max sum maxdec=2 nolabels;
by &v1;
var accel_followup_yrs;
run;
%mend;
%ptime(v1=rank_se);
%ptime(v1=waso_m_cat3);
%ptime(v1=waso_m_cat);
%ptime(v1=wasocountv7);


/*TABLE 2-3 EXPOSURE-OUTCOME FREQUENCIES*/ 
%macro fq (v1);
proc freq data=a5; 
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq(v1=rank_se);
%fq(v1=waso_m_cat3);
%fq(v1=waso_m_cat);
%fq(v1=wasocountv7);

/*Unadjusted analyses-categorical exposures-no sleep variables, age, or outcome variables are missing and thus valid inference is made prior to imputation*/
/*Results are equivalent when associations are evaluated in a multiply imputed model, though SAS will return a warning that no data is missing*/
%macro uns(v0, v1);
proc phreg data=a5; 
class &v1(ref=&v0);
model accel_followup_days*pca(0)=&v1/rl;
run;
%mend;
%uns(v0='0', v1=rank_se);
%uns(v0='0', v1=WASO_M_CAT3);
%uns(v0='0', v1=WASO_M_CAT);
%uns(v0='1', v1=wasocountv7);

/*Unadjusted analyses-linear exposures-no sleep variables, age, or outcome variables are missing and thus valid inference is made prior to imputation*/
/*Results are equivalent when associations are evaluated in a multiply imputed model, though SAS will return a warning that no data is missing*/
%macro unsl(v1);
proc phreg data=a5; 
model accel_followup_days*pca(0)=&v1/rl;
run;
%mend;
%unsl(v1=WASO_H_AVG);
%unsl(v1=WASOC30M);

/*Read in imputed dataset for covariate-adjusted analyses*/
data dat_impute; 
set i.X;
run; 

/*Rename variables due to issues with variable length or characters in subsequent analyses*/
data a6mi; 
rename nswork_v7=nswork; 
set dat_impute;
run;

/*Sort data prior to conducting imputed analyses*/
proc sort data=a6mi; by n_eid _imputation_; run;
proc sort data=a6mi; by _imputation_; run;

/*Sleep Efficiency (categorical)*/
proc phreg data=a6mi; 
  by _imputation_;
  class rank_se(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=rank_se nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=q.PEsta_rank_se CovB=cov_E_rankse; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_rank_se CovB=cov_E_rankse;
class rank_se nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects rank_se nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_rank_se; 
run;

data q.pa_rank_se; set q.va_rank_se; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO (categorical; 4 level)*/
proc phreg data=a6mi; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=q.PEsta_WASO_M_CAT3 CovB=q.cov_E_WASO_M_CAT3; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_WASO_M_CAT3 CovB=q.cov_E_WASO_M_CAT3;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_WASO_M_CAT3; 
run;

data q.pa_WASO_M_CAT3; set q.va_WASO_M_CAT3; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO (binary)*/
proc phreg data=a6mi; 
  by _imputation_;
  class WASO_M_CAT(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/ covb;
  ods output ParameterEstimates=q.PEsta_WASO_M_CAT CovB=q.cov_E_WASO_M_CAT; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_WASO_M_CAT CovB=q.cov_E_WASO_M_CAT;
class WASO_M_CAT nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_WASO_M_CAT; 
run;

data q.pa_WASO_M_CAT; set q.va_WASO_M_CAT; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Frequency of WASO (categorical)*/
proc phreg data=a6mi; 
  by _imputation_;
  class wasocountv7(ref='1') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=wasocountv7 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=q.PEsta_wasocountv7 CovB=q.cov_E_wasocountv7; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_wasocountv7 CovB=q.cov_E_wasocountv7;
class wasocountv7 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects wasocountv7 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_wasocountv7; 
run;

data q.pa_wasocountv7; set q.va_wasocountv7; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO (continuous)*/
proc phreg data=a6mi; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=q.PEsta_WASO_H_AVG CovB=q.cov_WASO_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_WASO_H_AVG CovB=q.cov_WASO_H_AVG;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_WASO_H_AVG; 
run;

data q.pa_WASO_H_AVG; set q.va_WASO_H_AVG; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Frequency of WASO (continuous)*/
proc phreg data=a6mi; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASOC30M nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=q.PEsta_WASOC30M CovB=q.cov_WASOC30M; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_WASOC30M CovB=q.cov_WASOC30M;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASOC30M nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=q.va_WASOC30M;
run;

data q.pa_WASOC30M; set q.va_WASOC30M; 
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

/*WASO*/
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASO_H_AVG 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=q.PEsta_LT_WASO_H_AVG CovB=q.cov_LT_WASO_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_LT_WASO_H_AVG CovB=q.cov_LT_WASO_H_AVG;
modeleffects WASO_H_AVG
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASO_H_AVG/mult;
ods output TestMultStat=tt.test_1LAD_WASO_H_AVG TestParameterEstimates=tt.test_1LAD_p_WASO_H_AVG;
run;

/*Frequency of WASO*/
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASOC30M 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=q.PEsta_LT_WASOC30M CovB=q.cov_LT_WASOC30M; 
run;

proc mianalyze parms(classvar=ClassVal)=q.PEsta_LT_WASOC30M CovB=q.cov_LT_WASOC30M;
modeleffects WASOC30M 
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASOC30M/mult;
ods output TestMultStat=tt.test_1LAD_WASOC30M TestParameterEstimates=tt.test_1LAD_p_WASOC30M;
run;
