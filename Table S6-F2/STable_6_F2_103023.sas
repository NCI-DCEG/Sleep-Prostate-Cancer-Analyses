/********************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: STable_6_F2_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s) and Figure(s): Supplementary Table 6 & Figure 2.

 PROGRAM FUNCTION: Generate summary statistics of characteristics presented in Supplementary Table 6 and Figure 2 of published manuscript.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Data are pre-imputed and imputed in these analyses.
********************************************************************************************************************************************/

******************************************************************************;
*     SYSTEM OPTIONS                                                         *;
******************************************************************************;
options nofmterr;
******************************************************************************;
*   INPUT FILES                                                              *;
******************************************************************************;
libname i '';/*Read in library for datasets*/
******************************************************************************;
*   OUTPUT FILES                                                             *;
******************************************************************************;
libname o '';/*Set library for output*/
****************************************************************************;
*  PROGRAM                                                                  *;
****************************************************************************;
/*Read in Non-imputed dataset to get frequencies for Supplementary Table 6*/
data a5; set i.X; run;

/*SUPPLEMENTAL TABLE 6 EXPOSURE-OUTCOME FREQUENCIES*/ 
%macro fq0 (v1);
proc freq data=a5; 
where psatest=0;
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq0(v1=waso_m_cat3);
%fq0(v1=waso_m_cat);

%macro fq1 (v1);
proc freq data=a5; 
where psatest=1;
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq1(v1=waso_m_cat3);
%fq1(v1=waso_m_cat);

/*READ IN IMPUTED ANALYTIC DATASET; n=34,260*/
data dat_impute; set  i.X; run;

/*RENAME EMPLOYMENT/NIGHT SHIFT WORK VARIABLE SO IT WILL READ INTO MODEL*/
data a6mi; rename nswork_v7=nswork; set dat_impute;
run;

proc sort data=a6mi; by n_eid _imputation_; run;
proc sort data=a6mi; by _imputation_; run;

/*CREATING DUMMY VARIABLES FOR INTERACTION ANALYSES*/
data a7a; set a6mi;

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

WASO_M_CAT3_0=0;
if WASO_M_CAT3=0 then WASO_M_CAT3_0=1;
WASO_M_CAT3_1=0;
if WASO_M_CAT3=1 then WASO_M_CAT3_1=1;
WASO_M_CAT3_2=0;
if WASO_M_CAT3=2 then WASO_M_CAT3_2=1;
WASO_M_CAT3_3=0;
if WASO_M_CAT3=3 then WASO_M_CAT3_3=1;
run;


/*CREATING INTERACTION TERMS FOR MODELS*/
data a7;
set a7a;
WASOAge=WASO_H_AVG*accel_entry_age;
WASObmi=WASO_H_AVG*bmi;
WASOpsa=WASO_H_AVG*PSATEST_1;
WASOwake=WASO_H_AVG*WAKE_ACTIVE_H_AVG;
run;

/*SUPPLEMENTARY TABLE 6 RESULTS*/

/*CALCULATING P-VALUE FOR INTERACTION TERMS*/
/*WASOHAVG-PSA (Supplemental Table 6)*/
proc sort data=a7; by n_eid _imputation_; run;
proc sort data=a7; by _imputation_; run;
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASO_H_AVG WASOpsa
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_WASOPSA CovB=o.cov_LT_WASOPSA; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_WASOPSA CovB=o.cov_LT_WASOPSA;
modeleffects WASO_H_AVG WASOpsa
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4 accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASOpsa/mult;
ods output TestMultStat=o.test_1LAD_WASOPSA TestParameterEstimates=o.test_1LAD_p_WASOPSA;
run;

/*SUPPLEMENTARY TABLE 6 RESULTS*/
/*Effect Modification-PSATest-No PSA Test History*/
data a6_psa1; set a6mi;
where psatest=0;
run;

/*WASO_M_CAT3*/
proc sort data=a6_psa1; by n_eid _imputation_; run;
proc sort data=a6_psa1; by _imputation_; run;
proc phreg data=a6_psa1; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3psa1 CovB=o.cov_E_WASO_M_CAT3psa1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3psa1 CovB=o.cov_E_WASO_M_CAT3psa1;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3psa1; 
run;

data o.pa_WASO_M_CAT3psa1; set o.va_WASO_M_CAT3psa1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc sort data=a6_psa1; by n_eid _imputation_; run;
proc sort data=a6_psa1; by _imputation_; run;
proc phreg data=a6_psa1; 
  by _imputation_;
  class WASO_M_CAT(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CATpsa1 CovB=o.cov_E_WASO_M_CATpsa1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CATpsa1 CovB=o.cov_E_WASO_M_CATpsa1;
class WASO_M_CAT nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CATpsa1; 
run;

data o.pa_WASO_M_CATpsa1; set o.va_WASO_M_CATpsa1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc sort data=a6_psa1; by n_eid _imputation_; run;
proc sort data=a6_psa1; by _imputation_; run;
proc phreg data=a6_psa1; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_AVGpsa1 CovB=o.cov_E_WASO_AVGpsa1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_AVGpsa1 CovB=o.cov_E_WASO_AVGpsa1;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_AVGpsa1; 
run;

data o.pa_WASO_AVGpsa1; set o.va_WASO_AVGpsa1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*SUPPLEMENTARY TABLE 6 RESULTS*/
/*Effect Modification-PSATest-PSA Test History*/
data a6_psa2; set a6mi;
where psatest=1;
run;

/*WASO_M_CAT3*/
proc sort data=a6_psa2; by n_eid _imputation_; run;
proc sort data=a6_psa2; by _imputation_; run;
proc phreg data=a6_psa2; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3psa2 CovB=o.cov_E_WASO_M_CAT3psa2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3psa2 CovB=o.cov_E_WASO_M_CAT3psa2;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3psa2; 
run;

data o.pa_WASO_M_CAT3psa2; set o.va_WASO_M_CAT3psa2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc sort data=a6_psa2; by n_eid _imputation_; run;
proc sort data=a6_psa2; by _imputation_; run;
proc phreg data=a6_psa2; 
  by _imputation_;
  class WASO_M_CAT(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CATpsa2 CovB=o.cov_E_WASO_M_CATpsa2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CATpsa2 CovB=o.cov_E_WASO_M_CATpsa2;
class WASO_M_CAT nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CATpsa2; 
run;

data o.pa_WASO_M_CATpsa2; set o.va_WASO_M_CATpsa2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc sort data=a6_psa2; by n_eid _imputation_; run;
proc sort data=a6_psa2; by _imputation_; run;
proc phreg data=a6_psa2; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_AVGpsa2 CovB=o.cov_E_WASO_AVGpsa2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_AVGpsa2 CovB=o.cov_E_WASO_AVGpsa2;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
modeleffects WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_AVGpsa2; 
run;

data o.pa_WASO_AVGpsa2; set o.va_WASO_AVGpsa2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/**************************************************************************************************************************************************************************************/
/*FIGURE 2 Results*/

/*CALCULATING P-VALUE FOR INTERACTION TERMS*/
/*WASOHAVG-Age*/
proc sort data=a7; by n_eid _imputation_; run;
proc sort data=a7; by _imputation_; run;

proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASO_H_AVG WASOAge
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_WASOAge CovB=o.cov_LT_WASOAge; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_WASOAge CovB=o.cov_LT_WASOAge;
modeleffects WASO_H_AVG WASOAge
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASOAge/mult;
ods output TestMultStat=o.test_1LAD_WASOAge TestParameterEstimates=o.test_1LAD_p_WASOAge;
run;

/*CALCULATING P-VALUE FOR INTERACTION TERMS*/
/*WASOHAVG-BMI*/

proc sort data=a7; by n_eid _imputation_; run;
proc sort data=a7; by _imputation_; run;
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASO_H_AVG WASObmi
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
bmi EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_WASOBMI CovB=o.cov_LT_WASOBMI; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_WASOBMI CovB=o.cov_LT_WASOBMI;
modeleffects WASO_H_AVG WASObmi
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4 accel_entry_age 
bmi EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASObmi/mult;
ods output TestMultStat=o.test_1LAD_WASOBMI TestParameterEstimates=o.test_1LAD_p_WASOBMI;
run;

/*CALCULATING P-VALUE FOR INTERACTION TERMS*/
/*WASOHAVG-WAKE*/

proc sort data=a7; by n_eid _imputation_; run;
proc sort data=a7; by _imputation_; run;
proc phreg data=a7; 
by _imputation_;
model accel_followup_days*pca(0)=WASO_H_AVG WASOwake WAKE_ACTIVE_H_AVG
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
/rl covb;
  ods output ParameterEstimates=o.PEsta_LT_WASOWAKE CovB=o.cov_LT_WASOWAKE; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_LT_WASOWAKE CovB=o.cov_LT_WASOWAKE;
modeleffects WASO_H_AVG WASOwake WAKE_ACTIVE_H_AVG
NSWORK_1 NSWORK_2 NSWORK_3 NSWORK_4  accel_entry_age 
BMICAT2_2 BMICAT2_3 BMICAT2_4 BMICAT2_5 EVER_SMOKE_1 EVER_SMOKE_2
HEALTH_2 HEALTH_3 HEALTH_4 ALC3_1 ALC3_2 ALC3_3 ALC3_4
educ_2 educ_3 INCOME_2 INCOME_3 INCOME_4 INCOME_5
RANK_ACC_PA_1 RANK_ACC_PA_2 RANK_ACC_PA_3 RACE_BI_1
RANK_TOWNSEND_1 RANK_TOWNSEND_2 RANK_TOWNSEND_3
PSATEST_1 PROSTATE_FAM_1 COFFEE_1 COFFEE_2 COFFEE_3 COFFEE_4
TEA_1 TEA_2 TEA_3 TEA_4 DIABETES_1
;
test WASOwake/mult;
ods output TestMultStat=o.test_1LAD_WASOWAKE TestParameterEstimates=o.test_1LAD_p_WASOWAKE;
run;


/*Conducting Effect modification analyses*/
/*Effect Modification-BMI <30kg/m^2*/
data a6_bmi1; set a6mi;
where bmicat2<3;
run;

proc sort data=a6_bmi1; by n_eid _imputation_; run;
proc sort data=a6_bmi1; by _imputation_; run;

/*WASO_M_CAT3*/
proc phreg data=a6_bmi1; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3bmi1 CovB=o.cov_E_WASO_M_CAT3bmi1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3bmi1 CovB=o.cov_E_WASO_M_CAT3bmi1;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3bmi1; 
run;

data o.pa_WASO_M_CAT3bmi1; set o.va_WASO_M_CAT3bmi1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Effect Modification-BMI >=30kg/m^2*/
data a6_bmi2; set a6mi;
where bmicat2>=3;
run;

/*WASO_M_CAT3*/

proc sort data=a6_bmi2; by n_eid _imputation_; run;
proc sort data=a6_bmi2; by _imputation_; run;
proc phreg data=a6_bmi2; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3bmi2 CovB=o.cov_E_WASO_M_CAT3bmi2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3bmi2 CovB=o.cov_E_WASO_M_CAT3bmi2;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3bmi2; 
run;

data o.pa_WASO_M_CAT3bmi2; set o.va_WASO_M_CAT3bmi2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Effect Modification-Age <65 years*/
data a6_age1; set a6mi;
where accel_entry_age<65;
run;

/*WASO_M_CAT3*/
proc sort data=a6_age1; by n_eid _imputation_; run;
proc sort data=a6_age1; by _imputation_; run;

proc phreg data=a6_age1; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3age1 CovB=o.cov_E_WASO_M_CAT3age1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3age1 CovB=o.cov_E_WASO_M_CAT3age1;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3age1; 
run;

data o.pa_WASO_M_CAT3age1; set o.va_WASO_M_CAT3age1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Effect Modification-Age >=65 years*/
data a6_age2; set a6mi;
where accel_entry_age>=65;
run;

/*WASO_M_CAT3*/
proc sort data=a6_age2; by n_eid _imputation_; run;
proc sort data=a6_age2; by _imputation_; run;
proc phreg data=a6_age2; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3age2 CovB=o.cov_E_WASO_M_CAT3age2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3age2 CovB=o.cov_E_WASO_M_CAT3age2;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3age2; 
run;

data o.pa_WASO_M_CAT3age2; set o.va_WASO_M_CAT3age2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Effect Modification-ActiveWASO-Below Median*/
data a6_awaso1; set a6mi;
where rank_active<=1;
run;

/*WASO_M_CAT3*/
proc sort data=a6_awaso1; by n_eid _imputation_; run;
proc sort data=a6_awaso1; by _imputation_; run;
proc phreg data=a6_awaso1; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3awaso1 CovB=o.cov_E_WASO_M_CAT3awaso1; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3awaso1 CovB=o.cov_E_WASO_M_CAT3awaso1;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3awaso1; 
run;

data o.pa_WASO_M_CAT3awaso1; set o.va_WASO_M_CAT3awaso1; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Effect Modification-ActiveWASO-Above Median*/
data a6_awaso2; set a6mi;
where rank_active>1;
run;

/*WASO_M_CAT3*/
proc sort data=a6_awaso2; by n_eid _imputation_; run;
proc sort data=a6_awaso2; by _imputation_; run;
proc phreg data=a6_awaso2; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3awaso2 CovB=o.cov_E_WASO_M_CAT3awaso2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3awaso2 CovB=o.cov_E_WASO_M_CAT3awaso2;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3awaso2; 
run;

data o.pa_WASO_M_CAT3awaso2; set o.va_WASO_M_CAT3awaso2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Adjust-ActiveWASO*/
/*WASO_M_CAT3*/
proc sort data=a6mi; by n_eid _imputation_; run;
proc sort data=a6mi; by _imputation_; run;
proc phreg data=a6mi; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 wake_active_h_avg nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3awasoadj CovB=o.cov_E_WASO_M_CAT3awasoadj; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3awasoadj CovB=o.cov_E_WASO_M_CAT3awasoadj;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 wake_active_h_avg nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3awasoadj; 
run;

data o.pa_WASO_M_CAT3awasoadj; set o.va_WASO_M_CAT3awasoadj; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude Shift Workers*/
data a6mi_NS; set a6mi;
if nswork in (1, 2) then delete;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_NS; by n_eid _imputation_; run;
proc sort data=a6mi_NS; by _imputation_; run;
proc phreg data=a6mi_NS; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3NS CovB=o.cov_E_WASO_M_CAT3NS; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3NS CovB=o.cov_E_WASO_M_CAT3NS;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3NS; 
run;

data o.pa_WASO_M_CAT3NS; set o.va_WASO_M_CAT3NS; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude Diabetes*/
data a6mi_DM; set a6mi;
where diabetes=0;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_DM; by n_eid _imputation_; run;
proc sort data=a6mi_DM; by _imputation_; run;
proc phreg data=a6mi_DM; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3DM CovB=o.cov_E_WASO_M_CAT3DM; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3DM CovB=o.cov_E_WASO_M_CAT3DM;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea;
ods output parameterestimates=o.va_WASO_M_CAT3DM; 
run;

data o.pa_WASO_M_CAT3DM; set o.va_WASO_M_CAT3DM; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude Prostate Issues*/
data a6mi_PI; set a6mi;
if p_hyperplasia=1 then delete; 
if p_disorders=1 then delete;
if p_inflammation=1 then delete;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_PI; by n_eid _imputation_; run;
proc sort data=a6mi_PI; by _imputation_; run;
proc phreg data=a6mi_PI; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3PI CovB=o.cov_E_WASO_M_CAT3PI; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3PI CovB=o.cov_E_WASO_M_CAT3PI;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3PI; 
run;

data o.pa_WASO_M_CAT3PI; set o.va_WASO_M_CAT3PI; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude first 2 years*/
data a6mi_2; set a6mi;
where accel_followup_days>730.5;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_2; by n_eid _imputation_; run;
proc sort data=a6mi_2; by _imputation_; run;
proc phreg data=a6mi_2; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT32yr CovB=o.cov_E_WASO_M_CAT32yr; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT32yr CovB=o.cov_E_WASO_M_CAT32yr;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT32yr; 
run;

data o.pa_WASO_M_CAT32yr; set o.va_WASO_M_CAT32yr; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude first 3 years*/
data a6mi_3; set a6mi;
where accel_followup_days>1095.75;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_3; by n_eid _imputation_; run;
proc sort data=a6mi_3; by _imputation_; run;
proc phreg data=a6mi_3; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT33yr CovB=o.cov_E_WASO_M_CAT33yr; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT33yr CovB=o.cov_E_WASO_M_CAT33yr;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT33yr; 
run;

data o.pa_WASO_M_CAT33yr; set o.va_WASO_M_CAT33yr; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*Exclude first 4 years*/
data a6mi_4; set a6mi;
where accel_followup_days>1461;
run;

/*WASO_M_CAT3*/
proc sort data=a6mi_4; by n_eid _imputation_; run;
proc sort data=a6mi_4; by _imputation_; run;
proc phreg data=a6mi_4; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT34yr CovB=o.cov_E_WASO_M_CAT34yr; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT34yr CovB=o.cov_E_WASO_M_CAT34yr;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT34yr; 
run;

data o.pa_WASO_M_CAT34yr; set o.va_WASO_M_CAT34yr; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;
