/***************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: STable_7_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s) and Figure(s): Supplementary Table 7.

 PROGRAM FUNCTION: Generate summary statistics and HR, 95%CI presented in Supplementary Table 7 of published manuscript.
				   Analyses based on derived analytic sample of 34,260 men for prostate cancer analyses.
				   Data are generated from pre-imputed and imputed data.
***************************************************************************************************************************************/

******************************************************************************;
*     SYSTEM OPTIONS                                                         *;
******************************************************************************;
options nofmterr;
******************************************************************************;
*   INPUT FILES                                                              *;
******************************************************************************;
libname i '';/*Read in input library*/
******************************************************************************;
*   OUTPUT FILES                                                             *;
******************************************************************************;
libname o ''; /*Set output library*/

****************************************************************************;
*  PROGRAM                                                                  *;
****************************************************************************;
/*READ IN ANALYTIC DATASET; n=34,260*/
data a5; set i.X; run;

/*SUPPLEMENTAL TABLE 7 EXPOSURE-OUTCOME FREQUENCIES*/
/*Prostate Cancer Incidence*/
%macro fq (v1);
proc freq data=a5; 
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq(v1=waso_m_cat3);
%fq(v1=waso_m_cat);
%fq(v1=rank_WASO);

/*Non-fatal Prostate Cancer Incidence*/
%macro fq (v1);
proc freq data=a5; 
where pca_mort2=0;
table &v1*pca/nopercent norow nocol nocum;
run;
%mend;
%fq(v1=waso_m_cat3);
%fq(v1=waso_m_cat);
%fq(v1=rank_WASO);

/*Prostate Cancer Mortality*/
%macro fq (v1);
proc freq data=a5; 
table &v1*pca_mort2/nopercent norow nocol nocum;
run;
%mend;
%fq(v1=waso_m_cat3);
%fq(v1=waso_m_cat);
%fq(v1=rank_WASO);



/*Table 7 MODEL*/
/*REDUCED PRIMARY ANALYSES MODELS*/
data a6mi; 
/*SET IMPUTED DATASET FOR ANALYSES*/
set i.X;
run;

proc sort data=a6mi; by n_eid _imputation_; run;
proc sort data=a6mi; by _imputation_; run;

/*WASO_M_CAT3*/
proc phreg data=a6mi; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') bmicat2;
   model accel_followup_days*pca(0)=WASO_M_CAT3 accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT37R CovB=o.cov_E_WASO_M_CAT37R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT37R CovB=o.cov_E_WASO_M_CAT37R;
class WASO_M_CAT3 bmicat2;
modeleffects WASO_M_CAT3 accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_M_CAT37R; 
run;

data o.pa_WASO_M_CAT37R; set o.va_WASO_M_CAT37R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc phreg data=a6mi; 
  by _imputation_;
  class WASO_M_CAT(ref='0') bmicat2;
   model accel_followup_days*pca(0)=WASO_M_CAT accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT7R CovB=o.cov_E_WASO_M_CAT7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT7R CovB=o.cov_E_WASO_M_CAT7R;
class WASO_M_CAT bmicat2;
modeleffects WASO_M_CAT accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_M_CAT7R; 
run;

data o.pa_WASO_M_CAT7R; set o.va_WASO_M_CAT7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*RANK_WASO*/
proc phreg data=a6mi; 
  by _imputation_;
  class rank_WASO(ref='0') bmicat2;
   model accel_followup_days*pca(0)=rank_WASO accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_WASO7R CovB=o.cov_E_rank_WASO7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_WASO7R CovB=o.cov_E_rank_WASO7R;
class rank_WASO bmicat2;
modeleffects rank_WASO accel_entry_age bmicat2;
ods output parameterestimates=o.va_rank_WASO7R; 
run;

data o.pa_rank_WASO7R; set o.va_rank_WASO7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc phreg data=a6mi; 
  by _imputation_;
  class bmicat2;
   model accel_followup_days*pca(0)=WASO_H_AVG accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_H_AVG7R CovB=o.cov_WASO_H_AVG7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_H_AVG7R CovB=o.cov_WASO_H_AVG7R;
class bmicat2;
modeleffects WASO_H_AVG accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_H_AVG7R; 
run;

data o.pa_WASO_H_AVG7R; set o.va_WASO_H_AVG7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*NON-FATAL PROSTATE CANCER ANALYSES*/
data dat_impute1; 
/*SET IMPUTED DATASET FOR ANALYSES*/
set i.X;
/*Exclude prostate cancer-specific fatalities among those diagnosed with prostate cancer over followup for Non-fatal analyses*/
if pca_mort2=1 then delete;
run;

data a6minf; 
rename nswork_v7=nswork; 
set dat_impute1;
run;

proc sort data=a6minf; by n_eid _imputation_; run;
proc sort data=a6minf; by _imputation_; run;

/*WASO_M_CAT3*/
proc phreg data=a6minf; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3 CovB=o.cov_E_WASO_M_CAT3; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3 CovB=o.cov_E_WASO_M_CAT3;
class WASO_M_CAT3 nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT3 nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3; 
run;

data o.pa_WASO_M_CAT3; set o.va_WASO_M_CAT3; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc phreg data=a6minf; 
  by _imputation_;
  class WASO_M_CAT(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT CovB=o.cov_E_WASO_M_CAT; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT CovB=o.cov_E_WASO_M_CAT;
class WASO_M_CAT nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_M_CAT nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_M_CAT; 
run;

data o.pa_WASO_M_CAT; set o.va_WASO_M_CAT; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*RANK_WASO*/
proc phreg data=a6minf; 
  by _imputation_;
  class rank_WASO(ref='0') nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=rank_WASO nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_WASO CovB=o.cov_E_rank_WASO; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_WASO CovB=o.cov_E_rank_WASO;
class rank_WASO nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects rank_WASO nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_rank_WASO; 
run;

data o.pa_rank_WASO; set o.va_rank_WASO; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc phreg data=a6minf; 
  by _imputation_;
  class nswork bmicat2 ever_smoke HEALTH alc3 educ income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
   model accel_followup_days*pca(0)=WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_H_AVG CovB=o.cov_WASO_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_H_AVG CovB=o.cov_WASO_H_AVG;
class nswork bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
modeleffects WASO_H_AVG nswork accel_entry_age bmicat2 ever_smoke HEALTH alc3 educ 
income rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee tea diabetes;
ods output parameterestimates=o.va_WASO_H_AVG; 
run;

data o.pa_WASO_H_AVG; set o.va_WASO_H_AVG; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*reduced models*/
/*WASO_M_CAT3*/
proc phreg data=a6minf; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') bmicat2;
   model accel_followup_days*pca(0)=WASO_M_CAT3 accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT37R CovB=o.cov_E_WASO_M_CAT37R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT37R CovB=o.cov_E_WASO_M_CAT37R;
class WASO_M_CAT3 bmicat2;
modeleffects WASO_M_CAT3 accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_M_CAT37R; 
run;

data o.pa_WASO_M_CAT37R; set o.va_WASO_M_CAT37R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc phreg data=a6minf; 
  by _imputation_;
  class WASO_M_CAT(ref='0') bmicat2;
   model accel_followup_days*pca(0)=WASO_M_CAT accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT7R CovB=o.cov_E_WASO_M_CAT7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT7R CovB=o.cov_E_WASO_M_CAT7R;
class WASO_M_CAT bmicat2;
modeleffects WASO_M_CAT accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_M_CAT7R; 
run;

data o.pa_WASO_M_CAT7R; set o.va_WASO_M_CAT7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*RANK_WASO*/
proc phreg data=a6minf; 
  by _imputation_;
  class rank_WASO(ref='0') bmicat2;
   model accel_followup_days*pca(0)=rank_WASO accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_WASO7R CovB=o.cov_E_rank_WASO7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_WASO7R CovB=o.cov_E_rank_WASO7R;
class rank_WASO bmicat2;
modeleffects rank_WASO accel_entry_age bmicat2;
ods output parameterestimates=o.va_rank_WASO7R; 
run;

data o.pa_rank_WASO7R; set o.va_rank_WASO7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc phreg data=a6minf; 
  by _imputation_;
  class bmicat2;
   model accel_followup_days*pca(0)=WASO_H_AVG accel_entry_age bmicat2/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_H_AVG7R CovB=o.cov_WASO_H_AVG7R; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_H_AVG7R CovB=o.cov_WASO_H_AVG7R;
class bmicat2;
modeleffects WASO_H_AVG accel_entry_age bmicat2;
ods output parameterestimates=o.va_WASO_H_AVG7R; 
run;

data o.pa_WASO_H_AVG7R; set o.va_WASO_H_AVG7R; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*PROSTATE CANCER MORTALITY ANALYSES*/
data dat_impute2; 
/*SET IMPUTED DATASET FOR ANALYSES*/
set i.X;
/*Create lower dimensionality variables for mortality analysis*/ 
tea2=tea;
if tea=4 then tea2=3;
coffee2=coffee;
if coffee=4 then coffee2=3;
nswork2=nswork_v7;
if nswork_v7=0 then nswork2=0;
if nswork_v7 in (1,2) then nswork2=1;
if nswork_v7 in (3,4) then nswork2=2;
run;

data a6mif; rename nswork_v7=nswork; set dat_impute2;run;

proc sort data=a6mif; by n_eid _imputation_; run;
proc sort data=a6mif; by _imputation_; run;

/*WASO_M_CAT3*/
proc phreg data=a6mif; 
  by _imputation_;
  class WASO_M_CAT3(ref='0') nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
   model accel_mort_followup_days*pca_mort2(0)=WASO_M_CAT3 nswork2 accel_entry_age bmi ever_smoke alc3 
rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3 CovB=o.cov_E_WASO_M_CAT3; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3 CovB=o.cov_E_WASO_M_CAT3;
class WASO_M_CAT3 nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
modeleffects WASO_M_CAT3 nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
ods output parameterestimates=o.va_WASO_M_CAT3; 
run;

data o.pa_WASO_M_CAT3; set o.va_WASO_M_CAT3; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc phreg data=a6mif; 
  by _imputation_;
  class WASO_M_CAT(ref='0') nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
   model accel_mort_followup_days*pca_mort2(0)=WASO_M_CAT nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT CovB=o.cov_E_WASO_M_CAT; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT CovB=o.cov_E_WASO_M_CAT;
class WASO_M_CAT nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
modeleffects WASO_M_CAT nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
ods output parameterestimates=o.va_WASO_M_CAT; 
run;

data o.pa_WASO_M_CAT; set o.va_WASO_M_CAT; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*RANK_WASO*/
proc phreg data=a6mif; 
  by _imputation_;
  class rank_WASO(ref='0') nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
   model accel_mort_followup_days*pca_mort2(0)=rank_WASO nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_WASO CovB=o.cov_E_rank_WASO; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_WASO CovB=o.cov_E_rank_WASO;
class rank_WASO nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
modeleffects rank_WASO nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
ods output parameterestimates=o.va_rank_WASO; 
run;

data o.pa_rank_WASO; set o.va_rank_WASO; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc phreg data=a6mif; 
  by _imputation_;
  class nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
   model accel_mort_followup_days*pca_mort2(0)=WASO_H_AVG nswork2 accel_entry_age bmi ever_smoke alc3
rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_H_AVG CovB=o.cov_WASO_H_AVG; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_H_AVG CovB=o.cov_WASO_H_AVG;
class nswork2 ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
modeleffects WASO_H_AVG nswork2 accel_entry_age bmi ever_smoke alc3 rank_acc_pa race_bi rank_townsend psatest prostate_fam coffee2 tea2 diabetes;
ods output parameterestimates=o.va_WASO_H_AVG; 
run;

data o.pa_WASO_H_AVG; set o.va_WASO_H_AVG; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*MV-Reduced*/
proc sort data=a6mif; by n_eid _imputation_; run;
proc sort data=a6mif; by _imputation_; run;

/*WASO_M_CAT3*/
proc phreg data=a6mif; 
  by _imputation_;
  class WASO_M_CAT3(ref='0');
   model accel_mort_followup_days*pca_mort2(0)=WASO_M_CAT3 accel_entry_age bmi/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CAT3MVR2 CovB=o.cov_E_WASO_M_CAT3MVR2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CAT3MVR2 CovB=o.cov_E_WASO_M_CAT3MVR2;
class WASO_M_CAT3;
modeleffects WASO_M_CAT3 accel_entry_age bmi;
ods output parameterestimates=o.va_WASO_M_CAT3MVR2; 
run;

data o.pa_WASO_M_CAT3MVR2; set o.va_WASO_M_CAT3MVR2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_M_CAT*/
proc phreg data=a6mif; 
  by _imputation_;
  class WASO_M_CAT(ref='0');
   model accel_mort_followup_days*pca_mort2(0)=WASO_M_CAT accel_entry_age bmi/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_M_CATMVR2 CovB=o.cov_E_WASO_M_CATMVR2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_M_CATMVR2 CovB=o.cov_E_WASO_M_CATMVR2;
class WASO_M_CAT;
modeleffects WASO_M_CAT accel_entry_age bmi;
ods output parameterestimates=o.va_WASO_M_CATMVR2; 
run;

data o.pa_WASO_M_CATMVR2; set o.va_WASO_M_CATMVR2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*RANK_WASO*/
proc phreg data=a6mif; 
  by _imputation_;
  class rank_WASO(ref='0');
   model accel_mort_followup_days*pca_mort2(0)=rank_WASO accel_entry_age bmi/rl covb;
  ods output ParameterEstimates=o.PEsta_rank_WASOMVR2 CovB=o.cov_E_rank_WASOMVR2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_rank_WASOMVR2 CovB=o.cov_E_rank_WASOMVR2;
class rank_WASO;
modeleffects rank_WASO accel_entry_age bmi;
ods output parameterestimates=o.va_rank_WASOMVR2; 
run;

data o.pa_rank_WASOMVR2; set o.va_rank_WASOMVR2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;

/*WASO_H_AVG*/
proc phreg data=a6mif; 
  by _imputation_;
   model accel_mort_followup_days*pca_mort2(0)=WASO_H_AVG accel_entry_age bmi/rl covb;
  ods output ParameterEstimates=o.PEsta_WASO_H_AVGMVR2 CovB=o.cov_WASO_H_AVGMVR2; 
run;

proc mianalyze parms(classvar=ClassVal)=o.PEsta_WASO_H_AVGMVR2 CovB=o.cov_WASO_H_AVGMVR2;
modeleffects WASO_H_AVG accel_entry_age bmi;
ods output parameterestimates=o.va_WASO_H_AVGMVR2; 
run;

data o.pa_WASO_H_AVGMVR2; set o.va_WASO_H_AVGMVR2; 
HR=round(exp(Estimate),.01);
format ci $20.;

ci=" (" ||put(exp(estimate-(1.96*(stderr))),6.2)|| ", " ||put(exp(estimate+(1.96*(stderr))),6.2)|| ")";
ci=compress(ci," ");

hr_ci=" "||put(exp(estimate),6.2)||" "||put(ci,$20.)||" ";
hr_ci=left(hr_ci);
label hr_ci="HR (95% CI)";
run;
