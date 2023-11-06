/***************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: PrCA_Imputation_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 

 PROGRAM FUNCTION: Imputes missing data for covariates to be used in modelling.
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
libname o '';/*Set output library*/

****************************************************************************;
*  PROGRAM                                                                  *;
****************************************************************************;

/*READ IN ANALYTIC DATASET; n=34,260*/
data a5a; 
set i.X; 
/*RECODING SOCIALJETLAG_H_CAT=3 TO MISSING GIVEN IT IS MISSING INDICATOR CATEGORY*/
if socialjetlag_h_cat=3 then socialjetlag_h_cat=.;
run;

/*MISSING DATA-IMPUTATION*/;
proc sort data=a5; by n_eid;
run;

data impute;
set a5;
keep n_eid
SLEEPDURATION_H_AVG 
sleepduration_acc_65
midsleep_avg_h_cat
SLEEPMIDPOINT_T_IAVG
rank_onset
SLEEPONSET_H_AVG
rank_wakeup
WAKEUP_H_AVG
SE_AVG
rank_se
SOCIALJETLAG_H_CAT
SOCIALJETLAG_H 
WASO_M_CAT
WASO_M_CAT3
rank_WASO
wasocountv7
WASOC30M
WASO_H_AVG
rank_active
WAKE_ACTIVE_H_AVG
accel_followup_days
pca
accel_entry_age
accel_end_age
pca_mort2
accel_mort_followup_days
bmi

nswork_v7
employed
rank_acc_pa 
bmicat2 
HEALTH 
alc3 
ever_smoke 
educ 
income
race_bi 
rank_townsend 
psatest 
prostate_fam 
coffee 
tea 
diabetes

p_disorders
p_hyperplasia
p_inflammation
; 
run;

proc sort data=impute;
by n_eid;
run;

proc mi data=impute nimpute=10 out=dat_impute seed=12345;
class nswork_v7 bmicat2 alc3 health ever_smoke educ income employed race_bi psatest coffee tea diabetes rank_townsend;
var nswork_v7 bmi bmicat2 alc3 health ever_smoke educ income employed race_bi psatest coffee tea diabetes rank_townsend; 
fcs discrim (nswork_v7 bmicat2 alc3 health ever_smoke educ income employed race_bi psatest coffee tea diabetes rank_townsend/classeffects=include);
run;

proc sort data=dat_impute; by n_eid _imputation_;
run;

/*Set output dataset that is multiply imputed*/
data o.X; set dat_impute;run;



