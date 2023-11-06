/**************************************************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: PrCA_Setup_ST5_F1_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman
 Table(s) & Figure(s): Supplemental Table 5; Partial Figure 1.
 
 PROGRAM FUNCTION: Generates complete dataset for analyses in published manuscript and provides variables listed in Supplemental Table 5.
				   Data come from UK Biobank baseline data, and derived files based on cancer registry, death linkage data, and accelerometry data (Data Field 90001)
				   Derives analytic sample of 34,260 men for prostate cancer analyses illustrating some exclusions presented in Figure 1.
				   Generates non-imputed dataset and has annotations for variables.
				   Data are pre-imputation.
***************************************************************************************************************************************************************************/

******************************************************************************;
* OPTIONS                                                                    *;
******************************************************************************;
options nofmterr;
******************************************************************************;
* INPUT FILES                                                                *;
******************************************************************************;
libname i ''; /*Read in input library*/
******************************************************************************;
* OUTPUT FILES                                                               *;
******************************************************************************;
libname o ''; /*Set output library*/
****************************************************************************;
*  PROGRAM                                                                 *;
****************************************************************************;

/*Read in dataset and keep variables of interest*/
data prca;  
set i.X;
keep 
n_eid/*Participant ID */
id2keep /*Flag variable for men with valid data in accelerometer cohort*/
accel_followup_start /*Derived end date of participant follow-up from accel protocol*/
current_analysis /*Current analysis flag based on exclusion criteria released by UK Biobank*/
DOB /*Participant date of birth based on N_34_0_0 (Year of birth) & N_52_0_0 (Month of Birth) */
INC_ENTRY_DATE /*Date of entry into incidence analyses based on assessment date (S_53_0_0) */
INC_ENTRY_AGE /*Calculated variable based on DOB and assessment date (S_53_0_0)*/
INC_EXIT_STATUS/*Binary variable to indicate incident cancer (0,1)*/
/*Based on S_40005, S_40006, S_40011, S_40012, S_40013, S_53_0_0, Center censor date, death date from UK Biobank death file, and date of lost to follow-up, whichever came first; Missing if prevalent cancer */
INC_EXIT_DATE /*Date of incidence analysis exit */
/*Based on S_40005, S_40006, S_40011, S_40012, S_40013, S_53_0_0, Center censor date, death date from UK Biobank death file, and date of lost to follow-up, whichever came first; Missing if prevalent cancer */
INC_EXIT_ICDCODE /*Cancer ICD Code in incidence analyses based on S_40005, S_40006*/
MORT_ENTRY_DATE /*Date of entry into mortality analyses based on assessment date (S_53_0_0) */
MORT_EXIT_DATE /*Mortality analysis exit date based on assigned date of death based on maximum ascertainment date for study centers if not dead, ltf_date if lost to follow-up, or death date if dead based on death files */
MORT_EXIT_ICDCODE /*Based on ICD Code for priary case of death from UK Biobank Death Linkage file*/
MORT_EXIT_STATUS /*Binary variable for whether a participant died during follow-up*/
PREVALENT_CANCER_REG /*Prevalent cancer flag (0,1) if any registry-reported cancer before assessment date(S_53_0_0) based on S_40005, S_40006, S_40011, S_40012, S_40013*/
IS_LTF /*Flag for lost to follow-up based on S_191_0_0*/
IS_DEAD /*Flag for death (0,1) based on UK Biobank death files data*/
prostate_fam /*Family history of prostate cancer (N_20107, N_20110, & N_20111) */
N_31_0_0 /*Self-reported sex*/
N_1558_0_0 /*Alcohol intake frequency*/
N_20116_0_0 /*Smoking status-ever/never/current*/
N_189_0_0 /*Townsend Deprivation Index-New version-22189 on UK Biobank showcase*/
N_738_0_0 /*Average household income before tax*/
N_21000_0_0/*Ethnic background*/
N_2178_0_0 /*Overall health rating*/
N_22001_0_0 /*Genetic Sex*/
N_2365_0_0 /*Ever had PSA Test-Baseline*/
N_3426_0_0 /*Job involves shift work-Baseline-Only among those employed and not endorsing never/rarely in shift work*/
N_826_0_0 /*Job involves shift work-Baseline-Only among those endorsing ever employed*/
N_2443_0_0/*Self-reported Diabetes Status*/
BMI /*BMI (continuous) based on N_23104_0_0 & N_21001_0_0 */
BMIWHO /*BMI (catgegorical) based on N_23104_0_0 & N_21001_0_0 */
educ /*Arrayed education/qualifications variable based on N_6138*/
employed /*Arrayed employment variable based on N_6142*/

N_1488_0_0/*Tea intake*/
N_1498_0_0 /*Coffee Intake*/

N_132073_0_0 /*Source of report of N40 (hyperplasia of prostate)*/
N_132075_0_0 /*Source of report of N41 (inflammatory diseases of prostate)*/
N_132077_0_0 /*Source of report of N42 (other disorders of prostate)*/
S_132072_0_0/*Date N40 first reported (hyperplasia of prostate)*/
S_132074_0_0 /*Date N41 first reported (inflammatory diseases of prostate)*/
S_132076_0_0 /*Date N42 first reported (other disorders of prostate)*/

/*Sleep and physical activity variables derived from UK Biobank data field 90001*/
sleepduration_acc_6 /*Sleep Duration, categorical (original)*/
SLEEPDURATION_H_AVG /*Sleep Duration, continuous*/
midsleep_avg_h_cat /*Sleep Midpoint, categorical*/
SLEEPMIDPOINT_T_IAVG /*Sleep Midpoint, continuous*/
SLEEPONSET_H_AVG /*Sleep onset, continuous*/
WAKEUP_H_AVG /*Wakeup time, continuous*/
sleep_efficiency_per_avg /*Sleep Efficiency, continuous (original)*/
SOCIALJETLAG_H /*Social jetlag, continuous*/
WASO_M_CAT /*WASO, categorical (binary)*/
WASO_M_CAT3 /*WASO, categorical (4-level)*/
WASO_H_AVG /*WASO, continuous*/
WASOC30M /*Frequency of WASO, continuous*/
WAKE_ACTIVE_H_AVG /*Active WASO, continuous*/
MVPA_HRS_0M_AVG/*Moderate-Vigorous Physical Activity (hours/day)*/
;
run;


/*Subset to current_analysis (n=502,366). Note this has changed to n=502,356 based on 10 additional exclusions that occurred after completion of analysis*/
data prca; 
set prca;
where current_analysis=1;
sex=N_31_0_0;
run;

/*SUBSET TO ACCEL SUBCOHORT USING ID2KEEP n=37142 UNIQUE IDs-8 missing due to subset from current_analysis excluding those no longer in the study*/
/*RENAME WASOCOUNT VARIABLE AND ROUND FOR CATEGORICAL VARIABLE CODING*/
data prca1; set prca;
where id2keep=1;
WASOCOUNT_30M_STANDARD2=round(WASOCOUNT_30M_STANDARD, 1);
run;

/*EXCLUDE PARTICIPANTS WITH EXCLUSION_CURRENT=1 TO ENSURE DATA IS SUBSET TO POPULATION CONSENTING TO REMAIN IN THE UK BIOBANK; n=37142 unique IDs*/
/*CODE PROSTATE CANCER OUTCOMES--INCIDENCE (pca), ALL PROSTATE CANCER MORTALITY (pca_mort), and PROSTATE CANCER MORTALITY BASED ON INCIDENT PROSTATE CANCER (pca_mort2)*/
data prca2; set prca1;
if exclusion_current=1 then delete;
/*Prostate cancer*/
pca=0;
if INC_EXIT_ICDCODE='C61' then pca=1;
/*All prostate cancer mortality*/
pca_mort=0; 
if MORT_EXIT_ICDCODE='C61' then pca_mort=1;
/*Prostate cancer mortality among incident cases*/
pca_mort2=0;
if MORT_EXIT_ICDCODE='C61' & pca=1 then pca_mort2=1;

run;

/*EXCLUDES PARTICIPANTS WITH MISSING SEX, WOMEN, AND THOSE WHOSE GENETIC SEX WERE FEMALE; n=37117*/
data prca3; set prca2;
where sex=1; /*Selects for men*/
if N_22001_0_0=0 then delete; /*n_22001_0_0 (Genetic female)=0*/
run;

/*EXCLUDES THOSE WITH PREVALENT CANCER BEFORE ENTRY INTO UK BIOBANK BASELINE AND CREATES FOLLOW-UP TIME BASED ON ACCEL ENTRY START; n=35677*/
data prca4; set prca3;
where PREVALENT_CANCER_REG=0; /*PREVALENT CANCER EXCLUSION*/
/*Update follow-up and age calculcation based on accel subcohort entry start rather than baseline*/
accel_followup_days=inc_exit_date-accel_followup_start; /*Time on study for incidence analyses*/
accel_mort_followup_days=mort_exit_date-accel_followup_start; /*Time on study for mortality analyses*/
accel_entry_age=(accel_followup_start-DOB)/365.25;/*Age at accel subcohort entry*/
accel_end_age=(inc_exit_date-DOB)/365.25;/*Age at end of incidence follow-up*/
accel_mort_end_age=(mort_exit_date-DOB)/365.25;/*Age at end of mortality follow-up*/
accel_followup_yrs=accel_followup_days/365.25; /*Years of accel follow-up*/
RUN;

/*FLAGS AND EXCLUDES THOSE WHO DEVELOPED CANCER, WERE LTF, OR DIED PRIOR TO ACCEL COHORT ENTRY; n=34,260-final analytic sample*/
data prca5; set prca4;
if inc_exit_status in (0,1) & inc_exit_date > accel_followup_start;
run;

/*DATA CHECK FOR THOSE WHO DIED, WERE LOST TO FOLLOW-UP, OR DEVELOPED CANCER BEFORE ACCEL COHORT ENTRY*/
data prca_check; set prca4;
if inc_exit_status in (0,1) & inc_exit_date <= accel_followup_start;
run;

proc freq data=PRCA_CHECK;
table mort_exit_status inc_exit_status is_ltf mort_exit_status*is_ltf*inc_exit_status/missing;
run;

/*SUBSETS SAMPLE TO MEN WITH ACCEL_FOLLOWUP TIME >0; n=34,260*/
data prca6; set prca5;
where accel_followup_days>0;
run;

/*CREATING ADDITIONAL SLEEP VARIABLES OF INTEREST*/
data prca7; set prca6;

/*CATEGORICAL SLEEP DURATION*/
sleepduration_acc_65=.;
if sleepduration_acc_6=0 then sleepduration_acc_65=0;/*<5 hours*/
if sleepduration_acc_6=1 then sleepduration_acc_65=1;/*<6 hours*/
if sleepduration_acc_6=2 then sleepduration_acc_65=2;/*6-7 hours*/
if sleepduration_acc_6=3 then sleepduration_acc_65=3;/*7-8 hours*/
if sleepduration_acc_6=4 then sleepduration_acc_65=4;/*8-9 hours*/
if sleepduration_acc_6=5 then sleepduration_acc_65=4;/*>9 hours*/
if sleepduration_acc_6=6 then sleepduration_acc_65=5; /*Missing*/

/*WASO >=30 MINUTE COUNTS OVER PROTOCOL STANDARDIZED TO 7 DAYS*/
wasocountv7=.;
if WASOCOUNT_30M_STANDARD2=0 then wasocountv7=1;
if WASOCOUNT_30M_STANDARD2=1 then wasocountv7=1;
if WASOCOUNT_30M_STANDARD2=2 then wasocountv7=2;
if WASOCOUNT_30M_STANDARD2=3 then wasocountv7=3;
if WASOCOUNT_30M_STANDARD2=4 then wasocountv7=4;
if WASOCOUNT_30M_STANDARD2=5 then wasocountv7=5;
if WASOCOUNT_30M_STANDARD2=6 then wasocountv7=5;
if WASOCOUNT_30M_STANDARD2=7 then wasocountv7=5;
if WASOCOUNT_30M_STANDARD2=. then wasocountv7=.;

/*SOCIAL JETLAG CATEGORIES BASED ON HOURLY CUTOFFS*/
SOCIALJETLAG_H_CAT=.;
if 0<=socialjetlag_h<1 then socialjetlag_h_cat=0;
if 1<=socialjetlag_h<2 then socialjetlag_h_cat=1;
if socialjetlag_h>=2 then socialjetlag_h_cat=2;
if socialjetlag_h=. then socialjetlag_h_cat=3;

/*RENAMED SLEEP EFFICIENCY VARIABLE DUE TO ISSUE WITH VARIABLE NAME LENGTH*/
SE_AVG=.;
SE_AVG=sleep_efficiency_per_avg;

run;

/*CREATING COVARIATES OF INTEREST-Assigning missing to "Don't Know" & "Prefer Not to Answer" responses as appropriate*/
data prca8; set prca7;

/*Overall health*/
HEALTH=N_2178_0_0;
if N_2178_0_0<0 then HEALTH=.;

/*Self-report of diabetes*/
diabetes=N_2443_0_0;
if N_2443_0_0<0 then diabetes=.;

/*Binary race-ethnicity ehtnicity variable for models-due to high proportion of men who self-reported white race in the data*/
/*Note: Categories are based on races and were concatenated from participant report of race and ethnicity
See https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000 for additional details*/
race_bi=0;/*Black, Mixed, Asian, Other*/
if N_21000_0_0 in (1, 1001, 1002, 1003) then race_bi=1; /*White*/
if N_21000_0_0=. then race_bi=.;/*Missing*/
if N_21000_0_0<0 then race_bi=.;/*Missing*/

/*Race-ethnicity variable based on top categories on UK Biobank Data Showcase: https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000 */
/*Note: Categories are based on races and were concatenated from participant report of race and ethnicity*/
race_v1=0;
if N_21000_0_0 in (1, 1001, 1002, 1003) then race_v1=1; /*White*/
if N_21000_0_0 in (4, 4001, 4002, 4003) then race_v1=2; /*Black*/
if N_21000_0_0 in (2, 2001, 2002, 2003, 2004) then race_v1=3; /*Mixed*/
if N_21000_0_0 in (3, 3001, 3002, 3003, 3004, 5) then race_v1=4; /*Asian*/
if N_21000_0_0=6 then race_v1=5; /*Other*/
if N_21000_0_0=. then race_v1=.;/*Missing*/
if N_21000_0_0<0 then race_v1=.;/*Missing*/

/*Average total household income before tax*/
income=N_738_0_0;
if N_738_0_0<0 then income=.;

/*Smoking status-never, former, current*/
ever_smoke=N_20116_0_0 ;
if N_20116_0_0 <0 then ever_smoke=.;

/*Alcohol intake frequency*/
alc2=N_1558_0_0;
if  N_1558_0_0<0 then alc2=.;

/*Concatenated, alcohol intake frequency due to data sparsity*/
alc3=.;
if alc2=1 then alc3=1;/*Daily-almost daily*/
if alc2 in (2, 3) then alc3=2;/*3-4x/week & 1-2x/week--Weekly*/
if alc2=4 then alc3=3;/*1-3x/month--Monthly*/
if alc2=5 then alc3=4;/*Special occasions only*/
if alc2=6 then alc3=0;/*Never*/

/*NOTES FOR BMI*/
/*BMI based on both n_23104_0_0 (BMI-Impedance) & n_21001_0_0 (BMI-Measurement)*/
/*Categories were presented below and were pre-coded in a derived dataset*/
/*Categories:
0 = <18.5
1 = 18.5-24.9
2 = 25.0-29.9
3 = 30.0-34.9
4 = 35.0-39.9
5 = >=40.0
99 = Unknown
*/

/*Modified BMI variable to address missingness in imputation and recode <18.5 and 18.5-24.9 groups together due to data sparsity*/
bmicat2=bmiwho;
if bmiwho=99 then bmicat2=.;/*Unknown/Missing*/
if bmiwho=0 then bmicat2=1; /*Recoded <18.5 as <=24.9*/

/*Prostate Issues-recoded based on occurrence date prior to accel follow-up start*/
/*NOTE numbers match equivalent N_ variables, which indicate source of report*/ 
/*No apparent flags indicating missing, placeholder, or potentially erroneous dates*/

/*Date of report of prostate disorders (before accel cohort entry)*/
p_disorders=.;
if accel_followup_start>S_132076_0_0 then p_disorders=1;
else p_disorders=0;
if S_132076_0_0=. then p_disorders=.;

/*Date of report of prostate hyperplasia (before accel cohort entry)*/
p_hyperplasia=.;
if accel_followup_start>S_132072_0_0 then p_hyperplasia=1;
else p_hyperplasia=0;
if S_132072_0_0=. then p_hyperplasia=.;

/*Date of report of prostate inflammation (before accel cohort entry)*/
p_inflammation=.;
if accel_followup_start>S_132074_0_0 then p_inflammation=1;
else p_inflammation=0;
if S_132074_0_0=. then p_inflammation=.;

/*History of Ever Having a PSA Test*/
psatest=.;
if n_2365_0_0>=0 then psatest=n_2365_0_0;
if n_2365_0_0<0 then psatest=.;/*Missing*/

/*Coffee Intake*/
coffee=.;
if N_1498_0_0 = 0 or N_1498_0_0= -10 then coffee=0;
if N_1498_0_0=1 or N_1498_0_0=2 then coffee=1;
if N_1498_0_0=3 or N_1498_0_0=4 then coffee=2;
if N_1498_0_0=5 or N_1498_0_0=6 then coffee=3;
if N_1498_0_0 gt 6 then coffee=4;
if N_1498_0_0=. or N_1498_0_0=-1 or N_1498_0_0=-3 then coffee=.;

/*Tea intake*/
tea=.;
if N_1488_0_0 = 0 or N_1488_0_0=-10 then tea=0;
if N_1488_0_0=1 or N_1488_0_0=2 then tea=1;
if N_1488_0_0=3 or N_1488_0_0=4 then tea=2;
if N_1488_0_0=5 or N_1488_0_0=6 then tea=3;
if N_1488_0_0 gt 6 then tea=4;
if N_1488_0_0=. or N_1488_0_0=-1 or N_1488_0_0=-3 then tea=.;

run;

/*Create quartile versions of variables based on proc rank*/

/*Townsend Deprivation Index Quartiles*/
proc sort data=prca8; by n_189_0_0;run;
proc rank data=prca8 out=prca5t groups=4;
var n_189_0_0;
ranks rank_townsend;
run;

/*Accelerometer-derived MVPA Quartiles*/
proc sort data=prca5t; by MVPA_HRS_0M_AVG;run;
proc rank data=prca5t out=a5u groups=4;
var MVPA_HRS_0M_AVG;
ranks rank_acc_pa;
run;

/*Sleep Onset Time Quartiles*/
proc sort data=a5u; by SLEEPONSET_H_AVG;run;
proc rank data=a5u out=a5v groups=4;
var SLEEPONSET_H_AVG;
ranks rank_onset;
run;

/*Wakeup Time Quartiles*/
proc sort data=a5v; by WAKEUP_H_AVG;run;
proc rank data=a5v out=a5w groups=4;
var WAKEUP_H_AVG;
ranks rank_wakeup;
run;

/*Active Wake Time Quartiles*/
proc sort data=a5w; by WAKE_ACTIVE_H_AVG;run;
proc rank data=a5w out=a5x groups=4;
var WAKE_ACTIVE_H_AVG;
ranks rank_active;
run;

/*WASO Quartiles*/
proc sort data=a5x; by WASO_H_AVG;run;
proc rank data=a5x out=a5xx groups=4;
var WASO_H_AVG;
ranks rank_WASO;
run;

/*Sleep Efficiency Quartiles*/
proc sort data=a5xx; by SLEEP_EFFICIENCY_PER_AVG;run;
proc rank data=a5xx out=a5y groups=4;
var SLEEP_EFFICIENCY_PER_AVG;
ranks rank_se;
run;

data prca9; set a5y;
/*Employment/Shift Work Combined Variable*/
nswork_v7=.;
if employed=1 then nswork_v7=0; /*Daytime Work Only*/
if employed=1 & n_826_0_0=1 then nswork_v7=0; /*Daytime Work Only*/
if employed=1 & n_3426_0_0=1 then nswork_v7=0; /*Daytime Work Only*/
if employed=1 & n_826_0_0=. then nswork_v7=0; /*Daytime Work Only*/
if employed=1 & n_3426_0_0=. then nswork_v7=0; /*Daytime Work Only*/
if employed=1 & n_826_0_0>=2 & n_3426_0_0=1 then nswork_v7=1; /*Mixed Shift Work*/
if employed=1 & n_826_0_0>=2 & n_3426_0_0 in (-1, -3) then nswork_v7=1; /*Mixed Shift Work - Missing night shift work information, but has shift work information*/
if employed=1 & n_826_0_0>=2 & n_3426_0_0>=2 then nswork_v7=2; /*Night Shift Work*/
if employed=1 & n_826_0_0 in (-1, -3) & n_3426_0_0=1 then nswork_v7=1;/*Rotating Shift Work - based on missing shift work, but answered night shift work*/
if employed=1 & n_826_0_0 in (-1, -3) & n_3426_0_0>=2 then nswork_v7=2;/*Night Shift Work - based on missing shift work, but answered night shift work*/
if employed=1 & n_826_0_0 in (-1, -3) & n_3426_0_0 in (-1, -3) then nswork_v7=0;/*Daytime Work based on missing both shift work variables*/
if employed=2 then nswork_v7=3;/*Retired*/
if employed=3 then nswork_v7=4;/*Student/Volunteer/Disabled/Looking after family/Unemployed/Do not know*/
if employed=. then nswork_v7=.; /*Prefer Not to Answer + Missing Employment Status*/
run;

/*Renamed WASO Count variables to run in SAS models-Issues due to name length*/
data a77; 
rename WASOCOUNT_30M_STANDARD2=WASOC30M2 /*Renamed WASO Count, rounded for categorical variable creation*/
WASOCOUNT_30M_STANDARD=WASOC30M;/*Renamed WASO Count, not rounded*/
set prca9;
run;

proc sort data=a77; by n_eid;run;

/*GENERATING PROSTATE INCIDENCE AND MORTALITY DATASET-Numbers for denominator remain the same; n=34,260 IDs*/
data a77_m; set a77;
if mort_exit_status in (0,1) & mort_exit_date > accel_followup_start;
run;

/*Ensuring follow-up time is >0 for analyses; denominator numbers remain the same, n=34,260*/
data a77_m2; set a77_m;
where accel_mort_followup_days>0;
run;

/*OUTPUT CURATED INCIDENCE & MORTALITY DATASET FOR MAIN COX PH ANALYSIS, & FOR MISSING DATA IMPUTATION; n=34,260 IDs*/
data o.X; 
set a77_m2; 
keep
n_eid /*Participant ID*/
SLEEPDURATION_H_AVG /*Sleep Duration, continuous*/
sleepduration_acc_65 /*Sleep Duration, categorical*/
midsleep_avg_h_cat /*Sleep Midpoint, categorical*/
SLEEPMIDPOINT_T_IAVG /*Sleep Midpoint, continuous*/
rank_onset /*Sleep onset, categorical*/
SLEEPONSET_H_AVG /*Sleep onset, continuous*/
rank_wakeup /*Wakeup time, categorical*/
WAKEUP_H_AVG /*Wakeup time, continuous*/
SE_AVG /*Sleep Efficiency, continuous*/
rank_se /*Sleep Efficiency, categorical*/
SOCIALJETLAG_H_CAT /*Social jetlag, categorical*/
SOCIALJETLAG_H /*Social jetlag, continuous*/
WASO_M_CAT /*WASO, categorical (binary)*/
WASO_M_CAT3 /*WASO, categorical (4-level)*/
rank_WASO /*WASO, cateogrical (quartiles)*/
WASO_H_AVG /*WASO, continuous*/
wasocountv7 /*Frequency of WASO, categorical*/
WASOC30M /*Frequency of WASO, continuous*/
rank_active /*Active WASO, categorical*/
WAKE_ACTIVE_H_AVG /*Active WASO, continuous*/
MVPA_HRS_0M_AVG /*Moderate-Vigorous Physical Activity, continuous*/
accel_followup_days /*Time on study for incidence analyses*/
pca /*Prostate Cancer Incidence, binary*/
accel_entry_age /*Age at entry into accel subcohort*/
accel_end_age /*Age at end of incidence analysis for accel subcohort*/
pca_mort2 /*Prostate cancer-speicifc mortality among incident prostate cancer cases*/
accel_mort_followup_days /*Time on study for mortality analyses*/
bmi /*BMI, continuous*/

nswork_v7 /*Derived employment/shift work variable, categorical*/
employed /*Recoded employment variable based on array, categorical*/
rank_acc_pa /*Moderate-Vigorous physical, categorical*/
bmicat2 /*Recoded BMI, categorical*/
HEALTH /*Recoded health status variable, categorical*/
alc3 /*Recoded alcohol variable, categorical*/
ever_smoke /*Recoded smoking status variable, categorical*/
educ /*Recoded education variable, categorical*/
income /*Recoded income variable, categorical*/
race_bi /*Reocded race-ethnicity, binary*/
RACE_V1 /*Recoded race-ethnicity, categorical*/
rank_townsend /*Townsend Deprivation Index, categorical*/
psatest /*History of PSA Testing, categorical*/
prostate_fam /*Family history of prostate cancer, categorical*/
coffee /*Coffee intake, categorical*/
tea /*Tea intake, categorical*/
diabetes /*Diabetes status, binary*/

p_disorders /*Prostate disorders before entry into accel subcohort, indicator for sensitivity analyses*/
p_hyperplasia /*Prostate hyperplasia before entry into accel subcohort, indicator for sensitivity analyses*/
p_inflammation /*Prostate inflammation before entry into accel subcohort, indicator for sensitivity analyses*/

day_count_men /*Count of the number of valid days men provided over accelerometer follow-up*/
N_189_0_0 /*Original Townsend Deprivation Index variable*/
accel_followup_yrs /*Years of follow-up in incidence analyses*/
;

run;
