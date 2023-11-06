/*********************************************************************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PUBLICATION NAME: "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank"
 PROGRAM NAME: STable_1_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 
 Publication Table(s) and Figure(s): Supplementary Table 1.

 PROGRAM FUNCTION: Generate summary statistics of characteristics presented in Supplementary Table 1 of published manuscript.
				   Analyses based on full UK Biobank cohort (502,366) and full actigraphy cohort (103,696).
				   Data are derived from full UK Biobank cohort dataset and are subset to actigraphy cohort using S_90003_0_0=. (The request device start date).
				   Please note that additional participants will have been removed from the UK Biobank since this publication.
				   Therefore, distributions of summary statistics will have changed, but should be similar to the published manuscript.
				   Finally, bmiwho, education (qualifications), and employed are based on arrays of baseline variables using an array function run in a separate script that will be made
				   available upon request.
********************************************************************************************************************************************************************************************/

******************************************************************************;
*     SYSTEM OPTIONS                                                         *;
******************************************************************************;
options nofmterr;
******************************************************************************;
*   INPUT FILES                                                              *;
******************************************************************************;
libname i ''; /*Read in input library*/
****************************************************************************;
*  PROGRAM                                                                  *;
****************************************************************************;

/*Read in combined datasets and keep covariates for analyses*/
data X;  set i.X;
keep 
n_eid
current_analysis /*Exclusion code to subset to current analytic population based on UK Biobank provided exclusion documents*/
exclusion_current /*Exclusion code to subset to current analytic population based on UK Biobank provided exclusion documents*/
S_90003_0_0 /*Accelerometer start date*/
N_21022_0_0 /*Age at recruitment*/
N_1160_0_0 /*Baseline, self-reported sleep duration*/
N_189_0_0 /*Townsend Deprivation Index*/
employed /*Arrayed, baseline employment variable*/
N_3426_0_0 /*Job involves shift work*/
N_826_0_0 /*Job involves shift work*/
n_31_0_0 /*SEX*/
bmiwho /*Based on N_23104_0_0; N_21001_0_0*/
N_2178_0_0 /*Overall health rating*/
N_1558_0_0 /*Alcohol intake frequency*/
N_20116_0_0 /*Smoking status-ever/never/current*/
educ /*Arrayed, baseline education variable*/
N_738_0_0 /*Average household income before tax*/
N_1498_0_0 /*Coffee Intake*/
N_1488_0_0/*Tea intake*/
N_21000_0_0/*Ethnic background*/
N_2443_0_0/*Diabetes status*/
;
run;

data prca; set X;
where current_analysis=1;
run;

data prca1; set prca;
if exclusion_current=1 then delete;
run;

/*CREATING COVARIATES OF INTEREST*/
data prca2; set prca1;

/*Sleep Duration, continuous*/
sleepdur=.;
if N_1160_0_0>=0 then sleepdur=N_1160_0_0;
if N_1160_0_0<0 then sleepdur=.;

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

/*Self-reported sex*/
sex=n_31_0_0;

/*Recode of bmiwho variable*/
bmicat2=bmiwho;
if bmiwho=99 then bmicat2=.;
if bmiwho=0 then bmicat2=1;

/*Overall health status*/
HEALTH=N_2178_0_0;
if N_2178_0_0<0 then HEALTH=.;

/*Alcohol intake frequency recode*/
alc2=N_1558_0_0;
if  N_1558_0_0<0 then alc2=.;
alc3=.;
if alc2=1 then alc3=1;
if alc2 in (2, 3) then alc3=2;
if alc2=4 then alc3=3;
if alc2=5 then alc3=4;
if alc2=6 then alc3=5;

/*Self-reported smoking status*/
ever_smoke=N_20116_0_0 ;
if N_20116_0_0 <0 then ever_smoke=.;

/*Self-reported income*/
income=N_738_0_0;
if N_738_0_0<0 then income=.;

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

/*Self-reported race-ethnicity*/
race_v1=0;
if N_21000_0_0 in (1, 1001, 1002, 1003) then race_v1=1; /*white*/
if N_21000_0_0 in (4, 4001, 4002, 4003) then race_v1=2; /*black*/
if N_21000_0_0 in (2, 2001, 2002, 2003, 2004) then race_v1=3; /*Mixed*/
if N_21000_0_0 in (3, 3001, 3002, 3003, 3004, 5) then race_v1=4; /*Asian*/
if N_21000_0_0=6 then race_v1=5; /*Other*/
if N_21000_0_0=. then race_v1=.;/*Missing*/
if N_21000_0_0<0 then race_v1=.;/*Missing*/

/*Diabetes status*/
diabetes=N_2443_0_0;
if N_2443_0_0<0 then diabetes=.;

run;

data prca3; set prca2;
if S_90003_0_0=. then delete;*keeps accelerometer sample only;
run;

/*Supplementary Table 1S*/

/*Baseline Cohort*/
proc means data=prca2 n nmiss mean std min q1 median q3 max maxdec=2; 
var 
N_21022_0_0/*Age at baseline*/ 
sleepdur
n_189_0_0/*Townsend Deprivation Index*/ 
;
run;

proc freq data=prca2;
table
nswork_v7
sex
bmicat2 
HEALTH 
alc3
ever_smoke
educ
income
coffee 
tea
race_v1
diabetes
;
run;

/*Actigraphy subcohort*/
proc means data=prca3 n nmiss mean std min q1 median q3 max maxdec=2; 
var 
N_21022_0_0/*Age at baseline*/ 
sleepdur /*Sleep Duration (self-report, baseline)*/
n_189_0_0/*Townsend Deprivation Index*/ 
;
run;

proc freq data=prca3;
table
nswork_v7
sex
bmicat2 
HEALTH 
alc3
ever_smoke
educ
income
coffee 
tea
race_v1
diabetes
;
run;
