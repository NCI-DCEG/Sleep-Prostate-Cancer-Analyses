/***************************************************************************************************************************************
 STUDY NAME: Actigraphy sleep and prostate cancer analysis
 PROGRAM NAME: Array_Var_Code_103023.sas 
 PROGRAM LOCATION: NCI DCEG GitHub Repository
 PROGRAMMER: Josh Freeman 

 PROGRAM FUNCTION: Create derived education and employment variables based on array data from UK Biobank.
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

/*READ IN ANALYTIC DATASET FROM BASELINE (DATA IS SUBSET OF ONLY ID, EMPLOYMENT, AND QUALIFICATIONS VARIABLES FROM UK Biobank*/
data a; set i.X; run;

/*QUALIFICATIONS RECODE*/
/*Recode using array call*/
data b; set a;
array qual(6) N_6138_0_0-N_6138_0_5;
qualv2=.;
do i=1 to 6;
if qual(i)>. & qual(i)=-3 then qualv2=.; /*Prefer not to answer-Recoded as missing*/
if qual(i)>. & qual(i)=-7 then qualv2=3; /*Do not know*/
if qual(i)>. & qual(i) in (1, 5, 6) then qualv2=1; /*College or University Degree; Other qualifications; NVQ/HND/HNC/Equivalent*/
else if qual(i)>. & qual(i) NOT in (-3, -7) & qualv2^=1 then qualv2=2; /*A/AS; O/GCSE; CSE or equivalents*/
end;
drop i;
run;

/*EMPLOYMENT RECODE*/
/*Recode using array call*/
data c; set b;
array employ(7) N_6142_0_0-N_6142_0_6;
do i=1 to 7;
if employ{i}=1 then employ_d=1; /*Employed/Self-employed*/
else if employ{i}=2 & employ_d^=1 then employ_d=2; /*Retired*/
else if employ{i} in (-7, 3, 4, 5, 6, 7) & employ_d^=1 & employ_d^=2 then employ_d=3; /*None of the above; Unemployed; Student; Volunteer work; Looking after family; Unable to work*/
else if employ{i}=-3 then employ_d=.; /*Prefer not to answer-recoded as missing*/
end;
drop i;
run;

/*Create derived output dataset*/
data o.array_employ_qual; set c;
LABEL     
employ_d="Employment-Derived Array Var"
qualv2="Qualifications-Derived Array Var";
keep n_eid employ_d qualv2;
run;
