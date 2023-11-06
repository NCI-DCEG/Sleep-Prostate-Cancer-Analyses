# Sleep-Prostate-Cancer-Analyses

Contained in this directory are annotated, analytic scripts created in SAS v9.4 that were used to create results output 
presented in the mansucript, "Actigraphy-derived measures of sleep and risk of prostate cancer in the UK Biobank."

Each subdirectory contains a script specific to a table or figure presented in the manuscript, except in some cases in which code was aggregated into a single script for processing efficiency.
Additionally included are the setup file (to create derived variables used in the analyses), a separate array file (to array  variables), and the imputation script (to impute the data).

Not included in the directory are scripts to derive the analytic sample in figure 1 as this was completed across separate scripts, some of which are being returned to the UK Biobank. 
The setup file partially shows some of the sample derivation from the number of men with valid accelerometer data through exclusion of prevalent cancer to achieve the analytic sample (n=34,260).

Data used in these analyses are not included given the data are provided by the UK Biobank and we do not have permission to share it as per UK Biobank Application Number 43456.
UK Biobank data are globally available to approved researchers through the UK Biobank research portal (https://www.ukbiobank.ac.uk/). Coding and data for the accelerometer-derived sleep 
and physical activity variables are not included given these files are being returned to the UK Biobank, which will provide a digital object identifier for the data and code, and 
will be released at a later date.

If there are any questions or comments, please reach out to Josh Freeman (josh.freeman@nih.gov) or Charles Matthews (charles.matthews2@nih.gov).

List of directories and brief description of the file within:

1. Array_Vars: Array_Var_Code_103023.sas--Creates derived education and employment variables based on array data from UK Biobank.
2. Imputation: PrCA_Imputation_103023.sas--Imputes missing covariate data.
3. Setup-F1-T5: PrCA_Setup_ST5_F1_103023.sas--Sets up pre-imputed dataset and recodes variables to be used in analyses and imputed. Also provides list of variables from Supplementary Table 5 & partial derivation of the analytic sample in Figure 1.
4. Table 1: Table_1_103023.sas--Generates summary statistics of characteristics presented in Table 1 of published manuscript and in first paragraph of results.
5. Table 2: Table_2_103023.sas--Generates results presented in Table 2 of manuscript.
6. Table 3: Table_3_103023.sas--Generates results presented in Table 3 of manuscript.
7. Table S1: STable_1_103023.sas--Generates results presented in Supplementary Table 1 of Supplementary Material.
8. Table S6-F2: STable_6_F2_103023.sas--Generates results presented in Supplementary Table 6 and Figure 2 of Supplementary Material and manuscript.
9. Table S7: STable_7_103023.sas--Generates results presented in Supplementary Table 7 of Supplementary Material.
10. Table S234: STable_234_103023.sas--Generates results presented in Supplementary Tables 2, 3, & 4 of Supplementary Material.



