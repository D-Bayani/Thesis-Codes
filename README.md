# Thesis-Codes

Thesis Title: Managing Uncertainty in Oncology Drug Assessments: Applications in Multiple Myeloma

Author: Diana Bayani, PhD Candidate at NUS (Singapore) | dbayani@u.nus.edu


_Instructions per chapter below:_


**Chapter 2**
The IPD data set should first be set up to contain the following variables: USUBJID (subject codes), ARM (intervention or control), ISS, ECOG, AGE, SEX, EVENT (progression or outcome), and TIME. A separate spreadsheet (.csv) containing the summary measures from the trial with aggregate data must be prepared. This will be uploaded into R and used in the analysis. The spreadsheet should contain the following information: N, AGE, AGE.SD, SEX (prop.male), ECOG (prop.ecog1), and ISS (prop.iss3). We refer to the trial with IPD as the “Main trial” and the trial with aggregate data as the “Anchor trial”.


**Chapter 3**
The data set should be set up to contain patient-level anonymized data containing separate sheets for demographic data, billing, case movement, service utilization, diagnosis, and pharmacy prescribing data. 


**Chapter 4**
Patient-level data or pseudo patient-level data must be available for each trial. To generate pseudo patient-level data, a digitizer software must be used together with the validated algorithm by Guyot, et al. This can also be constructed via this link: https://www.trialdesign.org/one-page-shell.html#IPDfromKM


**Chapter 7**
The code for running a probabilistic sensitivity analysis on R using the partitioned survival model is presented below. For each scenario or scheme implemented, a different .csv file must be saved for parameter values, costs, and outcomes. The outputs may be uploaded onto https://savi.shef.ac.uk/SAVI/ to generate the EVPI/P-SUB outputs.


History:
10/25/2023: Codes uploaded
11/15/2023: Code for PSM updated
