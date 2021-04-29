Read me file for Major Trauma Triage Study Helath Economic Model

Please read Pollard et al for more details.


Below are details of the files contained within this repository:
README.txt - this file
LICENCE - licence file 

Patient Characteristics files
means.csv - this file that contains the mean values for patient's age (years), gender (1= male, 0=female), Injury severity score (ISS), Glasgow Coma Scale (GCS), Blunt trauma injury status (1= blunt trauma, 0=penetrating trauma) 
covaraiance_dutch_v2.csv - this file is a csv file that contains the variance-coviarance matrix to correlate all of the above patient characteristics when sampling individual patients within the model. 

age_tab_dutch_v2.csv - this file is a csv file that contains the age groups (years) for patients with suspected major trauma, numbers of people in each group and values used to simulate the age disribution of patients in the model
blunt_tab_dutch_v2.csv - this file is a csv file that contains the calssifcation of patients into blunt or penetrating trauma (1= blunt, 0 = pentrating), numbers of people with each type of injury and values used to simulate whether patients have blunt or penetrating trauma
GCS_tab_dutch_v2.csv - this file is a csv file that contains the Galsogw Coma Scale (GCS) for patients with suspected major trauma, numbers of people with each GCS score and values used to simulate the disribution Glasgow Coma Scale for patients in the model
ISS_tab_dutch_v2.csv - this file is a csv file that contains the Injury Severity Scale (ISS) for patients with suspected major trauma, numbers of people with each ISSscore and values used to simulate the disribution ISS for patients in the model
male_tab_dutch_v2.csv - this file is a csv file that contains the classification of patients to either male or female (1=male, 0=female), numbers of people of either gender and values used within the model simulation to determine gender.

N.B. - for full public access, it was necessary to merge some categories in these files to prevent the reporting of small numbers. Access to the original files can be made available upon a resonable reques to: d.j.pollard@sheffield.ac.uk

Model parameter files
lifetime-health care costs.csv - this file is a csv file that lists the lifetime health cre costs for patients in the UK by age and gender.
ONSlifetables.csv - this file contains a copy of the ONS lifetables 20XX-20XX
parameters.csv -this file contains a copy of all the paramters used within the model

R files
All analyses in the paper were run using R version 4.0.2
ideal-winner.Rproj - this is the R project file that must be loaded prior to running the model
Master major trauma model file.R - this file contains all of the user defined model inputs and code required to run the model that is needed to produce a set of results with this model
Functions.R - this file containts all of the code that is used to run the model in R

Folders
These folders contain all the variations of the Master Majr Trauma model file.R that are needed to reporoduce our analysis
Deterministic - This folder contains all variants of the Master major trauma model file.R needed to produce our base case deterministic results (note not recommended as this model requires probabilistic analysis)
Base Case - New TARN - This folder contains all variants of the Master major trauma model file.R needed to produce our base case analysis using a probabilistic sensitivity analysis
Base Case - This folder is an old version, that contains all variants of the Master major trauma model file.R needed to produce an analysis using the 2005 TARN risk equation for within hospital mortality in UK for patients with suspected major trauma
25% benefit - ISS 9 to 15 - This folder contains all variants of the Master major trauma model file.R needed to produce our sensitivity analysis in which patients with an ISS between 9 and 15 (inclusive) recieve 25% of the benefit from MTC care compared to patients with an ISS of 16 or over
50% benefit - ISS 9 to 15 - This folder contains all variants of the Master major trauma model file.R needed to produce our sensitivity analysis in which patients with an ISS between 9 and 15 (inclusive) recieve 50% of the benefit from MTC care compared to patients with an ISS of 16 or over
75% benefit - ISS 9 to 15 - This folder contains all variants of the Master major trauma model file.R needed to produce our sensitivity analysis in which patients with an ISS between 9 and 15 (inclusive) recieve 75% of the benefit from MTC care compared to patients with an ISS of 16 or over
