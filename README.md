# 22q
Project which examines the effect of remote CNB testing on 22q research participants

# Processing Scripts

The project begins with three files - one CNB summary file which contains correct responses and reaction time on a variety of CNB tests and two itemwise data sets (one for in-person tests and the other for remote tests) which give question by question responses for each of the tests. Processing begins with the Add_auto_validation_flags.R script, which filters the data set to create a 22q cross-sectional data set for the project. Also, flags are added on a test by test basis to label tests which have poor quality (the auto-validation rules document is used to create the flags). 

Next, the data set which is output from the Add_auto_validation_flags.R script is input into the Add_SMVE.R script, which adds a second quality control flag procedure based on person-fit statistics. This script also removes individuals if their entire session is low quality. 

Now that quality control flags have been added to the data set, the scripts Remote&SexEffectsCNB.R and Remote_models_factors.R can be run on the output of Add_SMVE.R. The Remote&SexEffectsCNB.R script generates regression tables and model plots which examine the effect of test location on accuracy and reaction time on *individual* tests, while the Remote_models_factors.R script does the same thing on *composite* scores. 
