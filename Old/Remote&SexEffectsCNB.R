# Import the necessary functions and data
library(ggdist)
library(lme4)
library(ggpubr)
library(lmerTest)
library(visreg)
library(ggh4x)
library(gridExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(MatchIt)
library(gtsummary)
library(table1)
library(tidyverse)

CNB_cross_qc <- read_csv("/Users/hillmann/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")
codebook <- read_csv("~//Projects/22q/Data/Summary/bbl_cross_battery_codebook.csv")


# Turn poor quality data to NA

CNB_before_qc_tests_completed <- CNB_cross_qc %>% 
  select(test_sessions.bblid,test_sessions_v.gender,test_sessions_v.age,test_sessions_v.battery,remote,test_num,matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("_CR$"),matches("SCTAP_TOT$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("TP$"),matches("_flag$")) %>% 
  select(-CPF_B.CPF_W_RTCR,-SPCPTN90.SCPN90_TPRT,-SPCPTNL.SCPN_TPRT,-SPCPTN90.SCPN90_TP,-SPCPTNL.SCPN_TP) %>% 
  rowwise() %>% 
  mutate(Tests_completed_pre_qc = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SPCPTN90_TP)))) %>% 
  ungroup()

CNB_cross_clean <- CNB_cross_qc %>% 
  select(test_sessions.bblid,test_sessions_v.gender,test_sessions_v.age,test_sessions_v.battery,remote,test_num,matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("_CR$"),matches("SCTAP_TOT$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("TP$"),matches("_flag$")) %>% 
  select(-CPF_B.CPF_W_RTCR,-SPCPTN90.SCPN90_TPRT,-SPCPTNL.SCPN_TPRT,-SPCPTN90.SCPN90_TP,-SPCPTNL.SCPN_TP) %>% 
  rowwise() %>% 
  mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,ADT36_A.ADT36A_CR),.fns = ~ ifelse(sum(c(PFscores_ADT_flag,ADT_comment_flag,ADT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(CPF_B.CPF_RTCR,CPF_B.CPF_CR),.fns = ~ ifelse(sum(c(PFscores_CPF_flag,CPF_comment_flag,CPF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(ER40_D.ER40D_RTCR,ER40_D.ER40D_CR),.fns = ~ ifelse(sum(c(PFscores_ER40_flag,ER40_comment_flag,ER40_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MEDF36_A.MEDF36A_RTCR,MEDF36_A.MEDF36A_CR),.fns = ~ ifelse(sum(c(PFscores_MEDF_flag,MEDF_comment_flag,MEDF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MPRACT.MP2RTCR),.fns = ~ ifelse(sum(c(MPRACT_comment_flag,MPRACT_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(PCET_A.PCET_RTCR,PCET_A.PCET_ACC2),.fns = ~ ifelse(sum(c(PCET_comment_flag,PCET_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(PMAT24_A.PMAT24_A_RTCR,PMAT24_A.PMAT24_A_CR),.fns = ~ ifelse(sum(c(PFscores_PMAT_flag,PMAT_comment_flag,PMAT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(SVOLT_A.SVOLT_RTCR,SVOLT_A.SVOLT_CR),.fns = ~ ifelse(sum(c(PFscores_SVOLT_flag,SVOLT_comment_flag,SVOLT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(VSPLOT15.VSPLOT15_RTCR,VSPLOT15.VSPLOT15_CR),.fns = ~ ifelse(sum(c(PFscores_VSPLOT_flag,VSPLOT_comment_flag,VSPLOT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(SLNB2_90.SLNB2_MRTC,SLNB2_90.SLNB2_MCR),.fns = ~ ifelse(sum(c(SLNB_comment_flag,SLNB_AV_flag) == "F") >= 1,NA,.x))) %>%
  mutate(across(.cols = c(SPCPTN90_TPRT,SPCPTN90_TP),.fns = ~ ifelse(sum(c(SPCPTN_comment_flag,SPCPTN_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(SCTAP.SCTAP_TOT),.fns = ~ ifelse(sum(c(SCTAP_comment_flag,SCTAP_AV_flag) == "F") >= 1,NA,.x))) %>% 
  ungroup() %>% 
  #mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR),.fns = ~ log(.x))) %>% 
  select(!matches("_flag$"))

CNB_cross_max_tests <- CNB_cross_clean %>% 
  rowwise() %>% 
  mutate(Tests_completed = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SPCPTN90_TP)))) %>% 
  ungroup() %>% 
  group_by(test_sessions_v.battery) %>% 
  summarize(Number_of_tests = max(Tests_completed))

CNB_cross_tests_completed_under_10 <- CNB_cross_clean %>% 
  left_join(CNB_cross_max_tests) %>% 
  rowwise() %>% 
  mutate(Tests_completed = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SPCPTN90_TP)))) %>% 
  mutate(Tests_completed_percentage = 100*Tests_completed/Number_of_tests) %>% 
  ungroup() %>% 
  filter(Tests_completed_percentage < 10) 

Over_50_removed_qc <- CNB_cross_clean %>% 
  rowwise() %>% 
  mutate(Tests_completed_post_qc = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SPCPTN90_TP)))) %>% 
  ungroup() %>% 
  left_join(CNB_before_qc_tests_completed[,c('test_sessions.bblid','Tests_completed_pre_qc')]) %>% 
  mutate(percent_tests_dropped = 100*(Tests_completed_pre_qc-Tests_completed_post_qc)/Tests_completed_pre_qc) %>% 
  relocate(test_sessions.bblid,percent_tests_dropped) %>% 
  arrange(desc(percent_tests_dropped)) %>% 
  filter(percent_tests_dropped >= 50)
  
# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("er40","pvrt","volt","cpf","cpw","gng","mpract","pcet","pmat24","medf36","adt36","plot","tap","cpt","lnb","sctap","slnb2","SPCPTN90","SPCPTNL","svolt","vsplot15"),
                       "Test_name" = c("Penn Emotion Recognition Test","Penn Verbal Reasoning Test","Visual Object Learning Test","Penn Face Memory Test",'Penn Word Memory Test',"Go-No-Go Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test","Penn Computerized Finger Tapping Test",
                                       "Penn Continuous Performance Test","Letter-N-Back Test","Penn Computerized Finger Tapping Test","Letter-N-Back","Penn Continuous Performance Test - Numbers","Penn Continuous Performance Test - Letters","Visual Object Learning Test","Penn Line Orientation Test"))
Test_map$Prefix <- str_to_upper(Test_map$Prefix)

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Suffix" = c("_cr","_rtcr","_tot","_acc2","_tprt","_ptp","_mcr","_mrtc","_pc","_mp2rtcr","_cat","_tp"),
                         "Label" = c("Correct Responses","Median Reaction Time \n Correct Responses (ms)","Average Taps \n (Dominant and Non-dominant hand added together)","Accuracy",
                                     "Median Response Time \n True Positives (ms)","True Positive (%)","Total True Positive Responses","Median Response Time \n Correct Responses","Correct Responses (%)","Median Reaction Time \n Correct Responses (ms)","Categories Achieved","True Positive Responses"))
Metric_map$Suffix <- str_to_upper(Metric_map$Suffix)

# Run models on the cross-sectional data 

# Start by examining reaction time

theme_set(theme_minimal())
theme_update(axis.title.y = element_text(size = 10))

CNB_cross_clean <- CNB_cross_clean %>% 
  mutate(Age_centered = as.numeric(scale(test_sessions_v.age,center = T,scale = F))) %>% 
  rename(Test_Location = remote)

response_cols_RT <- CNB_cross_clean %>% 
  select(matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("TAP")) %>% 
  colnames()

RT_cross_mod_output <- list()
RT_cross_mod_output_full <- list()
ylabels_RT <- c()
plot_titles_RT <- c()
cntr <- 1
CNB_cross_clean_tmp_RT <- CNB_cross_clean %>% 
  rename(Sex = test_sessions_v.gender)



for(test in response_cols_RT){
  # Parse test to generate the correct labels for the model plots
  test_noPeriod <- str_replace_all(test,pattern = "\\.",replacement = "_")
  test_split <- str_split(test_noPeriod,pattern = "_")[[1]]
  test_prefix <- test_split[1]
  test_suffix <- paste0("_",test_split[length(test_split)])
  
  Plot_title <- Test_map %>% 
    filter(Prefix == test_prefix) %>% 
    pull(Test_name)
  
  ylabel <- Metric_map %>% 
    filter(Suffix == test_suffix) %>% 
    pull(Label)
  
  test_name <- paste(Plot_title,ylabel)
  test_name_noN <- str_replace_all(test_name,pattern = " \n","")
  test_name_no_space <- str_replace_all(test_name_noN,pattern = " ",replacement = "_")
  test_name_no_space <- str_replace_all(test_name_no_space,pattern = "\\(|\\)",replacement = "")
  test_name_no_space <- str_replace_all(test_name_no_space,pattern = "-",replacement = "_")
  test_name_no_space <- str_replace_all(test_name_no_space,pattern = "_ms$",replacement = "")

  colnames(CNB_cross_clean_tmp_RT) <- ifelse(test == colnames(CNB_cross_clean_tmp_RT),test_name_no_space,colnames(CNB_cross_clean_tmp_RT))

  f <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
  f_full <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + Sex:Age_centered + Sex:Test_Location + Test_Location:Age_centered + Sex:Test_Location:Age_centered"))

  mod <- lm(f,data = CNB_cross_clean_tmp_RT)
  mod_full <- lm(f_full,data = CNB_cross_clean_tmp_RT)
  RT_cross_mod_output[[cntr]] <- mod
  RT_cross_mod_output_full[[cntr]] <- mod_full
  ylabels_RT[cntr] <- ylabel
  plot_titles_RT[cntr] <- Plot_title
  cntr <- cntr + 1
}

# Similar analysis for Correct responses

CR_cross_mod_output <- list()
CR_cross_mod_output_full <- list()
CNB_CR_data <- list()
ylabels_CR <- c()
plot_titles_CR <- c()
cntr <- 1
CNB_cross_clean_tmp_CR <- CNB_cross_clean %>% 
  rename(Sex = test_sessions_v.gender)
response_cols_CR <- CNB_cross_clean %>% 
  select(matches("_CR$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("TP$")) %>% 
  colnames()

for(test in response_cols_CR){
  
  test_noPeriod <- str_replace_all(test,pattern = "\\.",replacement = "_")
  test_split <- str_split(test_noPeriod,pattern = "_")[[1]]
  test_prefix <- test_split[1]
  test_suffix <- paste0("_",test_split[length(test_split)])
  
  Plot_title <- Test_map %>% 
    filter(Prefix == test_prefix) %>% 
    pull(Test_name)
  
  ylabel <- Metric_map %>% 
    filter(Suffix == test_suffix) %>% 
    pull(Label)
  
  test_name <- paste(Plot_title,ylabel)
  test_name_noN <- str_replace_all(test_name,pattern = "\\(|\\)| \n","")
  test_name_no_space <- str_replace_all(test_name_noN,pattern = " |-",replacement = "_")
  test_name_no_space <- str_replace_all(test_name_no_space,pattern = "%",replacement = 'Percent')
  test_name_no_space <- str_replace_all(test_name_no_space,pattern = "_ms$",replacement = "")

  colnames(CNB_cross_clean_tmp_CR) <- ifelse(test == colnames(CNB_cross_clean_tmp_CR),test_name_no_space,colnames(CNB_cross_clean_tmp_CR))

  f <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
  f_full <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + Sex:Age_centered + Sex:Test_Location + Test_Location:Age_centered + Sex:Test_Location:Age_centered"))

  mod <- lm(f,data = CNB_cross_clean_tmp_CR)
  mod_full <- lm(f_full,data = CNB_cross_clean_tmp_CR)
  CR_cross_mod_output[[cntr]] <- mod
  CR_cross_mod_output_full[[cntr]] <- mod_full
  ylabels_CR[cntr] <- ylabel
  plot_titles_CR[cntr] <- Plot_title
  cntr <- cntr + 1
}

# Write RT and CR model output and plots to file 
tab_model(CR_cross_mod_output[[1]],CR_cross_mod_output[[2]],CR_cross_mod_output[[3]],CR_cross_mod_output[[4]],CR_cross_mod_output[[5]],CR_cross_mod_output[[6]],CR_cross_mod_output[[7]],CR_cross_mod_output[[8]],CR_cross_mod_output[[9]],CR_cross_mod_output[[10]],file = "/Users/hillmann/Projects/22q/Results/By_test/Correct_responses/Cross_sectional/With_qc/22qRemote_CR_cross_ModelOutput.html")
tab_model(RT_cross_mod_output[[1]],RT_cross_mod_output[[2]],RT_cross_mod_output[[3]],RT_cross_mod_output[[4]],RT_cross_mod_output[[5]],RT_cross_mod_output[[6]],RT_cross_mod_output[[7]],RT_cross_mod_output[[8]],RT_cross_mod_output[[9]],RT_cross_mod_output[[10]],RT_cross_mod_output[[11]],RT_cross_mod_output[[12]],file = "/Users/hillmann/Projects/22q/Results/By_test/Reaction_time/Cross_sectional/With_qc/22qRemote_RT_cross_ModelOutput.html")
tab_model(CR_cross_mod_output_full[[1]],CR_cross_mod_output_full[[2]],CR_cross_mod_output_full[[3]],CR_cross_mod_output_full[[4]],CR_cross_mod_output_full[[5]],CR_cross_mod_output_full[[6]],CR_cross_mod_output_full[[7]],CR_cross_mod_output_full[[8]],CR_cross_mod_output_full[[9]],CR_cross_mod_output_full[[10]],file = "/Users/hillmann/Projects/22q/Results/By_test/Correct_responses/Cross_sectional/With_qc/22qRemote_CR_cross_full_ModelOutput.html")
tab_model(RT_cross_mod_output_full[[1]],RT_cross_mod_output_full[[2]],RT_cross_mod_output_full[[3]],RT_cross_mod_output_full[[4]],RT_cross_mod_output_full[[5]],RT_cross_mod_output_full[[6]],RT_cross_mod_output_full[[7]],RT_cross_mod_output_full[[8]],RT_cross_mod_output_full[[9]],RT_cross_mod_output_full[[10]],RT_cross_mod_output_full[[11]],RT_cross_mod_output_full[[12]],file = "/Users/hillmann/Projects/22q/Results/By_test/Reaction_time/Cross_sectional/With_qc/22qRemote_RT_cross_full_ModelOutput.html")

# p_vals_RT_full <- c()
# p_vals_CR_full <- c()
# 
# for(i in 1:10){
#   p_vals_CR_full <- c(p_vals_CR_full,as.numeric(summary(CR_cross_mod_output_full[[i]])$coefficients[c(8,10,11,12),4]))
# }
# for(i in 1:12){
#   p_vals_RT_full <- c(p_vals_RT_full,as.numeric(summary(RT_cross_mod_output_full[[i]])$coefficients[c(8,10,11,12),4]))
# }

p_vals_RT <- c()

for(i in 1:12){
  p_vals_RT <- c(p_vals_RT,as.numeric(summary(RT_cross_mod_output[[i]])$coefficients[8,4]))
}

pdf(file = "/Users/hillmann/Projects/22q/Results/By_test/Reaction_time/Cross_sectional/With_qc/22qRemote_RT_cross_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(RT_cross_mod_output)){
  grid.arrange(visreg(RT_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Male")),visreg(RT_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/By_test/Correct_responses/Cross_sectional/With_qc/22qRemote_CR_cross_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(CR_cross_mod_output)){
  grid.arrange(visreg(CR_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Male")),visreg(CR_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/By_test/Reaction_time/Cross_sectional/With_qc/22qRemote_RT_cross_full_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(RT_cross_mod_output_full)){
  grid.arrange(visreg(RT_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Male")),visreg(RT_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/By_test/Correct_responses/Cross_sectional/With_qc/22qRemote_CR_cross_full_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(CR_cross_mod_output_full)){
  grid.arrange(visreg(CR_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Male")),visreg(CR_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Female")))
}
dev.off()

# # Match data set 
# 
# CNB_cross_clean <- CNB_cross_clean %>% 
#   mutate(remote = case_when(Test_Location == "Remote" ~ 1,Test_Location == "In-person" ~ 0,TRUE ~ NA_real_))
# mod <- matchit(remote ~ test_sessions_v.age + test_sessions_v.gender,data = CNB_cross_clean,ratio = 1)
# remote_indx <- as.numeric(rownames(mod$match.matrix))
# inPerson_indx <- as.numeric(mod$match.matrix[,1])
# 
# CNB_cross_clean_matched <- CNB_cross_clean %>% 
#   relocate(Test_Location) %>% 
#   slice(c(remote_indx,inPerson_indx)) 
# 
# # Run models with matched data, starting with correct responses
# 
# CR_cross_matched_mod_output <- list()
# CNB_cross_clean_matched_tmp_CR <- CNB_cross_clean_matched
# cntr <- 1
# 
# for(test in response_cols_CR){
#   test_noPeriod <- str_replace_all(test,pattern = "\\.",replacement = "_")
#   test_split <- str_split(test_noPeriod,pattern = "_")[[1]]
#   test_prefix <- test_split[1]
#   test_suffix <- paste0("_",test_split[length(test_split)])
#   
#   Plot_title <- Test_map %>% 
#     filter(Prefix == test_prefix) %>% 
#     pull(Test_name)
#   
#   ylabel <- Metric_map %>% 
#     filter(Suffix == test_suffix) %>% 
#     pull(Label)
#   
#   test_name <- paste(Plot_title,ylabel)
#   test_name_noN <- str_replace_all(test_name,pattern = "\\(|\\)| \n","")
#   test_name_no_space <- str_replace_all(test_name_noN,pattern = " |-",replacement = "_")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "%",replacement = 'Percent')
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "_ms$",replacement = "")
#   
#   colnames(CNB_cross_clean_matched_tmp_CR) <- ifelse(test == colnames(CNB_cross_clean_matched_tmp_CR),test_name_no_space,colnames(CNB_cross_clean_matched_tmp_CR))
#   
#   #age terms test_sessions_v.age + I(test_sessions_v.age^2) + I(test_sessions_v.age^3) +
#   f <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Test_Location"))
#   mod <- lm(f,data = CNB_cross_clean_matched_tmp_CR)
#   CR_cross_matched_mod_output[[cntr]] <- mod
#   cntr <- cntr + 1
# }
#   
# tab_model(CR_cross_matched_mod_output[[1]],CR_cross_matched_mod_output[[2]],CR_cross_matched_mod_output[[3]],CR_cross_matched_mod_output[[4]],CR_cross_matched_mod_output[[5]],CR_cross_matched_mod_output[[6]],CR_cross_matched_mod_output[[7]],CR_cross_matched_mod_output[[8]],CR_cross_matched_mod_output[[9]],CR_cross_matched_mod_output[[10]],CR_cross_matched_mod_output[[11]])
# 
# # p_vals_RT_match <- c()
# # for(i in 1:11){
# #   p_vals_RT_match[i] <- summary(RT_cross_matched_mod_output[[i]])$coefficients[4,4]
# # }
# 
# # Reaction Time models with matched data
# 
# RT_cross_matched_mod_output <- list()
# CNB_cross_clean_matched_tmp_RT <- CNB_cross_clean_matched
# cntr <- 1
# 
# for(test in response_cols_RT){
#   test_noPeriod <- str_replace_all(test,pattern = "\\.",replacement = "_")
#   test_split <- str_split(test_noPeriod,pattern = "_")[[1]]
#   test_prefix <- test_split[1]
#   test_suffix <- paste0("_",test_split[length(test_split)])
#   
#   Plot_title <- Test_map %>% 
#     filter(Prefix == test_prefix) %>% 
#     pull(Test_name)
#   
#   ylabel <- Metric_map %>% 
#     filter(Suffix == test_suffix) %>% 
#     pull(Label)
#   
#   test_name <- paste(Plot_title,ylabel)
#   test_name_noN <- str_replace_all(test_name,pattern = "\\(|\\)| \n","")
#   test_name_no_space <- str_replace_all(test_name_noN,pattern = " |-",replacement = "_")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "%",replacement = 'Percent')
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "_ms$",replacement = "")
#   
#   colnames(CNB_cross_clean_matched_tmp_RT) <- ifelse(test == colnames(CNB_cross_clean_matched_tmp_RT),test_name_no_space,colnames(CNB_cross_clean_matched_tmp_RT))
#   #age terms test_sessions_v.age + I(test_sessions_v.age^2) + I(test_sessions_v.age^3) +
#   f <- as.formula(paste(test_name_no_space,"~","test_num + I(test_num^2) + Test_Location"))
#   mod <- lm(f,data = CNB_cross_clean_matched_tmp_RT)
#   RT_cross_matched_mod_output[[cntr]] <- mod
#   cntr <- cntr + 1
# }
# 
# tab_model(RT_cross_matched_mod_output[[1]],RT_cross_matched_mod_output[[2]],RT_cross_matched_mod_output[[3]],RT_cross_matched_mod_output[[4]],RT_cross_matched_mod_output[[5]],RT_cross_matched_mod_output[[6]],RT_cross_matched_mod_output[[7]],RT_cross_matched_mod_output[[8]],RT_cross_matched_mod_output[[9]],RT_cross_matched_mod_output[[10]],RT_cross_matched_mod_output[[11]])
# 
# 
# 
