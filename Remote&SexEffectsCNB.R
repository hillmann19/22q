# Import the necessary functions and data

library(tidyverse)
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

CNB <- read_csv("~/Projects/22q/Data/cnb_all_202109.csv")
codebook <- read_csv("~/Downloads/bbl_cross_battery_codebook.csv")

# Filter to include only 22q subjects, select only necessary columns

CNB <- CNB %>% 
  select(test_sessions.bblid,test_sessions.datasetid,test_sessions.siteid, 
         test_sessions.famid, test_sessions.subid, test_sessions_v.age, test_sessions_v.battery, 
         test_sessions_v.dob, test_sessions_v.dotest, test_sessions_v.education, test_sessions_v.feducation,
         test_sessions_v.gender, test_sessions_v.handedness, test_sessions_v.meducation, deleted_sample, cnbagegrp, 
         platform, ADT36_A.valid_code, ADT36_A.ADT36A_CR, ADT36_A.ADT36A_PC, ADT36_A.ADT36A_RTCR, CPF_B.valid_code, 
         CPF_B.CPF_CR, CPF_B.CPF_RTCR, 
         CPF_B.CPF_W_RTCR, ER40_D.valid_code, ER40_D.ER40D_CR, ER40_D.ER40D_RTCR, MEDF36_A.valid_code, 
         MEDF36_A.MEDF36A_CR, MEDF36_A.MEDF36A_RTCR, MPRACT.valid_code, MPRACT.MP2RTCR, PCET_A.valid_code, PCET_A.PCET_RTCR, 
         PCET_A.PCET_CAT, PCET_A.PCET_ACC2, PMAT24_A.valid_code, PMAT24_A.PMAT24_A_CR, 
         PMAT24_A.PMAT24_A_RTCR, SCTAP.valid_code, SCTAP.SCTAP_TOT, SLNB2_90.valid_code, 
         SLNB2_90.SLNB2_MCR, SLNB2_90.SLNB2_MRTC, SPCPTN90.valid_code, SPCPTN90.SCPN90_TP,
         SPCPTN90.SCPN90_TPRT, SPCPTNL.valid_code, SPCPTNL.SCPN_TPRT, SPCPTNL.SCPN_TP,
         SVOLT_A.SVOLT_RTCR, VSPLOT15.valid_code, VSPLOT15.VSPLOT15_CR, VSPLOT15.VSPLOT15_RTCR) %>% 
  filter(deleted_sample == 1) %>% 
  mutate(test_sessions_v.gender = case_when(test_sessions_v.gender == "F" ~ "Female",test_sessions_v.gender == "M" ~ "Male",TRUE ~ NA_character_)) %>% 
  mutate(remote = ifelse(platform == "webcnp","In-person","Remote")) %>% 
  mutate(gender_remote = paste(test_sessions_v.gender,remote))

#write_csv(CNB,file = "/Users/hillmann/Projects/22q/Data/cnb_22q_202109.csv")

# Remove timepoints where subjects are over 35; create number of tests variable

CNB_under35 <- CNB %>% 
  filter(test_sessions_v.age <= 35) %>%
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "^([[:digit:]])/",replacement = "0\\1/")) %>% # pad months with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "/([[:digit:]])/",replacement = "/0\\1/")) %>% #pad days with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% # Add 20 to year (15 becomes 2015)
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions_v.dotest) %>% 
  mutate(test_num = row_number()) %>% 
  ungroup()

# First, create cross-sectional data set by taking only the last test from repeat subjects

CNB_repeats_list <- CNB_under35 %>% 
  group_by(test_sessions.bblid) %>% 
  filter(n() > 1) %>% 
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% 
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  arrange(test_sessions_v.dotest) %>% 
  ungroup() %>% 
  group_split(test_sessions.bblid)

find_last_test <- function(df_by_bblid){
  if(any(df_by_bblid$remote == "Remote")){
   last_test <- df_by_bblid %>% 
      arrange(test_sessions_v.dotest) %>% 
      filter(remote == "Remote") %>% 
      slice_tail(n = 1)
  } else{
   last_test <- df_by_bblid %>% 
      arrange(test_sessions_v.dotest) %>% 
      slice_tail(n = 1)
  }
  return(last_test)
}
  
Last_tests <- map_dfr(CNB_repeats_list,find_last_test)

CNB_cross <- CNB_under35 %>% 
  group_by(test_sessions.bblid) %>% 
  filter(n() == 1) %>% 
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% 
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  arrange(test_sessions_v.dotest) %>% 
  ungroup() %>% 
  bind_rows(Last_tests) 

# Cap values at 6 sd 

tests <- CNB_cross %>% 
     select(!(matches("^test") | matches("valid_code") | matches("_AR$") | "remote" | "gender_remote" | "deleted_sample" | "cnbagegrp" | "platform")) %>% 
     colnames()
for(test in tests){
  CNB_cross[[test]] <- ifelse(CNB_cross[[test]] > mean(CNB_cross[[test]],na.rm = TRUE) + 6*sd(CNB_cross[[test]],na.rm = TRUE),mean(CNB_cross[[test]],na.rm = TRUE) + 6*sd(CNB_cross[[test]],na.rm = TRUE),CNB_cross[[test]])
  CNB_cross[[test]] <- ifelse(CNB_cross[[test]] < mean(CNB_cross[[test]],na.rm = TRUE) - 6*sd(CNB_cross[[test]],na.rm = TRUE),mean(CNB_cross[[test]],na.rm = TRUE) - 6*sd(CNB_cross[[test]],na.rm = TRUE),CNB_cross[[test]])
}


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

response_cols_RT <- CNB_cross %>% 
  select(matches("RTCR$"),matches("MRTC"),matches("TPRT")) %>% 
  select(!(matches("_AR$")|CPF_B.CPF_W_RTCR|SPCPTNL.SCPN_TPRT)) %>% 
  colnames()

CNB_cross <- CNB_cross %>% 
  mutate(Age_centered = as.numeric(scale(test_sessions_v.age,center = T,scale = F))) %>% 
  rename(Test_Location = remote)

RT_cross_mod_output <- list()
RT_cross_mod_output_full <- list()
ylabels_RT <- c()
plot_titles_RT <- c()
test_names_full <- c()
cntr <- 1
CNB_cross_tmp <- CNB_cross 
for(test in response_cols_RT){
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
  
  colnames(CNB_cross_tmp) <- ifelse(test == colnames(CNB_cross_tmp),test_name_no_space,colnames(CNB_cross_tmp))
  
  f <- as.formula(paste(test_name_no_space,"~","test_num + test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + test_sessions_v.gender:Age_centered + test_sessions_v.gender:Test_Location + Test_Location:Age_centered + test_sessions_v.gender:Test_Location:Age_centered"))
  f_full <- as.formula(paste(test_name_no_space,"~","test_num + test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + test_sessions_v.gender:Age_centered + test_sessions_v.gender:Test_Location + Test_Location:Age_centered + test_sessions_v.gender:Test_Location:Age_centered"))
  mod <- lm(f,data = CNB_cross_tmp)
  mod_full <- lm(f_full,data = CNB_cross_tmp)
  RT_cross_mod_output[[cntr]] <- mod
  RT_cross_mod_output_full[[cntr]] <- mod_full
  ylabels_RT[cntr] <- ylabel
  plot_titles_RT[cntr] <- Plot_title
  test_names_full[cntr] <- test_name_no_space
  cntr <- cntr + 1
}

# Similar analysis for Correct responses

response_cols_CR <- CNB_cross %>% 
  select(matches("_CR$"),matches("_TP$"),matches("_MCR"),matches("_PC$"),matches("ACC2$"),matches("_CAT$")) %>% 
  select(-SPCPTNL.SCPN_TP) %>% 
  colnames()

CR_cross_mod_output <- list()
CR_cross_mod_output_full <- list()
ylabels_CR <- c()
plot_titles_CR <- c()
test_names_full <- c()
cntr <- 1

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
  
  
  colnames(CNB_cross_tmp) <- ifelse(test == colnames(CNB_cross_tmp),test_name_no_space,colnames(CNB_cross_tmp))
  
  f <- as.formula(paste(test_name_no_space,"~","test_num + test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
  f_full <- as.formula(paste(test_name_no_space,"~","test_num + test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + test_sessions_v.gender:Age_centered + test_sessions_v.gender:Test_Location + Test_Location:Age_centered + test_sessions_v.gender:Test_Location:Age_centered"))
  mod <- lm(f,data = CNB_cross_tmp)
  mod_full <- lm(f_full,data = CNB_cross_tmp)
  CR_cross_mod_output[[cntr]] <- mod
  CR_cross_mod_output_full[[cntr]] <- mod_full
  ylabels_CR[cntr] <- ylabel
  plot_titles_CR[cntr] <- Plot_title
  test_names_full[cntr] <- test_name_no_space
  cntr <- cntr + 1
}

# Write RT and CR model output and plots to file 
tab_model(CR_cross_mod_output[[1]],CR_cross_mod_output[[2]],CR_cross_mod_output[[3]],CR_cross_mod_output[[4]],CR_cross_mod_output[[5]],CR_cross_mod_output[[6]],CR_cross_mod_output[[7]],CR_cross_mod_output[[8]],CR_cross_mod_output[[9]],CR_cross_mod_output[[10]],CR_cross_mod_output[[11]],file = "/Users/hillmann/Projects/22q/Results/Correct_responses/Cross_sectional/22qRemote_CR_cross_ModelOutput.html")
tab_model(RT_cross_mod_output[[1]],RT_cross_mod_output[[2]],RT_cross_mod_output[[3]],RT_cross_mod_output[[4]],RT_cross_mod_output[[5]],RT_cross_mod_output[[6]],RT_cross_mod_output[[7]],RT_cross_mod_output[[8]],RT_cross_mod_output[[9]],RT_cross_mod_output[[10]],RT_cross_mod_output[[11]],file = "/Users/hillmann/Projects/22q/Results/Reaction_time/Cross_sectional/22qRemote_RT_cross_ModelOutput.html")
tab_model(CR_cross_mod_output_full[[1]],CR_cross_mod_output_full[[2]],CR_cross_mod_output_full[[3]],CR_cross_mod_output_full[[4]],CR_cross_mod_output_full[[5]],CR_cross_mod_output_full[[6]],CR_cross_mod_output_full[[7]],CR_cross_mod_output_full[[8]],CR_cross_mod_output_full[[9]],CR_cross_mod_output_full[[10]],CR_cross_mod_output_full[[11]],file = "/Users/hillmann/Projects/22q/Results/Correct_responses/Cross_sectional/22qRemote_CR_cross_full_ModelOutput.html")
tab_model(RT_cross_mod_output_full[[1]],RT_cross_mod_output_full[[2]],RT_cross_mod_output_full[[3]],RT_cross_mod_output_full[[4]],RT_cross_mod_output_full[[5]],RT_cross_mod_output_full[[6]],RT_cross_mod_output_full[[7]],RT_cross_mod_output_full[[8]],RT_cross_mod_output_full[[9]],RT_cross_mod_output_full[[10]],RT_cross_mod_output_full[[11]],file = "/Users/hillmann/Projects/22q/Results/Reaction_time/Cross_sectional/22qRemote_RT_cross_full_ModelOutput.html")

pdf(file = "/Users/hillmann/Projects/22q/Results/Reaction_time/Cross_sectional/22qRemote_RT_cross_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(RT_cross_mod_output)){
  grid.arrange(visreg(RT_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(test_sessions_v.gender = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Male")),visreg(RT_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(test_sessions_v.gender = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/Correct_responses/Cross_sectional/22qRemote_CR_cross_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(CR_cross_mod_output)){
  grid.arrange(visreg(CR_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(test_sessions_v.gender = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Male")),visreg(CR_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(test_sessions_v.gender = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/Reaction_time/Cross_sectional/22qRemote_RT_cross_full_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(RT_cross_mod_output_full)){
  grid.arrange(visreg(RT_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(test_sessions_v.gender = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Male")),visreg(RT_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_RT[i],cond = list(test_sessions_v.gender = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_RT[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/Correct_responses/Cross_sectional/22qRemote_CR_cross_full_ModelPlots.pdf",width = 14,height = 8)
for(i in 1:length(CR_cross_mod_output_full)){
  grid.arrange(visreg(CR_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(test_sessions_v.gender = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Male")),visreg(CR_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels_CR[i],cond = list(test_sessions_v.gender = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles_CR[i],": Female")))
}
dev.off()


# Models for longitudinal data
# 
# # Starting with Reaction Time
# 
# RT_long_mod_output <- list()
# RT_long_dfs <- list()
# f_list <- list()
# ylabels <- c()
# plot_titles <- c()
# test_names_full <- c()
# cntr <- 1
# CNB_long <- CNB %>% 
#   group_by(test_sessions.bblid) %>% 
#   filter(sum(remote == "Remote") > 0,sum(remote == "In-person") > 0) %>% 
#   ungroup() %>% 
#   mutate(Age_centered = as.numeric(scale(test_sessions_v.age,center = T,scale = F)))
# 
# CNB_long_tmp <- CNB_long %>% 
#   rename(Test_Location = remote)
# 
# for(test in response_cols_RT){
#   CNB_long_tmp_rows <- CNB_long_tmp %>% 
#     filter(!is.na(.data[[test]])) %>% 
#     group_by(test_sessions.bblid) %>% 
#     filter(sum(Test_Location == "Remote") > 0,sum(Test_Location == "In-person") > 0) %>% 
#     ungroup() %>% 
#     nrow()
#   
#   if(CNB_long_tmp_rows < 20){
#     print(test)
#     next
#   }
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
#   test_name_noN <- str_replace_all(test_name,pattern = " \n","")
#   test_name_no_space <- str_replace_all(test_name_noN,pattern = " ",replacement = "_")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "\\(|\\)",replacement = "")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "-",replacement = "_")
#   
#   colnames(CNB_long_tmp) <- ifelse(test == colnames(CNB_long_tmp),test_name_no_space,colnames(CNB_long_tmp))
#   f <- as.formula(paste(test_name_no_space,"~","test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + test_sessions_v.gender:Test_Location + Test_Location:Age_centered + test_sessions_v.gender:Age_centered + test_sessions_v.gender:Age_centered:Test_Location + (1|test_sessions.bblid)"))
#   CNB_long_tmp_clean <- CNB_long_tmp %>% 
#     filter(!is.na(.data[[test_name_no_space]])) %>% 
#     group_by(test_sessions.bblid) %>% 
#     filter(sum(Test_Location == "Remote") > 0,sum(Test_Location == "In-person") > 0) %>% 
#     ungroup()
#   RT_long_dfs[[cntr]] <- CNB_long_tmp_clean
#   ylabels[cntr] <- ylabel
#   plot_titles[cntr] <- Plot_title
#   test_names_full[cntr] <- test_name_no_space
#   f_list[[cntr]] <- f
#   cntr <- cntr + 1
# }
# 
# df1 <- RT_long_dfs[[1]]
# df2 <- RT_long_dfs[[2]]
# df3 <- RT_long_dfs[[3]]
# df4 <- RT_long_dfs[[4]]
# df5 <- RT_long_dfs[[5]]
# df6 <- RT_long_dfs[[6]]
# df7 <- RT_long_dfs[[7]]
# df8 <- RT_long_dfs[[8]]
# RT_long_mod_output[[1]] <- lmer(f_list[[1]],data = df1)
# RT_long_mod_output[[2]] <- lmer(f_list[[2]],data = df2)
# RT_long_mod_output[[3]] <- lmer(f_list[[3]],data = df3)
# RT_long_mod_output[[4]] <- lmer(f_list[[4]],data = df4)
# RT_long_mod_output[[5]] <- lmer(f_list[[5]],data = df5)
# RT_long_mod_output[[6]] <- lmer(f_list[[6]],data = df6)
# RT_long_mod_output[[7]] <- lmer(f_list[[7]],data = df7)
# RT_long_mod_output[[8]] <- lmer(f_list[[8]],data = df8)
# 
# tab_model(RT_long_mod_output[[1]],RT_long_mod_output[[2]],RT_long_mod_output[[3]],RT_long_mod_output[[4]],RT_long_mod_output[[5]],RT_long_mod_output[[6]],RT_long_mod_output[[7]],RT_long_mod_output[[8]],file = "/Users/hillmann/Projects/22q/Results/22qRemote_RT_long_ModelOutput.html")
# 
# pdf(file = "/Users/hillmann/Projects/22q/Results/22qRemote_RT_long_ModelPlots.pdf",width = 14,height = 8)
# for(i in 1:length(RT_long_mod_output)){
#   grid.arrange(visreg(RT_long_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(test_sessions_v.gender = "Male"), gg = TRUE,overlay = TRUE,partial = F,rug = F) + 
#                  geom_point(data = RT_long_dfs[[i]][RT_long_dfs[[i]]$test_sessions_v.gender == "Male",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location")) + 
#                  geom_line(data = RT_long_dfs[[i]][RT_long_dfs[[i]]$test_sessions_v.gender == "Male",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location",group = "test_sessions.bblid"),color = "grey",alpha = .75) + 
#                  labs(color = "",fill = "",title = paste0(plot_titles[i],": Male")), 
#                visreg(RT_long_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(test_sessions_v.gender = "Female"), gg = T,overlay = T,partial = F,rug = F) + 
#                  geom_point(data = RT_long_dfs[[i]][RT_long_dfs[[i]]$test_sessions_v.gender == "Female",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location")) + 
#                  geom_line(data = RT_long_dfs[[i]][RT_long_dfs[[i]]$test_sessions_v.gender == "Female",],aes_string(x = "Age_centered",y = test_names_full[i],color ="Test_Location",group = "test_sessions.bblid"),color = "grey",alpha = .75) + 
#                  labs(color = "",fill = "",title = paste0(plot_titles[i],": Female")))
# }
# dev.off()
# 
# 
# # Same analysis for Correct Resposes 
# 
# CR_long_mod_output <- list()
# CR_long_dfs <- list()
# f_list <- list()
# ylabels <- c()
# plot_titles <- c()
# test_names_full <- c()
# cntr <- 1
# 
# for(test in response_cols_CR){
#   CNB_long_tmp_rows <- CNB_long_tmp %>% 
#     filter(!is.na(.data[[test]])) %>% 
#     group_by(test_sessions.bblid) %>% 
#     filter(sum(Test_Location == "Remote") > 0,sum(Test_Location == "In-person") > 0) %>% 
#     ungroup() %>% 
#     nrow()
#   
#   if(CNB_long_tmp_rows < 20){
#     print(test)
#     next
#   }
# 
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
#   test_name_noN <- str_replace_all(test_name,pattern = " \n","")
#   test_name_no_space <- str_replace_all(test_name_noN,pattern = " ",replacement = "_")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "\\(|\\)",replacement = "")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "-",replacement = "_")
#   test_name_no_space <- str_replace_all(test_name_no_space,pattern = "%",replacement = "Percent")
# 
#   colnames(CNB_long_tmp) <- ifelse(test == colnames(CNB_long_tmp),test_name_no_space,colnames(CNB_long_tmp))
#   CNB_long_tmp_clean <- CNB_long_tmp %>% 
#     filter(!is.na(.data[[test_name_no_space]])) %>% 
#     group_by(test_sessions.bblid) %>% 
#     filter(sum(Test_Location == "Remote") > 0,sum(Test_Location == "In-person") > 0) %>% 
#     ungroup()
#   
#   CR_long_dfs[[cntr]] <- CNB_long_tmp_clean
#   f_list[[cntr]] <- as.formula(paste(test_name_no_space,"~","test_sessions_v.gender + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + test_sessions_v.gender:Test_Location + Test_Location:Age_centered + test_sessions_v.gender:Age_centered + test_sessions_v.gender:Age_centered:Test_Location + (1|test_sessions.bblid)"))
#   ylabels[cntr] <- ylabel
#   plot_titles[cntr] <- Plot_title
#   test_names_full[cntr] <- test_name_no_space
#   cntr <- cntr + 1
# }
# 
# df1 <- CR_long_dfs[[1]]
# df2 <- CR_long_dfs[[2]]
# df3 <- CR_long_dfs[[3]]
# df4 <- CR_long_dfs[[4]]
# df5 <- CR_long_dfs[[5]]
# df6 <- CR_long_dfs[[6]]
# df7 <- CR_long_dfs[[7]]
# df8 <- CR_long_dfs[[8]]
# CR_long_mod_output[[1]] <- lmer(f_list[[1]],data = df1)
# CR_long_mod_output[[2]] <- lmer(f_list[[2]],data = df2)
# CR_long_mod_output[[3]] <- lmer(f_list[[3]],data = df3)
# CR_long_mod_output[[4]] <- lmer(f_list[[4]],data = df4)
# CR_long_mod_output[[5]] <- lmer(f_list[[5]],data = df5)
# CR_long_mod_output[[6]] <- lmer(f_list[[6]],data = df6)
# CR_long_mod_output[[7]] <- lmer(f_list[[7]],data = df7)
# CR_long_mod_output[[8]] <- lmer(f_list[[8]],data = df8)
# 
# tab_model(CR_long_mod_output[[1]],CR_long_mod_output[[2]],CR_long_mod_output[[3]],CR_long_mod_output[[4]],CR_long_mod_output[[5]],CR_long_mod_output[[6]],CR_long_mod_output[[7]],CR_long_mod_output[[8]],file = "/Users/hillmann/Projects/22q/Results/22qRemote_CR_long_ModelOutput.html")
# 
# pdf(file = "/Users/hillmann/Projects/22q/Results/22qRemote_CR_long_ModelPlots.pdf",width = 14,height = 8)
# for(i in 1:length(CR_long_mod_output)){
#   grid.arrange(visreg(CR_long_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(test_sessions_v.gender = "Male"), gg = TRUE,overlay = TRUE,partial = F,rug = F) + 
#                  geom_point(data = CR_long_dfs[[i]][CR_long_dfs[[i]]$test_sessions_v.gender == "Male",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location")) + 
#                  geom_line(data = CR_long_dfs[[i]][CR_long_dfs[[i]]$test_sessions_v.gender == "Male",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location",group = "test_sessions.bblid"),color = "grey",alpha = .75) + 
#                  labs(color = "",fill = "",title = paste0(plot_titles[i],": Male")), 
#                visreg(CR_long_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(test_sessions_v.gender = "Female"),gg = T,overlay = T,partial = F,rug = F) + 
#                  geom_point(data=CR_long_dfs[[i]][CR_long_dfs[[i]]$test_sessions_v.gender == "Female",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location")) + 
#                  geom_line(data = CR_long_dfs[[i]][CR_long_dfs[[i]]$test_sessions_v.gender == "Female",],aes_string(x = "Age_centered",y = test_names_full[i],color = "Test_Location",group = "test_sessions.bblid"),color = "grey",alpha = .75) + 
#                  labs(color = "",fill = "",title = paste0(plot_titles[i],": Female")))
# }
# dev.off()

# #Extra: 22q remote plots 
# 
# # Create plots of cognitive test performance by sex and remote
# LongitudinalPlots <- list()
# response_vars <- CNB_cross %>% 
#   select(!(contains("valid_code"))) %>% 
#   select(ADT36_A.ADT36A_CR:VSPLOT15.VSPLOT15_RTCR) %>% 
#   select(-CPF_B.CPF_W_RTCR) %>% 
#   colnames()
# 
# cntr <- 1
# theme_set(theme_minimal())
# 
# for(test in response_vars){
#   
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
#   N.df <- CNB_cross[,c(test,"gender_remote","test_sessions_v.age")] %>% 
#     filter(if_all(everything(), ~ !is.na(.))) %>% 
#     group_by(gender_remote) %>% 
#     summarize(n = n()) %>% 
#     mutate(gender_remote_N = factor(paste0(gender_remote,": ","N = ",n))) %>% 
#     arrange(gender_remote_N)
#   
#   if(nrow(N.df) == 2){
#     LongitudinalPlots[[cntr]] <- CNB_cross %>% 
#       left_join(N.df) %>% 
#       filter(!is.na(gender_remote_N)) %>% 
#       ggplot(aes_string(x = "test_sessions_v.age",y = test,color = "gender_remote_N")) + geom_point(size = .6) + geom_smooth(se = FALSE) + labs(x = "Age",y = ylabel,title = Plot_title,color = "") + scale_color_manual(values = c("#ca0020","#0571b0")) + theme(legend.position = "bottom")
#     cntr <- cntr + 1
#   } else{
#     LongitudinalPlots[[cntr]] <- CNB_cross %>% 
#       left_join(N.df) %>% 
#       filter(!is.na(gender_remote_N)) %>% 
#       ggplot(aes_string(x = "test_sessions_v.age",y = test,color = "gender_remote_N")) + geom_point(size = .6) + geom_smooth(se = FALSE) + labs(x = "Age",y = ylabel,title = Plot_title,color = "") + scale_color_manual(values = c("#ca0020","#f4a582","#0571b0","#92c5de")) + theme(legend.position = "bottom")
#     cntr <- cntr + 1
#   }
# }
# 
# # Create plots of cognitive test performance by siteid
# SiteIDPlots <- list()
# 
# response_vars <- CNB_cross %>% 
#   select(!(contains("valid_code"))) %>% 
#   select(ADT36_A.ADT36A_CR:VSPLOT15.VSPLOT15_RTCR) %>% 
#   colnames()
# 
# cntr <- 1
# theme_set(theme_minimal())
# for(test in response_vars){
#   
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
#   SiteIDPlots[[cntr]] <- CNB_cross %>% 
#     ggplot(aes(x = test_sessions.siteid,y = .data[[test]],color = test_sessions.siteid,fill = test_sessions.siteid)) + geom_point(size = 1.3,alpha = .4,position = position_jitter(seed = 1,width = .1)) + 
#     geom_boxplot(alpha = .1,width = .3,outlier.shape = NA) + labs(x = "Site",y = ylabel,title = Plot_title) + 
#     scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") + theme(legend.position = "none") + 
#     ggdist::stat_halfeye(justification = 1.25,.width = 0,point_colour = NA,side = "left",adjust = .6) 
#   cntr <- cntr + 1
# }
# 
# pdf(file = "/Users/hillmann/Projects/22q/Results/Sex&remoteEffects22qplots.pdf",width = 14,height = 8)  
# LongitudinalPlots[[1]]
# LongitudinalPlots[[2]]
# LongitudinalPlots[[3]]
# LongitudinalPlots[[4]]
# LongitudinalPlots[[5]]
# LongitudinalPlots[[6]]
# LongitudinalPlots[[7]]
# LongitudinalPlots[[8]]
# LongitudinalPlots[[9]]
# LongitudinalPlots[[10]]
# LongitudinalPlots[[11]]
# LongitudinalPlots[[12]]
# LongitudinalPlots[[13]]
# LongitudinalPlots[[14]]
# LongitudinalPlots[[15]]
# LongitudinalPlots[[16]]
# LongitudinalPlots[[17]]
# LongitudinalPlots[[18]]
# LongitudinalPlots[[19]]
# LongitudinalPlots[[20]]
# LongitudinalPlots[[21]]
# LongitudinalPlots[[22]]
# LongitudinalPlots[[23]]
# LongitudinalPlots[[24]]
# LongitudinalPlots[[25]]
# LongitudinalPlots[[26]]
# 
# SiteIDPlots[[1]]
# SiteIDPlots[[2]]
# SiteIDPlots[[3]]
# SiteIDPlots[[4]]
# SiteIDPlots[[5]]
# SiteIDPlots[[6]]
# SiteIDPlots[[7]]
# SiteIDPlots[[8]]
# SiteIDPlots[[9]]
# SiteIDPlots[[10]]
# SiteIDPlots[[11]]
# SiteIDPlots[[12]]
# SiteIDPlots[[13]]
# SiteIDPlots[[14]]
# SiteIDPlots[[15]]
# SiteIDPlots[[16]]
# SiteIDPlots[[17]]
# SiteIDPlots[[18]]
# SiteIDPlots[[19]]
# SiteIDPlots[[20]]
# SiteIDPlots[[21]]
# SiteIDPlots[[22]]
# SiteIDPlots[[23]]
# SiteIDPlots[[24]]
# SiteIDPlots[[25]]
# SiteIDPlots[[26]]
# dev.off()



