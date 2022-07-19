# Import the necessary functions and data
library(tidyverse)
library(missForest)
library(parallel)

CNB_cross_qc <- read_csv("/Users/hillmann/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")
codebook <- read_csv("~//Projects/22q/Data/Summary/bbl_cross_battery_codebook.csv")

# Turn poor quality data to NA

CNB_cross_clean <- CNB_cross_qc %>% 
  select(test_sessions.bblid,test_sessions_v.gender,test_sessions_v.age,remote,test_num,matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("_CR$"),matches("SCTAP_TOT$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("CAT$"),matches("TP$"),matches("_flag$")) %>% 
  select(-CPF_B.CPF_W_RTCR,-SPCPTN90.SCPN90_TPRT,-SPCPTNL.SCPN_TPRT,-SPCPTN90.SCPN90_TP,-SPCPTNL.SCPN_TP) %>% 
  rowwise() %>% 
  mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,ADT36_A.ADT36A_CR),.fns = ~ ifelse(sum(c(PFscores_ADT_flag,ADT_comment_flag,ADT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(CPF_B.CPF_RTCR,CPF_B.CPF_CR),.fns = ~ ifelse(sum(c(PFscores_CPF_flag,CPF_comment_flag,CPF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(ER40_D.ER40D_RTCR,ER40_D.ER40D_CR),.fns = ~ ifelse(sum(c(PFscores_ER40_flag,ER40_comment_flag,ER40_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MEDF36_A.MEDF36A_RTCR,MEDF36_A.MEDF36A_CR),.fns = ~ ifelse(sum(c(PFscores_MEDF_flag,MEDF_comment_flag,MEDF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MPRACT.MP2RTCR),.fns = ~ ifelse(sum(c(MPRACT_comment_flag,MPRACT_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(PCET_A.PCET_RTCR,PCET_A.PCET_ACC2,PCET_A.PCET_CAT),.fns = ~ ifelse(sum(c(PCET_comment_flag,PCET_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(PMAT24_A.PMAT24_A_RTCR,PMAT24_A.PMAT24_A_CR),.fns = ~ ifelse(sum(c(PFscores_PMAT_flag,PMAT_comment_flag,PMAT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(SVOLT_A.SVOLT_CR,SVOLT_A.SVOLT_RTCR),.fns = ~ ifelse(sum(c(PFscores_SVOLT_flag,SVOLT_comment_flag,SVOLT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(VSPLOT15.VSPLOT15_RTCR,VSPLOT15.VSPLOT15_CR),.fns = ~ ifelse(sum(c(PFscores_VSPLOT_flag,VSPLOT_comment_flag,VSPLOT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(SLNB2_90.SLNB2_MRTC,SLNB2_90.SLNB2_MCR),.fns = ~ ifelse(sum(c(SLNB_comment_flag,SLNB_AV_flag) == "F") >= 1,NA,.x))) %>%
  mutate(across(.cols = c(SPCPTN90_TPRT,SPCPTN90_TP),.fns = ~ ifelse(sum(c(SPCPTN_comment_flag,SPCPTN_AV_flag) == "F") >= 1,NA,.x))) %>% 
  mutate(across(.cols = c(SCTAP.SCTAP_TOT),.fns = ~ ifelse(sum(c(SCTAP_comment_flag,SCTAP_AV_flag) == "F") >= 1,NA,.x))) %>% 
  ungroup() %>% 
  select(!matches("_flag$"))

# Create composite scores for each neurocognitive domain

CNB_scaled <- CNB_cross_clean %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,remote,test_num,ADT36_A.ADT36A_CR,ADT36_A.ADT36A_RTCR,CPF_B.CPF_CR,CPF_B.CPF_RTCR,ER40_D.ER40D_CR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_CR,MEDF36_A.MEDF36A_RTCR,MPRACT.MP2RTCR,PCET_A.PCET_RTCR,PCET_A.PCET_ACC2,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_CR,PMAT24_A.PMAT24_A_RTCR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TP,SPCPTN90_TPRT,SVOLT_A.SVOLT_CR,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR,VSPLOT15.VSPLOT15_RTCR) %>% 
  mutate(across(.cols = ADT36_A.ADT36A_CR:last_col(),.fns = ~ as.numeric(scale(.x)))) %>%
  mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR),.fns = ~ -1*.x)) %>%
  mutate(test_sessions_v.gender = factor(test_sessions_v.gender))

CNB_z <- missForest(CNB_scaled %>% select(test_sessions_v.age,test_sessions_v.gender,test_num,ADT36_A.ADT36A_CR:last_col()) %>% as.data.frame())$ximp %>% 
  tibble() %>% 
  cbind(CNB_scaled[,c("test_sessions.bblid","remote")]) %>% 
  relocate(remote,.after = test_num) %>% 
  relocate(test_sessions.bblid) 
  
CNB_z <- CNB_z %>%  
  mutate(ADT_eff = (ADT36_A.ADT36A_CR+ADT36_A.ADT36A_RTCR)/2,CPF_eff = (CPF_B.CPF_CR+CPF_B.CPF_RTCR)/2,ER40_eff = (ER40_D.ER40D_CR+ER40_D.ER40D_RTCR)/2,MEDF_eff = (MEDF36_A.MEDF36A_CR+MEDF36_A.MEDF36A_RTCR)/2,PCET_eff = (PCET_A.PCET_RTCR+PCET_A.PCET_ACC2)/2,PMAT_eff = (PMAT24_A.PMAT24_A_CR+PMAT24_A.PMAT24_A_RTCR)/2,SLNB_eff = (SLNB2_90.SLNB2_MCR+SLNB2_90.SLNB2_MRTC)/2,SPCPTN_eff = (SPCPTN90_TP+SPCPTN90_TPRT)/2,SVOLT_eff = (SVOLT_A.SVOLT_CR+SVOLT_A.SVOLT_RTCR)/2,VSPLOT_eff = (VSPLOT15.VSPLOT15_CR+VSPLOT15.VSPLOT15_RTCR)/2) %>% 
  rowwise() %>% 
  mutate(overall_acc_z = mean(c_across(cols = c(ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SLNB2_90.SLNB2_MCR,SPCPTN90_TP,SVOLT_A.SVOLT_CR,VSPLOT15.VSPLOT15_CR)),na.rm = T)) %>% 
  mutate(overall_speed_z = mean(c_across(cols = c(ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR)),na.rm = T)) %>% 
  mutate(overall_efficiency_z = mean(c_across(cols = matches("_eff$")))) %>% 
  mutate(exec_acc_z = mean(c_across(cols = c(PCET_A.PCET_ACC2,SLNB2_90.SLNB2_MCR,SPCPTN90_TP)),na.rm = T)) %>% 
  mutate(mem_acc_z = mean(c_across(cols = c(CPF_B.CPF_CR,SVOLT_A.SVOLT_CR)),na.rm = T)) %>% 
  mutate(complex_acc_z = mean(c_across(cols = c(VSPLOT15.VSPLOT15_CR,PMAT24_A.PMAT24_A_CR)),na.rm = T)) %>% 
  mutate(social_acc_z = mean(c_across(cols = c(ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,ADT36_A.ADT36A_CR)),na.rm = T)) %>% 
  mutate(exec_speed_z = mean(c_across(cols = c(PCET_A.PCET_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT)),na.rm = T)) %>% 
  mutate(mem_speed_z = mean(c_across(cols = c(CPF_B.CPF_RTCR,SVOLT_A.SVOLT_RTCR)),na.rm = T)) %>% 
  mutate(complex_speed_z = mean(c_across(cols = c(VSPLOT15.VSPLOT15_RTCR,PMAT24_A.PMAT24_A_RTCR)),na.rm = T)) %>% 
  mutate(social_speed_z = mean(c_across(cols = c(ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,ADT36_A.ADT36A_RTCR)),na.rm = T)) %>% 
  mutate(motor_speed_z = mean(c_across(cols = c(MPRACT.MP2RTCR,SCTAP.SCTAP_TOT)),na.rm = T)) %>% 
  ungroup() %>% 
  mutate(across(.cols = matches("_z$"),.fns = ~ as.numeric(scale(.x)))) %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,remote,test_num,matches("_z$")) 

# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("overall","exec","mem","complex","social","motor"),
                       "Test_name" = c("Overall","Executive Functioning","Memory","Complex Cognition","Social Cognition","Motor"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Suffix" = c("speed","acc","efficiency"),
                         "Label" = c("Speed","Accuracy","Efficiency"))
# Run models on the cross-sectional data 

CNB_z <- CNB_z %>% 
  mutate(Age_centered = as.numeric(scale(test_sessions_v.age,center = T,scale = F))) %>% 
  rename(Test_Location = remote) %>% 
  rename(Sex = test_sessions_v.gender)

response_cols <- CNB_z %>% 
  select(matches("speed|acc|efficiency")) %>% 
  colnames()

domain_cross_mod_output <- list()
domain_cross_mod_output_full <- list()
ylabels<- c()
plot_titles<- c()
cntr <- 1
CNB_z_tmp <- CNB_z

for(domain in response_cols){
  # Parse test to generate the correct labels for the model plots
  domain_split <- str_split(domain,pattern = "_")[[1]]
  domain_prefix <- domain_split[1]
  domain_suffix <- domain_split[2]
  
  Plot_title <- Test_map %>% 
    filter(Prefix == domain_prefix) %>% 
    pull(Test_name)
  
  ylabel <- Metric_map %>% 
    filter(Suffix == domain_suffix) %>% 
    pull(Label)
  
  ylabel_full <- paste(ylabel,"(Z-score)")
  
  domain_name <- paste0(Plot_title,ylabel)
  domain_name_no_space <- str_replace_all(domain_name,pattern = " ",replacement = "_")
  
  colnames(CNB_z_tmp) <- ifelse(domain == colnames(CNB_z_tmp),domain_name_no_space,colnames(CNB_z_tmp))
  
  f <- as.formula(paste(domain_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
  f_full <- as.formula(paste(domain_name_no_space,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location + Sex:Age_centered + Sex:Test_Location + Test_Location:Age_centered + Sex:Test_Location:Age_centered"))
  
  mod <- lm(f,data = CNB_z_tmp)
  mod_full <- lm(f_full,data = CNB_z_tmp)
  domain_cross_mod_output[[cntr]] <- mod
  domain_cross_mod_output_full[[cntr]] <- mod_full
  ylabels[cntr] <- ylabel_full
  plot_titles[cntr] <- Plot_title
  cntr <- cntr + 1
}

# Write RT and CR model output and plots to file 
tab_model(domain_cross_mod_output[[1]],domain_cross_mod_output[[2]],domain_cross_mod_output[[3]],domain_cross_mod_output[[4]],domain_cross_mod_output[[5]],domain_cross_mod_output[[6]],domain_cross_mod_output[[7]],domain_cross_mod_output[[8]],domain_cross_mod_output[[9]],domain_cross_mod_output[[10]],domain_cross_mod_output[[11]],domain_cross_mod_output[[12]],file = "/Users/hillmann/Projects/22q/Results/By_domain/Regression_tables/No_interaction.html")
tab_model(domain_cross_mod_output_full[[1]],domain_cross_mod_output_full[[2]],domain_cross_mod_output_full[[3]],domain_cross_mod_output_full[[4]],domain_cross_mod_output_full[[5]],domain_cross_mod_output_full[[6]],domain_cross_mod_output_full[[7]],domain_cross_mod_output_full[[8]],domain_cross_mod_output_full[[9]],domain_cross_mod_output_full[[10]],domain_cross_mod_output_full[[11]],domain_cross_mod_output_full[[12]],file = "/Users/hillmann/Projects/22q/Results/By_domain/Regression_tables/Full_model.html")


# p_vals_no_int <- c()
# p_vals_full <- c()
# for(i in 1:12){
#   p_vals_no_int[i] <- summary(domain_cross_mod_output[[i]])$coefficients[8,4]
#   p_vals_full[i] <- summary(domain_cross_mod_output_full[[i]])$coefficients[8,4]
# }


pdf(file = "/Users/hillmann/Projects/22q/Results/By_domain/Plots/No_interaction.pdf",width = 14,height = 8)
for(i in 1:length(domain_cross_mod_output)){
  grid.arrange(visreg(domain_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles[i],": Male")),visreg(domain_cross_mod_output[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles[i],": Female")))
}
dev.off()

pdf(file = "/Users/hillmann/Projects/22q/Results/By_domain/Plots/Full_Model.pdf",width = 14,height = 8)
for(i in 1:length(domain_cross_mod_output_full)){
  grid.arrange(visreg(domain_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(Sex = "Male"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles[i],": Male")),visreg(domain_cross_mod_output_full[[i]],xvar = "Age_centered",by = "Test_Location",xlab = "Age (centered at mean)",ylab = ylabels[i],cond = list(Sex = "Female"),gg = TRUE,overlay = T) + labs(color = "Test Location",fill = "Test Location",title = paste0(plot_titles[i],": Female")))
}
dev.off()


# Match data set 

# CNB_z <- CNB_z %>% 
#   mutate(remote = case_when(Test_Location == "Remote" ~ 1,Test_Location == "In-person" ~ 0,TRUE ~ NA_real_))
# mod <- matchit(remote ~ test_sessions_v.age + test_sessions_v.gender,data = CNB_z,ratio = 1)
# remote_indx <- as.numeric(rownames(mod$match.matrix))
# inPerson_indx <- as.numeric(mod$match.matrix[,1])
# 
# CNB_z_matched <- CNB_z %>% 
#   relocate(Test_Location) %>% 
#   slice(c(remote_indx,inPerson_indx)) 

# Run models with matched data

# cross_matched_mod_output <- list()
# CNB_z_matched_tmp <- CNB_z_matched
# cntr <- 1
# 
# for(domain in response_cols){
#   # Parse test to generate the correct labels for the model plots
#   domain_split <- str_split(domain,pattern = "_")[[1]]
#   domain_prefix <- domain_split[1]
#   domain_suffix <- domain_split[2]
#   
#   Plot_title <- Test_map %>% 
#     filter(Prefix == domain_prefix) %>% 
#     pull(Test_name)
#   
#   ylabel <- Metric_map %>% 
#     filter(Suffix == domain_suffix) %>% 
#     pull(Label)
#   
#   ylabel_full <- paste(ylabel,"(Z-score)")
#   
#   domain_name <- paste0(Plot_title,ylabel)
#   domain_name_no_space <- str_replace_all(domain_name,pattern = " ",replacement = "_")
#   
#   colnames(CNB_z_matched_tmp) <- ifelse(domain == colnames(CNB_z_matched_tmp),domain_name_no_space,colnames(CNB_z_matched_tmp))
#   
#   f <- as.formula(paste(domain_name_no_space,"~","test_num + I(test_num^2) + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
#   
#   mod <- lm(f,data = CNB_z_matched_tmp)
#   cross_matched_mod_output[[cntr]] <- mod
#   cntr <- cntr + 1
# }
# 
# tab_model(cross_matched_mod_output[[1]],cross_matched_mod_output[[2]],cross_matched_mod_output[[3]],cross_matched_mod_output[[4]],cross_matched_mod_output[[5]],cross_matched_mod_output[[6]],cross_matched_mod_output[[7]],cross_matched_mod_output[[8]],cross_matched_mod_output[[9]],cross_matched_mod_output[[10]],cross_matched_mod_output[[11]])
