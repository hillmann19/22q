library(visreg)
library(table1)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(missForest)
library(lmerTest)

CNB <- read_csv("~/Projects/22q_Longitudinal/Data/22qlong_dx_sips_rawcnb_stdcnb_merged.csv")

# Get cross-sectional data set

CNB_under35 <- CNB %>% 
  mutate(Test_Location = case_when(platform == "webcnp" ~ "In-person",platform == "webcnp-surveys" ~ "Remote",TRUE ~ NA_character_)) %>% 
  select(bblid,test_sessions.datasetid,test_sessions_v.age,test_sessions_v.dotest,Test_Location,gaf_c,matches("_valid$"),matches("_asr$")) %>% 
  filter(test_sessions_v.age <= 35) %>%
  filter(if_any(.cols = matches("_asr"),.fns = ~ !is.na(.x))) %>% 
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "^([[:digit:]])/",replacement = "0\\1/")) %>% # pad months with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "/([[:digit:]])/",replacement = "/0\\1/")) %>% #pad days with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% # Add 20 to year (15 becomes 2015)
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  group_by(bblid) %>% 
  arrange(bblid,test_sessions_v.dotest) %>% 
  mutate(test_num = row_number()) %>% 
  ungroup() 


# First, create cross-sectional data set by taking only the last test from repeat subjects

CNB_repeats_list <- CNB_under35 %>% 
  group_by(bblid) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  group_split(bblid)

find_last_test <- function(df_by_bblid){
  if(any(df_by_bblid$Test_Location == "Remote")){
    last_test <- df_by_bblid %>% 
      arrange(test_sessions_v.dotest) %>% 
      filter(Test_Location == "Remote") %>% 
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
  group_by(bblid) %>%
  filter(n() == 1) %>% 
  arrange(test_sessions_v.dotest) %>% 
  ungroup() %>% 
  bind_rows(Last_tests) 


# Turn poor quality data to NA

CNB_cross_clean <- CNB_cross %>% 
  mutate(across(.cols = c(abf_az_asr,abf_sz_asr),.fns = ~ ifelse(pcet_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(att_az_asr,att_sz_asr),.fns = ~ ifelse(cpt_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(wm_az_asr,wm_sz_asr),.fns = ~ ifelse(lnb_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(fmem_az_asr,fmem_sz_asr),.fns = ~ ifelse(cpf_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(smem_az_asr,smem_sz_asr),.fns = ~ ifelse(volt_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(nvr_az_asr,nvr_sz_asr),.fns = ~ ifelse(pmat_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ ifelse(volt_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ ifelse(er40_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ ifelse(medf_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ ifelse(adt_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ case_when(.x < -6 ~ -6,.x > 6 ~ 6,TRUE ~ .x))) %>% 
  mutate(gaf_c = ifelse(gaf_c == 999,NA,gaf_c))
  
  
CNB_iq <- missForest(as.data.frame(CNB_cross_clean %>% select(matches("_asr"))))$ximp %>% 
  tibble() %>% 
  cbind(CNB_cross_clean[,c("bblid","Test_Location","test_num","gaf_c")]) %>% 
  relocate(bblid) %>% 
  relocate(Test_Location,.after = bblid) %>% 
  relocate(test_num,.after = Test_Location) %>% 
  rowwise() %>% 
  mutate(cnb_iq = mean(c_across(cols = matches("az_asr")))) %>% 
  ungroup() %>% 
  mutate(cnb_iq = as.numeric(scale(cnb_iq)))


# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model on the CNB IQ model 

iq_mod <- lm(cnb_iq ~ gaf_c + test_num + I(test_num^2) + Test_Location,data = CNB_iq)
summary(iq_mod)

# Run mixed-effect models 

acc_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_split(Test,pattern = "_")[[1]][2]) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_name = str_replace_all(Test_name,pattern = " ",replacement = "_")) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Accuracy")

speed_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_split(Test,pattern = "_")[[1]][2]) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_name = str_replace_all(Test_name,pattern = " ",replacement = "_")) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Speed")

acc_model_no_interaction <- lmer(Z_score ~ gaf_c + test_num + I(test_num^2) + Test + Test_Location + (1|bblid),data = acc_long)
acc_model_interaction <- lmer(Z_score ~ gaf_c + test_num + I(test_num^2) + Test + Test_Location + Test:Test_Location + (1|bblid),data = acc_long)

summary(acc_model_no_interaction)

speed_model_no_interaction <- lmer(Z_score ~ gaf_c + test_num + I(test_num^2) + Test + Test_Location + (1|bblid),data = speed_long)
speed_model_interaction <- lmer(Z_score ~ gaf_c + test_num + I(test_num^2) + Test + Test_Location + Test:Test_Location + (1|bblid),data = speed_long)

summary(speed_model_no_interaction)
summary(speed_model_interaction)


# # Running plots on data that hasn't been age and sex regressed 
# 
# CNB_raw_data <- CNB_cross_clean %>% 
#   left_join(CNB[,c(1,2,213,219,227:282)])
# 
# response_cols_CR <- CNB_raw_data %>% 
#   select(matches("_CR$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("TP$"),matches("_pc$")) %>% 
#   select(-pvrt_cr,-cpw_cr,-pvrt_pc,-cpt_ptp) %>% 
#   colnames()
# 
# CR_cross_mod_output <- list()
# 
# CNB_cross_clean_tmp_CR <- CNB_raw_data %>% 
#   rename(Sex = test_sessions_v.gender) %>% 
#   mutate(Age_centered = as.numeric(scale(test_sessions_v.age,scale = F)))
# 
# cntr <- 1
# for(test in response_cols_CR){
#   f <- as.formula(paste(test,"~","test_num + I(test_num^2) + Sex + Age_centered + I(Age_centered^2) + I(Age_centered^3) + Test_Location"))
#   
#   mod <- lm(f,data = CNB_cross_clean_tmp_CR)
#   CR_cross_mod_output[[cntr]] <- mod
#   cntr <- cntr + 1
# }
# 
# tab_model(CR_cross_mod_output[[1]],CR_cross_mod_output[[2]],CR_cross_mod_output[[3]],CR_cross_mod_output[[4]],CR_cross_mod_output[[5]],CR_cross_mod_output[[6]],CR_cross_mod_output[[7]],CR_cross_mod_output[[8]],CR_cross_mod_output[[9]])
# 
# plot(CR_cross_mod_output[[1]])




