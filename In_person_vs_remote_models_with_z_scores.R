library(visreg)
library(table1)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(missForest)
library(lmerTest)
library(VGAM)
library(psych)
library(ggpubr)

CNB <- read_csv("~/Projects/22q_Longitudinal/Data/22qlong_dx_sips_rawcnb_stdcnb_merged.csv")

# Get cross-sectional data set

CNB_under35 <- CNB %>% 
  mutate(Test_Location = case_when(platform == "webcnp" ~ "In-person",platform == "webcnp-surveys" ~ "Remote",TRUE ~ NA_character_)) %>% 
  select(bblid,test_sessions.datasetid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.dotest,Test_Location,matches("_valid$"),matches("_asr$")) %>% 
  select(!(matches("lan..._asr|^vmem.*"))) %>% 
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
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ ifelse(plot_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ ifelse(er40_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ ifelse(medf_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ ifelse(adt_valid %in% c("F","N"),NA,.x))) %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ as.numeric(scale(.x)))) %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ case_when(.x > 6 ~ 6,.x < -6 ~ -6,TRUE ~ .x))) 

# CNB_for_table <- CNB_cross_clean %>% 
#   rename(`Test Location` = "Test_Location",Age = test_sessions_v.age,Sex = test_sessions_v.gender,`# of CNB tests` = test_num) %>% 
#   mutate(Sex = case_when(Sex == "M" ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>% 
#   mutate(`# of CNB tests` = factor(`# of CNB tests`)) %>% 
#   mutate(`# of CNB tests` = fct_collapse(`# of CNB tests`,`4+` = c("4","5","6","7")))
# 
# table1(~ `Test Location` + Age + Sex + `# of CNB tests`, data = CNB_for_table)

CNB_iq <- missForest(as.data.frame(CNB_cross_clean %>% select(test_num,matches("_asr"))))$ximp %>% 
  tibble() %>% 
  cbind(CNB_cross_clean[,c("bblid","Test_Location","test_sessions_v.age")]) %>% 
  relocate(bblid) %>% 
  relocate(Test_Location,.after = bblid) %>% 
  relocate(test_num,.after = Test_Location) %>% 
  rowwise() %>% 
  mutate(CNB_IQ = mean(c_across(cols = matches("az_asr")))) %>%
  ungroup() %>% 
  mutate(CNB_IQ = as.numeric(scale(CNB_IQ))) %>%
  rename("#_of_CNB_tests" = test_num)


# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model on the CNB IQ model 

iq_mod <- lm(CNB_IQ ~ `#_of_CNB_tests` + I(`#_of_CNB_tests`^2) + Test_Location,data = CNB_iq)
  
# Run mixed-effect models 

acc_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_name = str_replace_all(Test_name,pattern = " ",replacement = "_")) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Accuracy") %>% 
  rename(Z_score_accuracy = Z_score) %>% 
  rename("#_of_CNB_tests" = test_num) 

speed_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_name = str_replace_all(Test_name,pattern = " ",replacement = "_")) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Speed") %>% 
  rename(Z_score_speed = Z_score) %>% 
  rename("#_of_CNB_tests" = test_num)

acc_model_no_interaction <- lmer(Z_score_accuracy ~ `#_of_CNB_tests` + I(`#_of_CNB_tests`^2) + Test_name + Test_Location + (1|bblid),data = acc_long)
acc_model_interaction <- lmer(Z_score_accuracy ~  `#_of_CNB_tests` + I(`#_of_CNB_tests`^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = acc_long)

speed_model_no_interaction <- lmer(Z_score_speed ~ `#_of_CNB_tests` + I(`#_of_CNB_tests`^2) + Test_name + Test_Location + (1|bblid),data = speed_long)
speed_model_interaction <- lmer(Z_score_speed ~ `#_of_CNB_tests` + I(`#_of_CNB_tests`^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = speed_long)


# Find basic boxplots comparing in-person to remote 

# theme_set(theme_minimal())
# theme_update(axis.title.x = element_text(size = 16,vjust = -.5),axis.title.y = element_text(size = 16))

# acc_long %>%
#   mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
#   mutate(Test = str_to_upper(Test)) %>% 
#   ggplot(aes(x = Test,y = Z_score_accuracy,color = Test_Location)) + geom_boxplot() + stat_compare_means(aes(label = paste("Wilcoxon p-value: \n",..p.format..)),show.legend = F) + labs(color = "Test Location",y = "Accuracy (Z-score)")
# 
# speed_long %>%
#   mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
#   mutate(Test = str_to_upper(Test)) %>% 
#   ggplot(aes(x = Test,y = Z_score_speed,color = Test_Location)) + geom_boxplot() + stat_compare_means(aes(label = paste("Wilcoxon p-value: \n",..p.format..)),show.legend = F)  + labs(color = "Test Location",y = "Speed (Z-score)")
# 
# CNB_iq %>% 
#   ggplot(aes(x = Test_Location, y = CNB_IQ,color = Test_Location)) + ggdist::stat_halfeye(aes(fill = Test_Location),adjust = .5, width = .6, .width = 0, justification = -.3) + 
#   geom_boxplot(width = .25, outlier.shape = NA) +
#   geom_point(size = 1.3,alpha = .3,position = position_jitter(seed = 1, width = .1)) + guides(color = F,fill = F) + labs(x = "Test Location",y = "CNB IQ") 

