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

CNB <- read_csv("/Users/hillmann/Projects/22q_Longitudinal/Data/22qlong_dx_sips_rawcnb_stdcnb_merged_updatedvalidcodes.csv")
CNB_all <- read_csv("~/Projects/22q_Longitudinal/Data/cnb_merged_webcnp_surveys_allbbl.csv")

# Get cross-sectional data set 

CNB_under35 <- CNB %>% 
  mutate(Test_Location = case_when(platform == "webcnp" ~ "In-person",platform == "webcnp-surveys" ~ "Remote",TRUE ~ NA_character_)) %>% 
  mutate(Test_Location = factor(Test_Location,levels = c("In-person","Remote"))) %>% 
  select(bblid,test_sessions.datasetid,test_sessions_v.age,test_sessions_v.gender,race,test_sessions_v.dotest,test_sessions_v.starttime,test_sessions_v.endtime,test_sessions_v.valid_code,Test_Location,matches("_valid$"),matches("_asr$")) %>% 
  select(!(matches("lan_.._asr|^vmem.*|pvrt_valid|cpw_valid"))) %>% 
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

# Look at QC breakdown by test location
# Valid_tests <- CNB_cross %>%
#   select(bblid,abf_az_asr:sm_sz_asr) %>%
#   mutate(abf_keep = ifelse(!is.na(abf_az_asr)|!is.na(abf_sz_asr),1,0)) %>%
#   mutate(att_keep = ifelse(!is.na(att_az_asr)|!is.na(att_sz_asr),1,0)) %>%
#   mutate(wm_keep = ifelse(!is.na(wm_az_asr)|!is.na(wm_sz_asr),1,0)) %>%
#   mutate(fmem_keep = ifelse(!is.na(fmem_az_asr)|!is.na(fmem_sz_asr),1,0)) %>%
#   mutate(smem_keep = ifelse(!is.na(smem_az_asr)|!is.na(smem_sz_asr),1,0)) %>%
#   mutate(nvr_keep = ifelse(!is.na(nvr_az_asr)|!is.na(nvr_sz_asr),1,0)) %>%
#   mutate(spa_keep = ifelse(!is.na(spa_az_asr)|!is.na(spa_sz_asr),1,0)) %>%
#   mutate(eid_keep = ifelse(!is.na(eid_az_asr)|!is.na(eid_sz_asr),1,0)) %>%
#   mutate(edi_keep = ifelse(!is.na(edi_az_asr)|!is.na(edi_sz_asr),1,0)) %>%
#   mutate(adi_keep = ifelse(!is.na(adi_az_asr)|!is.na(adi_sz_asr),1,0)) %>%
#   mutate(mot_keep = ifelse(!is.na(mot_sz_asr),1,0)) %>%
#   mutate(sm_keep = ifelse(!is.na(sm_sz_asr),1,0)) %>%
#   select(bblid,matches("_keep$")) %>%
#   pivot_longer(cols = matches("_keep$"),names_to = "Test",values_to = "Keep") %>%
#   mutate(Test = str_replace_all(Test,pattern = "_keep",replacement = "")) %>%
#   mutate(Test = case_when(Test == "abf" ~ "pcet",Test == "att" ~ "cpt",Test == "wm" ~ "lnb",Test == "fmem" ~ "cpf",Test == "smem" ~ "volt",
#                                Test == "nvr" ~ "pmat",Test== "spa" ~ "plot",Test == "eid" ~ "er40",Test == "edi" ~ "medf",Test == "adi" ~ "adt",Test == "mot" ~ "tap",Test == "sm" ~ "mpraxis")) %>%
#   select(bblid,Test,Keep)
# 
# df_for_table <- CNB_cross %>%
#   select(bblid,Test_Location,matches("_valid")) %>%
#   pivot_longer(cols = matches("_valid$"),names_to = "Test",values_to = "Valid_code") %>%
#   mutate(Invalid = ifelse(Valid_code == "N","Invalid","Valid")) %>%
#   mutate(Invalid = factor(Invalid,levels = c("Valid","Invalid"))) %>%
#   mutate(Test = str_replace_all(Test,pattern = "_valid",replacement = "")) %>%
#   left_join(Valid_tests) %>%
#   filter(Keep == 1) %>%
#   rename(`Quality control` = Invalid)
# 
# table1(~ Valid_code | Test_Location,data = df_for_table)

# Turn poor quality data to NA

CNB_cross_clean <- CNB_cross %>% 
  mutate(across(.cols = c(abf_az_asr,abf_sz_asr),.fns = ~ case_when(pcet_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(att_az_asr,att_sz_asr),.fns = ~ case_when(cpt_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(wm_az_asr,wm_sz_asr),.fns = ~ case_when(lnb_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(fmem_az_asr,fmem_sz_asr),.fns = ~ case_when(cpf_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(smem_az_asr,smem_sz_asr),.fns = ~ case_when(volt_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(nvr_az_asr,nvr_sz_asr),.fns = ~ case_when(pmat_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ case_when(plot_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ case_when(er40_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ case_when(medf_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ case_when(adt_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(mot_sz_asr),.fns = ~ case_when(tap_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(sm_sz_asr),.fns = ~ case_when(mpraxis_valid == "N" ~ NA_real_,TRUE ~ .x))) %>% 
  ungroup() %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ as.numeric(scale(.x)))) %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ case_when(.x > 6 ~ 6,.x < -6 ~ -6,TRUE ~ .x)))

CNB_iq <- missForest(as.data.frame(CNB_cross_clean %>% select(test_sessions_v.age,test_num,Test_Location,matches("_asr"))))$ximp %>% 
  tibble() %>% 
  cbind(CNB_cross_clean[,c("bblid")]) %>% 
  relocate(bblid) %>% 
  relocate(Test_Location,.after = bblid) %>% 
  relocate(test_num,.after = Test_Location) %>% 
  rowwise() %>% 
  mutate(CNB_IQ = mean(c_across(cols = matches("az_asr")))) %>%
  ungroup() %>% 
  mutate(CNB_IQ = as.numeric(scale(CNB_IQ))) %>%
  mutate(`Prior CNB tests` = test_num-1) %>% 
  select(-test_num)
  
# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model on the CNB IQ model 

iq_mod <- lm(CNB_IQ ~ `Prior CNB tests` + I(`Prior CNB tests`^2) + Test_Location,data = CNB_iq)

# Run mixed-effect models 

acc_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Accuracy") %>% 
  rename(Z_score_accuracy = Z_score) %>% 
  mutate(Prior_CNB_tests = test_num - 1) %>% 
  select(-test_num)

speed_long <- CNB_cross_clean %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Speed") %>% 
  rename(Z_score_speed = Z_score) %>% 
  mutate(Prior_CNB_tests = test_num - 1) %>% 
  select(-test_num)

acc_model_no_interaction <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = acc_long)
acc_model_interaction <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = acc_long)

speed_model_no_interaction <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = speed_long)
speed_model_interaction <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = speed_long)

# Get model plots
theme_set(theme_minimal())
theme_update(text = element_text(size = 15))

acc_no_int_labels <- c("Prior CNB tests",expression(paste("Prior CNB tests" ^2)),"Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(acc_no_int_labels) <- c("Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

speed_no_int_labels <- c("Prior CNB tests",expression(paste("Prior CNB tests" ^2)),"Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(speed_no_int_labels) <- c("Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

acc_int_labels <- c(acc_no_int_labels,str_c(c("Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(acc_int_labels) <- c(names(acc_no_int_labels),str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_namePenn Conditional Exclusion Test","Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test","Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

speed_int_labels <- c(speed_no_int_labels,str_c(c("Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(speed_int_labels) <- c(names(speed_no_int_labels),str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

plot_model(acc_model_no_interaction,title = "CNB Accuracy: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = acc_no_int_labels,axis.title = "Coefficient Estimate")
plot_model(speed_model_no_interaction,title = "CNB Speed: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = speed_no_int_labels,axis.title = "Coefficient Estimate")

theme_update(axis.text.y = element_text(size = 12))

plot_model(acc_model_interaction,title = "CNB Accuracy: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = acc_int_labels,axis.title = "Coefficient Estimate")
plot_model(speed_model_interaction,title = "CNB Speed: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = speed_int_labels,axis.title = "Coefficient Estimate")

remote_effects_acc <- as.numeric(summary(acc_model_interaction)$coefficients[13:22,5])
coefficient_names_acc <- rownames(summary(acc_model_interaction)$coefficients[13:22,])

p_val_df_acc <- data.frame(Test = coefficient_names_acc,P_val = remote_effects_acc)

p_val_df_acc$p_val_corrected <- p.adjust(p_val_df_acc$P_val,method = "fdr")

remote_effects_speed <- as.numeric(summary(speed_model_interaction)$coefficients[15:26,5])
coefficient_names_speed <- rownames(summary(speed_model_interaction)$coefficients[15:26,])

p_val_df_speed <- data.frame(Test = coefficient_names_speed,P_val = remote_effects_speed)

p_val_df_speed$p_val_corrected <- p.adjust(p_val_df_speed$P_val,method = "fdr")

# Create demographic tables for presentation

CNB_for_table <- CNB_cross_clean %>%
  mutate(`Prior CNB tests` = test_num - 1) %>% 
  rename(`Test Location` = "Test_Location",Age = test_sessions_v.age,Sex = test_sessions_v.gender) %>%
  mutate(Sex = case_when(Sex == "M" ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>%
  mutate(`Prior CNB tests` = factor(`Prior CNB tests`)) %>%
  mutate(`Prior CNB tests` = fct_collapse(`Prior CNB tests`,`3+` = c("3","4","5","6"))) %>% 
  mutate(race = case_when(race == 1 ~ "White",race == 2 ~ "Black/African American",race == 3 ~ "Native American",race == 4 ~ "Asian",race == 5 ~ "More than one race",race == 6 ~ "Hawaiian/Pacific Islander",race == 9 ~ NA_character_,TRUE ~ NA_character_)) %>% 
  rename(Race = race)

table1(~ Age + Sex + Race + `Prior CNB tests`|`Test Location`, data = CNB_for_table)

#Plot for CPT

theme_set(theme_minimal())

CNB_CPT <- CNB_cross_clean %>%
  mutate(`Prior CNB tests` = test_num - 1)

# att_az_asr_residuals <- summary(lm(att_az_asr ~ `Prior CNB tests` + I(`Prior CNB tests`^2),data = CNB_CPT))$residuals
# CNB_CPT$att_az_asr <- ifelse(is.na(CNB_CPT$att_az_asr),NA, att_az_asr_residuals)

CNB_CPT %>%
  ggplot(aes(x = Test_Location, y = att_az_asr)) +
  ggdist::stat_halfeye(aes(fill = Test_Location),
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA) +
  geom_boxplot(aes(color = Test_Location),
    width = .25,
    outlier.shape = NA
  ) +
  geom_point(aes(color = Test_Location),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") + labs(x = "Test Location",y = "PCPT Accuracy") + guides(color = "none",fill = "none") + stat_compare_means(label.x = 1.55)


# RE-DO ANALYSIS USING ONLY V AND VC DATA

vsPlot <- CNB_all %>% 
  filter(!is.na(test_sessions.bblid.clean),!is.na(test_sessions_v.dotest)) %>% 
  rename(bblid = test_sessions.bblid.clean) %>% 
  select(bblid,test_sessions_v.dotest,VSPLOT15.VSPLOT15_SUM_DEG_OFF,VSPLOT24.VSPLOT15_SUM_DEG_OFF,VSPLOT24.VSPLOT15_RT,VSPLOT15.VSPLOT15_RT) %>% 
  filter(if_any(.cols = matches("VSPLOT"),.fns = ~ !is.na(.x)))

CNB_cross_v <- CNB_cross %>% 
  left_join(vsPlot) %>% 
  mutate(VSPLOT_SUM_DEG_OFF = case_when(is.na(VSPLOT15.VSPLOT15_SUM_DEG_OFF) ~ VSPLOT24.VSPLOT15_SUM_DEG_OFF,TRUE ~ VSPLOT15.VSPLOT15_SUM_DEG_OFF)) %>% 
  mutate(VSPLOT_RT = case_when(is.na(VSPLOT24.VSPLOT15_RT) ~ VSPLOT15.VSPLOT15_RT,TRUE ~ VSPLOT24.VSPLOT15_RT)) %>% 
  mutate(plot_valid = case_when(Test_Location == "Remote" ~ case_when(VSPLOT_SUM_DEG_OFF > 500|VSPLOT_RT < 1000 ~ "0",TRUE ~ "V"),TRUE ~ plot_valid)) %>% 
  mutate(across(.cols = c(abf_az_asr,abf_sz_asr),.fns = ~ case_when(pcet_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(att_az_asr,att_sz_asr),.fns = ~ case_when(cpt_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(wm_az_asr,wm_sz_asr),.fns = ~ case_when(lnb_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(fmem_az_asr,fmem_sz_asr),.fns = ~ case_when(cpf_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(smem_az_asr,smem_sz_asr),.fns = ~ case_when(volt_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(nvr_az_asr,nvr_sz_asr),.fns = ~ case_when(pmat_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ case_when(plot_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ case_when(er40_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ case_when(medf_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ case_when(adt_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(mot_sz_asr),.fns = ~ case_when(tap_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  mutate(across(.cols = c(sm_sz_asr),.fns = ~ case_when(mpraxis_valid %in% c("V","VC") ~ .x,TRUE ~ NA_real_))) %>% 
  ungroup() %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ as.numeric(scale(.x)))) %>% 
  mutate(across(.cols = matches("_asr$"),.fns = ~ case_when(.x > 6 ~ 6,.x < -6 ~ -6,TRUE ~ .x)))

CNB_iq <- missForest(as.data.frame(CNB_cross_v %>% select(test_sessions_v.age,test_num,Test_Location,matches("_asr"))))$ximp %>% 
  tibble() %>% 
  cbind(CNB_cross_v[,c("bblid")]) %>% 
  relocate(bblid) %>% 
  relocate(Test_Location,.after = bblid) %>% 
  relocate(test_num,.after = Test_Location) %>% 
  rowwise() %>% 
  mutate(CNB_IQ = mean(c_across(cols = matches("az_asr")))) %>%
  ungroup() %>% 
  mutate(CNB_IQ = as.numeric(scale(CNB_IQ))) %>%
  mutate(`Prior CNB tests` = test_num-1) %>% 
  select(-test_num)

# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model on the CNB IQ model 

iq_mod <- lm(CNB_IQ ~ `Prior CNB tests` + I(`Prior CNB tests`^2) + Test_Location,data = CNB_iq)

# Run mixed-effect models 

acc_long_v <- CNB_cross_v %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Accuracy") %>% 
  rename(Z_score_accuracy = Z_score) %>% 
  mutate(Prior_CNB_tests = test_num - 1) %>% 
  select(-test_num)

speed_long_v <- CNB_cross_v %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Speed") %>% 
  rename(Z_score_speed = Z_score) %>% 
  mutate(Prior_CNB_tests = test_num - 1) %>% 
  select(-test_num)

acc_model_no_interaction_v <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = acc_long_v)
acc_model_interaction_v <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = acc_long_v)

speed_model_no_interaction_v <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = speed_long_v)
speed_model_interaction_v <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = speed_long_v)

# Get model plots
theme_set(theme_minimal())
theme_update(text = element_text(size = 15))

plot_model(acc_model_no_interaction_v,title = "CNB Accuracy: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = acc_no_int_labels,axis.title = "Coefficient Estimate")
plot_model(speed_model_no_interaction_v,title = "CNB Speed: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = speed_no_int_labels,axis.title = "Coefficient Estimate")

theme_update(axis.text.y = element_text(size = 12))

plot_model(acc_model_interaction_v,title = "CNB Accuracy: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = acc_int_labels,axis.title = "Coefficient Estimate")
plot_model(speed_model_interaction_v,title = "CNB Speed: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = speed_int_labels,axis.title = "Coefficient Estimate")

remote_effects_acc <- as.numeric(summary(acc_model_interaction)$coefficients[13:22,5])
coefficient_names_acc <- rownames(summary(acc_model_interaction)$coefficients[13:22,])

p_val_df_acc <- data.frame(Test = coefficient_names_acc,P_val = remote_effects_acc)

p_val_df_acc$p_val_corrected <- p.adjust(p_val_df_acc$P_val,method = "fdr")

remote_effects_speed<- as.numeric(summary(speed_model_interaction)$coefficients[15:26,5])
coefficient_names_speed <- rownames(summary(speed_model_interaction)$coefficients[15:26,])

p_val_df_speed <- data.frame(Test = coefficient_names_speed,P_val = remote_effects_speed)

p_val_df_speed$p_val_corrected <- p.adjust(p_val_df_speed$P_val,method = "fdr")

#Plot for CPT -- V and VC only

theme_set(theme_minimal())

CNB_CPT_v <- CNB_cross_v %>%
  mutate(`Prior CNB tests` = test_num - 1)

# att_az_asr_residuals <- summary(lm(att_az_asr ~ `Prior CNB tests` + I(`Prior CNB tests`^2),data = CNB_CPT))$residuals
# CNB_CPT$att_az_asr <- ifelse(is.na(CNB_CPT$att_az_asr),NA, att_az_asr_residuals)

CNB_CPT_v %>%
  ggplot(aes(x = Test_Location, y = att_az_asr)) +
  ggdist::stat_halfeye(aes(fill = Test_Location),
                       adjust = .5,
                       width = .6,
                       .width = 0,
                       justification = -.3,
                       point_colour = NA) +
  geom_boxplot(aes(color = Test_Location),
               width = .25,
               outlier.shape = NA
  ) +
  geom_point(aes(color = Test_Location),
             size = 1.3,
             alpha = .3,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") + labs(x = "Test Location",y = "PCPT Accuracy") + guides(color = "none",fill = "none") + stat_compare_means(label.x = 1.55)
