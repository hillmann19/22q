library(tidyverse)

CNB_for_Study <- read_csv("~/Projects/22q/Data/Summary/cnb_all_202109.csv")
athena_3360_2096 <- read_csv("~/Projects/22q/Data/itemwise/athena_3360_2096.csv")
athena_254_360 <- read_csv("~/Projects/22q/Data/itemwise/athena_254_360.csv")

# Filter to include only 22q subjects, select only necessary columns


CNB_22q <- CNB_for_Study %>% 
  select(test_sessions.bblid,test_sessions.datasetid,test_sessions.siteid, 
         test_sessions.famid, test_sessions.subid, test_sessions_v.age, test_sessions_v.battery, 
         test_sessions_v.dob, test_sessions_v.dotest, test_sessions_v.education, test_sessions_v.feducation,
         test_sessions_v.gender, test_sessions_v.handedness, test_sessions_v.meducation, deleted_sample, cnbagegrp, 
         platform,ADT36_A.ADT36A_CR, ADT36_A.ADT36A_PC, ADT36_A.ADT36A_RTCR,ADT36_A.ADT36A_LRR_750,ADT36_A.ADT36A_CNT_200,ADT36_A.ADT36A_LRSR,CPF_B.CPF_CNT_200,CPF_B.CPF_LRSR, 
         CPF_B.CPF_CR, CPF_B.CPF_RT,CPF_B.CPF_RTCR,CPF_B.CPF_W_RTCR, ER40_D.ER40D_RT,ER40_D.ER40D_HAP,ER40_D.ER40D_CNT_250,ER40_D.ER40D_LRSR,ER40_D.ER40D_CR, ER40_D.ER40D_RTCR, MEDF36_A.MEDF36A_LRR_750,MEDF36_A.MEDF36A_CNT_200,MEDF36_A.MEDF36A_LRSR,
         MEDF36_A.MEDF36A_CR, MEDF36_A.MEDF36A_RTCR,MPRACT.MP2,MPRACT.MP2RTCR, PCET_A.PCET_LRSR,PCET_A.PCET_RTCR, 
         PCET_A.PCET_CAT, PCET_A.PCET_ACC2, PMAT24_A.PMAT24_A_CR, 
         PMAT24_A.PMAT24_A_RTTO,PMAT24_A.PMAT24_A_RTCR, SCTAP.SCTAP_UNF_TOT,SCTAP.SCTAP_EXCESS_TOT,SCTAP.SCTAP_UNF_STD,SCTAP.SCTAP_TOT, SLNB2_90.SLNB2_TP0,SLNB2_90.SLNB2_FP0,SLNB2_90.SLNB2_LRR,SLNB2_90.SLNB2_LRNR,SLNB2_90.SLNB2_CNT_200, 
         SLNB2_90.SLNB2_MCR, SLNB2_90.SLNB2_MRTC,SPCPTNL.SCPT_LRNR,SPCPTN90.SCPN90_TP,
         SPCPTN90.SCPN90_TPRT,SPCPTNL.SCPN_TPRT, SPCPTNL.SCPN_TP,SVOLT_A.SVOLT_CR,SVOLT_A.SVOLT_RT,SVOLT_A.SVOLT_NR,
         SVOLT_A.SVOLT_CNT_200,SVOLT_A.SVOLT_LRR_200,SVOLT_A.SVOLT_LRSR,SVOLT_A.SVOLT_RTCR,SVOLT_A.SVOLT_CR, VSPLOT15.VSPLOT15_TOT_RT,VSPLOT15.VSPLOT15_SUM_DEG_OFF,VSPLOT15.VSPLOT15_SUM_EXCESS,VSPLOT15.VSPLOT15_SUM_DEFICIT,VSPLOT15.VSPLOT15_CR, VSPLOT15.VSPLOT15_RTCR) %>% 
  mutate(SPCPTN90_TP = ifelse(is.na(SPCPTN90.SCPN90_TP),SPCPTNL.SCPN_TP,SPCPTN90.SCPN90_TP)) %>% 
  mutate(SPCPTN90_TPRT = ifelse(is.na(SPCPTN90.SCPN90_TPRT),SPCPTNL.SCPN_TPRT,SPCPTN90.SCPN90_TPRT)) %>% 
  filter(deleted_sample == 1) %>% 
  mutate(test_sessions_v.gender = case_when(test_sessions_v.gender == "F" ~ "Female",test_sessions_v.gender == "M" ~ "Male",TRUE ~ NA_character_)) %>% 
  mutate(remote = case_when(platform == "webcnp" ~ "In-person",platform == "webcnp-surveys" ~ "Remote",platform == "crowdsource-surveys" ~ "Remote",TRUE ~ NA_character_)) %>% 
  mutate(gender_remote = paste(test_sessions_v.gender,remote)) %>% 
  filter(if_any(.cols = c(matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("_CR$"),matches("SCTAP_TOT$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("CAT$"),matches("TP$")),.fns = ~ !is.na(.x))) 

# Remove timepoints where subjects are over 35; create number of tests variable

CNB_under35 <- CNB_22q %>% 
  filter(test_sessions_v.age <= 35) %>%
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "^([[:digit:]])/",replacement = "0\\1/")) %>% # pad months with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "/([[:digit:]])/",replacement = "/0\\1/")) %>% #pad days with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% # Add 20 to year (15 becomes 2015)
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions.bblid,test_sessions_v.dotest) %>% 
  mutate(test_num = row_number()) %>% 
  ungroup() 


# First, create cross-sectional data set by taking only the last test from repeat subjects

CNB_repeats_list <- CNB_under35 %>% 
  group_by(test_sessions.bblid) %>% 
  filter(n() > 1) %>% 
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
  arrange(test_sessions_v.dotest) %>% 
  ungroup() %>% 
  bind_rows(Last_tests) 

# Create SPCPT LRNR variable - separate calculations for those who received numbers and letters vs. numbers only

athena_3360_2096_resp_long_both <- athena_3360_2096 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTNL.SCPN_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTNL")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("RESP")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTNL.SCPT_QID000001_0_RESP:last_col(),names_to = "Question",values_to = "RESP") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTNL.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_RESP",replacement = "")) 

SPCPTNL_trial_ordered_both <- athena_3360_2096 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTNL.SCPN_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTNL")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("TRIAL")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTNL.SCPT_QID000001_0_TRIAL:last_col(),names_to = "Question",values_to = "Trial") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTNL.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_TRIAL",replacement = "")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions.bblid,Trial) %>% 
  ungroup() %>% 
  left_join(athena_3360_2096_resp_long_both) %>% 
  group_split(test_sessions.bblid)

athena_3360_2096_resp_long_N_only <- athena_3360_2096 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTN90.SCPN90_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTN90")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("RESP")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTN90.SCPT_QID000001_0_RESP:last_col(),names_to = "Question",values_to = "RESP") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTN90.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_RESP",replacement = "")) 

SPCPTNL_trial_ordered_N_only <- athena_3360_2096 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTN90.SCPN90_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTN90")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("TRIAL")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTN90.SCPT_QID000001_0_TRIAL:last_col(),names_to = "Question",values_to = "Trial") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTN90.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_TRIAL",replacement = "")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions.bblid,Trial) %>% 
  ungroup() %>% 
  left_join(athena_3360_2096_resp_long_N_only) %>% 
  group_split(test_sessions.bblid)

athena_254_360_resp_long_N_only <- athena_254_360 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTN90.SCPN90_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTN90")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("RESP")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTN90.SCPT_QID000001_0_RESP:last_col(),names_to = "Question",values_to = "RESP") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTN90.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_RESP",replacement = "")) 

SPCPTNL_remote_trial_ordered_N_only <- athena_254_360 %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(SPCPTN90.SCPN90_TP)) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,contains("SPCPTN90")) %>% 
  select(test_sessions.bblid,test_sessions.datasetid,matches("TRIAL")) %>% 
  semi_join(CNB_cross) %>% 
  pivot_longer(cols = SPCPTN90.SCPT_QID000001_0_TRIAL:last_col(),names_to = "Question",values_to = "Trial") %>% 
  mutate(Question = str_replace_all(Question,pattern = "SPCPTN90.SCPT_",replacement = "")) %>% 
  mutate(Question = str_replace_all(Question,pattern = "_TRIAL",replacement = "")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions.bblid,Trial) %>% 
  ungroup() %>% 
  left_join(athena_254_360_resp_long_N_only) %>% 
  group_split(test_sessions.bblid)

find_LRNR <- function(df){
  
  # Find number trials 
  df <- df %>% 
    arrange(Trial) %>% 
    slice_head(n = 90)
  
  LRNR <- max(rle(df$RESP)$lengths[rle(df$RESP)$values == 0])
  df <- df %>% 
    distinct(test_sessions.bblid,.keep_all = T) %>% 
    mutate(SCPTN_LRNR = LRNR) %>% 
    select(test_sessions.bblid,test_sessions.datasetid,SCPTN_LRNR)
  
  return(df)
}

# Find LRNR for subjects that completed both tests and N only for in-person

LRNR_athena_3360_2096_both_tests <- map_dfr(SPCPTNL_trial_ordered_both,find_LRNR)

LRNR_athena_3360_2096_N_only <- map_dfr(SPCPTNL_trial_ordered_N_only,find_LRNR)

LRNR_athena_254_360_N_only <- map_dfr(SPCPTNL_remote_trial_ordered_N_only,find_LRNR)


# Add in LRNR for athena_254_360 later when I have the data 

LRNR_all <- rbind(LRNR_athena_3360_2096_N_only,LRNR_athena_3360_2096_both_tests,LRNR_athena_254_360_N_only)

CNB_cross <- CNB_cross %>% 
  left_join(LRNR_all)

# Create flags for individuals with poor data quality 

CNB_cross_qc <- CNB_cross %>% 
  mutate(ADT_flag = case_when(ADT36_A.ADT36A_LRR_750 > 3 | ADT36_A.ADT36A_CNT_200 > 3 | ADT36_A.ADT36A_LRSR > 11 ~ "F",TRUE ~ "V")) %>% 
  mutate(CPF_flag = case_when(CPF_B.CPF_CNT_200 >= 21 | CPF_B.CPF_LRSR >= 21 ~ "F",TRUE ~ "V")) %>%
  mutate(ER40_flag = case_when(ER40_D.ER40D_CR < 12 | ER40_D.ER40D_HAP < 6 | ER40_D.ER40D_CNT_250 > 2 | ER40_D.ER40D_LRSR > 11 ~ "F",TRUE ~ "V")) %>%
  mutate(MEDF_flag = case_when(MEDF36_A.MEDF36A_LRR_750 > 3 | MEDF36_A.MEDF36A_CNT_200 > 3 | MEDF36_A.MEDF36A_LRSR > 11 ~ "F",TRUE ~ "V")) %>%
  mutate(MPRACT_flag = case_when(MPRACT.MP2 < 16 ~ "F",TRUE ~ "V")) %>%
  mutate(PCET_flag = case_when(PCET_A.PCET_LRSR > 10 ~ "F",TRUE ~ "V")) %>%
  mutate(PMAT_flag = case_when(PMAT24_A.PMAT24_A_CR < 5 && PMAT24_A.PMAT24_A_RTTO < 2000  ~ "F",TRUE ~ "V")) %>%
  mutate(SCTAP_flag = case_when(SCTAP.SCTAP_UNF_TOT < 70 | SCTAP.SCTAP_EXCESS_TOT > 38 | SCTAP.SCTAP_UNF_STD > 17  ~ "F",TRUE ~ "V")) %>%
  mutate(SLNB_flag = case_when(SLNB2_90.SLNB2_TP0 < 6 | SLNB2_90.SLNB2_FP0 > 4 | SLNB2_90.SLNB2_LRR > 8 | SLNB2_90.SLNB2_LRNR > 26 | SLNB2_90.SLNB2_CNT_200 > 2 ~ "F",TRUE ~ "V")) %>%
  mutate(SPCPTN_flag = case_when(SCPTN_LRNR > 36 ~ "F",TRUE ~ "V")) %>% 
  mutate(SVOLT_flag = case_when(SVOLT_A.SVOLT_CNT_200 >= 11 | SVOLT_A.SVOLT_LRSR >= 11 ~ "F",TRUE ~ "V")) %>% 
  mutate(VSPLOT_flag = case_when(VSPLOT15.VSPLOT15_SUM_DEG_OFF > 400 && (VSPLOT15.VSPLOT15_SUM_EXCESS > 200 | VSPLOT15.VSPLOT15_SUM_DEFICIT > 45)  ~ "F",TRUE ~ "V")) 


write_csv(CNB_cross_qc,file = "/Users/hillmann/Projects/22q/Data/Summary/cnb_22q_cross_qc_05_12_2022.csv")

