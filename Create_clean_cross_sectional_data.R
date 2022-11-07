# Read in necessary data and packages
library(tidyverse)
library(psych)
library(PerFit)
library(table1)
select <- dplyr::select
CNB <- read_csv('/Users/hillmann/Projects/22q/Data/22qlong_dx_sips_rawcnb_stdcnb_merged_updatedvalidcodes_nov2022.csv')
athena_3360_2096 <- read_csv("~/Projects/22q/Data/itemwise/athena_3360_2096.csv")
athena_3360_2107 <- read_csv("~/Projects/22q/Data/itemwise/athena_3360_2107.csv")
athena_195_360 <- read_csv("~/Projects/22q/Data/itemwise/athena_195_360.csv")
athena_195_369 <- read_csv("~/Projects/22q/Data/itemwise/athena_195_369.csv")
reviewer_comments <- read_csv("Projects/22q/Data/Assessor_flags_to_update_10_31_2022_EM.csv", 
                              col_types = cols(...17 = col_skip()), 
                              skip = 1)

reviewer_comments <- reviewer_comments %>% 
  rename(test_sessions.datasetid = test_sessions_datasetid)

# Get cross-sectional data set 

CNB_under35 <- CNB %>% 
  rename(bblid = test_sessions.bblid.clean) %>% 
  mutate(Test_Location = case_when(platform == "webcnp" ~ "In-person",platform == "webcnp-surveys" ~ "Remote",TRUE ~ NA_character_)) %>% 
  mutate(Test_Location = factor(Test_Location,levels = c("In-person","Remote"))) %>% 
  select(bblid,test_sessions.datasetid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.dotest,Test_Location,matches('_genus$'),matches("_valid$"),matches("_asr$")) %>% 
  select(!(matches("lan_.._asr|vmem_.._asr|pvrt|^cpw|^gng|^aim|^digsym"))) %>% 
  filter(test_sessions_v.age <= 35) %>%
  filter(if_any(.cols = matches("_asr"),.fns = ~ !is.na(.x))) %>% 
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
  bind_rows(Last_tests) %>% 
  mutate(bblid = as.character(bblid))

# Create VSPLOT validation flags

test_grep_15 <- "^VSPLOT15.VSPLOT15"
test_grep_24 <- "^VSPLOT24.VSPLOT15"

athena_3360_2096_qc_15 <- athena_3360_2096 %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_15)) %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches("VSPLOT.*TOT_RT$"),matches('VSPLOT.*SUM_DEG_OFF$'))

athena_195_360_qc_15 <- athena_195_360 %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_15)) %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches("VSPLOT.*TOT_RT$"),matches('VSPLOT.*SUM_DEG_OFF$'))

athena_3360_2096_qc_24 <- athena_3360_2096 %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_24)) %>%
  select(test_sessions.bblid,test_sessions.datasetid,matches("VSPLOT.*TOT_RT$"),matches('VSPLOT.*SUM_DEG_OFF$'))

colnames(athena_3360_2096_qc_15) <- str_replace_all(colnames(athena_3360_2096_qc_15),pattern = "^VSPLOT15\\.",replacement = "")
colnames(athena_195_360_qc_15) <- str_replace_all(colnames(athena_195_360_qc_15),pattern = "^VSPLOT15\\.",replacement = "")
colnames(athena_3360_2096_qc_24) <- str_replace_all(colnames(athena_3360_2096_qc_24),pattern = "^VSPLOT24\\.",replacement = "")

vsplot_qc_data <- rbind(athena_3360_2096_qc_15,athena_195_360_qc_15,athena_3360_2096_qc_24) %>%
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
  rename(bblid = test_sessions.bblid) %>%
  filter(if_any(.cols = matches("VSPLOT"),.fns = ~ !is.na(.x))) %>%
  semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid"))

CNB_cross <- CNB_cross %>%
  left_join(vsplot_qc_data) %>%
  mutate(plot_valid = case_when(plot_valid == 'N' ~ 'N',VSPLOT15_TOT_RT < 100000 | VSPLOT15_SUM_DEG_OFF > 500 ~ 'F',TRUE ~ 'V'))

add_PFscores <- function(test,test_type,output_df){
  
  if(test == "PMAT"){
  
    test_grep_a <- "^PMAT24_A.PMAT24_A_QID"
    test_grep_b <- "^PMAT24_B.PMAT24_B_QID"
    
    athena_3360_2096_narrow_a <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_a <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_369_narrow_a <- athena_195_369 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_369_narrow_b <- athena_195_369 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_3360_2107_narrow_b <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    colnames(athena_3360_2096_narrow_a) <- str_replace_all(colnames(athena_3360_2096_narrow_a),pattern = "PMAT24_A","PMAT24")
    colnames(athena_195_360_narrow_a) <- str_replace_all(colnames(athena_195_360_narrow_a),pattern = "PMAT24_A","PMAT24")
    colnames(athena_195_369_narrow_a) <- str_replace_all(colnames(athena_195_369_narrow_a),pattern = "PMAT24_A","PMAT24")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = "PMAT24_B","PMAT24")
    colnames(athena_3360_2107_narrow_b) <- str_replace_all(colnames(athena_3360_2107_narrow_b),pattern = "PMAT24_B","PMAT24")
    
    itemwise_data <- rbind(athena_3360_2096_narrow_a,athena_195_360_narrow_a,athena_195_369_narrow_b,athena_3360_2107_narrow_b) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("PMAT24"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      arrange(bblid) %>% 
      select(!matches("^PMAT24.PMAT24_QID001")) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) %>% 
      mutate(across(.cols = matches("TTR$"),.fns = ~ case_when(.x < 0 ~ NA_real_,TRUE ~ .x))) 
    
  } else if(test == "SLNB"){ 
    
    test_grep = "^SLNB2_90.SLNB2_QID"
    
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR"),matches("_TTR")) %>% 
      select(!matches("8[[:digit:]][[:digit:]][[:digit:]]$")) 
    
    colnames(athena_3360_2096_narrow) <- str_replace_all(colnames(athena_3360_2096_narrow),pattern = "\\.\\.\\.[[:digit:]]{1,}$",replacement = "")
    
    athena_195_360_narrow <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR"),matches("_TTR")) 
    
    itemwise_data <- rbind(athena_3360_2096_narrow,athena_195_360_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("SLNB"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid" = "test_sessions.datasetid")) %>% 
      arrange(bblid) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) 
    
    
  }else if(test == "VSPLOT"){
    test_grep_15 <- "^VSPLOT15.VSPLOT15"
    test_grep_24 <- "^VSPLOT24.VSPLOT15"
    
    
    athena_3360_2096_narrow_15 <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_15)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_RT_[[:digit:]]")) 
    
    athena_195_360_narrow_15 <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_15)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_RT_[[:digit:]]")) 
    
    athena_3360_2096_narrow_24 <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_24)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_RT_[[:digit:]]")) 

    colnames(athena_3360_2096_narrow_15) <- str_replace_all(colnames(athena_3360_2096_narrow_15),pattern = "^VSPLOT15\\.",replacement = "")
    colnames(athena_195_360_narrow_15) <- str_replace_all(colnames(athena_195_360_narrow_15),pattern = "^VSPLOT15\\.",replacement = "")
    colnames(athena_3360_2096_narrow_24) <- str_replace_all(colnames(athena_3360_2096_narrow_24),pattern = "^VSPLOT24\\.",replacement = "")
     
    itemwise_data <- rbind(athena_3360_2096_narrow_15,athena_195_360_narrow_15,athena_3360_2096_narrow_24) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("VSPLOT"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x))) 
      
  } else if(test == "SPCPTN"){
    
    test_grep_nl <- "SPCPTNL"
    test_grep_n <- "SPCPTN90"
    
    athena_3360_2096_narrow_both <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_nl)) %>% 
      select(!matches("19|2[[:digit:]]|3[[:digit:]]")) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_3360_2096_narrow_N_only <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_n)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_both <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_nl)) %>% 
      select(!matches("19|2[[:digit:]]|3[[:digit:]]")) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_N_only <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_n)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    colnames(athena_3360_2096_narrow_both) <- str_replace_all(colnames(athena_3360_2096_narrow_both),pattern = "^SPCPTNL",replacement = "SPCPTN90")
    colnames(athena_195_360_narrow_both) <- str_replace_all(colnames(athena_195_360_narrow_both),pattern = "^SPCPTNL",replacement = "SPCPTN90")

    itemwise_data <- rbind(athena_3360_2096_narrow_both,athena_3360_2096_narrow_N_only,athena_195_360_narrow_both,athena_195_360_narrow_N_only) %>% 
      mutate(across(.cols = matches("TTR$"),.fns = ~ case_when(.x < 0 ~ NA_real_,TRUE ~ .x))) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("SPCPTN90"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x)))
    
  }
  else if(test == "CPF"){
    test_grep_a <- "^CPF_A.CPF_TRIAL"
    test_grep_b <- "^CPF_B.CPF_TRIAL"
    
    athena_3360_2107_narrow_a <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_369_narrow_a <- athena_195_369 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_b <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_3360_2096_narrow_b <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    colnames(athena_3360_2107_narrow_a) <- str_replace_all(colnames(athena_3360_2107_narrow_a),pattern = "^CPF_A",replacement = "CPF")
    colnames(athena_195_369_narrow_a) <- str_replace_all(colnames(athena_195_369_narrow_a),pattern = "^CPF_A",replacement = "CPF")
    colnames(athena_3360_2096_narrow_b) <- str_replace_all(colnames(athena_3360_2096_narrow_b),pattern = "^CPF_B",replacement = "CPF")
    colnames(athena_195_360_narrow_b) <- str_replace_all(colnames(athena_195_360_narrow_b),pattern = "^CPF_B",replacement = "CPF")
    
    itemwise_data <- rbind(athena_3360_2107_narrow_a,athena_195_369_narrow_a,athena_3360_2096_narrow_b,athena_195_360_narrow_b) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("CPF"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) %>% 
      mutate(across(.cols = matches("_TTR$"),.fns = ~ case_when(.x < 0 ~ NA_real_,TRUE ~ .x))) 

  }
  else if(test == "ER40"){
    test_grep_A <- "^ER40_A.ER40A_QID"
    test_grep_C <- "^ER40_C.ER40C_QID"
    test_grep_D <- "^ER40_D.ER40D_QID"
    
    colnames(athena_3360_2107) <- str_replace_all(colnames(athena_3360_2107),pattern = '(ER40_A\\.ER40)[^A]',replacement = '\\1A_')
    colnames(athena_195_369) <- str_replace_all(colnames(athena_195_369),pattern = '(ER40_A\\.ER40)[^A]',replacement = '\\1A_')
    
    athena_3360_2096_narrow_d <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_D)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 

    athena_3360_2107_narrow_a <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_A)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_3360_2107_narrow_c <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_C)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_d <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_D)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_369_narrow_a <- athena_195_369 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_A)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$"))
    
    colnames(athena_3360_2096_narrow_d) <- str_replace_all(colnames(athena_3360_2096_narrow_d),pattern = "^ER40_D\\.ER40D",replacement = "ER40.ER40")
    colnames(athena_3360_2107_narrow_a) <- str_replace_all(colnames(athena_3360_2107_narrow_a),pattern = "^ER40_A\\.ER40A",replacement = "ER40.ER40")
    colnames(athena_3360_2107_narrow_c) <- str_replace_all(colnames(athena_3360_2107_narrow_c),pattern = "^ER40_C\\.ER40C",replacement = "ER40.ER40")
    colnames(athena_195_360_narrow_d) <- str_replace_all(colnames(athena_195_360_narrow_d),pattern = "^ER40_D\\.ER40D",replacement = "ER40.ER40")
    colnames(athena_195_369_narrow_a) <- str_replace_all(colnames(athena_195_369_narrow_a),pattern = "^ER40_A\\.ER40A",replacement = "ER40.ER40")
    
    itemwise_data <- rbind(athena_3360_2096_narrow_d,athena_3360_2107_narrow_a,athena_3360_2107_narrow_c,athena_195_360_narrow_d,athena_195_369_narrow_a) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("ER40"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) 
    
  }
  
  else if(test == "SVOLT"){
    
    test_grep = "^SVOLT_A.SVOLT_TRIAL"
    
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    itemwise_data <- rbind(athena_3360_2096_narrow,athena_195_360_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("SVOLT"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) 
    
  }
  
  else if(test == "MEDF"){
    test_grep_a <-  "^MEDF36_A.MEDF36A_QID"
    test_grep_b <- "^MEDF36_B.MEDF36B_QID"
    
    athena_3360_2096_narrow_a <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_3360_2107_narrow_b <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_a <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    colnames(athena_3360_2096_narrow_a) <- str_replace_all(colnames(athena_3360_2096_narrow_a),pattern = 'MEDF36_A.MEDF36A',replacement = "MEDF36.MEDF36")
    colnames(athena_3360_2107_narrow_b) <- str_replace_all(colnames(athena_3360_2107_narrow_b),pattern = 'MEDF36_B.MEDF36B',replacement = "MEDF36.MEDF36")
    colnames(athena_195_360_narrow_a) <- str_replace_all(colnames(athena_195_360_narrow_a),pattern = 'MEDF36_A.MEDF36A',replacement = "MEDF36.MEDF36")
    
    itemwise_data <- rbind(athena_3360_2096_narrow_a,athena_3360_2107_narrow_b,athena_195_360_narrow_a) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("MEDF"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) 
    
  }
  
  else if(test == "ADT"){
    test_grep_a <-  "^ADT36_A.ADT36A_QID"
    test_grep_b <- "^ADT36_B\\.q"
    test_grep_b_full <- "^ADT36_B.ADT36B_QID"
    
    athena_3360_2096_narrow_a <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_360_narrow_a <- athena_195_360 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_a)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_CORR$"),matches("_TTR$")) 
    
    athena_195_369_narrow_b <- athena_195_369 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_corr$"),matches("_ttr$")) 
    
    athena_3360_2107_narrow_b <- athena_3360_2107 %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches(test_grep_b_full)) %>% 
      select(test_sessions.bblid,test_sessions.datasetid,matches("_corr$"),matches("_ttr$")) 
    
    colnames(athena_3360_2096_narrow_a) <- str_replace_all(colnames(athena_3360_2096_narrow_a),pattern = 'ADT36_A.ADT36A',replacement = "ADT36.ADT36")
    colnames(athena_195_360_narrow_a) <- str_replace_all(colnames(athena_195_360_narrow_a),pattern = 'ADT36_A.ADT36A',replacement = "ADT36.ADT36")
    colnames(athena_3360_2107_narrow_b) <- str_replace_all(colnames(athena_3360_2107_narrow_b),pattern = 'ADT36_B.ADT36B',replacement = "ADT36.ADT36")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = 'ADT36_B\\.',replacement = "ADT36.ADT36_")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = '_corr',replacement = "_CORR")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = '_ttr',replacement = "_TTR")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = 'q([0-9][0-9])',replacement = "QID0000\\1")
    colnames(athena_195_369_narrow_b) <- str_replace_all(colnames(athena_195_369_narrow_b),pattern = 'q([0-9])',replacement = "QID00000\\1")
    
    itemwise_data <- rbind(athena_3360_2096_narrow_a,athena_195_360_narrow_a,athena_195_369_narrow_b,athena_3360_2107_narrow_b) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      rename(bblid = test_sessions.bblid) %>% 
      filter(if_any(.cols = matches("ADT"),.fns = ~ !is.na(.x))) %>% 
      semi_join(CNB_cross,by = c("bblid","test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("_CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ case_when(is.na(.x) ~ 0,TRUE ~ .x))) %>% 
      mutate(across(.cols = matches("_TTR$"),.fns = ~ case_when(.x < 0 ~ NA_real_,TRUE ~ .x))) 
  }
  
  if(test == "VSPLOT"){
    dat <- itemwise_data[,c(grep("_CORR",colnames(itemwise_data)),grep("_RT_[[:digit:]]",colnames(itemwise_data)))]
  }else{
    dat <- itemwise_data[,c(grep("_CORR",colnames(itemwise_data)),grep("_TTR",colnames(itemwise_data)))]
  }
  
  dat <- as.data.frame(dat)
  items <- ncol(dat)/2
  dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
  res <- matrix(NA,dim(dat)[1],items)
  for(j in 1:items){
    if(all(is.na(dat[,(j+items)]))|all(is.na(dat[,j]))){
      res[,j] <- NA
    }else{
      mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
      res[,j] <- scale(residuals(mod,na.action=na.exclude))
    }
  }
  res2 <- res
  res2[abs(res2) < 2] <- 0
  res2[abs(res2) > 2] <- 1
  outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
  dat2 <- dat[,1:items]
  pfit1 <- r.pbis(dat2)$PFscores
  pfit1_vec <- ifelse(rowSums(dat2,na.rm = TRUE) == ncol(dat2),1,pfit1 %>% pull())
  pfit1_final <- ifelse(rowSums(dat2,na.rm = TRUE) == 0,0,pfit1_vec)
  pfit2 <- E.KB(dat2)$PFscores
  pfit2_vec <- ifelse(rowSums(dat2,na.rm = TRUE) == ncol(dat2),1,pfit2 %>% pull())
  pfit2_final <- ifelse(rowSums(dat2,na.rm = TRUE) == 0,0,pfit2_vec)
  acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))][,1:3],na.rm = T)
  if(test_type == "memory"){
    sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1_final) + (0.50*pfit2_final)
  } else if(test_type == "non-memory"){
    sc <- (0.34*outlier_score_2cut) + (0.22*pfit1_final) + (0.44*pfit2_final)
  }
  
  col_name <- paste0("PFscores_",test)
  itemwise_data[,col_name] <- sc
  
  output_df <- output_df %>% 
    left_join(itemwise_data %>% select(bblid,test_sessions.datasetid,matches("^PFscores_"),matches("_complete")), by = c("bblid","test_sessions.datasetid"))
  
  return(output_df)
}

tests <- c("CPF","SVOLT","PMAT","ER40","MEDF","SLNB","ADT","VSPLOT","SPCPTN")
test_type <- c("memory","memory","non-memory","non-memory","non-memory","memory","non-memory","non-memory","non-memory")

CNB_with_SMVE <- CNB_cross

for(i in 1:length(tests)){
  CNB_with_SMVE <- add_PFscores(test = tests[i],test_type = test_type[i],output_df = CNB_with_SMVE)
}

CNB_with_SMVE <- CNB_with_SMVE %>% 
  mutate(across(.cols = c(PFscores_CPF,PFscores_SVOLT,PFscores_PMAT,PFscores_ER40,PFscores_MEDF,PFscores_SLNB,PFscores_ADT,PFscores_VSPLOT,PFscores_SPCPTN),.fns = ~ case_when(.x <= quantile(.x,.05,na.rm = T) ~ "F",is.na(.x) ~ NA_character_,TRUE ~ "V"),.names = "{.col}_flag")) 

reviewer_comments <- reviewer_comments %>% 
  mutate(Overall_valid_comments =  case_when(is.na(Overall_valid_comments) ~ 'V',TRUE ~ Overall_valid_comments)) 

CNB_with_completed <- CNB_with_SMVE %>% 
  mutate(ADT_completed = ifelse(!is.na(adi_az_asr) | !is.na(adi_sz_asr),'Completed','Not Completed')) %>% 
  mutate(CPF_completed = ifelse(!is.na(fmem_az_asr) | !is.na(fmem_sz_asr),'Completed','Not Completed')) %>% 
  mutate(ER40_completed = ifelse(!is.na(eid_az_asr) | !is.na(eid_sz_asr),'Completed','Not Completed')) %>% 
  mutate(MEDF_completed = ifelse(!is.na(edi_az_asr) | !is.na(edi_sz_asr),'Completed','Not Completed')) %>% 
  mutate(MPRACT_completed = ifelse(!is.na(sm_sz_asr),'Completed','Not Completed')) %>% 
  mutate(PCET_completed = ifelse(!is.na(abf_az_asr) | !is.na(abf_sz_asr),'Completed','Not Completed')) %>% 
  mutate(PMAT_completed = ifelse(!is.na(nvr_az_asr) | !is.na(nvr_sz_asr),'Completed','Not Completed')) %>% 
  mutate(SCTAP_completed = ifelse(!is.na(mot_sz_asr),'Completed','Not Completed')) %>% 
  mutate(ADT_completed = ifelse(!is.na(adi_az_asr) | !is.na(adi_sz_asr),'Completed','Not Completed')) %>% 
  mutate(SLNB_completed = ifelse(!is.na(wm_az_asr) | !is.na(wm_sz_asr),'Completed','Not Completed')) %>% 
  mutate(CPT_completed = ifelse(!is.na(att_az_asr) | !is.na(att_sz_asr),'Completed','Not Completed')) %>% 
  mutate(VOLT_completed = ifelse(!is.na(smem_az_asr) | !is.na(smem_sz_asr),'Completed','Not Completed')) %>% 
  mutate(PLOT_completed = ifelse(!is.na(spa_az_asr) | !is.na(spa_sz_asr),'Completed','Not Completed')) %>% 
  mutate(ADT_completed = ifelse(!is.na(adi_az_asr) | !is.na(adi_sz_asr),'Completed','Not Completed')) 
  
# Remove Ns on valid codes, Fs on assessor comments
CNB_all_QC <- CNB_with_completed %>% 
  mutate(bblid = as.integer(bblid)) %>% 
  left_join(reviewer_comments) %>%
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ case_when(adt_valid == 'N' | ADT_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(fmem_az_asr,fmem_sz_asr),.fns = ~ case_when(cpf_valid == 'N' | CPF_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ case_when(er40_valid == 'N' | ER40_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ case_when(medf_valid == 'N' | MEDF_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(sm_sz_asr),.fns = ~ case_when(mpraxis_valid == 'N' | MPRACT_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(abf_az_asr,abf_sz_asr),.fns = ~ case_when(pcet_valid == 'N' | PCET_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(nvr_az_asr,nvr_sz_asr),.fns = ~ case_when(pmat_valid == 'N' | PMAT_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(mot_sz_asr),.fns = ~ case_when(tap_valid == 'N' | SCTAP_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(wm_az_asr,wm_sz_asr),.fns = ~ case_when(lnb_valid == 'N' | SLNB_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(att_az_asr,att_sz_asr),.fns = ~ case_when(cpt_valid == 'N' | SPCPTN_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(smem_az_asr,smem_sz_asr),.fns = ~ case_when(volt_valid == 'N' | SVOLT_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ case_when(plot_valid == 'N' | VSPLOT_valid_comments == 'F' ~ NA_real_,TRUE ~ .x))) %>% 
  rowwise() %>% 
  mutate(ADT_flags = sum(c_across(c(adt_valid,PFscores_ADT_flag)) == 'F',na.rm = T)) %>% 
  mutate(CPF_flags = sum(c_across(c(cpf_valid,PFscores_CPF_flag)) == 'F',na.rm = T)) %>%
  mutate(ER40_flags = sum(c_across(c(er40_valid,PFscores_ER40_flag)) == 'F',na.rm = T)) %>%
  mutate(MEDF_flags = sum(c_across(c(medf_valid,PFscores_MEDF_flag)) == 'F',na.rm = T)) %>%
  mutate(PMAT_flags = sum(c_across(c(pmat_valid,PFscores_PMAT_flag)) == 'F',na.rm = T)) %>%
  mutate(SLNB_flags = sum(c_across(c(lnb_valid,PFscores_SLNB_flag)) == 'F',na.rm = T)) %>%
  mutate(CPT_flags = sum(c_across(c(cpt_valid,PFscores_SPCPTN_flag)) == 'F',na.rm = T)) %>%
  mutate(VOLT_flags = sum(c_across(c(volt_valid,PFscores_SVOLT_flag)) == 'F',na.rm = T)) %>%
  mutate(PLOT_flags = sum(c_across(c(plot_valid,PFscores_VSPLOT_flag)) == 'F',na.rm = T)) %>%
  ungroup()
  
# Remove tests flagged on a majority of available QC metrics
CNB_all_QC <- CNB_all_QC %>% 
  filter(!(bblid %in% c(21223,123456))) %>% # Remove bblids which were from tests or parents of 22q
  mutate(across(.cols = c(adi_az_asr,adi_sz_asr),.fns = ~ case_when(ADT_flags == 2 ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(across(.cols = c(fmem_az_asr,fmem_sz_asr),.fns = ~ case_when(CPF_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(eid_az_asr,eid_sz_asr),.fns = ~ case_when(ER40_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(edi_az_asr,edi_sz_asr),.fns = ~ case_when(MEDF_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(nvr_az_asr,nvr_sz_asr),.fns = ~ case_when(PMAT_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(wm_az_asr,wm_sz_asr),.fns = ~ case_when(SLNB_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(att_az_asr,att_sz_asr),.fns = ~ case_when(CPT_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(smem_az_asr,smem_sz_asr),.fns = ~ case_when(VOLT_flags == 2 ~ NA_real_,TRUE ~ .x))) %>%
  mutate(across(.cols = c(spa_az_asr,spa_sz_asr),.fns = ~ case_when(PLOT_flags == 2 ~ NA_real_,TRUE ~ .x))) %>% 
  mutate(bblid = as.character(bblid))


# Write clean data set to file 

CNB_final <- CNB_all_QC %>%
  filter(Overall_valid_comments != 'F') %>% 
  mutate(ADT_removed = case_when(ADT_completed == 'Not Completed' ~ NA_character_,adt_valid == 'N'|ADT_valid_comments == 'F'|ADT_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(CPF_removed = case_when(CPF_completed == 'Not Completed' ~ NA_character_,cpf_valid == 'N'|CPF_valid_comments == 'F'|CPF_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(ER40_removed = case_when(ER40_completed == 'Not Completed' ~ NA_character_,er40_valid == 'N'|ER40_valid_comments == 'F'|ER40_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(MEDF_removed = case_when(MEDF_completed == 'Not Completed' ~ NA_character_,medf_valid == 'N'|MEDF_valid_comments == 'F'|MEDF_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(PMAT_removed = case_when(PMAT_completed == 'Not Completed' ~ NA_character_,pmat_valid == 'N'|PMAT_valid_comments == 'F'|PMAT_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(SLNB_removed = case_when(SLNB_completed == 'Not Completed' ~ NA_character_,lnb_valid == 'N'|SLNB_valid_comments == 'F'|SLNB_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(CPT_removed = case_when(CPT_completed == 'Not Completed' ~ NA_character_,cpt_valid == 'N'|SPCPTN_valid_comments == 'F'|CPT_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(VOLT_removed = case_when(VOLT_completed == 'Not Completed' ~ NA_character_,volt_valid == 'N'|SVOLT_valid_comments == 'F'|VOLT_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(PLOT_removed = case_when(PLOT_completed == 'Not Completed' ~ NA_character_,plot_valid == 'N'|VSPLOT_valid_comments == 'F'|PLOT_flags == 2 ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(SCTAP_removed = case_when(SCTAP_completed == 'Not Completed' ~ NA_character_,tap_valid == 'N'|SCTAP_valid_comments == 'F' ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(MPRACT_removed = case_when(MPRACT_completed == 'Not Completed' ~ NA_character_,mpraxis_valid == 'N'|MPRACT_valid_comments == 'F' ~ 'Removed',TRUE ~ 'Valid')) %>% 
  mutate(PCET_removed = case_when(PCET_completed == 'Not Completed' ~ NA_character_,pcet_valid == 'N'|PCET_valid_comments == 'F' ~ 'Removed',TRUE ~ 'Valid')) %>% 
  select(!matches('_flags'))


write_csv(CNB_final,'/Users/hillmann/Projects/22q/Data/22q_in_person_vs_remote_CNB_11_04_2022_NH.csv')

  
# Merge old reviewer comments so that work doesn't have to be repeated

# reviewer_flags <- old_reviewer_comments %>% 
#   select(test_sessions.bblid,test_sessions.datasetid,contains("(1 = yes,"),Notes) %>% 
#   mutate(across(.cols = matches("[vV]alid"),.fns = ~ case_when(.x == 1 ~ 'V',.x == 2 ~ 'F',TRUE ~ NA_character_))) %>% 
#   rename(bblid = test_sessions.bblid) %>% 
#   mutate(bblid = as.character(bblid))
#   
# colnames(reviewer_flags) <- str_replace_all(colnames(reviewer_flags),pattern = ' .*',replacement = '')
# colnames(reviewer_flags) <- str_replace_all(colnames(reviewer_flags),pattern = '_[vV]alid$',replacement = '_valid_comments')
# 
# reviewer_flags_to_update <- CNB_with_SMVE %>% 
#   left_join(reviewer_flags) %>% 
#   select(bblid,test_sessions.datasetid,matches("valid_comments"),Notes) 
# 
# write_csv(reviewer_flags_to_update,file = '/Users/hillmann/Projects/22q/Data/Assessor_flags_to_update_10_06_2022.csv')

