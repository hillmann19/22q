# Read in data and packages

library(tidyverse)
library(table1)
CNB <- read_csv("~/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")

CNB_completed <- CNB %>%
  select(test_sessions.bblid,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR) %>%
  pivot_longer(cols = ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Response") %>%
  mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "36|24|15|90|2",replacement = "")) %>%
  mutate(Completed = case_when(!is.na(Response) ~ 1,is.na(Response) ~ 0)) %>% 
  select(-Response) %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) 

CNB$test_sessions.bblid <- as.character(CNB$test_sessions.bblid)

PF_score_NA <- CNB %>% 
  select(!matches("acc$|res$|pfit")) %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,remote,matches("PFscores")) %>% 
  pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Test",values_to = "PFscore") %>% 
  mutate(Test = str_replace_all(Test,pattern = "PFscores_",replacement = "")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  filter(is.na(PFscore)) 

CNB %>% 
  select(test_sessions.bblid,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR) %>%
  pivot_longer(cols = ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Response") %>%
  mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "36|24|15|90|2",replacement = "")) %>% 
  semi_join(PF_score_NA) %>% 
  View()

intersection_df_SMVE <- CNB %>% 
  select(!matches("acc$|res$|pfit")) %>% 
  mutate(across(.cols = contains("PFscores"),~ case_when(.x < quantile(.x,.05,na.rm = T) ~ "F",TRUE ~ "V"),.names = "{.col}_flag")) %>%
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,remote,matches("flag$")) %>% 
  mutate(across(.cols = contains("comment_flag"),.fns = ~ case_when(.x == 1 ~ "V",.x == 2 ~ "F",TRUE ~ "V"))) %>%
  pivot_longer(cols = ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
  mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment",str_detect(Test,pattern = "^PFscores") ~ "SMVE")) %>% 
  mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB",str_detect(Test,pattern = "SCTAP") ~ "SCTAP",str_detect(Test,pattern = "PCET") ~ "PCET",str_detect(Test,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Test,pattern = "MPRACT") ~ "MPRACT",str_detect(Test,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
  arrange(test_sessions.bblid,Test) %>% 
  filter(Test %in% c("ADT","CPF","ER40","MEDF","PMAT","SLNB","SVOLT")) %>% 
  mutate(QC_combined = case_when(AV == "V" & Comment == "V" & SMVE == "V" ~ "All valid",AV == "F" & Comment == "V" & SMVE == "V" ~ "Only AV flagged",AV == "V" & Comment == "F" & SMVE == "V" ~ "Only reviewer flagged",AV == "V" & Comment == "V" & SMVE == "F" ~ "Only SMVE flagged",AV == "F" & Comment == "F" & SMVE == "V" ~ "AV and reviewer flagged",AV == "V" & Comment == "F" & SMVE == "F" ~ "Reviewer and SMVE flagged",AV == "F" & Comment == "V" & SMVE == "F" ~ "AV and SMVE flagged",AV == "F" & Comment == "F" & SMVE == "F" ~ "All flagged")) %>% 
  mutate(QC_combined = factor(QC_combined,levels = c("All valid","Only AV flagged","Only reviewer flagged","Only SMVE flagged","AV and SMVE flagged","AV and reviewer flagged","Reviewer and SMVE flagged","All flagged"))) %>% 
  rename(Age = test_sessions_v.age,Sex = test_sessions_v.gender,Remote = remote,"QC combined" = QC_combined)

table1(~ Age + Sex + `QC combined`|Remote,data = intersection_df_SMVE)

intersection_df_no_SMVE <- CNB %>% 
  select(!matches("PFscores")) %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,remote,matches("flag$")) %>% 
  mutate(across(.cols = contains("comment_flag"),.fns = ~ case_when(.x == 1 ~ "V",.x == 2 ~ "F",TRUE ~ "V"))) %>%
  pivot_longer(cols = ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
  mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment")) %>% 
  mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB",str_detect(Test,pattern = "SCTAP") ~ "SCTAP",str_detect(Test,pattern = "PCET") ~ "PCET",str_detect(Test,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Test,pattern = "MPRACT") ~ "MPRACT",str_detect(Test,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
  arrange(test_sessions.bblid,Test) %>% 
  filter(!(Test %in% c("ADT","CPF","ER40","MEDF","PMAT","SLNB","SVOLT"))) %>% 
  mutate(QC_combined = case_when(AV == "V" & Comment == "V" ~ "Both valid",AV == "V" & Comment == "F" ~ "Only reviewer flagged",AV == "F" & Comment == "V" ~ "Only auto-validation flagged",AV == "F" & Comment == "F" ~ "Both flagged")) %>% 
  mutate(QC_combined = factor(QC_combined,levels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"))) %>% 
  rename(Age = test_sessions_v.age,Sex = test_sessions_v.gender,Remote = remote,"QC combined" = QC_combined)

table1(~ Age + Sex + `QC combined`|Remote,data = intersection_df_no_SMVE)


  