# Read in data and packages
library(table1)
library(tidyverse)
rename <- dplyr::rename
select <- dplyr::select
CNB <- read_csv("~/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")

CNB_completed <- CNB %>%
  select(test_sessions.bblid,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90_TP,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR) %>%
  pivot_longer(cols = ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Response") %>%
  mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "36|24|15|90|2",replacement = "")) %>%
  mutate(Completed = case_when(!is.na(Response) ~ 1,is.na(Response) ~ 0)) %>% 
  select(-Response) %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) 

CNB$test_sessions.bblid <- as.character(CNB$test_sessions.bblid)

intersection_df_SMVE <- CNB %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,Notes,remote,matches("flag$")) %>% 
  pivot_longer(cols = ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
  mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment",str_detect(Test,pattern = "^PFscores") ~ "SMVE")) %>% 
  mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB",str_detect(Test,pattern = "SCTAP") ~ "SCTAP",str_detect(Test,pattern = "PCET") ~ "PCET",str_detect(Test,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Test,pattern = "MPRACT") ~ "MPRACT",str_detect(Test,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
  arrange(test_sessions.bblid,Test) %>% 
  filter(Test %in% c("ADT","CPF","ER40","MEDF","PMAT","SLNB","SVOLT","VSPLOT","SPCPTN")) %>% 
  mutate(QC_combined = case_when(AV == "V" & Comment == "V" & SMVE == "V" ~ "All valid",AV == "F" & Comment == "V" & SMVE == "V" ~ "Only auto-validation flagged",AV == "V" & Comment == "F" & SMVE == "V" ~ "Only reviewer flagged",AV == "V" & Comment == "V" & SMVE == "F" ~ "Only SMVE flagged",AV == "F" & Comment == "F" & SMVE == "V" ~ "Auto-validation and reviewer flagged",AV == "V" & Comment == "F" & SMVE == "F" ~ "Reviewer and SMVE flagged",AV == "F" & Comment == "V" & SMVE == "F" ~ "Auto-validation and SMVE flagged",AV == "F" & Comment == "F" & SMVE == "F" ~ "All flagged")) %>% 
  mutate(QC_combined = factor(QC_combined,levels = c("All valid","Only auto-validation flagged","Only reviewer flagged","Only SMVE flagged","Auto-validation and SMVE flagged","Auto-validation and reviewer flagged","Reviewer and SMVE flagged","All flagged"))) %>% 
  rename(Age = test_sessions_v.age) %>% 
  rename(Sex = test_sessions_v.gender) %>% 
  rename(Remote = remote) %>% 
  colnames()
  rename("QC combined" = QC_combined) 

table1(~ Age + Sex + `QC combined`|Remote,data = intersection_df_SMVE)

intersection_df_no_SMVE <- CNB %>% 
  select(!matches("PFscores")) %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,Notes,remote,matches("flag$")) %>% 
  pivot_longer(cols = ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
  mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment")) %>% 
  mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB",str_detect(Test,pattern = "SCTAP") ~ "SCTAP",str_detect(Test,pattern = "PCET") ~ "PCET",str_detect(Test,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Test,pattern = "MPRACT") ~ "MPRACT",str_detect(Test,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
  arrange(test_sessions.bblid,Test) %>% 
  filter(!(Test %in% c("ADT","CPF","ER40","MEDF","PMAT","SVOLT","VSPLOT"))) %>% 
  mutate(QC_combined = case_when(AV == "V" & Comment == "V" ~ "Both valid",AV == "V" & Comment == "F" ~ "Only reviewer flagged",AV == "F" & Comment == "V" ~ "Only auto-validation flagged",AV == "F" & Comment == "F" ~ "Both flagged")) %>% 
  mutate(QC_combined = factor(QC_combined,levels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"))) %>% 
  rename(Age = test_sessions_v.age,Sex = test_sessions_v.gender,Remote = remote,"QC combined" = QC_combined)

table1(~ `QC combined`|Test,data = intersection_df_no_SMVE)

# Compare pre and post qc - what percent of tests are completed

CNB_post <- CNB %>% 
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
  mutate(across(.cols = c(SVOLT_A.SVOLT_RTCR),.fns = ~ ifelse(sum(c(PFscores_SVOLT_flag,SVOLT_comment_flag,SVOLT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(VSPLOT15.VSPLOT15_RTCR,VSPLOT15.VSPLOT15_CR),.fns = ~ ifelse(sum(c(PFscores_VSPLOT_flag,VSPLOT_comment_flag,VSPLOT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(SLNB2_90.SLNB2_MRTC,SLNB2_90.SLNB2_MCR),.fns = ~ ifelse(sum(c(SLNB_comment_flag,SLNB_AV_flag) == "F") >= 1,NA,.x))) %>%
  mutate(across(.cols = c(SPCPTN90_TPRT,SPCPTN90_TP),.fns = ~ ifelse(sum(c(SPCPTN_comment_flag,SPCPTN_AV_flag) == "F") >= 1,NA,.x))) %>%
  mutate(across(.cols = c(SCTAP.SCTAP_TOT),.fns = ~ ifelse(sum(c(SCTAP_comment_flag,SCTAP_AV_flag) == "F") >= 1,NA,.x))) %>%
  ungroup() %>%
  select(!matches("_flag$")) %>% 
  select(test_sessions.bblid,ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,MPRACT.MP2RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SCTAP.SCTAP_TOT) %>% 
  rowwise() %>% 
  mutate(Number_of_tests_post = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SCTAP.SCTAP_TOT)))) %>% 
  pivot_longer(cols = ADT36_A.ADT36A_RTCR:SCTAP.SCTAP_TOT,names_to = "Test",values_to = "Value") %>% 
  group_by(test_sessions.bblid) %>% 
  mutate(Tests_post = str_c(Test[!is.na(Value)],collapse = ", ")) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Test,values_from = Value) %>% 
  filter(Number_of_tests_post <= 2) %>% 
  select(test_sessions.bblid,Tests_post,Number_of_tests_post)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")

CNB %>% 
  select(test_sessions.bblid,matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("SCTAP_TOT$")) %>% 
  select(-CPF_B.CPF_W_RTCR,-SPCPTN90.SCPN90_TPRT,-SPCPTNL.SCPN_TPRT) %>% 
  semi_join(CNB_post[,c("test_sessions.bblid")]) %>%
  left_join(CNB_post) %>% 
  rowwise() %>% 
  mutate(Number_of_tests_pre = sum(!is.na(c_across(cols = ADT36_A.ADT36A_RTCR:SCTAP.SCTAP_TOT)))) %>% 
  ungroup() %>% 
  pivot_longer(cols = ADT36_A.ADT36A_RTCR:SCTAP.SCTAP_TOT,names_to = "Test",values_to = "Value") %>% 
  group_by(test_sessions.bblid) %>% 
  mutate(Tests_pre = str_c(Test[!is.na(Value)],collapse = ", ")) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Test,values_from = Value) %>% 
  select(test_sessions.bblid,Number_of_tests_pre,Tests_pre,Number_of_tests_post,Tests_post) %>% 
  mutate(test_sessions.bblid = factor(test_sessions.bblid)) %>% 
  pivot_longer(cols = Number_of_tests_pre:last_col(),names_to = c(".value","Time"),names_pattern = "(.*)_(.*)$") %>% 
  mutate(Time = case_when(Time == "post" ~ "After QC",Time == "pre" ~ "Before QC")) %>% 
  mutate(Time = factor(Time,levels = c("Before QC","After QC"))) %>% 
  mutate(Tests = str_replace_all(Tests,pattern = "MPRACT.MP2RTCR","MPRACT")) %>% 
  mutate(Tests = str_replace_all(Tests,pattern = "MEDF36_A.MEDF36A_RTCR",replacement = "MEDF")) %>% 
  mutate(Tests = str_replace_all(Tests,pattern = "SPCPTN90_TPRT",replacement = "CPT")) %>% 
  mutate(Tests = str_replace_all(Tests,pattern = "PCET_A.PCET_RTCR",replacement = "PCET")) %>% 
  mutate(Tests = str_replace_all(Tests,pattern = "SCTAP.SCTAP_TOT",replacement = "SCTAP")) %>% 
  ggplot(aes(x = Number_of_tests,y = fct_reorder(test_sessions.bblid,Number_of_tests,.fun = mean),color = Time,fill = Time)) + geom_bar(stat = "identity",position = position_dodge(width = .9)) + labs(x = "# Tests Completed",y = "BBLID",color = "",fill = "") + scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + geom_text(aes(label = Tests),size = 3,hjust = -.05,position = position_dodge(width = .9)) + xlim(c(0,3.3))

# Generate data set with those flagged by at least one (but not all) qc methods

Performance_df_long <- CNB %>% 
  select(test_sessions.bblid,ADT36_A.ADT36A_CR:VSPLOT15.VSPLOT15_RTCR) %>% 
  select(!(matches("LRR|CNT|CPF_B.CPF_W_RTCR|LRSR|HAP|MP2$|RTTO|UNF|EXCESS|TP0|FP0|LRNR|SPCPTNL|NR$|TOT_RT$|SUM_|PC$|_RT$|CAT$"))) %>% 
  pivot_longer(cols = ADT36_A.ADT36A_CR:last_col(),names_to = "Variable",values_to = "Value") %>% 
  mutate(Test = case_when(str_detect(Variable,pattern = "ER40") ~ "ER40",str_detect(Variable,pattern = "PMAT") ~ "PMAT",str_detect(Variable,pattern = "CPF") ~ "CPF",str_detect(Variable,pattern = "SVOLT") ~ "SVOLT",str_detect(Variable,pattern = "ADT") ~ "ADT",str_detect(Variable,pattern = "MEDF") ~ "MEDF",str_detect(Variable,pattern = "SLNB") ~ "SLNB",str_detect(Variable,pattern = "SCTAP") ~ "SCTAP",str_detect(Variable,pattern = "PCET") ~ "PCET",str_detect(Variable,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Variable,pattern = "MPRACT") ~ "MPRACT",str_detect(Variable,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
  mutate(Metric = case_when(str_detect(Variable,pattern = "_CR$") ~ "Correct responses",str_detect(Variable,pattern = "RTCR$") ~ "Reaction time correct responses",str_detect(Variable,pattern = "_ACC2$") ~ "Accuracy",str_detect(Variable,pattern = "_TOT$") ~ "Total taps",str_detect(Variable,pattern = "_MCR$") ~ "Total True Positives (1-back and 2-back)",str_detect(Variable,pattern = "_MRTC$") ~ "Median response time (1-back and 2-back correct)",str_detect(Variable,pattern = "_TP$") ~ "True Positives",str_detect(Variable,pattern = "_TPRT$") ~ "Reaction time true positives")) %>% 
  mutate(Metric_type = case_when(str_detect(Variable,pattern = "_CR$|_ACC2$|_MCR$|_TP$") ~ "Accuracy",TRUE ~ "Reaction time")) %>% 
  select(-Variable) %>% 
  pivot_wider(names_from = c("Metric_type"),values_from = c("Value","Metric")) %>% 
  rename(Label_Accuracy = Metric_Accuracy,Label_reaction_time = `Metric_Reaction time`,Accuracy = Value_Accuracy,Reaction_time = `Value_Reaction time`) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  group_by(Test) %>% 
  mutate(Accuracy_Z = as.numeric(scale(Accuracy))) %>% 
  mutate(Reaction_time_Z = as.numeric(scale(Reaction_time))) %>% 
  ungroup() %>% 
  relocate(Accuracy_Z,.after = Accuracy) %>% 
  relocate(Reaction_time_Z,.after = Reaction_time) %>% 
  relocate(Label_Accuracy,.before = Accuracy) %>% 
  relocate(Label_reaction_time,.before = Reaction_time)

AV_vars_flagged <- CNB %>% 
  select(test_sessions.bblid,ADT36_A.ADT36A_LRR_750,ADT36_A.ADT36A_CNT_200,ADT36_A.ADT36A_LRSR,CPF_B.CPF_RT,CPF_B.CPF_CNT_200,CPF_B.CPF_LRSR,
         ER40_D.ER40D_CR,ER40_D.ER40D_HAP,ER40_D.ER40D_RT,ER40_D.ER40D_CNT_250,ER40_D.ER40D_LRSR,
         MEDF36_A.MEDF36A_LRR_750,MEDF36_A.MEDF36A_CNT_200,MEDF36_A.MEDF36A_LRSR,MPRACT.MP2,PCET_A.PCET_LRSR,
         PMAT24_A.PMAT24_A_CR,PMAT24_A.PMAT24_A_RTTO,SCTAP.SCTAP_UNF_TOT,SCTAP.SCTAP_EXCESS_TOT,SCTAP.SCTAP_UNF_STD,
         SLNB2_90.SLNB2_TP0,SLNB2_90.SLNB2_FP0,SLNB2_90.SLNB2_LRR,SLNB2_90.SLNB2_LRNR,SLNB2_90.SLNB2_CNT_200,
         SCPTN_LRNR,SVOLT_A.SVOLT_RT,SVOLT_A.SVOLT_CNT_200,SVOLT_A.SVOLT_LRR_200,SVOLT_A.SVOLT_LRSR,SVOLT_A.SVOLT_NR,
         VSPLOT15.VSPLOT15_TOT_RT,VSPLOT15.VSPLOT15_SUM_DEG_OFF,VSPLOT15.VSPLOT15_SUM_EXCESS,VSPLOT15.VSPLOT15_SUM_DEFICIT) %>% 
  pivot_longer(ADT36_A.ADT36A_LRR_750:last_col(),names_to = "Variable",values_to = "Value") %>% 
  mutate(Test = case_when(str_detect(Variable,pattern = "ER40") ~ "ER40",str_detect(Variable,pattern = "PMAT") ~ "PMAT",str_detect(Variable,pattern = "CPF") ~ "CPF",str_detect(Variable,pattern = "SVOLT") ~ "SVOLT",str_detect(Variable,pattern = "ADT") ~ "ADT",str_detect(Variable,pattern = "MEDF") ~ "MEDF",str_detect(Variable,pattern = "SLNB") ~ "SLNB",str_detect(Variable,pattern = "SCTAP") ~ "SCTAP",str_detect(Variable,pattern = "PCET") ~ "PCET",str_detect(Variable,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Variable,pattern = "MPRACT") ~ "MPRACT",str_detect(Variable,pattern = "VSPLOT") ~ "VSPLOT")) %>%
  mutate(Var_flagged = case_when(Variable == "ADT36_A.ADT36A_LRR_750" & Value > 3 ~ "F",Variable == "ADT36_A.ADT36A_CNT_200" & Value > 3 ~ "F",Variable == "ADT36_A.ADT36A_LRSR" & Value > 11 ~ "F",
         Variable == "CPF_B.CPF_CNT_200" & Value >= 21 ~ "F",Variable == "CPF_B.CPF_LRSR" & Value >= 21 ~ "F",
         Variable == "ER40_D.ER40D_CR" & Value < 12 ~ "F",Variable == "ER40_D.ER40D_HAP" & Value < 6 ~ "F",Variable == "ER40_D.ER40D_CNT_250" & Value > 2 ~ "F",Variable == "ER40_D.ER40D_LRSR" & Value > 11 ~ "F",
         Variable == "MEDF36_A.MEDF36A_LRR_750 " & Value > 3 ~ "F",Variable == "MEDF36_A.MEDF36A_CNT_200 " & Value > 3 ~ "F",Variable == "MEDF36_A.MEDF36A_LRSR " & Value > 11 ~ "F",
         Variable == "MPRACT.MP2" & Value < 16 ~ "F",Variable == "PCET_A.PCET_LRSR" & Value > 10 ~ "F",Variable == "PMAT24_A.PMAT24_A_CR" & Value < 5 ~ "F",
         Variable == "PMAT24_A.PMAT24_A_RTTO" & Value < 2000 ~ "F",Variable == "SCTAP.SCTAP_UNF_TOT" & Value < 70 ~ "F",Variable == "SCTAP.SCTAP_EXCESS_TOT" & Value > 38 ~ "F",Variable == "SCTAP.SCTAP_UNF_STD" & Value > 17 ~ "F",
         Variable == "SLNB2_90.SLNB2_TP0" & Value < 6 ~ "F", Variable == "SLNB2_90.SLNB2_FP0" & Value > 4 ~ "F", Variable == "SLNB2_90.SLNB2_LRR" & Value > 8 ~ "F", Variable == "SLNB2_90.SLNB2_LRNR" & Value > 26 ~ "F", Variable == "SLNB2_90.SLNB2_CNT_200" & Value > 2 ~ "F",
         Variable == "SCPTN_LRNR" & Value > 36 ~ "F",Variable == "SVOLT_A.SVOLT_CNT_200" & Value >= 11 ~ "F",Variable == "SVOLT_A.SVOLT_LRSR" & Value >= 11 ~ "F",
        Variable == "VSPLOT15.VSPLOT15_SUM_DEG_OFF" & Value > 500 ~ "F",Variable == "VSPLOT15.VSPLOT15_SUM_EXCESS" & Value > 200 ~ "F",Variable == "VSPLOT15.VSPLOT15_SUM_DEFICIT" & Value > 45 ~ "F",TRUE ~ "V")) %>% 
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  select(-Completed) %>% 
  group_by(test_sessions.bblid,Test) %>% 
  mutate(Vars_flagged = case_when(n() == sum(Var_flagged == "V") ~ "No variables flagged",TRUE ~ str_c(Variable[Var_flagged == "F"],Value[Var_flagged == 'F'],sep = ":",collapse = ","))) %>% 
  ungroup() %>% 
  rename(Vars_flagged_AV = Vars_flagged) %>% 
  distinct(test_sessions.bblid,Test,.keep_all = T) %>% 
  select(test_sessions.bblid,Test,Vars_flagged_AV) 
  
SMVE_subscores <- CNB %>% 
  select(test_sessions.bblid,contains("PFscores")) %>% 
  select(!matches("flag$")) %>% 
  pivot_longer(PFscores_CPF:last_col(),names_to = "Variable",values_to = "Value") %>% 
  mutate(Test = case_when(str_detect(Variable,pattern = "ER40") ~ "ER40",str_detect(Variable,pattern = "PMAT") ~ "PMAT",str_detect(Variable,pattern = "CPF") ~ "CPF",str_detect(Variable,pattern = "SVOLT") ~ "SVOLT",str_detect(Variable,pattern = "ADT") ~ "ADT",str_detect(Variable,pattern = "MEDF") ~ "MEDF",str_detect(Variable,pattern = "SLNB") ~ "SLNB",str_detect(Variable,pattern = "SCTAP") ~ "SCTAP",str_detect(Variable,pattern = "PCET") ~ "PCET",str_detect(Variable,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Variable,pattern = "MPRACT") ~ "MPRACT",str_detect(Variable,pattern = "VSPLOT") ~ "VSPLOT")) %>%
  left_join(CNB_completed) %>% 
  filter(Completed == 1) %>% 
  mutate(PFscore_type = case_when(str_detect(Variable,pattern = "res$") ~ "Outlier_reaction_time",str_detect(Variable,pattern = "acc$") ~ "Accuracy_easy",str_detect(Variable,pattern = "pfit1$") ~ "Person_fit_1",str_detect(Variable,pattern = "pfit2$") ~ "Person_fit_2",TRUE ~ "SMVE")) %>% 
  select(-Completed,-Variable) %>% 
  pivot_wider(names_from = "PFscore_type",values_from = "Value") %>% 
  group_by(Test) %>% 
  mutate(Outlier_reaction_time_z = as.numeric(scale(Outlier_reaction_time))) %>% 
  #mutate(Accuracy_easy_z = as.numeric(scale(Accuracy_easy))) %>% 
  mutate(Person_fit_1_z = as.numeric(scale(Person_fit_1))) %>% 
  mutate(Person_fit_2_z = as.numeric(scale(Person_fit_2))) %>% 
  ungroup() %>% 
  relocate(Outlier_reaction_time_z,.after = Outlier_reaction_time) %>% 
  relocate(Person_fit_1_z,.after = Person_fit_1) %>% 
  relocate(Person_fit_2_z,.after = Person_fit_2) 

qc_disagreement <- intersection_df_no_SMVE %>% 
  mutate(SMVE = NA_character_) %>% 
  bind_rows(intersection_df_SMVE) %>% 
  rename(SMVE_flag = SMVE,AV_flag = AV,Comment_flag = Comment,QC_combined = "QC combined") %>% 
  left_join(Performance_df_long) %>% 
  left_join(AV_vars_flagged) %>% 
  left_join(SMVE_subscores) %>% 
  relocate(SMVE_flag,.before = QC_combined) %>% 
  arrange(test_sessions.bblid,Test) %>% 
  filter(!(QC_combined %in% c("All valid","Both valid","All flagged","Both flagged"))) 


write_csv(qc_disagreement,file = "/Users/hillmann/Projects/22q/Data/QA/CNB/qc_disagreement.csv")


