# Load in necessary data and packages
library(ggpubr)
library(tidyverse)

CNB <- read_csv("~/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")

CNB_clean <- CNB %>% 
  select(test_sessions.bblid,test_sessions_v.gender,test_sessions_v.age,remote,test_num,SPCPTN90_test,matches("RTCR$"),matches("MRTC$"),matches("TPRT$"),matches("_CR$"),matches("SCTAP_TOT$"),matches("ACC2$"),matches("PTP$"),matches("MCR$"),matches("CAT$"),matches("TP$"),matches("_flag$")) %>% 
  select(-CPF_B.CPF_W_RTCR,-SPCPTN90.SCPN90_TPRT,-SPCPTNL.SCPN_TPRT,-SPCPTN90.SCPN90_TP,-SPCPTNL.SCPN_TP) %>% 
  rowwise() %>% 
  mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,ADT36_A.ADT36A_CR),.fns = ~ ifelse(sum(c(PFscores_ADT_flag,ADT_comment_flag,ADT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(CPF_B.CPF_RTCR,CPF_B.CPF_CR),.fns = ~ ifelse(sum(c(PFscores_CPF_flag,CPF_comment_flag,CPF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(ER40_D.ER40D_RTCR,ER40_D.ER40D_CR),.fns = ~ ifelse(sum(c(PFscores_ER40_flag,ER40_comment_flag,ER40_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MEDF36_A.MEDF36A_RTCR,MEDF36_A.MEDF36A_CR),.fns = ~ ifelse(sum(c(PFscores_MEDF_flag,MEDF_comment_flag,MEDF_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(MPRACT.MP2RTCR),.fns = ~ ifelse(sum(c(MPRACT_comment_flag,MPRACT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(PCET_A.PCET_RTCR,PCET_A.PCET_ACC2,PCET_A.PCET_CAT),.fns = ~ ifelse(sum(c(PCET_comment_flag,PCET_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(PMAT24_A.PMAT24_A_RTCR,PMAT24_A.PMAT24_A_CR),.fns = ~ ifelse(sum(c(PFscores_PMAT_flag,PMAT_comment_flag,PMAT_AV_flag) == "F") >= 2,NA,.x))) %>% 
  mutate(across(.cols = c(SVOLT_A.SVOLT_RTCR),.fns = ~ ifelse(sum(c(PFscores_SVOLT_flag,SVOLT_comment_flag,SVOLT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(VSPLOT15.VSPLOT15_RTCR,VSPLOT15.VSPLOT15_CR),.fns = ~ ifelse(sum(c(PFscores_VSPLOT_flag,VSPLOT_comment_flag,VSPLOT_AV_flag) == "F") >= 2,NA,.x))) %>%
  mutate(across(.cols = c(SLNB2_90.SLNB2_MRTC,SLNB2_90.SLNB2_MCR),.fns = ~ ifelse(sum(c(PFscores_SLNB_flag,SLNB_comment_flag,SLNB_AV_flag) == "F") >= 2,NA,.x))) %>%
  #mutate(across(.cols = c(SPCPTN90_TPRT,SPCPTN90_TP),.fns = ~ ifelse(sum(c(SPCPTN_comment_flag,SPCPTN_AV_flag,PFscores_SPCPTN_flag) == "F") >= 2,NA,.x))) %>% remove # once PFscores qc is complete
  mutate(across(.cols = c(SCTAP.SCTAP_TOT),.fns = ~ ifelse(sum(c(SCTAP_comment_flag,SCTAP_AV_flag) == "F") >= 2,NA,.x))) %>% 
  ungroup() %>% 
  select(!matches("_flag$"))


CNB_z <- CNB_clean %>% 
  mutate(across(.cols = c(ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR),.fns = ~ -1*.x)) %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,remote,test_num,ADT36_A.ADT36A_CR,ADT36_A.ADT36A_RTCR,CPF_B.CPF_CR,CPF_B.CPF_RTCR,ER40_D.ER40D_CR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_CR,MEDF36_A.MEDF36A_RTCR,MPRACT.MP2RTCR,PCET_A.PCET_RTCR,PCET_A.PCET_ACC2,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_CR,PMAT24_A.PMAT24_A_RTCR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TP,SPCPTN90_TPRT,SVOLT_A.SVOLT_CR,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR,VSPLOT15.VSPLOT15_RTCR) %>% 
  mutate(across(.cols = ADT36_A.ADT36A_CR:last_col(),.fns = ~ as.numeric(scale(.x)))) %>% 
  rowwise() %>% 
  mutate(overall_acc_z = mean(c_across(cols = c(ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SLNB2_90.SLNB2_MCR,SPCPTN90_TP,SVOLT_A.SVOLT_CR,VSPLOT15.VSPLOT15_CR)),na.rm = T)) %>% 
  mutate(overall_speed_z = mean(c_across(cols = c(ADT36_A.ADT36A_RTCR,CPF_B.CPF_RTCR,ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,PCET_A.PCET_RTCR,PMAT24_A.PMAT24_A_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_RTCR)),na.rm = T)) %>% 
  mutate(exec_acc_z = mean(c_across(cols = c(PCET_A.PCET_ACC2,SLNB2_90.SLNB2_MCR,SPCPTN90_TP)),na.rm = T)) %>% 
  mutate(mem_acc_z = mean(c_across(cols = c(CPF_B.CPF_CR,SVOLT_A.SVOLT_CR)),na.rm = T)) %>% 
  mutate(complex_acc_z = mean(c_across(cols = c(VSPLOT15.VSPLOT15_CR,PMAT24_A.PMAT24_A_CR)),na.rm = T)) %>% 
  mutate(social_acc_z = mean(c_across(cols = c(ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,ADT36_A.ADT36A_CR)),na.rm = T)) %>% 
  mutate(overall_acc_z = mean(c_across(cols = c(ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SLNB2_90.SLNB2_MCR,SPCPTN90_TP,SVOLT_A.SVOLT_CR,VSPLOT15.VSPLOT15_CR)),na.rm = T)) %>% 
  mutate(exec_speed_z = mean(c_across(cols = c(PCET_A.PCET_RTCR,SLNB2_90.SLNB2_MRTC,SPCPTN90_TPRT)),na.rm = T)) %>% 
  mutate(mem_speed_z = mean(c_across(cols = c(CPF_B.CPF_RTCR,SVOLT_A.SVOLT_RTCR)),na.rm = T)) %>% 
  mutate(complex_speed_z = mean(c_across(cols = c(VSPLOT15.VSPLOT15_RTCR,PMAT24_A.PMAT24_A_RTCR)),na.rm = T)) %>% 
  mutate(social_speed_z = mean(c_across(cols = c(ER40_D.ER40D_RTCR,MEDF36_A.MEDF36A_RTCR,ADT36_A.ADT36A_RTCR)),na.rm = T)) %>% 
  mutate(motor_z = mean(c_across(cols = c(MPRACT.MP2RTCR,SCTAP.SCTAP_TOT)),na.rm = T)) %>% 
  ungroup() %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,remote,test_num,matches("_z$"))

theme_set(theme_minimal())
theme_update(text = element_text(size = 16),legend.position = "bottom")

acc_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = overall_acc_z,color = remote,fill = remote)) + geom_point(size = .75) + geom_smooth(method = "gam") + labs(x = "Age",y = "Overall Accuracy \n (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

speed_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = overall_speed_z,color = remote,fill = remote)) + geom_point(size = .75) + geom_smooth(method = "gam") + labs(x = "Age",y = "Overall Speed \n (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

motor_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = motor_z,color = remote,fill = remote)) + geom_point(size = .75) + geom_smooth(method = "gam") + labs(x = "Age",y = "Motor Speed \n (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

ggarrange(acc_plot,speed_plot,common.legend = T,labels = c("Accuracy","Speed"),legend = "bottom",vjust = 1)
ggarrange(acc_plot,motor_plot,common.legend = T,labels = c("Accuracy","Motor Speed"),legend = "bottom")

acc_alone_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = overall_acc_z,color = remote,fill = remote)) + geom_point(size = .75) + geom_smooth(method = "gam") + labs(x = "Age",y = "Accuracy \n (Z-score)",color = "Test Location",title = "Overall Accuracy") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")


EF_acc_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = exec_acc_z,color = remote,fill = remote)) + geom_point() + geom_smooth(method = "gam") + labs(x = "Age",y = "Accuracy (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

mem_acc_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = mem_acc_z,color = remote,fill = remote)) + geom_point() + geom_smooth(method = "gam") + labs(x = "Age",y = "Accuracy (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

complex_acc_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = complex_acc_z,color = remote,fill = remote)) + geom_point() + geom_smooth(method = "gam") + labs(x = "Age",y = "Accuracy (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

social_acc_plot <- CNB_z %>% 
  ggplot(aes(x = test_sessions_v.age,y = social_acc_z,color = remote,fill = remote)) + geom_point() + geom_smooth(method = "gam") + labs(x = "Age",y = "Accuracy (Z-score)",color = "Test Location") + scale_color_manual(values = c("#990000","#011F5B")) + scale_fill_manual(values = c("#990000","#011F5B")) + guides(fill = "none")

ggarrange(EF_acc_plot,mem_acc_plot,complex_acc_plot,social_acc_plot,common.legend = T,labels = c("Executive Functioning","Memory","Complex Cognition","Social Cognition"),legend = "bottom")


