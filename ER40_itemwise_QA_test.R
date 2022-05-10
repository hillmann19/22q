# Read in packages and Labels
library(PerFit)
library(gtsummary)
library(tidyverse)
athena_254_360 <- read_csv("~/Projects/22q/Data/itemwise/athena_254_360.csv")
athena_3360_2096 <- read_csv("~/Projects/22q/Data/itemwise/athena_3360_2096.csv")
CNB_22q <- read_csv("~/Projects/22q/Data/Summary/cnb_22q_202109.csv")
CNB_all <- read_csv("~/Projects/Evolution/Data/CNB/cnb_merged_webcnp_surveys_allbblprjcts_longform20220406.csv")
diag <- read_csv("~/Projects/22q/Data/Summary/diagnosis_wsmryvars_20220425.csv")


# CNB - get cross-sectional Data

CNB_under35 <- CNB_22q %>% 
  filter(test_sessions_v.age <= 35) %>%
  filter(!if_all(.cols = c("ADT36_A.ADT36A_CR","ADT36_A.ADT36A_PC","ADT36_A.ADT36A_RTCR","CPF_B.CPF_CR","CPF_B.CPF_RTCR","CPF_B.CPF_W_RTCR","ER40_D.ER40D_CR","ER40_D.ER40D_RTCR","MEDF36_A.MEDF36A_CR","MEDF36_A.MEDF36A_RTCR","MPRACT.MP2RTCR","PCET_A.PCET_RTCR","PCET_A.PCET_CAT","PCET_A.PCET_ACC2","PMAT24_A.PMAT24_A_CR","PMAT24_A.PMAT24_A_RTCR","SCTAP.SCTAP_TOT","SLNB2_90.SLNB2_MCR","SLNB2_90.SLNB2_MRTC","SPCPTN90.SCPN90_TP","SPCPTN90.SCPN90_TPRT","SPCPTNL.SCPN_TPRT","SPCPTNL.SCPN_TP","SVOLT_A.SVOLT_RTCR", "VSPLOT15.VSPLOT15_CR","VSPLOT15.VSPLOT15_RTCR"),.fns = ~ is.na(.x))) %>% 
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "^([[:digit:]])/",replacement = "0\\1/")) %>% # pad months with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "/([[:digit:]])/",replacement = "/0\\1/")) %>% #pad days with 0
  mutate(test_sessions_v.dotest = str_replace_all(test_sessions_v.dotest,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1")) %>% # Add 20 to year (15 becomes 2015)
  mutate(test_sessions_v.dotest = as.Date(test_sessions_v.dotest,format = "%m/%d/%Y")) %>% 
  group_by(test_sessions.bblid) %>% 
  arrange(test_sessions_v.dotest) %>% 
  mutate(test_num = row_number()) %>% 
  ungroup()

CNB_repeats_list <- CNB_under35 %>% 
  group_by(test_sessions.bblid) %>% 
  filter(n() > 1) %>% 
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
  ungroup() %>% 
  bind_rows(Last_tests)  %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid))

# Cap values at 6 sd 

# tests <- CNB_cross %>% 
#   dplyr::select(!(matches("^test") | matches("valid_code") | matches("_AR$") | "remote" | "gender_remote" | "deleted_sample" | "cnbagegrp" | "platform")) %>% 
#   colnames()
# for(test in tests){
#   CNB_cross[[test]] <- ifelse(CNB_cross[[test]] > mean(CNB_cross[[test]],na.rm = TRUE) + 6*sd(CNB_cross[[test]],na.rm = TRUE),mean(CNB_cross[[test]],na.rm = TRUE) + 6*sd(CNB_cross[[test]],na.rm = TRUE),CNB_cross[[test]])
#   CNB_cross[[test]] <- ifelse(CNB_cross[[test]] < mean(CNB_cross[[test]],na.rm = TRUE) - 6*sd(CNB_cross[[test]],na.rm = TRUE),mean(CNB_cross[[test]],na.rm = TRUE) - 6*sd(CNB_cross[[test]],na.rm = TRUE),CNB_cross[[test]])
# }

# Look into valid codes
# CNB_cross %>% 
#   select(contains("valid_code"),remote) %>% 
#   mutate(across(.cols = !contains("remote"),.fns = ~ as.character(.x))) %>% 
#   pivot_longer(cols = !contains("remote"),names_to = "Test",values_to = "Code") %>% 
#   mutate(Test = str_replace_all(Test,pattern = "_.*","")) %>% 
#   mutate(Test = str_replace_all(Test,pattern = "\\.valid$","")) %>% 
#   filter(!is.na(Code)) %>% 
#   tbl_summary(by = Test)

athena_3360_2096_er40 <- athena_3360_2096 %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,contains("ER40_D.ER40D_QID")) %>% 
  select(!matches("TRIAL$")) %>% 
  select(!matches("QID$")) %>% 
  filter(if_all(.cols = contains("ER40_D.ER40D_QID"),.fns = ~ !is.na(.x))) 

athena_254_360_er40 <- athena_254_360 %>% 
  select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,contains("ER40_D.ER40D_QID")) %>% 
  select(!matches("TRIAL$")) %>% 
  select(!matches("QID$")) %>% 
  filter(if_all(.cols = contains("ER40_D.ER40D_QID"),.fns = ~ !is.na(.x)))

itemwise_er40 <- rbind(athena_3360_2096_er40,athena_254_360_er40) %>% 
  filter(if_all(.cols = matches("TTR$"),.fns = ~ .x > 0)) 

# Run person statistic on all data

dat <- itemwise_er40[,c(grep("_CORR",colnames(itemwise_er40)),grep("_TTR",colnames(itemwise_er40)))]
dat <- as.data.frame(dat)
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])

res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
  mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
  res[,j] <- scale(residuals(mod,na.action=na.exclude))
}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)

itemwise_er40 <- data.frame(itemwise_er40,"PFscores_AllData" = as.numeric(sc$PFscores))

### Run person statistic on 22q data 

itemwise_er40_22q <- itemwise_er40 %>% 
  semi_join(CNB_cross,by = c("test_sessions.bblid","test_sessions_v.dotest"))

dat <- itemwise_er40_22q[,c(grep("_CORR",colnames(itemwise_er40_22q)),grep("_TTR",colnames(itemwise_er40_22q)))]
dat <- as.data.frame(dat)
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])

res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
  mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
  res[,j] <- scale(residuals(mod,na.action=na.exclude))
}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)

itemwise_er40_22q <- data.frame(itemwise_er40_22q,"PFscores_22q" = as.numeric(sc$PFscores))
  

itemwise_er40 <- itemwise_er40 %>% 
  left_join(itemwise_er40_22q[,c("test_sessions.bblid","test_sessions_v.dotest","PFscores_22q")])

# Add diagnosis to itemwise data 

diag <- diag %>% 
  select(BBLID,DODIAGNOSIS,dx_pscat) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "JAN",replacement = "01")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "FEB",replacement = "02")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "MAR",replacement = "03")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "APR",replacement = "04")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "MAY",replacement = "05")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "JUN",replacement = "06")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "JUL",replacement = "07")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "AUG",replacement = "08")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "SEP",replacement = "09")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "OCT",replacement = "10")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "NOV",replacement = "11")) %>% 
  mutate(DODIAGNOSIS = str_replace_all(DODIAGNOSIS,pattern = "DEC",replacement = "12")) %>% 
  mutate(Year = str_extract(DODIAGNOSIS,pattern = "[[:digit:]][[:digit:]]$")) %>% 
  mutate(DODIAGNOSIS = ifelse(as.numeric(Year) <= 22,str_replace_all(DODIAGNOSIS,pattern = "([[:digit:]][[:digit:]])$",replacement = "20\\1"),str_replace_all(DODIAGNOSIS,pattern = "([[:digit:]][[:digit:]])$",replacement = "19\\1"))) %>% 
  mutate(DODIAGNOSIS = as.Date(DODIAGNOSIS,format = "%d-%m-%Y")) %>% 
  mutate(BBLID = as.character(BBLID)) %>% 
  select(-Year)

itemwise_er40_l <- itemwise_er40 %>% 
  arrange(test_sessions.bblid,test_sessions_v.dotest) %>% 
  group_split(test_sessions.bblid)

Link_CNB_diag <- function(df){
  bblid <- unique(df$test_sessions.bblid)
  CNB_dates <- unique(df$test_sessions_v.dotest)
  diag_dates <- diag %>% 
    filter(BBLID == bblid) %>% 
    pull(DODIAGNOSIS)
  final_df <- data.frame()
  
  if(length(CNB_dates) == 0|length(diag_dates) == 0){
    final_df <- df
    final_df$DODIAGNOSIS <- NA
    final_df$dx_pscat <- NA
  } else{
    for(i in 1:length(CNB_dates)){
      cnb_date <- CNB_dates[i]
      min_diff <- min(abs(as.numeric(cnb_date - diag_dates)))
      diag_min <- diag_dates[which(abs(as.numeric(cnb_date - diag_dates)) == min_diff)]
      
      new_df <- df %>% 
        filter(test_sessions_v.dotest == cnb_date)
      new_df$DODIAGNOSIS <- diag_min
      new_df$dx_pscat <- diag %>% 
        filter(BBLID == bblid) %>% 
        filter(DODIAGNOSIS == diag_min) %>% 
        pull(dx_pscat)
      
      final_df <- rbind(final_df,new_df)
    }
  }
  
return(final_df)
}

itemwise_er40 <- map_dfr(itemwise_er40_l,Link_CNB_diag)

itemwise_er40 %>% 
  select(test_sessions.bblid,test_sessions_v.dotest,PFscores_22q) %>% 
  filter(!is.na(PFscores_22q)) %>% 
  filter(PFscores_22q <= mean(PFscores_22q,na.rm = T) - 2*sd(PFscores_22q,na.rm = T))


# Create summary plots and tables comparing SMVM for 22q and all data

PF_summary <- itemwise_er40 %>% 
  select(PFscores_AllData,PFscores_22q) %>% 
  pivot_longer(cols = everything(),names_to = "Sample",values_to = "Scott-Moore Validity Metric") %>% 
  mutate(Sample = str_replace_all(Sample,pattern = ".*_",replacement = "")) %>% 
  mutate(Sample = ifelse(Sample == "AllData","All Data",Sample)) %>% 
  group_by(Sample) %>% 
  summarize(Avg = mean(`Scott-Moore Validity Metric`,na.rm = T),SD = sd(`Scott-Moore Validity Metric`,na.rm = T))

theme_set(theme_minimal(base_size = 18))

itemwise_er40 %>% 
  select(PFscores_AllData,PFscores_22q) %>% 
  pivot_longer(cols = everything(),names_to = "Sample",values_to = "Scott-Moore Validity Metric") %>% 
  mutate(Sample = str_replace_all(Sample,pattern = ".*_",replacement = "")) %>% 
  mutate(Sample = ifelse(Sample == "AllData","All Data",Sample)) %>% 
  ggplot(aes(x = `Scott-Moore Validity Metric`,fill = Sample)) + geom_density(alpha = .5) + scale_fill_manual(values = c("#011F5b","#990000")) + scale_color_manual(values = c("#011F5b","#990000"))  + geom_vline(data = PF_summary,aes(xintercept = Avg,color = Sample)) + geom_vline(data = PF_summary,aes(xintercept = Avg - 2*SD,color = Sample),lty = 2) + guides(color = "none")

itemwise_er40 %>% 
  select(PFscores_AllData,PFscores_22q) %>% 
  pivot_longer(cols = everything(),names_to = "Sample",values_to = "Scott-Moore Validity Metric") %>% 
  mutate(Sample = str_replace_all(Sample,pattern = ".*_",replacement = "")) %>% 
  mutate(Sample = ifelse(Sample == "AllData","All Data",Sample)) %>% 
  group_by(Sample) %>% 
  mutate(remote = ifelse(`Scott-Moore Validity Metric` < mean(`Scott-Moore Validity Metric`,na.rm = T) - 2*sd(`Scott-Moore Validity Metric`,na.rm = T),"Remove","Keep")) %>% 
  ungroup() %>% 
  with(prop.table(table(remote,Sample),margin = 2))


CNB_all_clean <- CNB_all %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  filter(!is.na(test_sessions.bblid),!is.na(test_sessions_v.dotest),!is.na(test_sessions_v.gender)) %>% 
  filter(!(test_sessions.bblid == "20427" & test_sessions_v.age == 18))

## 22q remove characteristics
itemwise_er40 %>% 
  left_join(CNB_all_clean[,c("test_sessions.bblid","test_sessions_v.gender","test_sessions_v.dotest")]) %>% 
  distinct(.keep_all = T) %>% 
  filter(!is.na(test_sessions.bblid),!is.na(test_sessions_v.gender),!is.na(dx_pscat)) %>%
  pivot_longer(cols = c(PFscores_AllData,PFscores_22q),names_to = "Sample",values_to = "Scott-Moore Validity Metric") %>% 
  mutate(Sample = str_replace_all(Sample,pattern = ".*_",replacement = "")) %>% 
  mutate(Sample = ifelse(Sample == "AllData","All Data",Sample)) %>% 
  group_by(Sample) %>% 
  mutate(remove = ifelse(`Scott-Moore Validity Metric` < mean(`Scott-Moore Validity Metric`,na.rm = T) - 2*sd(`Scott-Moore Validity Metric`,na.rm = T),"Remove","Keep")) %>% 
  ungroup() %>% 
  filter(Sample == "22q") %>% 
  select(test_sessions_v.gender,test_sessions_v.age,dx_pscat,remove) %>% 
  rename(Sex = test_sessions_v.gender,Age = test_sessions_v.age,Diagnosis = dx_pscat) %>% 
  mutate(Sex = case_when(Sex == 'M' ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>% 
  tbl_summary(by = remove)
  
# All data remove characteristics 

itemwise_er40 %>% 
  left_join(CNB_all_clean[,c("test_sessions.bblid","test_sessions_v.gender","test_sessions_v.dotest")]) %>% 
  distinct(.keep_all = T) %>% 
  filter(!is.na(test_sessions.bblid),!is.na(test_sessions_v.gender),!is.na(dx_pscat)) %>%
  pivot_longer(cols = c(PFscores_AllData,PFscores_22q),names_to = "Sample",values_to = "Scott-Moore Validity Metric") %>% 
  mutate(Sample = str_replace_all(Sample,pattern = ".*_",replacement = "")) %>% 
  mutate(Sample = ifelse(Sample == "AllData","All Data",Sample)) %>% 
  group_by(Sample) %>% 
  mutate(remove = ifelse(`Scott-Moore Validity Metric` < mean(`Scott-Moore Validity Metric`,na.rm = T) - 2*sd(`Scott-Moore Validity Metric`,na.rm = T),"Remove","Keep")) %>% 
  ungroup() %>% 
  filter(Sample == "All Data") %>% 
  select(test_sessions_v.gender,test_sessions_v.age,dx_pscat,remove) %>% 
  rename(Sex = test_sessions_v.gender,Age = test_sessions_v.age,Diagnosis = dx_pscat) %>% 
  mutate(Sex = case_when(Sex == 'M' ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>% 
  tbl_summary(by = remove)




# Remove bottom 5% from data 
# x <- ADT30_1[,grepl("_CORR",colnames(ADT30_1)) | grepl("PFscores",colnames(ADT30_1))]
# qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
# x <- x[which(x$PFscores > qu),]
# x <- x[,colnames(x) != "PFscores"]



