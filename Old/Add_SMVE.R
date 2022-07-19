library(tidyverse)
library(PerFit)
select <- dplyr::select

Comment_flags <- read_csv("~/Projects/22q/Data/QA/CNB/cnb_22q_cross_qc_05_12_2022_complete.csv", 
                skip = 1)
athena_3360_2096 <- read_csv("~/Projects/22q/Data/itemwise/athena_3360_2096.csv")
athena_254_360 <- read_csv("~/Projects/22q/Data/itemwise/athena_254_360.csv")
CNB_order <- read_csv("/Users/hillmann/Projects/22q/Data/QA/CNB/22q_cnb_batteries_order.csv")
CNB <- read_csv("/Users/hillmann/Projects/22q/Data/Summary/cnb_22q_cross_qc_05_12_2022.csv")

Comment_flags_narrow <- Comment_flags %>% 
  select(test_sessions.bblid,test_sessions.datasetid,Notes,contains("(1 = yes, no = 2)"),contains("Overall_Valid")) %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid))

CNB <- CNB %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
  left_join(Comment_flags_narrow) 

# Clean CNB comment flags
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = " \\(1 = yes, no = 2\\)",replacement = "")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = " \\(1 = yes, 2 = no\\)",replacement = "")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "_flag",replacement = "_AV_flag")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "_valid",replacement = "_comment_flag")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "Overall_Valid",replacement = "Overall_valid")

# Remove individuals who are flagged on the overall session indicator
CNB <- CNB %>% 
  mutate(across(.cols = c(contains("comment_flag"),Overall_valid),.fns = ~ case_when(.x == 1 ~ "V",.x == 2 ~ "F",TRUE ~ "V"))) %>% 
  filter(Overall_valid != "F")
    
# First, find SVME for each test in the data set
  
add_PFscores <- function(test,test_grep,test_type,output_df){
  
  if(test == "PMAT"){
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) 
    
    athena_254_360_narrow <- athena_254_360 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) 
    
    itemwise_data <- rbind(athena_3360_2096_narrow,athena_254_360_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      semi_join(CNB,by = c("test_sessions.bblid","test_sessions_v.datasetid" = "test_sessions.datasetid")) %>% 
      arrange(test_sessions.bblid) %>% 
      select(!matches("^PMAT24_A.PMAT24_A_QID001")) %>% 
      filter(!if_all(.cols = matches(test_grep),.fns = ~ is.na(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x))) %>% 
      mutate(across(.cols = matches("TTR$"),.fns = ~ ifelse(.x < 0,NA,.x)))
    
  } else if(test == "SLNB"){ 
    
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR"),matches("_TTR")) %>% 
      select(!matches("8[[:digit:]][[:digit:]][[:digit:]]$")) 
    
    colnames(athena_3360_2096_narrow) <- str_replace_all(colnames(athena_3360_2096_narrow),pattern = "\\.\\.\\.[[:digit:]]{1,}$",replacement = "")
    
    athena_254_360_narrow <- athena_254_360 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR"),matches("_TTR")) %>% 
      select(!matches("8[[:digit:]][[:digit:]][[:digit:]]$"))
    
    itemwise_data <- rbind(athena_254_360_narrow,athena_3360_2096_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      semi_join(CNB,by = c("test_sessions.bblid","test_sessions_v.datasetid" = "test_sessions.datasetid")) %>% 
      arrange(test_sessions.bblid) %>% 
      filter(!if_all(.cols = matches(test_grep),.fns = ~ is.na(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x))) 
       
    
  }else if(test == "VSPLOT"){
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_RT_[[:digit:]]")) %>% 
      filter(if_any(.cols = matches(test_grep),.fns = ~ !is.na(.x))) 
    
    athena_254_360_narrow <- athena_254_360 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_RT_[[:digit:]]")) %>% 
      filter(if_any(.cols = matches(test_grep),.fns = ~ !is.na(.x))) 
    
    itemwise_data <- rbind(athena_3360_2096_narrow,athena_254_360_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      semi_join(CNB,by = c("test_sessions.bblid","test_sessions_v.datasetid" = "test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x)))
    
  } else if(test == "SPCPTN"){
    
    athena_3360_2096_narrow_both <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("SPCPTNL")) %>% 
      select(!matches("19|2[[:digit:]]|3[[:digit:]]")) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) %>% 
      filter(if_any(.cols = matches("SPCPTNL"),.fns = ~ !is.na(.x))) 
    
    athena_3360_2096_narrow_N_only <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("SPCPTN90")) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) %>% 
      filter(if_any(.cols = matches("SPCPTN90"),.fns = ~ !is.na(.x))) 
    
    colnames(athena_3360_2096_narrow_both) <- str_replace_all(colnames(athena_3360_2096_narrow_both),pattern = "^SPCPTNL",replacement = "SPCPTN90")
    
    athena_254_360_narrow <- athena_254_360 %>%
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("SPCPTN90")) %>%
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) %>%
      filter(if_any(.cols = matches("SPCPTN90"),.fns = ~ !is.na(.x)))
    
    itemwise_data <- rbind(athena_3360_2096_narrow_both,athena_3360_2096_narrow_N_only,athena_254_360_narrow) %>% 
      mutate(across(.cols = matches("TTR$"),.fns = ~ ifelse(.x < 0,NA,.x))) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      semi_join(CNB,by = c("test_sessions.bblid","test_sessions_v.datasetid" = "test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x)))
    
  }
    else{
    
    athena_3360_2096_narrow <- athena_3360_2096 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) %>% 
      filter(if_any(.cols = matches(test_grep),.fns = ~ !is.na(.x))) 
    
    colnames(athena_3360_2096_narrow) <- str_replace_all(colnames(athena_3360_2096_narrow),pattern = "TRIAL000001_0",replacement = "TRIAL000001")
    
    athena_254_360_narrow <- athena_254_360 %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches(test_grep)) %>% 
      select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.datasetid,test_sessions_v.dotest,matches("_CORR$"),matches("_TTR$")) %>% 
      filter(if_any(.cols = matches(test_grep),.fns = ~ !is.na(.x))) 
    
    itemwise_data <- rbind(athena_3360_2096_narrow,athena_254_360_narrow) %>% 
      mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
      semi_join(CNB,by = c("test_sessions.bblid","test_sessions_v.datasetid" = "test_sessions.datasetid")) %>% 
      mutate(across(.cols = contains("CORR"),.fns = ~ as.numeric(.x))) %>% 
      mutate(across(.cols = matches("_CORR"),.fns = ~ ifelse(is.na(.x),0,.x)))
    
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
  pfit1_vec <- ifelse(rowSums(dat2) == ncol(dat2),1,pfit1 %>% pull())
  pfit1_final <- ifelse(rowSums(dat2) == 0,0,pfit1_vec)
  pfit2 <- E.KB(dat2)$PFscores
  pfit2_vec <- ifelse(rowSums(dat2) == ncol(dat2),1,pfit2 %>% pull())
  pfit2_final <- ifelse(rowSums(dat2) == 0,0,pfit2_vec)
  acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))][,1:3],na.rm = T)
  if(test_type == "memory"){
    sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1_final) + (0.50*pfit2_final)
  } else if(test_type == "non-memory"){
    sc <- (0.34*outlier_score_2cut) + (0.22*pfit1_final) + (0.44*pfit2_final)
  }
  
  col_name <- paste0("PFscores_",test)
  itemwise_data[,col_name] <- sc
  if(test_type == "memory"){
    itemwise_data[,paste0(col_name,"_res")] <- 0.42*outlier_score_2cut
    itemwise_data[,paste0(col_name,"_acc")] <- 0.02*acc3e
    itemwise_data[,paste0(col_name,"_pfit1")] <- 0.05*pfit1_final
    itemwise_data[,paste0(col_name,"_pfit2")] <- 0.50*pfit2_final
  } else{
    itemwise_data[,paste0(col_name,"_res")] <- 0.34*outlier_score_2cut
    itemwise_data[,paste0(col_name,"_acc")] <- 0
    itemwise_data[,paste0(col_name,"_pfit1")] <- 0.22*pfit1_final
    itemwise_data[,paste0(col_name,"_pfit2")] <- 0.44*pfit2_final
  }
  
  output_df <- output_df %>% 
    left_join(itemwise_data[,c("test_sessions.bblid","test_sessions_v.datasetid",paste(col_name),paste0(col_name,"_res"),paste0(col_name,"_acc"),paste0(col_name,"_pfit1"),paste0(col_name,"_pfit2"))], by = c("test_sessions.bblid","test_sessions.datasetid" = "test_sessions_v.datasetid"))
  
  return(output_df)
}

tests <- c("CPF","SVOLT","PMAT","ER40","MEDF","SLNB","ADT","VSPLOT","SPCPTN")
test_type <- c("memory","memory","non-memory","non-memory","non-memory","memory","non-memory","non-memory","non-memory")
test_grep <- c("^CPF_B.CPF_TRIAL","^SVOLT_A.SVOLT_TRIAL","^PMAT24_A.PMAT24_A_QID","^ER40_D.ER40D_QID","^MEDF36_A.MEDF36A_QID","^SLNB2_90.SLNB2_QID","^ADT36_A.ADT36A_QID","^VSPLOT15.VSPLOT15","^SPCPTN")

CNB_with_SMVE <- CNB

for(i in 1:length(tests)){
  CNB_with_SMVE <- add_PFscores(test = tests[i],test_grep = test_grep[i],test_type = test_type[i],output_df = CNB_with_SMVE)
}

CNB_with_SMVE <- CNB_with_SMVE %>% 
  mutate(across(.cols = c(PFscores_CPF,PFscores_SVOLT,PFscores_PMAT,PFscores_ER40,PFscores_MEDF,PFscores_SLNB,PFscores_ADT,PFscores_VSPLOT,PFscores_SPCPTN),.fns = ~ case_when(.x < quantile(.x,.05,na.rm = T) ~ "F",TRUE ~ "V"),.names = "{.col}_flag")) 

write_csv(CNB_with_SMVE,file = "/Users/hillmann/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")




