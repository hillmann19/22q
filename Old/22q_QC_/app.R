#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(shiny)
library(gtsummary)
library(gt)
library(tidytext)
library(ggpubr)
df <- read_csv("~/Projects/22q/Data/QA/CNB/CNB_with_SMVE.csv")
CNB <- read_csv("~/Projects/22q/Data/QA/CNB/cnb_22q_cross_qc_05_12_2022_complete.csv", 
                skip = 1)
CNB_order <- read_csv("/Users/hillmann/Projects/22q/Data/QA/CNB/22q_cnb_batteries_order.csv")


# Clean CNB comment flags
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = " \\(1 = yes, no = 2\\)",replacement = "")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = " \\(1 = yes, 2 = no\\)",replacement = "")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "_flag",replacement = "_AV_flag")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "_valid",replacement = "_comment_flag")
colnames(CNB) <- str_replace_all(colnames(CNB),pattern = "Overall_Valid",replacement = "Overall_valid")

CNB$test_sessions.bblid <- as.character(CNB$test_sessions.bblid)
df$test_sessions.bblid <- as.character(df$test_sessions.bblid)

CNB_completed <- CNB %>%
  select(test_sessions.bblid,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR) %>%
  pivot_longer(cols = ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Response") %>%
  mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>%
  mutate(Test = str_replace_all(Test,pattern = "36|24|15|90|2",replacement = "")) %>%
  mutate(Completed = case_when(!is.na(Response) ~ 1,is.na(Response) ~ 0)) %>% 
  select(-Response) %>% 
  mutate(test_sessions.bblid = as.character(test_sessions.bblid)) 

df <- df %>% 
  mutate(Overall_valid = ifelse(is.na(Overall_valid),1,Overall_valid))

CNB_order_clean <- CNB_order %>% 
  filter(BatteryName %in% unique(CNB$test_sessions_v.battery)) %>% 
  mutate(Test = case_when(str_detect(TestName,pattern = "adt") ~ "ADT",str_detect(TestName,pattern = "cpf") ~ "CPF",str_detect(TestName,pattern = "er40") ~ "ER40",str_detect(TestName,pattern = "medf") ~ "MEDF",str_detect(TestName,pattern = "mpraxis") ~ "MPRACT",str_detect(TestName,pattern = "pcet") ~ "PCET",str_detect(TestName,pattern = "pmat") ~ "PMAT",str_detect(TestName,pattern = "sctap") ~ "SCTAP",str_detect(TestName,pattern = "slnb") ~ "SLNB",str_detect(TestName,pattern = "spcptn90") ~ "SPCPTN",str_detect(TestName,pattern = "svolt") ~ "SVOLT",str_detect(TestName,pattern = "vsplot") ~ "VSPLOT",TRUE ~ NA_character_)) %>% 
  select(-TestName) %>% 
  filter(!is.na(Test)) %>% 
  group_by(BatteryName) %>% 
  arrange(BatteryName,TestOrder) %>% 
  mutate(Test_Order_relative = row_number()) %>% 
  ungroup() 

theme_set(theme_minimal())
theme_update(axis.text.x = element_text(size = 6.5,angle = 45),legend.position = "none")

# Define UI for application that draws a histogram
ui <- navbarPage(title = "22q Performance Validity QC",

        # Show a plot of the generated distribution
        tabPanel("Distribution Plots",sidebarPanel(width = 2,position = "left",
          sliderInput("cutoff",
                      "% Removed",
                      min = 1,
                      max = 25,
                      value = 5),checkboxInput("removeOverallbad","Remove individuals with entire sessions of poor quality")),sidebarPanel(width = 10,position = "right",plotOutput("distPlot",height = "300px")),mainPanel(width = 12,gt::gt_output("RemovedTable"))),
          tabPanel("Itemwise Flags: A Closer Look",sidebarPanel(width = 2,position = "left",
                                                                selectInput("test",
                                                                            label = "Test",
                                                                            choices = c("All","ADT","CPF","ER40","MEDF","PMAT","SLNB","SVOLT"),selected = "All"),checkboxInput("cumulative","Cumulative plots"),uiOutput("dynamic")),
                                                  fillPage(tags$style(type = "text/css", "#itemwiseFlags {height: calc(82vh - 75px) !important;}"),plotOutput("itemwiseFlags",hover = "plot_hover",width = "95%",height = "95%"))),
        tabPanel("Test Fatigue",sidebarPanel(width = 2,position = "left",sliderInput("dropoff","Dropoff criteria",min = 0,max = 5,value = 2,step = .1),checkboxInput("byTest","Show flagged by Test")),fillPage(tags$style(type = "text/css", "#Fatigue {height: calc(80vh - 75px) !important;}"),plotOutput("Fatigue",width = "95%",height = "80%"))),
        tabPanel("Comparing QC Methods",sidebarPanel(width = 2,position = "left",selectInput("method","QC method to nest by",choices = c("Auto-validation","Reviewer","Multivariate CNB Validity Estimate","Score Dropoff"))), fillPage(tags$style(type = "text/css", "#all_table {height: calc(90vh - 75px) !important;}"),gt::gt_output("all_table")))
        )
    


# Define server logic required to draw a histogram
server <- function(input, output) {
  
    output$distPlot <- renderPlot({
      theme_set(theme_minimal())
      theme_update(text = element_text(size = 16),legend.position = "none")
      
      if(input$removeOverallbad){
        df <- df %>% 
          filter(Overall_valid == 1)
      }
  
        df %>% 
          select(test_sessions.bblid,contains("PFscores")) %>% 
          select(!matches("acc$|res$|pfit")) %>%  
          pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Test",values_to = "SVME") %>% 
          mutate(Test = str_replace_all(Test,pattern = ".*_",replacement = "")) %>% 
          left_join(CNB_completed) %>% 
          filter(Completed == 1) %>% 
          select(-Completed) %>% 
          group_by(Test) %>% 
          mutate(Flagged = ifelse(SVME < quantile(SVME,input$cutoff/100,na.rm = T),"F","V")) %>% 
          mutate(Cutoff_value = quantile(SVME,input$cutoff/100,na.rm = T)) %>% 
          ungroup() %>% 
          ggplot(aes(x = SVME,color = Test,fill = Test)) + geom_vline(aes(xintercept = Cutoff_value)) + geom_density(alpha = .5) + facet_wrap(~Test,scales = "free") + labs(x = "Multivariate CNB Validity Estimate")
        
    })
    
    output$RemovedTable = render_gt({
      
      if(input$removeOverallbad){
        df <- df %>% 
          filter(Overall_valid == 1)
      }
      
        df_for_table <- df %>%
          select(!matches("acc$|res$|pfit")) %>% 
          mutate(across(.cols = contains("PFscores"),~ case_when(.x < quantile(.x,input$cutoff/100,na.rm = T) ~ "F",TRUE ~ "V"),.names = "{.col}_flag")) %>%
          select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,matches("flag$")) %>% 
          select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,contains("ER40"),contains("PMAT"),contains("CPF"),contains("SVOLT"),contains("ADT"),contains("SLNB"),contains("MEDF")) %>% 
          mutate(across(.cols = contains("comment_flag"),.fns = ~ case_when(.x == 1 ~ "V",.x == 2 ~ "F",TRUE ~ "V"))) %>%
          pivot_longer(cols = ER40_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
          mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment",str_detect(Test,pattern = "^PFscores") ~ "itemwise")) %>% 
          mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB")) %>% 
          left_join(CNB_completed) %>%
          group_by(Test,test_sessions_v.battery) %>% 
          filter(n() != sum(Completed == 0)) %>%
          ungroup() %>% 
          select(-Completed) %>% 
          pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
          mutate(itemwise = case_when(itemwise == "V" ~ "Valid itemwise",itemwise == "F" ~ "Flagged itemwise")) %>% 
          mutate(Quality = case_when(AV == "V" & Comment == "V" ~ "Both valid",AV == "F" & Comment == "V" ~ "Only auto-validation flagged",AV == "V" & Comment == "F" ~ "Only reviewer flagged",AV == "F" & Comment == "F" ~ "Both flagged")) %>% 
          mutate(Quality = factor(Quality,levels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"))) 
        
        df_for_table_list <- df_for_table %>% 
          group_split(Test)
        
        build_table <- function(df){
          
          table <- df %>% 
            filter(Test == unique(df$Test)) %>% 
            mutate(AV = ifelse(AV == "V","Valid","Flagged")) %>% 
            mutate(Comment = ifelse(Comment == "V","Valid","Flagged")) %>% 
            select(itemwise,test_sessions_v.age,test_sessions_v.gender) %>% 
            tbl_summary(by = itemwise,label = list(test_sessions_v.age ~ "Age",test_sessions_v.gender ~ "Sex"))
          
          return(table)
        }
        
        QA_tables <- map(df_for_table_list,build_table)
        
        Overall <- df_for_table %>% 
          mutate(AV = ifelse(AV == "V","Valid","Flagged")) %>% 
          mutate(Comment = ifelse(Comment == "V","Valid","Flagged")) %>% 
          select(itemwise,test_sessions_v.age,test_sessions_v.gender) %>%  
          tbl_summary(by = itemwise,label = list(test_sessions_v.age ~ "Age",test_sessions_v.gender ~ "Sex")) %>% 
          modify_spanning_header(all_stat_cols() ~ "**Overall**")
        
        all_table <- tbl_merge(list(Overall,QA_tables[[1]],QA_tables[[2]],QA_tables[[3]],QA_tables[[4]],QA_tables[[5]],QA_tables[[6]],QA_tables[[7]]),tab_spanner = c("**Overall**","**ADT**","**CPF**","**ER40**","**MEDF**","**PMAT**","**SLNB**","**SVOLT**")) %>% 
          as_gt()
        all_table
      
    })
    
    output$itemwiseFlags = renderPlot({
      
      if(input$removeOverallbad){
        df <- df %>% 
          filter(Overall_valid == 1)
      }
      
      itemwise_flagged <- df %>% 
        select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
        pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>% 
        mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>% 
        separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>% 
        mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>% 
        filter(Validity_score_type == "SMVE") %>% 
        mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
        group_by(Test) %>% 
        mutate(Validity_score = case_when(Validity_score < quantile(Validity_score,input$cutoff/100,na.rm = T) ~ "F",TRUE ~ "V")) %>% 
        ungroup() %>% 
        left_join(CNB_completed) %>%
        group_by(Test,test_sessions_v.battery) %>% 
        filter(n() != sum(Completed == 0)) %>%
        ungroup() %>% 
        select(-Completed) %>% 
        rename(Validity_score_overall = "Validity_score") 
        
      theme_set(theme_minimal())
      theme_update(text = element_text(size = 15),legend.position = "bottom")
      
      if(!input$cumulative){
        if(input$removeOverallbad){
          df <- df %>% 
            filter(Overall_valid == 1)
        }
        df %>% 
          select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
          pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>%
          mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>%
          separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>%
          mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>%
          mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
          left_join(CNB_completed) %>%
          group_by(Test,test_sessions_v.battery) %>% 
          filter(n() != sum(Completed == 0)) %>% 
          ungroup() %>% 
          select(-Completed) %>%
          mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
          left_join(itemwise_flagged) %>%
          mutate(Validity_score_type = factor(Validity_score_type,levels = c("SMVE","acc","res","pfit1","pfit2"))) %>%
          group_by(test_sessions.bblid,Test) %>%
          arrange(test_sessions.bblid,Test,Validity_score_type) %>%
          mutate(Validity_score_overall = Validity_score_overall[1]) %>%
          ungroup() %>%
          filter(Validity_score_type != "SMVE") %>%
          mutate(Validity_score_type = case_when(Validity_score_type == "acc" ~ "Accuracy",Validity_score_type == "res" ~ "Outlier \n reaction time",Validity_score_type == "pfit1" ~ "Person Fit 1",Validity_score_type == "pfit2" ~ "Person Fit 2")) %>%
          mutate(Validity_score_overall = case_when(Validity_score_overall == "V" ~ "Valid",Validity_score_overall == "F" ~ "Flagged")) %>%
          ggplot(aes(x = Validity_score,fill = Validity_score_overall)) + geom_density(alpha = .5) + facet_grid(Validity_score_type ~ Test,scales = "free") + labs(x = "Score",fill = "Quality by Multivariate CNB Validity Estimate")
      } else{
        if(input$test != "All"){
          
          if(input$removeOverallbad){
            df <- df %>% 
              filter(Overall_valid == 1)
          }
          
          df %>% 
            select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
            pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>% 
            mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>% 
            separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>% 
            mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>% 
            mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
            left_join(CNB_completed) %>%
            group_by(Test,test_sessions_v.battery) %>% 
            filter(n() != sum(Completed == 0)) %>%
            ungroup() %>% 
            select(-Completed) %>% 
            mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
            left_join(itemwise_flagged) %>% 
            mutate(Validity_score_type = factor(Validity_score_type,levels = c("acc","res","pfit1","pfit2","SMVE"))) %>% 
            group_by(test_sessions.bblid,Test) %>% 
            arrange(test_sessions.bblid,Test,Validity_score_type) %>% 
            mutate(Validity_score_overall = Validity_score_overall[5]) %>% 
            mutate(Validity_score_cumulative = cumsum(Validity_score)) %>% 
            filter(sum(is.na(Validity_score)) == 0) %>% 
            ungroup() %>% 
            group_by(Test,Validity_score_type) %>%
            mutate(Validity_score_cumulative = as.numeric(scale(Validity_score_cumulative,scale = F))) %>%
            ungroup() %>%
            filter(Validity_score_type != "SMVE") %>% 
            filter(Test == input$test) %>% 
            mutate(Validity_score_type = case_when(Validity_score_type == "acc" ~ "Accuracy",Validity_score_type == "res" ~ "Accuracy + Outlier reaction time",Validity_score_type == "pfit1" ~ "Accuracy + Outlier reaction time + Person Fit 1",Validity_score_type == "pfit2" ~ "Accuracy + Outlier reaction time + Person Fit 1 + Person Fit 2")) %>% 
            mutate(Validity_score_overall = case_when(Validity_score_overall == "V" ~ "Valid",Validity_score_overall == "F" ~ "Flagged")) %>% 
            ggplot(aes(x = Validity_score_type,y = Validity_score_cumulative,color = Validity_score_overall)) + geom_line(aes(group = test_sessions.bblid,alpha = Validity_score_overall)) + scale_alpha_manual(values = c(1,.5)) + guides(alpha = "none") + labs(x = "",y = "Cumulative Validity Score \n (centered at mean)", color = "Quality by Multivariate CNB Validity Estimate",title = input$test) 
          
        } else{
          if(input$removeOverallbad){
            df <- df %>% 
              filter(Overall_valid == 1)
          }
          
          df %>% 
            select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
            pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>% 
            mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>% 
            separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>% 
            mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>% 
            mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
            left_join(CNB_completed) %>%
            group_by(Test,test_sessions_v.battery) %>% 
            filter(n() != sum(Completed == 0)) %>%
            ungroup() %>% 
            select(-Completed) %>% 
            mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
            left_join(itemwise_flagged) %>% 
            mutate(Validity_score_type = factor(Validity_score_type,levels = c("acc","res","pfit1","pfit2","SMVE"))) %>% 
            group_by(test_sessions.bblid,Test) %>% 
            arrange(test_sessions.bblid,Test,Validity_score_type) %>% 
            mutate(Validity_score_overall = Validity_score_overall[5]) %>% 
            mutate(Validity_score_cumulative = cumsum(Validity_score)) %>% 
            filter(sum(is.na(Validity_score)) == 0) %>% 
            ungroup() %>% 
            group_by(Test,Validity_score_type) %>%
            mutate(Validity_score_cumulative = as.numeric(scale(Validity_score_cumulative,scale = F))) %>%
            ungroup() %>%
            filter(Validity_score_type != "SMVE") %>% 
            mutate(Validity_score_type = case_when(Validity_score_type == "acc" ~ "Acc",Validity_score_type == "res" ~ "Acc + RT",Validity_score_type == "pfit1" ~ "Acc + RT + \n PF1",Validity_score_type == "pfit2" ~ "Acc + RT + \n PF1 + PF2")) %>% 
            mutate(Validity_score_overall = case_when(Validity_score_overall == "V" ~ "Valid",Validity_score_overall == "F" ~ "Flagged")) %>% 
            ggplot(aes(x = Validity_score_type,y = Validity_score_cumulative,color = Validity_score_overall)) + geom_line(aes(group = test_sessions.bblid,alpha = Validity_score_overall)) + scale_alpha_manual(values = c(1,.5)) + guides(alpha = "none") + facet_wrap(~Test) + labs(x = "",y = "Cumulative Validity Score \n (centered at mean)", color = "Quality by Multivariate CNB Validity Estimate") 
          
        }
      }
})
    
    output$Fatigue <- renderPlot({
      
      theme_set(theme_minimal())
      theme_update(axis.text.x = element_text(size = 12,angle = 45),legend.position = "none")
      
      if(input$removeOverallbad){
        df <- df %>% 
          filter(Overall_valid == 1)
      }
      
      df_for_fatigue_list <- df %>% 
        select(test_sessions.bblid,test_sessions_v.battery,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_CR,VSPLOT15.VSPLOT15_CR) %>% 
        pivot_longer(ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Score") %>% 
        mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
        mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>% 
        mutate(Test = str_replace_all(Test,pattern = "36|90|24|2|15",replacement = "")) %>% 
        left_join(CNB_completed) %>% 
        group_by(test_sessions_v.battery,Test) %>% 
        filter(sum(Completed == 0) != n()) %>% 
        ungroup() %>% 
        left_join(CNB_order_clean,by = c("Test","test_sessions_v.battery" = "BatteryName")) %>% 
        mutate(Score = case_when(Test == "MPRACT" ~ -Score,TRUE ~ Score)) %>% 
        group_by(Test) %>% 
        mutate(Score_z = as.numeric(scale(Score))) %>% 
        ungroup() %>% 
        group_split(test_sessions.bblid)
      
      # df %>% 
      #   select(ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_RTCR,VSPLOT15.VSPLOT15_CR) %>% 
      #   cor(use = "pairwise.complete.obs") %>% 
      #   View()
      
      find_fatigue <- function(data){
        data <- data %>% 
          arrange(TestOrder)
        
        Test_invalid_one <- c()
        Rolling_avg <- c()
        Test_invalid_one[1] <- 0
        Rolling_avg[1] <- NA

      
        for(i in 2:(nrow(data))){
          first_scores <- data$Score_z[1:(i-1)]
          next_score <- data$Score_z[i]
          first_scores_clean <- first_scores[which(Test_invalid_one == 0)]
          first_scores_avg <- mean(first_scores_clean,na.rm = T)
          if(is.na(first_scores_avg)|is.na(next_score)){
            Test_invalid_one[i] <- 0
            Rolling_avg[i] <- NA
          } else if(first_scores_avg > (next_score + input$dropoff)){
            Test_invalid_one[i] <- 1
            Rolling_avg[i] <- first_scores_avg

          }else{
            Test_invalid_one[i] <- 0
            Rolling_avg[i] <- NA

          }
        }
        
        data[,"Test_fatigue"] <- as.numeric(Test_invalid_one)
        data[,"Rolling_avg"] <- Rolling_avg
        
        return(data)
        
      }
      
      df_fatigue_codes <- map_dfr(df_for_fatigue_list,find_fatigue)
      
      if(input$byTest == FALSE){
        
        theme_set(theme_minimal())
        theme_update(text = element_text(size = 12),axis.text.x = element_text(angle = 45),legend.position = "none")
        
        prop_flagged <- df_fatigue_codes %>% 
          group_by(test_sessions.bblid) %>% 
          mutate(Ind_flagged = ifelse(sum(Test_fatigue) > 0,1,0)) %>% 
          ungroup() %>% 
          distinct(test_sessions.bblid,.keep_all = T) %>% 
          group_by(test_sessions_v.battery) %>% 
          summarize(Num_flagged_battery = paste0(sum(Ind_flagged),"/",n())) %>% 
          mutate(Num_flagged_battery = paste0(test_sessions_v.battery,", ",Num_flagged_battery," flagged"))
        
        
        plot1 <- df_fatigue_codes %>% 
          group_by(test_sessions.bblid) %>% 
          mutate(Ind_flagged = factor(ifelse(sum(Test_fatigue) > 0,1,0))) %>% 
          ungroup() %>% 
          mutate(Test_fatigue = factor(Test_fatigue)) %>% 
          filter(test_sessions_v.battery != "RareCNV") %>% 
          left_join(prop_flagged) %>% 
          filter(!is.na(Score_z)) %>% 
          ggplot(aes(x = reorder_within(Test,by = TestOrder,within = test_sessions_v.battery),y = Score_z,group = test_sessions.bblid,color = Ind_flagged,alpha = Ind_flagged)) + scale_color_manual(values = c("grey","black")) + scale_alpha_manual(values = c(.4,2)) + geom_point(aes(size = Test_fatigue,shape = Test_fatigue),fill = "black") + scale_size_manual(values = c(1,3)) + scale_shape_manual(values = c(21,23)) + geom_line() + facet_wrap(~Num_flagged_battery,scales = "free") + scale_x_reordered() + labs(x = "",y = "Z-score") 
        
        plot2 <- df_fatigue_codes %>% 
          group_by(test_sessions.bblid) %>% 
          mutate(Ind_flagged = factor(ifelse(sum(Test_fatigue) > 0,1,0))) %>% 
          ungroup() %>% 
          mutate(Test_fatigue = factor(Test_fatigue)) %>% 
          filter(test_sessions_v.battery == "RareCNV") %>% 
          left_join(prop_flagged) %>% 
          filter(!is.na(Score_z)) %>% 
          ggplot(aes(x = reorder_within(Test,by = TestOrder,within = test_sessions_v.battery),y = Score_z,group = test_sessions.bblid,color = Ind_flagged,alpha = Ind_flagged)) + scale_color_manual(values = c("grey","black")) + scale_alpha_manual(values = c(.4,2)) + geom_point(aes(size = Test_fatigue,shape = Test_fatigue),fill = "black") + scale_size_manual(values = c(1,3)) + scale_shape_manual(values = c(21,23)) + geom_line() + facet_wrap(~Num_flagged_battery,scales = "free") + scale_x_reordered() + labs(x = "",y = "Z-score") 
        
        ggarrange(plot1,plot2,labels = c("",""))
        
      } else{
        theme_set(theme_minimal())
        theme_update(text = element_text(size = 12),legend.position = "none")
        
        prop_flagged <- df_fatigue_codes %>% 
          group_by(Test) %>% 
          summarize(Num_flagged_test = paste0(sum(Test_fatigue),"/",n())) %>% 
          mutate(Test_num_flagged = paste0(Test,", ",Num_flagged_test," flagged"))
        
        
        df_fatigue_codes %>% 
          filter(Test_fatigue == 1) %>% 
          select(test_sessions.bblid,Test,Score_z,Rolling_avg) %>% 
          pivot_longer(cols = c(Score_z,Rolling_avg),names_to = "Type",values_to = "Score") %>% 
          mutate(Type = case_when(Type == "Score_z" ~ "Flagged \n Score",Type == "Rolling_avg" ~ "Rolling \n Average",TRUE ~ NA_character_)) %>% 
          left_join(prop_flagged) %>% 
          ggplot(aes(x = factor(Type,levels = c("Rolling \n Average","Flagged \n Score")),y = Score,color = test_sessions.bblid)) + geom_point() + geom_line(aes(group = test_sessions.bblid)) + facet_wrap(~Test_num_flagged,scales = "free") + labs(x = "Score Type",y = "Performance (Z-score)")
          
      }
      
      
      
    })
    
    output$all_table <- render_gt({
      if(input$removeOverallbad){
        df <- df %>% 
          filter(Overall_valid == 1)
      }
      
      df_for_fatigue_list <- df %>% 
        select(test_sessions.bblid,test_sessions_v.battery,ADT36_A.ADT36A_CR,CPF_B.CPF_CR,ER40_D.ER40D_CR,MEDF36_A.MEDF36A_CR,MPRACT.MP2RTCR,PCET_A.PCET_ACC2,PMAT24_A.PMAT24_A_CR,SCTAP.SCTAP_TOT,SLNB2_90.SLNB2_MCR,SPCPTN90.SCPN90_TP,SVOLT_A.SVOLT_CR,VSPLOT15.VSPLOT15_CR) %>% 
        pivot_longer(ADT36_A.ADT36A_CR:last_col(),names_to = "Test",values_to = "Score") %>% 
        mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
        mutate(Test = str_replace_all(Test,pattern = "\\..*",replacement = "")) %>% 
        mutate(Test = str_replace_all(Test,pattern = "36|90|24|2|15",replacement = "")) %>% 
        left_join(CNB_completed) %>% 
        group_by(test_sessions_v.battery,Test) %>% 
        filter(sum(Completed == 0) != n()) %>% 
        ungroup() %>% 
        left_join(CNB_order_clean,by = c("Test","test_sessions_v.battery" = "BatteryName")) %>% 
        mutate(Score = case_when(Test == "MPRACT" ~ -Score,TRUE ~ Score)) %>% 
        group_by(Test) %>% 
        mutate(Score_z = as.numeric(scale(Score))) %>% 
        ungroup() %>% 
        group_split(test_sessions.bblid)
      
      
      find_fatigue <- function(data){
        data <- data %>% 
          arrange(TestOrder)
        
        Test_invalid_one <- c()
        Rolling_avg <- c()
        Test_invalid_one[1] <- 0
        Rolling_avg[1] <- NA
        
        for(i in 2:(nrow(data))){
          first_scores <- data$Score_z[1:(i-1)]
          next_score <- data$Score_z[i]
          first_scores_clean <- first_scores[which(Test_invalid_one == 0)]
          first_scores_avg <- mean(first_scores_clean,na.rm = T)
          if(is.na(first_scores_avg)|is.na(next_score)){
            Test_invalid_one[i] <- 0
            Rolling_avg[i] <- NA
          }else if(first_scores_avg > (next_score + input$dropoff)){
            Test_invalid_one[i] <- 1
            Rolling_avg[i] <- first_scores_avg
            
          }else{
            Test_invalid_one[i] <- 0
            Rolling_avg[i] <- NA
            
          }
        }
        
        data[,"Test_fatigue"] <- as.numeric(Test_invalid_one)
        data[,"Rolling_avg"] <- Rolling_avg
        
        return(data)
      }
      
    df_fatigue_codes <- map_dfr(df_for_fatigue_list,find_fatigue)
        
     df_for_mega_table <- df %>% 
       select(!matches("acc$|res$|pfit")) %>% 
       mutate(across(.cols = contains("PFscores"),~ case_when(.x < quantile(.x,input$cutoff/100,na.rm = T) ~ "F",TRUE ~ "V"),.names = "{.col}_flag")) %>%
       select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,test_sessions_v.battery,matches("flag$")) %>% 
       mutate(across(.cols = contains("comment_flag"),.fns = ~ case_when(.x == 1 ~ "V",.x == 2 ~ "F",TRUE ~ "V"))) %>%
       pivot_longer(cols = ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>% 
       mutate(Flag_type = case_when(str_detect(Test,pattern = "_AV_") ~ "AV",str_detect(Test,pattern = "_comment_") ~ "Comment",str_detect(Test,pattern = "^PFscores") ~ "itemwise")) %>% 
       mutate(Test = case_when(str_detect(Test,pattern = "ER40") ~ "ER40",str_detect(Test,pattern = "PMAT") ~ "PMAT",str_detect(Test,pattern = "CPF") ~ "CPF",str_detect(Test,pattern = "SVOLT") ~ "SVOLT",str_detect(Test,pattern = "ADT") ~ "ADT",str_detect(Test,pattern = "MEDF") ~ "MEDF",str_detect(Test,pattern = "SLNB") ~ "SLNB",str_detect(Test,pattern = "SCTAP") ~ "SCTAP",str_detect(Test,pattern = "PCET") ~ "PCET",str_detect(Test,pattern = "SPCPTN") ~ "SPCPTN",str_detect(Test,pattern = "MPRACT") ~ "MPRACT",str_detect(Test,pattern = "VSPLOT") ~ "VSPLOT")) %>% 
       left_join(CNB_completed) %>%
       filter(Completed == 1) %>% 
       select(-Completed) %>% 
       pivot_wider(names_from = "Flag_type",values_from = "Valid_code") %>% 
       left_join(df_fatigue_codes[,c("test_sessions.bblid","Test","Test_fatigue")]) %>% 
       mutate(Test_fatigue = case_when(Test_fatigue == 0 ~ "V",Test_fatigue == 1 ~ "F")) %>% 
       pivot_longer(AV:last_col(),names_to = "QC_method",values_to = "Valid_code") %>% 
       mutate(QC_method = case_when(QC_method == "AV" ~ "Auto-validation",QC_method == "Comment" ~ "Reviewer",QC_method == "itemwise" ~ "Multivariate CNB Validity Estimate",QC_method == "Test_fatigue" ~ "Score Dropoff")) %>% 
       select(test_sessions.bblid,test_sessions_v.age,test_sessions_v.gender,Test,QC_method,Valid_code) %>% 
       mutate(Valid_code = case_when(Valid_code == "V" ~ "Valid",Valid_code == "F" ~ "Flagged")) 
    
     if(input$method == "Multivariate CNB Validity Estimate"){
       df_for_mega_table_list <- df_for_mega_table %>% 
         filter(Test %in% c("ADT","CPF","ER40","MEDF","PMAT","SLNB","SVOLT")) %>% 
         group_split(Test)
     }else{
       df_for_mega_table_list <- df_for_mega_table %>% 
         group_split(Test)
     }
     
     
     build_table_mega <- function(data){
       
       table <- data %>% 
         filter(Test == unique(data$Test)) %>% 
         pivot_wider(names_from = "QC_method",values_from = "Valid_code") %>% 
         select(-test_sessions.bblid,-Test) %>% 
         tbl_summary(by = input$method,label = list(test_sessions_v.age ~ "Age",test_sessions_v.gender ~ "Sex"))
       
       return(table)
     }
     
     QA_tables <- map(df_for_mega_table_list,build_table_mega)
     
     Overall <- df_for_mega_table %>% 
       pivot_wider(names_from = "QC_method",values_from = "Valid_code") %>% 
       select(-test_sessions.bblid,-Test) %>%  
       tbl_summary(by = input$method,label = list(test_sessions_v.age ~ "Age",test_sessions_v.gender ~ "Sex")) %>% 
       modify_spanning_header(all_stat_cols() ~ "**Overall**")
     
     if(input$method == "Multivariate CNB Validity Estimate"){
       all_table <- tbl_merge(list(Overall,QA_tables[[1]],QA_tables[[2]],QA_tables[[3]],QA_tables[[4]],QA_tables[[5]],QA_tables[[6]],QA_tables[[7]]),tab_spanner = c("**Overall**","**ADT**","**CPF**","**ER40**","**MEDF**","**PMAT**","**SLNB**","**SVOLT**")) %>% 
         as_gt()
       all_table
     } else{
       all_table <- tbl_merge(list(Overall,QA_tables[[1]],QA_tables[[2]],QA_tables[[3]],QA_tables[[4]],QA_tables[[5]],QA_tables[[6]],QA_tables[[7]],QA_tables[[8]],QA_tables[[9]],QA_tables[[10]],QA_tables[[11]],QA_tables[[12]]),tab_spanner = c("**Overall**","**ADT**","**CPF**","**ER40**","**MEDF**","**MPRACT**","**PCET**","**PMAT**","**SCTAP**","**SLNB**","**SPCPTN**","**SVOLT**","**VSPLOT**")) %>% 
         as_gt()
       all_table
     }
     
     
      
      
    }
    )
    
    output$dynamic <- renderUI({
      req(input$plot_hover) 
      if(input$cumulative && input$test != "All"){
        verbatimTextOutput("vals")
      }
    })
    
    output$vals <- renderPrint({
      if(input$test != "All" && input$cumulative){
        if(input$removeOverallbad){
          df <- df %>% 
            filter(Overall_valid == 1)
        }
        itemwise_flagged <- df %>% 
          select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
          pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>% 
          mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>% 
          separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>% 
          mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>% 
          filter(Validity_score_type == "SMVE") %>% 
          mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
          group_by(Test) %>% 
          mutate(Validity_score = case_when(Validity_score < quantile(Validity_score,input$cutoff/100,na.rm = T) ~ "F",TRUE ~ "V")) %>% 
          ungroup() %>% 
          left_join(CNB_completed) %>%
          group_by(Test,test_sessions_v.battery) %>% 
          filter(n() != sum(Completed == 0)) %>%
          ungroup() %>% 
          select(-Completed) %>% 
          rename(Validity_score_overall = "Validity_score")
        
        label_df_small <- df %>% 
          select(test_sessions.bblid,test_sessions_v.battery,matches("^PFscores")) %>%
          pivot_longer(cols = PFscores_CPF:last_col(),names_to = "Type",values_to = "Validity_score") %>% 
          mutate(Type = str_replace_all(Type,pattern = "^PFscores_",replacement = "")) %>% 
          separate(Type,into = c("Test","Validity_score_type"),sep = "_") %>% 
          mutate(Validity_score_type = ifelse(is.na(Validity_score_type),"SMVE",Validity_score_type)) %>% 
          mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>% 
          left_join(CNB_completed) %>%
          group_by(Test,test_sessions_v.battery) %>% 
          filter(n() != sum(Completed == 0)) %>%
          ungroup() %>% 
          select(-Completed) %>% 
          mutate(test_sessions.bblid = as.character(test_sessions.bblid)) %>%
          left_join(itemwise_flagged) %>% 
          mutate(Validity_score_type = factor(Validity_score_type,levels = c("acc","res","pfit1","pfit2","SMVE"))) %>% 
          group_by(test_sessions.bblid,Test) %>% 
          arrange(test_sessions.bblid,Test,Validity_score_type) %>% 
          mutate(Validity_score_overall = Validity_score_overall[5]) %>% 
          mutate(Validity_score_cumulative = cumsum(Validity_score)) %>% 
          filter(sum(is.na(Validity_score)) == 0) %>% 
          ungroup() %>% 
          group_by(Test,Validity_score_type) %>%
          mutate(Validity_score_cumulative = as.numeric(scale(Validity_score_cumulative,scale = F))) %>%
          ungroup() %>%
          filter(Validity_score_type != "SMVE") %>% 
          filter(Test == input$test) %>% 
          mutate(Validity_score_type = case_when(Validity_score_type == "acc" ~ "Accuracy",Validity_score_type == "res" ~ "Accuracy + Outlier reaction time",Validity_score_type == "pfit1" ~ "Accuracy + Outlier reaction time + Person Fit 1",Validity_score_type == "pfit2" ~ "Accuracy + Outlier reaction time + Person Fit 1 + Person Fit 2")) %>% 
          mutate(Validity_score_overall = case_when(Validity_score_overall == "V" ~ "Valid",Validity_score_overall == "F" ~ "Flagged")) %>% 
          filter(Test == input$test)
        small_df <- nearPoints(label_df_small,input$plot_hover,maxpoints = 1) 
        paste("BBL ID:",as.character(small_df$test_sessions.bblid))
      }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

# Extra completion plots 

# df %>%
#   select(test_sessions.bblid,test_sessions_v.battery,matches("AV_flag")) %>%
#   pivot_longer(ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>%
#   mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
#   left_join(CNB_completed) %>%
#   group_by(Test,test_sessions_v.battery) %>%
#   filter(n() != sum(Completed == 0)) %>%
#   summarize(Complete_percent = sum(Completed)/n(),N = n()) %>%
#   ungroup() %>%
#   left_join(CNB_order_clean,by = c("Test","test_sessions_v.battery" = "BatteryName")) %>%
#   group_by(Test,test_sessions_v.battery) %>%
#   arrange(test_sessions_v.battery,TestOrder) %>%
#   ungroup() %>%
#   ggplot(aes(x = reorder_within(Test,by = TestOrder,within = test_sessions_v.battery),y = Complete_percent,color = Test,fill = Test)) + geom_bar(stat = "identity") + facet_wrap(~test_sessions_v.battery,scales = "free") + scale_x_reordered() + labs(x = "",y = "% Completed")
# 
# df %>% 
#   select(test_sessions_v.battery) %>% 
#   tbl_summary(label = list(test_sessions_v.battery = "Battery Name"))

# df %>% 
#   select(test_sessions.bblid,test_sessions_v.battery,matches("AV_flag")) %>%
#   pivot_longer(ADT_AV_flag:last_col(),names_to = "Test",values_to = "Valid_code") %>%
#   mutate(Test = str_replace_all(Test,pattern = "_.*",replacement = "")) %>%
#   left_join(CNB_completed) %>%
#   group_by(Test,test_sessions_v.battery) %>%
#   filter(n() != sum(Completed == 0)) %>%
#   ungroup() %>%
#   left_join(CNB_order_clean,by = c("Test","test_sessions_v.battery" = "BatteryName")) %>%
#   group_by(test_sessions.bblid) %>% 
#   mutate(Incomplete_once = ifelse(sum(Completed) == n(),"All complete","Incomplete at least once")) %>% 
#   ungroup() %>% 
#   filter(Incomplete_once == "Incomplete at least once") %>%
#   group_by(test_sessions.bblid) %>% 
#   mutate(Num = cur_group_id()) %>% 
#   #mutate(Completed = factor(ifelse(Completed == 1,"Complete","Incomplete"),levels = c("Incomplete","Complete"))) %>% 
#   mutate(Completed = Completed + seq(-.2,.2,length = 99)[unique(Num)]) %>% 
#   ungroup() %>%
#   ggplot(aes(x = reorder_within(Test,by = TestOrder,within = test_sessions_v.battery),y = Completed,group = test_sessions.bblid,color = test_sessions_v.battery)) + geom_point() + geom_line() + facet_wrap(~test_sessions_v.battery,scales = "free") + scale_x_reordered() + labs(x = "",y = "Completed")



