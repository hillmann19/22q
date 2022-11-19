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
library(gtsummary)
library(ggh4x)
select <- dplyr::select
demo <- read_csv('/Users/hillmann/Projects/22q/Data/subject.csv')
CNB <- read_csv('/Users/hillmann/Projects/22q/Data/22q_in_person_vs_remote_no_relatives_CNB_11_11_2022_NH.csv') 

demo_race <- demo %>% 
  rename(bblid = BBLID) %>% 
  rename(Race = RACE) %>% 
  select(bblid,Race)

CNB$Test_Location <- factor(CNB$Test_Location)

#set.seed(19104)

# CNB_iq <- missForest(as.data.frame(CNB %>% select(test_sessions_v.age,Prior_CNB_test,Test_Location,matches("_asr"))))$ximp %>%
#   tibble() %>%
#   cbind(CNB[,c("bblid","Prior_CNB_tests")]) %>%
#   relocate(bblid) %>%
#   relocate(Test_Location,.after = bblid) %>%
#   rowwise() %>%
#   mutate(Overall_Accuracy = mean(c_across(cols = matches("az_asr")))) %>%
#   mutate(Overall_Speed = mean(c_across(cols = matches("sz_asr")))) %>%
#   ungroup()

#write_csv(CNB_iq,'/Users/hillmann/Projects/22q/Data/22q_in_person_vs_remote_CNB_imputed_11_15_2022_NH.csv')

CNB_iq <- read_csv('/Users/hillmann/Projects/22q/Data/22q_in_person_vs_remote_CNB_imputed_11_15_2022_NH.csv')

# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model for CNB IQ 

iq_mod <- lm(Overall_Accuracy ~ Prior_CNB_test + Test_Location,data = CNB_iq)

iq_mod_labels <- c("Intercept","Prior CNB Test","Remote")
names(iq_mod_labels) <- c("(Intercept)","Prior_CNB_test","Test_LocationRemote")
tab_model(iq_mod,pred.labels = iq_mod_labels)

speed_g_mod <- lm(Overall_Speed ~ Prior_CNB_test + Test_Location,data = CNB_iq)

speed_g_mod_labels <- c("Intercept","Prior CNB Test","Remote")
names(speed_g_mod_labels) <- c("(Intercept)","Prior_CNB_test","Test_LocationRemote")
tab_model(speed_g_mod,pred.labels = speed_g_mod_labels)


# Run mixed-effect models 

acc_long <- CNB_iq %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Accuracy") %>% 
  rename(Z_score_accuracy = Z_score) 
  

speed_long <- CNB_iq %>% 
  select(!matches("_valid$")) %>% 
  pivot_longer(cols = matches("asr$"),names_to = "Test",values_to = "Z_score") %>% 
  mutate(Test_prefix = str_replace_all(Test,pattern = "_.*",replacement = "")) %>% 
  mutate(Test_middle = str_replace_all(Test,pattern = "[[:alpha:]]+_(..)_.*",replacement = "\\1")) %>% 
  rowwise() %>% 
  mutate(Test_name = Test_map$Test_name[which(Test_map$Prefix == Test_prefix)]) %>% 
  mutate(Test_Type = case_when(Test_middle == "az" ~ "Accuracy",Test_middle == "sz" ~ "Speed",TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  filter(Test_Type == "Speed") %>% 
  rename(Z_score_speed = Z_score) 
  
acc_model_interaction <- lmer(Z_score_accuracy ~ 0 + Prior_CNB_test  + Test_name*Test_Location + (1|bblid),data = acc_long)
speed_model_interaction <- lmer(Z_score_speed ~ 0 + Prior_CNB_test + Test_name*Test_Location + (1|bblid),data = speed_long)

# Conduct F-tests for accuracy and speed models

anova(acc_model_interaction)
anova(speed_model_interaction)

# Get model plots
theme_set(theme_minimal())
theme_update(text = element_text(size = 12))

acc_no_int_labels <- c("Prior CNB test","Age Differentiation Test","Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(acc_no_int_labels) <- c("Prior_CNB_test","Test_nameAge Differentiation Test","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

speed_no_int_labels <- c("Prior CNB test","Age Differentiation Test","Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(speed_no_int_labels) <- c("Prior_CNB_test","Test_nameAge Differentiation Test","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

acc_int_labels <- c(acc_no_int_labels,str_c(c("Age Differentiation Test","Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(acc_int_labels) <- c(names(acc_no_int_labels),"Test_LocationRemote",str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_namePenn Conditional Exclusion Test","Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test","Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

speed_int_labels <- c(speed_no_int_labels,str_c(c("Age Differentiation Test","Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(speed_int_labels) <- c(names(speed_no_int_labels),"Test_LocationRemote",str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

plot_model(acc_model_interaction,title = "CNB Accuracy: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = acc_int_labels,axis.title = "Coefficient Estimate")
plot_model(speed_model_interaction,title = "CNB Speed: Mixed-Effect Model with Interaction",auto.label = FALSE,axis.labels = speed_int_labels,axis.title = "Coefficient Estimate")

remote_effects_acc <- as.numeric(summary(acc_model_interaction)$coefficients[12:21,5])
coefficient_names_acc <- rownames(summary(acc_model_interaction)$coefficients[12:21,])

p_val_df_acc <- data.frame(Test = coefficient_names_acc,P_val = remote_effects_acc)

p_val_df_acc$p_val_corrected <- p.adjust(p_val_df_acc$P_val,method = "fdr")

remote_effects_speed <- as.numeric(summary(speed_model_interaction)$coefficients[14:25,5])
coefficient_names_speed <- rownames(summary(speed_model_interaction)$coefficients[14:25,])

p_val_df_speed <- data.frame(Test = coefficient_names_speed,P_val = remote_effects_speed)

p_val_df_speed$p_val_corrected <- p.adjust(p_val_df_speed$P_val,method = "fdr")

# Create demographic tables for presentation

CNB_for_table <- CNB %>%
  left_join(demo_race) %>% 
  mutate(`Prior CNB tests` = test_num - 1) %>% 
  rename(`Test Location` = "Test_Location",Age = test_sessions_v.age,Sex = test_sessions_v.gender) %>%
  mutate(Sex = case_when(Sex == "M" ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>%
  mutate(`Prior CNB tests` = factor(`Prior CNB tests`)) %>%
  mutate(`Prior CNB tests` = fct_collapse(`Prior CNB tests`,`3+` = c("3","4","5","6","8"))) %>% 
  mutate(Race = case_when(Race == 1 ~ "White",Race == 2 ~ "Black/African American",Race == 3 ~ "Native American",Race == 4 ~ "Asian",Race == 5 ~ "More than one race",Race == 6 ~ "Hawaiian/Pacific Islander",Race == 9 ~ NA_character_,TRUE ~ NA_character_)) %>% 
  mutate(Seen_CNB = ifelse(`Prior CNB tests` == '0',"First CNB","Repeat CNB"))

CNB_for_table %>% 
  select(Age,Sex,Race,`Prior CNB tests`,`Test Location`) %>% 
  tbl_summary(by = 'Test Location') %>% 
  add_p(test = list(all_continuous() ~ "t.test", all_categorical() ~ "fisher.test"))

# Check for differences between remote vs in-person by number of tests seen

CNB_for_table %>% 
  select(Age,Sex,Race,Seen_CNB,`Test Location`) %>% 
  tbl_strata(
    strata = 'Seen_CNB',
    .tbl_fun = 
    ~ .x %>% 
      tbl_summary(by = 'Test Location') %>% 
      add_p(test = list(all_continuous() ~ "t.test", all_categorical() ~ "fisher.test")))

# Create table which examines tests removed by location 

CNB_completed <- CNB %>%
  select(bblid,test_sessions.datasetid,Test_Location,matches('_completed')) %>%
  pivot_longer(cols = matches('_completed$'),names_to = 'Test',values_to = 'Complete') %>%
  mutate(Test = str_replace_all(Test,pattern = '_completed',replacement = ''))

CNB_qc_table <- CNB %>%
  select(bblid,test_sessions.datasetid,Test_Location,matches('_removed')) %>%
  pivot_longer(cols = matches('_removed$'),names_to = 'Test',values_to = 'QC') %>%
  mutate(Test = str_replace_all(Test,pattern = '_removed',replacement = '')) %>%
  left_join(CNB_completed) %>%
  filter(Complete == 'Completed')

CNB_qc_table %>% 
  select(QC,Test_Location) %>% 
  tbl_summary(by = 'Test_Location') %>% 
  add_p()

#Plot for CPT

theme_set(theme_minimal())

CNB_CPT <- CNB %>%
  mutate(`Prior CNB tests` = test_num - 1)

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

# Compare In-person and Remote for prior CNB tests

CNB_iq %>% 
  mutate(Prior_CNB_tests = factor(Prior_CNB_tests)) %>%
  mutate(Prior_CNB_tests = fct_collapse(Prior_CNB_tests,`3+` = c("3","4","5","6","8"))) %>% 
  pivot_longer(cols = c(Overall_Accuracy,Overall_Speed),names_to = 'Metric',values_to = 'Score') %>% 
  mutate(Metric = str_replace_all(Metric,pattern = 'Overall_',replacement = '')) %>% 
  group_by(Test_Location,Prior_CNB_tests,Metric) %>% 
  summarize(Avg_score = mean(Score,na.rm = T),SE = sd(Score,na.rm = T)/sqrt(n()),DF = n() - 1) %>% 
  ggplot(aes(x = Prior_CNB_tests,y = Avg_score,fill = Test_Location)) + 
  geom_bar(stat = 'identity',position = 'dodge') + 
  facet_wrap(~Metric) + 
  geom_errorbar(aes(ymin = Avg_score - qt(.975,df = DF)*SE,ymax = Avg_score + qt(.975,df = DF)*SE),position = position_dodge(width = .9),width = .2) +
  geom_hline(yintercept = 0,color = 'black') + 
  labs(x = 'Prior CNB Tests',y = 'CNB Score',fill = 'Test Location')


# Create figures for paper

CNB_acc_summary <- CNB_iq %>% 
  select(bblid,Test_Location,matches('_az_asr$'),Overall_Accuracy) %>% 
  pivot_longer(cols = c(matches('_az_asr$'),Overall_Accuracy),names_to = 'Test',values_to = 'Score') %>% 
  mutate(Test = str_replace_all(Test,pattern = '_az_asr',replacement = '')) %>% 
  mutate(Test = str_to_upper(Test)) %>% 
  mutate(Test = case_when(Test == 'OVERALL_ACCURACY' ~ 'G',TRUE ~ Test)) %>% 
  mutate(Test = factor(Test,levels = c('G','ABF','ATT','WM','FMEM','SMEM','NVR','SPA','EID','EDI','ADI'))) %>% 
  group_by(Test,Test_Location) %>% 
  summarize(Score_avg = mean(Score),Score_SE = sd(Score)/sqrt(n()),DF = n() - 1) %>% 
  mutate(Domain = case_when(Test %in% c('ABF','ATT','WM') ~ 'Executive',
                            Test %in% c('FMEM','SMEM') ~ 'Memory',
                            Test %in% c('NVR','SPA') ~ 'Complex',
                            Test %in% c('EID','EDI','ADI') ~ 'Social',
                            Test == 'G' ~ 'Overall')) %>% 
  mutate(Domain = factor(Domain,levels = c('Overall','Executive','Memory','Complex','Social')))

theme_set(theme_minimal())
theme_update(ggh4x.axis.nestline = element_line(linetype = 1),
             ggh4x.axis.nesttext.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.text = element_text(size=14),
             legend.key.size = unit(1, 'cm'),
             legend.position = c(.8,1),
             title = element_text(size = 16))
 
acc_by_test_line_plot <- CNB_acc_summary %>% 
  ggplot(aes(x = interaction(Test,Domain),y = Score_avg,color = Test_Location,group = Test_Location)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin = Score_avg - qt(.975,df = DF)*Score_SE,
                    ymax = Score_avg + qt(.975,df = DF)*Score_SE,
                    color = Test_Location),
                width = .2) + 
  scale_x_discrete(guide = 'axis_nested',limit = c('G.Overall','','ABF.Executive','ATT.Executive','WM.Executive',
                                                   'FMEM.Memory','SMEM.Memory',
                                                   'NVR.Complex','SPA.Complex',
                                                   'EID.Social','EDI.Social','ADI.Social')) +
  scale_color_brewer(palette = 'Accent') +
  labs(x = '',y = 'Z-score',color = '',title = '') + 
  guides(color = guide_legend(nrow = 1)) +
  annotate("text", x = 4, y = -1.9, label = "**",size = 6) + 
  annotate("text", x = 10, y = -1.9, label = "***",size = 6)

theme_set(theme_minimal())
theme_update(ggh4x.axis.nestline = element_line(linetype = 1),
             ggh4x.axis.nesttext.x = element_text(size = 14),
             axis.title.y = element_text(size = 16),
             legend.text = element_text(size=14),
             legend.key.size = unit(.7, 'cm'),
             legend.position = c(.15,.2),
             title = element_text(size = 16))

acc_by_test_bar_plot <- CNB_acc_summary %>% 
  mutate(Test_Domain = interaction(Test,Domain)) %>% 
  ggplot() +
  geom_col(aes(x = Test_Domain,y = Score_avg,fill = Test_Location),
           position = position_dodge(.9),color = 'black') + 
  geom_errorbar(aes(x = interaction(Test,Domain),y = Score_avg,
                    ymin = Score_avg - Score_SE,
                    ymax = Score_avg + Score_SE,
                    fill = Test_Location),
                position = position_dodge(.9),
                width = .2) + 
  geom_hline(yintercept = 0) + 
  scale_x_discrete(guide = 'axis_nested',limit = c('G.Overall','','ABF.Executive','ATT.Executive','WM.Executive',
                                                   'FMEM.Memory','SMEM.Memory',
                                                   'NVR.Complex','SPA.Complex',
                                                   'EID.Social','EDI.Social','ADI.Social')) +
  scale_fill_brewer(palette = 'Accent') +
  labs(x = '',y = 'Z-score',fill = '',title = '') + 
  annotate("text", x = 4, y = .05, label = "**",size = 6) + 
  annotate("text", x = 10, y = .05, label = "***",size = 6) 


# Same plot for speed scores 

speed_cols <- CNB_iq %>% 
  select(matches('_sz_asr'),Overall_Speed) %>% 
  colnames()

pvals <- c()
for(col in speed_cols){
  pval <- t.test(CNB_iq[[col]] ~ CNB_iq$Test_Location)$p.value
  pvals <- c(pvals,pval)
}

speed_t_tests <- tibble('Test' = speed_cols,'P-value' = pvals)

CNB_speed_summary <- CNB_iq %>% 
  select(bblid,Test_Location,matches('_sz_asr$'),Overall_Speed) %>% 
  pivot_longer(cols = c(matches('_sz_asr$'),Overall_Speed),names_to = 'Test',values_to = 'Score') %>% 
  mutate(Test = str_replace_all(Test,pattern = '_sz_asr',replacement = '')) %>% 
  mutate(Test = str_to_upper(Test)) %>% 
  mutate(Test = case_when(Test == 'OVERALL_SPEED' ~ 'G',TRUE ~ Test)) %>% 
  mutate(Test = factor(Test,levels = c('G','ABF','ATT','WM','FMEM','SMEM','NVR','SPA','EID','EDI','ADI','MOT','SM'))) %>% 
  group_by(Test,Test_Location) %>% 
  summarize(Score_avg = mean(Score),Score_SE = sd(Score)/sqrt(n()),DF = n() - 1) %>% 
  mutate(Domain = case_when(Test %in% c('ABF','ATT','WM') ~ 'Executive',
                            Test %in% c('FMEM','SMEM') ~ 'Memory',
                            Test %in% c('NVR','SPA') ~ 'Complex',
                            Test %in% c('EID','EDI','ADI') ~ 'Social',
                            Test %in% c('MOT','SM') ~ 'Motor',
                            Test == 'G' ~ 'Overall')) %>% 
  mutate(Domain = factor(Domain,levels = c('Overall','Executive','Memory','Complex','Social','Motor')))

speed_by_test_bar_plot <- CNB_speed_summary %>% 
  ggplot() +
  geom_col(aes(x = interaction(Test,Domain),y = Score_avg,fill = Test_Location),
           position = position_dodge(.9),color = 'black') + 
  geom_errorbar(aes(x = interaction(Test,Domain),y = Score_avg,
                    ymin = Score_avg - Score_SE,
                    ymax = Score_avg + Score_SE,
                    fill = Test_Location),
                position = position_dodge(.9),
                width = .2) + 
  geom_hline(yintercept = 0) + 
  scale_x_discrete(guide = 'axis_nested',limit = c('G.Overall','',
                                                   'ABF.Executive','ATT.Executive','WM.Executive',
                                                   'FMEM.Memory','SMEM.Memory',
                                                   'NVR.Complex','SPA.Complex',
                                                   'EID.Social','EDI.Social','ADI.Social',
                                                   'MOT.Motor','SM.Motor')) +
  scale_fill_brewer(palette = 'Accent') +
  labs(x = '',y = 'Z-score',fill = '',title = '') + 
  annotate("text", x = 4, y = .5, label = "*",size = 6) + 
  annotate("text", x = 5, y = .5, label = "***",size = 6) + 
  annotate("text", x = 7, y = .5, label = "*",size = 6) +
  annotate("text", x = 8, y = .5, label = "***",size = 6) +
  annotate("text", x = 9, y = .5, label = "**",size = 6) +
  annotate("text", x = 14, y = .5, label = "*",size = 6)

theme_set(theme_minimal())
theme_update(ggh4x.axis.nestline = element_line(linetype = 1),
             ggh4x.axis.nesttext.x = element_text(size = 14),
             axis.title.y = element_text(size = 16),
             legend.text = element_text(size=14),
             legend.key.size = unit(1, 'cm'),
             legend.position = c(.8,1),
             title = element_text(size = 16))

speed_by_test_line_plot <- CNB_speed_summary %>% 
  ggplot(aes(x = interaction(Test,Domain),y = Score_avg,color = Test_Location,group = Test_Location)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin = Score_avg - Score_SE,
                    ymax = Score_avg + Score_SE,
                    color = Test_Location),
                width = .2) + 
  scale_x_discrete(guide = 'axis_nested',limit = c('G.Overall','',
                                                   'ABF.Executive','ATT.Executive','WM.Executive',
                                                   'FMEM.Memory','SMEM.Memory',
                                                   'NVR.Complex','SPA.Complex',
                                                   'EID.Social','EDI.Social','ADI.Social',
                                                   'MOT.Motor','SM.Motor')) +
  scale_color_brewer(palette = 'Accent') +
  labs(x = '',y = 'Z-score',color = '',title = '') + 
  guides(color = guide_legend(nrow = 1)) + 
  annotate("text", x = 4, y = -1.4, label = "*",size = 6) + 
  annotate("text", x = 5, y = -1.4, label = "***",size = 6) + 
  annotate("text", x = 7, y = -1.4, label = "*",size = 6) +
  annotate("text", x = 8, y = -1.4, label = "***",size = 6) +
  annotate("text", x = 9, y = -1.4, label = "**",size = 6) +
  annotate("text", x = 14, y = -1.4, label = "*",size = 6)

# Lineplots 
ggarrange(acc_by_test_line_plot,speed_by_test_line_plot,
          labels = c('Accuracy','Speed'),common.legend = T,legend = 'bottom',font.label = list(size = 18))

# Barplots 

ggarrange(acc_by_test_bar_plot,speed_by_test_bar_plot,
          labels = c('Accuracy','Speed'),common.legend = T,legend = 'bottom',font.label = list(size = 18))








  
  





