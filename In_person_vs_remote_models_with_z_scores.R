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
select <- dplyr::select

CNB <- read_csv("/Users/hillmann/Projects/22q/Data/22q_in_person_vs_remote_CNB_clean.csv") 

CNB$Test_Location <- factor(CNB$Test_Location)

CNB_iq <- missForest(as.data.frame(CNB %>% select(test_sessions_v.age,test_num,Test_Location,matches("_asr"))))$ximp %>% 
  tibble() %>% 
  cbind(CNB[,c("bblid")]) %>% 
  relocate(bblid) %>% 
  relocate(Test_Location,.after = bblid) %>% 
  relocate(test_num,.after = Test_Location) %>% 
  rowwise() %>% 
  mutate(Overall_Accuracy = mean(c_across(cols = matches("az_asr")))) %>%
  mutate(Overall_Speed = mean(c_across(cols = matches("sz_asr")))) %>% 
  ungroup() %>% 
  mutate(Overall_Accuracy = as.numeric(scale(Overall_Accuracy))) %>%
  mutate(Overall_Speed = as.numeric(scale(Overall_Speed))) %>%
  mutate(Prior_CNB_tests = test_num-1) %>% 
  select(-test_num)
  
# Use codebook to build data frame which maps test acronyms to test names
Test_map <- data.frame('Prefix' = c("eid","smem","fmem","sm","abf","nvr","edi","adi","spa","mot","att","wm"),
                       "Test_name" = c("Penn Emotion Recognition Test","Visual Object Learning Test","Penn Face Memory Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Matrix Analysis Test","Measured Emotion Differentiation Test","Age Differentiation Test","Penn Line Orientation Test",
                                       "Penn Finger Tapping Test","Penn Continuous Performance Test","Letter-N-Back Test"))

# Use codebook to map measurements to longer names
Metric_map <- data.frame("Middle" = c("az","sz"),
                         "Label" = c("Accuracy","Speed"))

# Run regression model for CNB IQ 

iq_mod <- lm(Overall_Accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_Location,data = CNB_iq)

iq_mod_labels <- c("Intercept","Prior CNB tests","Prior CNB tests<sup>2","Remote")
names(iq_mod_labels) <- c("(Intercept)","Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_LocationRemote")
tab_model(iq_mod,pred.labels = iq_mod_labels)

speed_g_mod <- lm(Overall_Speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_Location,data = CNB_iq)

speed_g_mod_labels <- c("Intercept","Prior CNB tests","Prior CNB tests<sup>2","Remote")
names(speed_g_mod_labels) <- c("(Intercept)","Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_LocationRemote")
tab_model(speed_g_mod,pred.labels = speed_g_mod_labels)


# Run mixed-effect models 

acc_long <- CNB %>% 
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

speed_long <- CNB %>% 
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

#acc_model_no_interaction <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = acc_long)
acc_model_interaction <- lmer(Z_score_accuracy ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = acc_long)

#speed_model_no_interaction <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + (1|bblid),data = speed_long)
speed_model_interaction <- lmer(Z_score_speed ~ Prior_CNB_tests + I(Prior_CNB_tests^2) + Test_name + Test_Location + Test_name:Test_Location + (1|bblid),data = speed_long)

# Conduct F-tests for accuracy and speed models

anova(acc_model_interaction)
anova(speed_model_interaction)

# Get model plots
theme_set(theme_minimal())
theme_update(text = element_text(size = 12))

acc_no_int_labels <- c("Prior CNB tests",expression(paste("Prior CNB tests" ^2)),"Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(acc_no_int_labels) <- c("Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

speed_no_int_labels <- c("Prior CNB tests",expression(paste("Prior CNB tests" ^2)),"Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test","Remote")
names(speed_no_int_labels) <- c("Prior_CNB_tests","I(Prior_CNB_tests^2)","Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test","Test_LocationRemote")

acc_int_labels <- c(acc_no_int_labels,str_c(c("Letter-N-Back Test","Measured Emotion Differentiation Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(acc_int_labels) <- c(names(acc_no_int_labels),str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_namePenn Conditional Exclusion Test","Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test","Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

speed_int_labels <- c(speed_no_int_labels,str_c(c("Letter-N-Back Test","Measured Emotion Differentiation Test","Motor Praxis Test","Penn Conditional Exclusion Test","Penn Continuous Performance Test","Penn Emotion Recognition Test","Penn Face Memory Test","Penn Finger Tapping Test","Penn Line Orientation Test","Penn Matrix Analysis Test","Visual Object Learning Test"),"Remote",sep = ":"))
names(speed_int_labels) <- c(names(speed_no_int_labels),str_c(c("Test_nameLetter-N-Back Test","Test_nameMeasured Emotion Differentiation Test","Test_nameMotor Praxis Test", "Test_namePenn Conditional Exclusion Test", "Test_namePenn Continuous Performance Test","Test_namePenn Emotion Recognition Test","Test_namePenn Face Memory Test","Test_namePenn Finger Tapping Test","Test_namePenn Line Orientation Test","Test_namePenn Matrix Analysis Test", "Test_nameVisual Object Learning Test"),"Test_LocationRemote",sep = ":"))

# plot_model(acc_model_no_interaction,title = "CNB Accuracy: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = acc_no_int_labels,axis.title = "Coefficient Estimate")
# plot_model(speed_model_no_interaction,title = "CNB Speed: No Interaction Mixed-Effect Model",auto.label = FALSE,axis.labels = speed_no_int_labels,axis.title = "Coefficient Estimate")

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

CNB_for_table <- CNB %>%
  mutate(`Prior CNB tests` = test_num - 1) %>% 
  rename(`Test Location` = "Test_Location",Age = test_sessions_v.age,Sex = test_sessions_v.gender) %>%
  mutate(Sex = case_when(Sex == "M" ~ "Male",Sex == "F" ~ "Female",TRUE ~ NA_character_)) %>%
  mutate(`Prior CNB tests` = factor(`Prior CNB tests`)) %>%
  mutate(`Prior CNB tests` = fct_collapse(`Prior CNB tests`,`3+` = c("3","4","5","7"))) %>% 
  mutate(race = case_when(race == 1 ~ "White",race == 2 ~ "Black/African American",race == 3 ~ "Native American",race == 4 ~ "Asian",race == 5 ~ "More than one race",race == 6 ~ "Hawaiian/Pacific Islander",race == 9 ~ NA_character_,TRUE ~ NA_character_)) %>% 
  rename(Race = race)

table1(~ Age + Sex + Race + `Prior CNB tests`|`Test Location`, data = CNB_for_table)

#Plot for CPT

theme_set(theme_minimal())
theme_update(text = element_text(size = 14))

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
  mutate(Prior_CNB_tests = fct_collapse(Prior_CNB_tests,`3+` = c("3","4","5","7"))) %>% 
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



