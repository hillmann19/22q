#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggrepel)

df <- read_csv(file = "/Users/hillmann/Projects/22q/Data/QA/CNB/CNB_22q_speed_by_acc.csv") %>% 
  #mutate(Valid_code = factor(Valid_code,levels = c("Valid","Flagged"))) 
  left_join(CNB_for_table) %>% 
  select(-Valid_code) %>% 
  rename("Valid_code" = `Data quality`) %>% 
  filter(Test != "SPCPTNL")

CNB_for_table <- read_csv(file = "/Users/hillmann/Projects/22q/Data/QA/CNB/CNB_comments_for_shiny")
  
  
ui <- fillPage(
  tags$style(type = "text/css", "#flaggedPlot {height: calc(95vh - 80px) !important;}"),
  plotOutput("flaggedPlot", height="95%", width="95%",hover = "plot_hover"),uiOutput("dynamic"),textInput("id","BBL ID"))
  
server <- function(input, output) {
  output$flaggedPlot <- renderPlot({
    theme_set(theme_minimal())
    theme_update(text = element_text(size = 18))
    if(input$id == ""){
      df %>%  
        ggplot(aes(x = Accuracy,y = Speed,color = Valid_code)) + geom_point(aes(size = Valid_code,shape = Valid_code)) + scale_color_manual(name = "Valid Code",labels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"),values = c("#b2182b","#fddbc7","#92c5de","#2166ac")) + scale_shape_manual(name = "Valid Code",labels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"),values = c(16,15,17,18)) + scale_size_manual(name = "Valid Code",labels = c("Both valid","Only auto-validation flagged","Only reviewer flagged","Both flagged"),values = c(1.5,4,4,6)) + facet_wrap(~Test,scales = "free") + labs(x = "Accuracy (Z-scored)",y = "Speed (Z-scored)",color = "Valid Code") 
     } 
      #else if(input$multiple == TRUE){
    #   two_plus_flags <- df %>% 
    #     group_by(test_sessions.bblid) %>% 
    #     filter(sum(Valid_code == "Flagged") >= 2) %>% 
    #     ungroup() %>% 
    #     with(unique(test_sessions.bblid))
    #   df %>%  
    #     ggplot(aes_string(x = "Accuracy",y = "Speed",color = "Valid_code")) + geom_point(aes_string(size = "Valid_code",shape = "Valid_code")) + scale_size_manual(name = "Valid Code",labels = c("Valid","Flagged"),values = c(3, 6)) + scale_shape_manual(name = "Valid Code",labels = c("Valid","Flagged"),values = c(20,18)) + scale_color_manual(name = "Valid Code",labels = c("Valid","Flagged"),values = c("#990000","#011F5b")) +
    #     geom_label_repel(data = df[df$test_sessions.bblid %in% two_plus_flags & df$Valid_code == "Flagged",],aes(label = test_sessions.bblid),size = 5,show.legend = F) + facet_wrap(~Test,scales = "free") + labs(x = "Accuracy (Z-scored)",y = "Speed (Z-scored)") + guides(size = "legend",shape = "legend",color = "legend")
    # } 
      else if(input$id != ""){
      df[,"ID"] <- factor(ifelse(df$test_sessions.bblid == input$id,input$id,paste("Not",input$id)),levels = c(input$id,paste("Not",input$id)))
      df %>%  
        ggplot(aes_string(x = "Accuracy",y = "Speed",color = "ID")) + geom_point(aes_string(size = "ID",shape = "ID")) + scale_size_manual(name = "ID",labels = c(input$id,paste("Not",input$id)),values = c(6, 3)) + scale_shape_manual(name = "ID",labels = c(input$id,paste("Not",input$id)),values = c(18,20)) + scale_color_manual(name = "ID",labels = c(input$id,paste("Not",input$id)),values = c("#990000","#011F5b")) +
        geom_label_repel(data = df[df$test_sessions.bblid == input$id,],aes(label = Valid_code),size = 4,show.legend = F,force_pull = .75) + facet_wrap(~Test,scales = "free") + labs(x = "Accuracy (Z-scored)",y = "Speed (Z-scored)") + guides(size = "legend",shape = "legend",color = "legend")
    }
    })

  
  output$dynamic <- renderUI({
    req(input$plot_hover) 
    verbatimTextOutput("vals")
  })
  
  output$vals <- renderPrint({
    hover <- input$plot_hover 
    # print(str(hover)) # list
    small_df <- nearPoints(df, input$plot_hover,maxpoints = 1)
    req(nrow(small_df) != 0)
    paste("BBL ID:",unique(small_df$test_sessions.bblid))
  })
  }


# Run the application 
shinyApp(ui = ui, server = server)
