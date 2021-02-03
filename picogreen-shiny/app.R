
library(shiny)
library(shinythemes)
library(tidyverse)
library(DT)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  
  # Application title
  titlePanel("Calculating DNA concentration using the picogreen assay:"),
  
  # Sidebar setup 
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "platename",
                label = "Plate Name:",
                value = "Who am I?"),
      textInput(inputId = "date",
                label = "Date:",
                value = Sys.Date()),
      fileInput(inputId = "fluor",
                label = "Fluorescence readings (csv only):",
                accept = ".csv"),
      fileInput(inputId = "platelayout",
                label = "Plate layout of sampleIDs (csv only):",
                accept = ".csv"),
      numericInput(inputId = "dilution",
                   label = "Volume (uL) of sample added to dye/TE solution:",
                   value = 1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("std.curve.plot"),
      DTOutput("calculations")
      
    )
  )
)

# Define server logic required to draw a histogram

server <- function(input, output, session) {
  
  pico.fluor <- reactive({
    req(input$fluor)
    ext <- tools::file_ext(input$fluor$name)
    switch(ext,
           csv = vroom::vroom(input$fluor$datapath, delim = ","),
           xlsx = readxl::read_excel(input$fluor$datapath, sheet = "Results"),
           validate("Invalid file; Please upload a .csv or .xlsx file")
    )
  })

  plate <-reactive({
    req(input$platelayout)
    ext <- tools::file_ext(input$platelayout$name)
    switch(ext,
           csv = vroom::vroom(input$platelayout$datapath, delim = ","),
           xlsx = readxl::read_excel(input$platelayout$datapath, sheet = "Results"),
           validate("Invalid file; Please upload a .csv or .xlsx file")
    )
  })

values <- reactive({
  pico.fluor() %>%
    as_tibble() %>%
    pivot_longer(-Row, names_to = "column", values_to = "Fluorescence") %>%
    mutate(Location = paste0(Row, column))
})

 all <- reactive({
   plate() %>%
     as_tibble() %>%
     pivot_longer(-Row, names_to = "column", values_to = "sample") %>%
     mutate(Location = paste0(Row, column)) %>%
     left_join(values()) %>%
     select(sample, Fluorescence, Location, Row, column)
 })
#This is the current dilution series for the standard curve - evidently would need adjusting if this changes. Could also make as a radiobutton
stand.val <- c(0.75,0.375,0.1875,0.09375,0.046875,0.0234375,0.01171875,0)

standards <- reactive({
  all() %>%
      filter(str_detect(sample, "S")) %>%
      group_by(sample) %>%
      summarize(m = mean(Fluorescence))
})

 background <- reactive({
   standards() %>%
      filter(sample=="S8") %>%
      select(m) %>%
      as.numeric()
})

  adj.fluor <- reactive({
    standards() %>%
      mutate(adj=m-background()) %>%
      mutate(val=stand.val)
  })

  stdcurve <- reactive({lm(val ~ adj, data = adj.fluor()) })
  intercept <- reactive({as.numeric(stdcurve()$coefficients[1])})
  slope <- reactive({as.numeric(stdcurve()$coefficients[2])})
  fit <-reactive({summary(stdcurve()) })

  calculated.conc <- reactive({
    all() %>%
      mutate(conc =(slope*all$adjFluor+intercept)*200/input$dilution) %>%
      dplyr::select(Location, sample, Fluorescence, conc) %>%
      datatable(extensions = "Buttons",
                options = list(dom = "frtipB",
                               buttons = list(
                                 list(
                                   extend = "collection",
                                   buttons = c("csv", "excel"),
                                   text = "Download table values"))),
                rownames = F,
                style = "bootstrap")
      #formatRound(columns = 2:8, digits=3)

  })

output$std.curve.plot <- plotOutput({
    adj.fluor() %>%
    filter(!sample=="S8") %>%
    ggplot(aes(x = adj, y = val))+
    geom_point()+
    geom_smooth(method="lm", se = F)+
    # annotate(geom="text", x=100, y=0.55, label=paste0("y = ", round(slope, digits = 8), "x + ", round(intercept, digits = 5)),
    #          color="black")+
    # annotate(geom="text", x=60000, y=0.65, label=paste0("r.squared = ", round(fit$r.squared, digits = 4)),
    #          color="black")+
    theme_linedraw()
})

  output$calculations <- renderDataTable({calculated.conc()})


}

shinyApp(ui, server)