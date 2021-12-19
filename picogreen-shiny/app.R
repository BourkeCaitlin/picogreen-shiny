
library(shiny)
library(shinythemes)
library(tidyverse)
library(DT)
library(bslib)

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4, "flatly"),
  
  # Application title

  h1("PICOGREEN CALCULATION OF DNA CONCENTRATION :", style="color:#214455;  padding:15px; font-weight: bold; font: rockwell; letter-spacing: 5px"),
  
  
  
  # Sidebar setup 
  sidebarLayout(
    sidebarPanel(
      tags$style(".well {background-color:#C2E3D3;}"),
      textInput(inputId = "platename",
                label = "Plate name:",
                value = "Who am I?"),
      textInput(inputId = "date",
                label = "Date:",
                value = Sys.Date()),
      fileInput(inputId = "fluor",
                label = "Fluorescence readings:",
                accept = ".csv"),
      fileInput(inputId = "platelayout",
                label = "Plate layout of sampleIDs:",
                accept = ".csv"),
      numericInput(inputId = "dilution",
                   label = "Volume (uL) of sample added to dye/TE solution:",
                   value = 1),
      numericInput(inputId = "totalPool",
                   label = "Total DNA (ng) required for pool:",
                   value = 120),
      numericInput(inputId = "normConc",
                   label = "Target concentration (ng/uL) for normalising:",
                   value = 15),
      downloadButton("report", "Generate full report")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("std.curve.plot"),
      DTOutput("calculations"),
      plotOutput("plot.factor.grid"),
      h3("Volume of water required to acheive normalised concentration:"),
      plotOutput("factorPlotWaterDil_grid"),
      h3("Plate DNA concentrations:"),
      plotOutput("factorPlotConc_grid")
  
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
           xlsx = readxl::read_excel(input$fluor$datapath),
           validate("Invalid file; Please upload a .csv or .xlsx file")
    )
  })

  plate <-reactive({
    req(input$platelayout)
    ext <- tools::file_ext(input$platelayout$name)
    switch(ext,
           csv = vroom::vroom(input$platelayout$datapath, delim = ","),
           xlsx = readxl::read_excel(input$platelayout$datapath),
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

  adj.fluor.stds <- reactive({
    standards() %>%
      mutate(adj=m-background()) %>%
      mutate(val=stand.val)
  })

  stdcurve <- reactive({lm(val ~ adj, data = adj.fluor.stds()) })
  intercept <- reactive({as.numeric(stdcurve()$coefficients[1])})
  slope <- reactive({as.numeric(stdcurve()$coefficients[2])})
  fit <-reactive({summary(stdcurve()) })
  
  stdcurve_plot <- reactive({
    adj.fluor.stds() %>%
      dplyr::filter(!sample=="S8") %>%
      ggplot(aes(x = adj, y = val))+
      geom_point()+
      geom_smooth(method="lm", se = F, colour = "#00A7A3")+
      annotate(geom="text", x=5e06, y=0.6, label=paste0("y = ", round(slope(), digits = 8), "x + ", round(intercept(), digits = 5)),
               color="black")+
      annotate(geom="text", x=5e06, y=0.5, label=paste0("r.squared = ", round(fit()$r.squared, digits = 4)),
               color="black")+
      labs(title = input$platename, subtitle = input$date)+
      theme_linedraw()
  })
  
  calculated.conc <- reactive({
    all() %>%
      mutate(adjFluor = Fluorescence - background()) %>% 
      mutate(`concDNA(ng/ul)` =(slope()*adjFluor+intercept())*200/input$dilution) %>%
      mutate(waterDil = `concDNA(ng/ul)` - input$normConc) %>% 
      mutate(DNAadd= input$normConc) %>% 
      mutate(volumePool= (input$totalPool)/(input$normConc)) %>% 
      dplyr::select(Location, sample, Fluorescence,adjFluor, `concDNA(ng/ul)`, Row, column, waterDil, DNAadd, volumePool)
  })
  
  
  DT.values <- reactive({
    calculated.conc() %>% 
      dplyr::select(Location, sample, Fluorescence,adjFluor, `concDNA(ng/ul)`, waterDil,  DNAadd, volumePool) %>% 
      filter(!sample=="B") %>%
      filter(!str_detect(sample, "S")) %>% 
      datatable() %>% 
      formatRound(columns  = 5, digits = 3)
      
  })
  
  factor.plot <- reactive({
    calculated.conc() %>% 
      mutate(conc_bin = cut(`concDNA(ng/ul)`, breaks = c(-Inf, 2, 5, 10, 20, 40, Inf))) %>% 
      mutate(fill=case_when(
        sample=="B" ~ "blank",
        conc_bin=="(40, Inf]" ~"40+",
        conc_bin=="(20,40]" ~"20 - 40",
        conc_bin=="(10,20]" ~"10 - 20",
        conc_bin=="(5,10]" ~"5 - 10",
        conc_bin=="(2,5]" ~ "2 - 5",
        conc_bin=="(-Inf,2]" ~"Less than 2"
      )) %>% 
      mutate(row = factor(Row, levels = c( "H", "G", "F", "E", "D", "C", "B", "A"))) %>%
      mutate(fill = factor(fill, levels = c("Less than 2", "2 - 5", "5 - 10", "10 - 20","20 - 40","40+", "blank"))) %>%
      ggplot(mapping = aes(x = as.numeric(column), y = factor(row), fill = factor(fill)))+
      geom_tile(colour = "white") +
      scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
      scale_fill_manual(values = c("#FBE0D3","#D1DDE6","#ACC3D1","#38808F","#661F3F", "#F56D65", "#E4E4E4"))+
      theme_linedraw()+
      labs(x = "columns", y = "rows" ,fill = "DNA conc (ng/uL)")
  })
  
  
  factorPlotWaterDil <- reactive({
   calculated.conc() %>% 
      mutate(conc_bin = cut(`concDNA(ng/ul)`, breaks = c(-Inf, 2, 5, 10, 20, 40, Inf))) %>% 
      mutate(fill=case_when(
        sample=="B" ~ "blank",
        str_detect(sample, "S")~"standard",
        waterDil<0~"too low to normalise",
        TRUE~"DNA"
      )) %>% 
       mutate(fill2=case_when(
        sample=="B" ~ "blank",
        conc_bin=="(40, Inf]" ~"40+",
        conc_bin=="(20,40]" ~"20 - 40",
        conc_bin=="(10,20]" ~"10 - 20",
        conc_bin=="(5,10]" ~"5 - 10",
        conc_bin=="(2,5]" ~ "2 - 5",
        conc_bin=="(-Inf,2]" ~"Less than 2"
      )) %>% 
      mutate(row = factor(Row, levels = c( "H", "G", "F", "E", "D", "C", "B", "A"))) %>%
      mutate(fill = factor(fill, levels = c("DNA", "standard","too low to normalise", "blank"))) %>% 
      mutate(fill2 = factor(fill2, levels = c("Less than 2", "2 - 5", "5 - 10", "10 - 20","20 - 40","40+", "blank"))) 
  })
  
  waterPlot <- reactive({
    ggplot(data = factorPlotWaterDil(), mapping = aes(x = as.numeric(column), y = as.factor(row), fill = factor(fill), label = waterDil))+
      geom_tile(colour = "white") +
      geom_text(data = factorPlotWaterDil()[!(factorPlotWaterDil()$sample=="B"|factorPlotWaterDil()$fill=="standard"),], aes(label = round(waterDil,1)))+
      scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
      scale_fill_manual(values = c("cornflowerblue","orange", "red","#E4E4E4"))+
      theme_linedraw()+
      labs(x = "columns", y = "rows" ,fill = "sample")
  })
  
 
 concPlot <- reactive({
   ggplot(data = factorPlotWaterDil(), mapping = aes(x = as.numeric(column), y = as.factor(row), fill = factor(fill2), label = `concDNA(ng/ul)`))+
     geom_tile(aes(fill = factor(fill2)),colour = "white") +
     geom_text(data = factorPlotWaterDil()[!factorPlotWaterDil()$sample=="B",], aes(label = round(`concDNA(ng/ul)`,2)))+
     scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
     scale_fill_manual(values = c("#FBE0D3","#D1DDE6","#ACC3D1","#38808F","#661F3F", "#F56D65", "#E4E4E4"))+
     theme_linedraw()+
     labs(x = "columns", y = "rows" ,fill = "sample")
 })
  
  
  
  
  
##### create output objects for display in the ui -----------------------------  

output$std.curve.plot <-renderPlot({
  stdcurve_plot()
})

output$calculations <- renderDataTable({DT.values()})


output$plot.factor.grid <-renderPlot({
  factor.plot()
})


output$factorPlotWaterDil_grid <- renderPlot({
  waterPlot()
})

output$factorPlotConc_grid <- renderPlot({
  concPlot()
})

output$report <- downloadHandler(
  filename = paste0(Sys.Date(),"_picogreen-report.html"),
  content = function(file) {
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)
    params <- list(plotstdcurve = stdcurve_plot(),
                   list_calculations = calculated.conc(),
                   platename = input$platename,
                   TEdil = input$dilution,
                   totalPool = input$totalPool,
                   waterDils = DT.values(),
                   normConc = input$normConc,
                  factor.plot = factor.plot(),
                  water.plot = waterPlot(),
                  conc.plot = concPlot(),
                  factorPlotWaterDil = factorPlotWaterDil()
                  )
    rmarkdown::render(tempReport, 
                      output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv()),
    )
  }
)

}

shinyApp(ui, server)