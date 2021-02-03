
library(shiny)
library(shinythemes)
library(tidyverse)
library(DT)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  
  # Application title

    
  h1("Calculating DNA concentration using the picogreen assay:", style="color:#81A78C;  padding:7px; font-weight: bold; "),
  
  # Sidebar setup 
  sidebarLayout(
    sidebarPanel(
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
      downloadButton("report", "Generate full report")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("std.curve.plot"),
      DTOutput("calculations"),
      plotOutput("plot.factor.grid")
  
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
      geom_smooth(method="lm", se = F, colour = "grey")+
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
      dplyr::select(Location, sample, Fluorescence,adjFluor, `concDNA(ng/ul)`, Row, column)
  })
  
  
  DT.values <- reactive({
    calculated.conc() %>% 
      dplyr::select(Location, sample, Fluorescence,adjFluor, `concDNA(ng/ul)`) %>% 
      filter(!sample=="B") %>%
      filter(!str_detect(sample, "S")) %>%
      datatable(extensions = "Buttons",
                options = list(dom = "frtipB",
                               buttons = list(
                                 list(
                                   extend = "collection",
                                   buttons = c("csv", "excel"),
                                   text = "Download table values"))),
                rownames = F,
                style = "bootstrap")
  })
  
  factor.plot <- reactive({
    calculated.conc() %>% 
      mutate(conc_bin = cut(`concDNA(ng/ul)`, breaks = c(-Inf, 2, 5, 10, 20, 40, Inf))) %>% 
      mutate(row = factor(Row, levels = c( "H", "G", "F", "E", "D", "C", "B", "A"))) %>%
      ggplot(mapping = aes(x = as.numeric(column), y = factor(Row), fill = conc_bin))+
      geom_tile(colour = "white") +
      scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))+
      scale_fill_manual(values = c("#FBE0D3","#D1DDE6","#ACC3D1",   "#38808F","#661F3F", "#F56D65"  ))+
      theme_linedraw()+
      labs(x = "columns", y = "rows" )
  })
  
  
  
  
  
##### create output objects for display in the ui -----------------------------  

output$std.curve.plot <-renderPlot({
  stdcurve_plot()
})

output$calculations <- renderDataTable({DT.values()})


output$plot.factor.grid <-renderPlot({
  factor.plot()
})

output$report <- downloadHandler(
  filename = paste0(Sys.Date(),"_picogreen-report.html"),
  content = function(file) {
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)
    params <- list(plotstdcurve = stdcurve_plot(),
                   list_calculations = calculated.conc(),
                   platename = input$platename,
                  factor.plot = factor.plot())
                   
    rmarkdown::render(tempReport, 
                      output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv()),
    )
  }
)

}

shinyApp(ui, server)