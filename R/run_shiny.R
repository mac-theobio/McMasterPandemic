#' run shiny app
#' @importFrom shiny fluidPage titlePanel mainPanel fluidRow radioButtons h3 conditionalPanel
#' @importFrom shiny column selectInput textInput tabsetPanel tabPanel checkboxInput sliderInput
#' @importFrom shiny actionButton tableOutput plotOutput reactive observeEvent renderPlot
#' @importFrom shiny renderTable shinyApp
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom scales log10_trans trans_breaks trans_format math_format
#' @importFrom ggplot2 element_text
#' @export
run_shiny <- function(){
#How I want stuff to look. Also collect input and things.
  ui <- fluidPage(
    titlePanel("McMasterPandemic Shiny"),
    mainPanel(
      fluidRow(
        radioButtons("useOwnParams", label = h3("Use file for parameters? Or input own"),
                     choices = list("Input parameters" = 1, "Use file" = 2),
                     selected = 2),
        #Only show the selector for the parameter file to upload if that's selected.
        conditionalPanel(
          condition = "input.useOwnParams == 2",
          column(2,
                 selectInput("fn",
                             label = "Chose file containing parameters",
                             choices = c("CI_base.csv",
                                         "CI_updApr1.csv",
                                         "ICU1.csv",
                                         "ICU_diffs.csv",
                                         "midas_estimates.csv",
                                         "stanford_estimates.csv"), selected = "ICU1.csv"))),
        fluidRow(
          column(2,
                 textInput("sd", "Simulation Start Date (lubridate mdy)", value = "01-01-2020")),
          column(2,
                 textInput("ed", "Simulation End Date  (lubridate mdy)", value = "06-01-2020")),
          column(2,
                 selectInput("timeChanges", "Include time-changing transmission rates", c("Yes", "No"), selected = "No"))
        ),
  #Only show the selector to input parameters if that's selected.
    tabsetPanel(
      selected = "plotControls",
      tabPanel(
        title = "Time changing transmission rates",
        value = "tcr",
          textInput("timeParsDates", label = "Enter dates of changes here, in lubridate ymd format, separated by commas", placeholder = "2020-02-20, 2020-05-20, 2020-07-02"),
          textInput("timeParsSymbols", label = "Enter the corresponding symbol for each date that you'd like to change here, separated by commas", placeholder = "beta0, beta0, alpha"),
          textInput("timeParsRelativeValues", label = "Enter relative value changes here, separated by commas", placeholder = "0.5, 0.1, 0.01")
      ),
      tabPanel(
        title = "Process and Observation error",
        value = "procObsErr",
        textInput("procError", label = "Enter the relevant parameter", value = 0),
        textInput("ObsError", label = "Enter the relevant parameter", value = 0)
      ),
      tabPanel(
        title = "Simulation Parameters",
        value = "manualParams",
          column(5,
                 textInput("beta0",
                            label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "beta0","meaning"]),
                 textInput("Ca",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Ca","meaning"]),
                 textInput("Cp",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cp","meaning"]),
                textInput("Cs",
                          label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cs","meaning"]),
                textInput("Cm",
                          label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cm","meaning"])),
          column(5,
                 textInput("alpha",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "alpha","meaning"]),
                 textInput("sigma",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "sigma","meaning"]),
                 textInput("gamma_a",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_a","meaning"]),
                 textInput("gamma_s",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_s","meaning"]),
                 textInput("gamma_m",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_m","meaning"])),
          column(5,
                 textInput("gamma_p",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_p","meaning"]),
                 textInput("rho",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "rho","meaning"]),
                 textInput("delta",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "delta","meaning"]),
                 textInput("mu",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "mu","meaning"]),
                 textInput("N",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "N","meaning"])),
          column(5,
                 textInput("E0",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "E0","meaning"]),
                 textInput("iso_m",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "iso_m","meaning"]),
                 textInput("iso_s",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "iso_s","meaning"]),
                 textInput("phi1",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "phi1","meaning"]),
                 textInput("phi2",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "phi2","meaning"])),
          column(5,
                 textInput("psi1",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi1","meaning"]),
                 textInput("psi2",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi2","meaning"]),
                 textInput("c_prop",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_prop","meaning"]),
                 textInput("c_delay_mean",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_delay_mean","meaning"]),
                 textInput("c_delay_cv",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_delay_cv","meaning"]),
                 textInput("psi3",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi3","meaning"]))),
    tabPanel(
      title = "Plot controls",
      value = "plotControls",
      fluidRow(
        column(2,
        sliderInput("titleSize", "Title size:",
                       min = 0, max = 25,
                       value = 20)),
      column(2,
           sliderInput("XtextSize", "X axis title size:",
                       min = 0, max = 25,
                       value = 10)),
      column(2,
           sliderInput("YtextSize", "Y axis title size:",
                       min = 0, max = 25,
                       value = 10)
           ),
      column(2,
        checkboxInput(inputId = "use_logYscale",
                    label = "Use log-y scale",
                    value = FALSE))
      ))),
  width = 12),
    fluidRow(
      actionButton("startButton", "Click to start the simulation")
    ),
  fluidRow(
    tableOutput("summary")
    ),
  fluidRow(
    plotOutput("plot")
  )
))

#Everything else.
  server <- function(input, output){
    #Take a string of a list of numbers and turn that into an actual list of numbers. Let's fix the input for the time-varying transmission option.
    justValues_f <- function(valuesString, mode){
      #Put a comma at the end to make life easier.
      valuesString <- paste(valuesString, ",", sep = "")
      justValueCounter <- 1
      justValueString <- ""
      justValues <- c()
      while (justValueCounter <= stringr::str_length(valuesString)){
        #Grab the current letter
        currentLetter <- substr(valuesString, start = justValueCounter, stop = justValueCounter)
        #If the current letter is not a comma or a space, add it to the justValueString string
        if (currentLetter != ',' && currentLetter != ' '){
          justValueString <- paste(justValueString,  currentLetter, sep = "")
        }
        else if (currentLetter == ',' || currentLetter != ' '){
          if (mode == "values"){
            justValue <- as.numeric(justValueString)
            justValues <- c(justValues, justValue)
          }
          if (mode == "dates"){
            justValue <- lubridate::ymd(justValueString)
            justValues <- c(justValues, justValue)
          }
          if (mode == "symbols"){
            justValue <- justValueString
            justValues <- c(justValues, justValue)
          }
          #Clear the working string after we add the value we want to the list.
          justValueString = ""
        }
        justValueCounter <- justValueCounter + 1
      }
      return(justValues)
    }
    get_factor_timePars <- reactive({
      relValues <- justValues_f(input$timeParsRelativeValues, mode = "values")
      dates <- justValues_f(input$timeParsDates, mode = "dates")
      symbols <- justValues_f(input$timeParsSymbols, mode = "symbols")
      currentPars <- data.frame("Date" = dates, "Symbol" = symbols, "Relative_value" = relValues)
      return(currentPars)
    })
    #Gather input params, if that option is selected.
    gatherValues <- reactive({
          params <- c(input$beta0,
          input$Ca,
          input$Cp,
          input$Cs,
          input$Cm,
          input$alpha,
          input$sigma,
          input$gamma_a,
          input$gamma_s,
          input$gamma_s,
          input$gamma_m,
          input$gamma_p,
          input$rho,
          input$delta,
          input$mu,
          input$N,
          input$E0,
          input$iso_m,
          input$iso_s,
          input$phi1,
          input$phi2,
          input$psi1,
          input$psi2,
          input$psi3,
          input$c_prop,
          input$c_delay_mean,
          input$c_delay_cv)
          return(params)
        })
#In conjunction with the reactive gatherValues, take input from the app and package it in a params_pansim object that run_sim and the like can read.
#Reactivity of the gatherValues makes this work.
#Fighting with namespaces, reactive expressions seem to complicate things a little bit.
      makeParams <- function(){
        params <- gatherValues()
        paramNames <- describe_params(read_params("ICU1.csv"))$symbol
        params <- gatherValues()
        #I don't want numbers to be strings
        params <- vapply(params, function(z) eval(parse(text=z)), numeric(1))
        names(params) <- paramNames
        class(params) <- "params_pansim"
        #Do this after because changing the numbers from strings removes their names.
        return(params)
      }
    observeEvent(input$startButton, {
      output$plot <- renderPlot({
        #If we're using the params we put in.
        if (input$useOwnParams == 1){
          params <- makeParams()
        }
        #If we're using params from a file.
        else{
          params <- read_params(system.file("params", input$fn, package="ShinySimulations"))
        }
        #Throw in proc and obs error as zero by default.
        params <- update(params, c(proc_disp = justValues_f(input$procError, mode = "values"), obs_disp = justValues_f(input$ObsError, mode = "values")))
        if (input$timeChanges == 'Yes'){
          time_pars <- get_factor_timePars()
          sim = run_sim(params, start_date = lubridate::mdy(input$sd), end_date = lubridate::mdy(input$ed), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"), params_timevar = time_pars)
      }
        else{
          sim = run_sim(params, start_date = lubridate::mdy(input$sd), end_date = lubridate::mdy(input$ed), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"))
        }
        #Allow for process and observation error, set to zero by default.
        p <- plot.pansim(sim) +
        labs(title = "Pandemic Simulation") +
        ggplot2::theme(
          plot.title = element_text(color = "black", size = input$titleSize, face = "bold"),
          axis.title.x = element_text(color = "black", size = input$Xsize, face = "bold"),
          axis.title.y = element_text(color = "black", size = input$Ysize, face = "bold"))
        ##Allow for log-y scaling, and adjust the tick-marks and labels accordingly.
        .x <- NULL ## undefined variable check
        if (input$use_logYscale == 1){
          p + scale_y_continuous(trans = log10_trans(),
                                 breaks = trans_breaks("log10", function(x) 10^x),
                                 labels = trans_format("log10", math_format(10^.x)))
        }
        else{
          p
        }
      })
      if (input$useOwnParams == 1){
        params <- makeParams()
      }
      #If we're using params from a file.
      else{
        params <- read_params(system.file("params", input$fn, package="ShinySimulations"))
      }
      output$summary <-renderTable({
        data.frame("Simulation Parameters" = names(summary(params)),"Value" = summary(params))
        })
      })
#Run.
  }
  shinyApp(ui, server)
}
