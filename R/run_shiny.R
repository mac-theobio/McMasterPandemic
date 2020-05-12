#' Run the McMasterPandemic Shiny
#'
#' @importFrom shiny fluidPage titlePanel mainPanel fluidRow radioButtons h3 conditionalPanel
#' @importFrom shiny column selectInput textInput tabsetPanel tabPanel checkboxInput sliderInput
#' @importFrom shiny actionButton tableOutput plotOutput reactive observeEvent renderPlot
#' @importFrom shiny renderTable shinyApp uiOutput
#' @importFrom anytime anytime
#' @importFrom ggplot2 scale_y_continuous theme_gray
#' @importFrom directlabels direct.label
#' @importFrom scales log10_trans trans_breaks trans_format math_format
#' @importFrom ggplot2 element_text
#' @return NULL
#' @export
run_shiny <- function(){
#How I want stuff to look. Also collect input and things.
  ui <- fluidPage(
    titlePanel("McMasterPandemic Shiny"),
    mainPanel(
      fluidRow(
        #Only show the selector for the parameter file to upload if that's selected.
          column(2,
                 selectInput("fn",
                             label = "Chose file containing parameters",
                             choices = c("CI_base.csv",
                                         "CI_updApr1.csv",
                                         "ICU1.csv",
                                         "ICU_diffs.csv"), selected = "ICU1.csv")),
        fluidRow(
          column(2,
                 textInput("sd", "Simulation Start Date (ymd)", value = "2020-01-01")),
          column(2,
                 textInput("ed", "Simulation End Date  (ymd)", value = "2020-08-01"))
        ),
  #Only show the selector to input parameters if that's selected.
    tabsetPanel(
      selected = "plotControls",
      tabPanel(
        title = "Time changing transmission rates",
        value = "tcr",
        column(8,
          textInput("timeParsDates", label = "Enter dates of changes here, in ymd format, separated by commas", placeholder = "2020-02-20, 2020-05-20, 2020-07-02", value = "2020-02-20, 2020-05-20, 2020-07-02"),
          textInput("timeParsSymbols", label = "Enter the corresponding symbol for each date that you'd like to change here, separated by commas", placeholder = "beta0, beta0, alpha", value = "beta0, beta0, alpha"),
          textInput("timeParsRelativeValues", label = "Enter relative value changes here, separated by commas", placeholder = "1, 1, 1", value = "1, 1, 1")
      ),
      column(8,
             plotOutput("paramsPlot"))),
      tabPanel(
        title = "Process and Observation error",
        value = "procObsErr",
        textInput("procError", label = "Enter the process error", value = 0),
        textInput("ObsError", label = "Enter the observation error", value = 0)
      ),
        tabPanel(
          title = "Simulation Parameters",
          value = "parametersPanel",
          column(5,
                 textInput("beta0",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "beta0","meaning"],
                           value = read_params(system.file("params", "ICU1.csv", package="ShinySimulations"))["beta0"]),
                 textInput("Ca",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Ca","meaning"],
                           value = read_params("ICU1.csv")["Ca"]),
                 textInput("Cp",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cp","meaning"],
                           value = read_params("ICU1.csv")["Cp"]),
                 textInput("Cs",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cs","meaning"],
                           value = read_params("ICU1.csv")["Cs"]),
                 textInput("Cm",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "Cm","meaning"],
                           value = read_params("ICU1.csv")["Cm"])),
          column(5,
                 textInput("alpha",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "alpha","meaning"],
                           value = read_params("ICU1.csv")["alpha"]),
                 textInput("sigma",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "sigma","meaning"],
                           value = read_params("ICU1.csv")["sigma"]),
                 textInput("gamma_a",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_a","meaning"],
                           value = read_params("ICU1.csv")["gamma_a"]),
                 textInput("gamma_s",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_s","meaning"],
                           value = read_params("ICU1.csv")["gamma_s"]),
                 textInput("gamma_m",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_m","meaning"],
                           value = read_params("ICU1.csv")["gamma_m"])),
          column(5,
                 textInput("gamma_p",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "gamma_p","meaning"],
                           value = read_params("ICU1.csv")["gamma_p"]),
                 textInput("rho",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "rho","meaning"],
                           value = read_params("ICU1.csv")["rho"]),
                 textInput("delta",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "delta","meaning"],
                           value = read_params("ICU1.csv")["delta"]),
                 textInput("mu",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "mu","meaning"],
                           value = read_params("ICU1.csv")["mu"]),
                 textInput("N",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "N","meaning"],
                           value = read_params("ICU1.csv")["N"])),
          column(5,
                 textInput("E0",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "E0","meaning"],
                           value = read_params("ICU1.csv")["E0"]),
                 textInput("nonhosp_mort",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "nonhosp_mort","meaning"],
                           value = read_params("ICU1.csv")["nonhosp_mort"]),
                 textInput("iso_m",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "iso_m","meaning"],
                           value = read_params("ICU1.csv")["iso_m"]),
                 textInput("iso_s",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "iso_s","meaning"],
                           value = read_params("ICU1.csv")["iso_s"]),
                 textInput("phi1",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "phi1","meaning"],
                           value = read_params("ICU1.csv")["phi1"]),
                 textInput("phi2",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "phi2","meaning"],
                           value = read_params("ICU1.csv")["phi2"])),
          column(5,
                 textInput("psi1",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi1","meaning"],
                           value = read_params("ICU1.csv")["psi1"]),
                 textInput("psi2",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi2","meaning"],
                           value = read_params("ICU1.csv")["psi2"]),
                 textInput("psi3",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "psi3","meaning"],
                           value = read_params("ICU1.csv")["psi3"]),
                 textInput("c_delay_mean",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_delay_mean","meaning"],
                           value = read_params("ICU1.csv")["c_delay_mean"]),
                 textInput("c_delay_cv",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_delay_cv","meaning"],
                           value = read_params("ICU1.csv")["c_delay_cv"]),
                 textInput("c_prop",
                           label = describe_params(read_params("ICU1.csv"))[describe_params(read_params("ICU1.csv"))$symbol == "c_prop","meaning"],
                           value = read_params("ICU1.csv")["c_prop"]))),
      tabPanel(
        title = "Plot controls",
        value = "plotControls",
        fluidRow(
          column(2,
        checkboxInput(inputId = "use_logYscale",
                    label = "Use log-y scale",
                    value = FALSE)),
        column(2,
        radioButtons(inputId = "use_directLabels",
                      label = ("Use direct labels or a legend"),
                      choices = list("Use direct labels" = 1, "Use a legend" = 2),
                      selected = 1))
      )),
    tabPanel(
      title = "Plot aesthetics",
      value = "plotaes",
      column(2,
      sliderInput("Globalsize", "Text size:",
                  min = 5, max = 45,
                  step = 0.25,
                  value = 12),
      sliderInput("lineThickness", "Line thickness:",
                  min = 0, max = 10,
                  step = 0.25,
                  value = 3)),
      column(2,
      radioButtons(inputId = "automaticSize",
                   label = ("Change individual text elements size"),
                   choices = list("No" = 1, "Yes" = 2),
                   selected = 1)),
      conditionalPanel(condition = "input.automaticSize == 2",
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
                         value = 10))
        ))),
  width = 12),
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
            justValue <- anytime::anydate(justValueString)
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
      #Render the parameter tab server-side, to enable reactivity and to let us edit things.

      #In conjunction with the reactive gatherValues, take input from the app and package it in a params_pansim object that run_sim and the like can read.
      makeParams <- function(){
        params <- c(input$beta0,
                    input$Ca,
                    input$Cp,
                    input$Cs,
                    input$Cm,
                    input$alpha,
                    input$sigma,
                    input$gamma_a,
                    input$gamma_s,
                    input$gamma_m,
                    input$gamma_p,
                    input$rho,
                    input$delta,
                    input$mu,
                    input$N,
                    input$E0,
                    input$nonhosp_mort,
                    input$iso_m,
                    input$iso_s,
                    input$phi1,
                    input$phi2,
                    input$psi1,
                    input$psi2,
                    input$psi3,
                    input$c_prop,
                    input$c_delay_mean,
                    input$c_delay_cv,
                    input$procError)
        paramNames <- describe_params(read_params("ICU1.csv"))$symbol
        #I don't want numbers to be strings
        params <- vapply(params, function(z) eval(parse(text=z)), numeric(1))
        names(params) <- paramNames
        class(params) <- "params_pansim"
        #Do this after because changing the numbers from strings removes their names.
        return(params)
      }
      output$plot <- renderPlot({
        #Make the params
        params <- makeParams()
        #Throw in proc and obs error as zero by default.
        params <- update(params, c(proc_disp = justValues_f(input$procError, mode = "values"), obs_disp = justValues_f(input$ObsError, mode = "values")))
        #Detect changes from default values for time-changing transmission rates, and apply these changes in the simulation.
        time_pars <- get_factor_timePars()
        defaultTCParams <- data.frame("Date" = justValues_f(c("2020-02-20, 2020-05-20, 2020-07-02"), mode = "dates"), "Symbol"  = justValues_f(c("beta0, beta0, alpha"), mode = "symbols"), "Relative_value"= justValues_f(c(1, 1, 1), mode = "values"))
        if (length(time_pars$Date) != length(defaultTCParams$Date)){
          #To avoid errors when the default parameters are a different length than the input ones, such as when we don't want three time changing parameter value changes, if the default values were shortened, that's still counted as a change.
          useTimeChanges <- TRUE
        }
        else{
          #Do these individually to account for different sizes of data.frames.
          #If the number of different elements between the not equal to zero, something has been changed.
          datesDifferent <- sum(time_pars$Date != defaultTCParams$Date) != 0
          symbolsDifferent <- sum(time_pars$Symbol != defaultTCParams$Symbol) != 0
          valuesDifferent <- sum(time_pars$Relative_value != defaultTCParams$Relative_value) != 0
          useTimeChanges <- datesDifferent || symbolsDifferent || valuesDifferent
        }
        if (useTimeChanges){
          sim = run_sim(params, start_date = anytime::anydate(input$sd), end_date = anytime::anydate(input$ed), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"), params_timevar = time_pars)
      }
        else{
          sim = run_sim(params, start_date = anytime::anydate(input$sd), end_date = anytime::anydate(input$ed), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"))
        }
        #Allow for process and observation error, set to zero by default.
        p <- plot.pansim(sim) + labs(title = "Pandemic Simulation")
        if (input$automaticSize == 2){
        p <- p + ggplot2::theme(
          plot.title = element_text(color = "black", size = input$titleSize),
          axis.title.x = element_text(color = "black", size = input$XtextSize),
          axis.title.y = element_text(color = "black", size = input$YtextSize))
        }
        else{
          p <- p + theme_gray(base_size = input$Globalsize)
        }
        #Line thickness is applied regardless.
        p <- p + geom_line(size = input$lineThickness)
        ##Allow for log-y scaling, and adjust the tick-marks and labels accordingly.
        .x <- NULL ## undefined variable check
        if (input$use_logYscale == 1){
          p <- p + scale_y_continuous(trans = log10_trans(),
                                 breaks = trans_breaks("log10", function(x) 10^x),
                                 labels = trans_format("log10", math_format(10^.x)))
        }
        else{
        }
        if (input$use_directLabels == 1){
          #Use a positioning metho to scale down the size for the direct labels, only label the last points for each graph, and add spacing beween the lines and the labels.
          p <- direct.label(p, list("last.points", cex = input$Globalsize/15, dl.trans(x = x + 0.05)))
          p
        }
        else{
          p
        }
      })
        output$summary <-renderTable({
            params <- makeParams()
            params <- update(params, c(proc_disp = justValues_f(input$procError, mode = "values"), obs_disp = justValues_f(input$ObsError, mode = "values")))
        data.frame("Symbol" = describe_params(summary(read_params("ICU1.csv")))$symbol, "Meaning" = describe_params(summary(read_params("ICU1.csv")))$meaning,"Value" = summary(params))
        })
        output$paramsPlot <- renderPlot({
          parameterChanges <- get_factor_timePars()
          #We want all the symbols to have a line starting from the begining of the graph regardless of whether that was actually specified or not.
          for (symbol in parameterChanges$Symbol){
            if (sum(parameterChanges$Symbol == symbol) < nrow(parameterChanges)){
              parameterChanges <- rbind(data.frame("Date" = parameterChanges[1, "Date"], "Symbol" = symbol, "Relative_value" = 1), parameterChanges)
            }
            else{
            }
          }
          p <- ggplot(parameterChanges,aes(anytime::anydate(Date), Relative_value, colour=Symbol)) + geom_line(size = 2)
          p <- p + geom_vline(xintercept=parameterChanges$Date,lty=2) + labs(title = "Parameter changes over time", x = "Date", y = "Relative value")
          p <- direct.label(p, list("last.points", cex = input$Globalsize/15, dl.trans(x = x + 0.05)))
          p
        })
#Run.
  }
  shinyApp(ui, server)
}
