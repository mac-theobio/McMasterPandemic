##' Simulate ggplot picking colours for a plot.
##'
##' Obtain a list of the colours that ggplot will use for a plot with n variable.
##'
##' @param n the number of variables in the simulated plot. Should match the number of variables in the simulation for this to be useful.
##' @return hues, a list of the colours (specified in hex) that ggplot2 will use.
##' @export
##'
color_list <- function(n) {
  hues = seq(15, 375, length = n + 1)
  colours <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  return(colours)
}

##' Colour checkboxes checkboxes.
##'
##' Generate the necessary HTML/CSS tags to change the colour of a checkbox. Colour curve toggle buttons to match corresponding curves on the graph.
##'
##' @param curves a list of input IDs of the checkboxes to modify the fill of.
##' @param colourList a list of hex obects specifiying the colour to change corresponding checkboxes in curves to.
##' @return a tagsList containing the tags
##' @export
##'
togglePanelColourManager <- function(curves, colourList){
  i <- 1
  listoftags <- list()
  while (i <= length(curves)){
    colour <- colourList[i]
    curve <- curves[i]
    ##Create the text for the HTML tags. The + selects the right divider after the input of the curve.
    theTag <- paste0("#", curve, "+ .state label:after {background-color: ", colour, " !important;}", sep = "")
    listoftags <- c(listoftags, theTag)
    i <- i + 1
  }
  return(listoftags)
}

##' Handle external browser options for the shiny.
##'
##'
##' @import shiny
##' @export
##' @param openinBrowser a logical indicating whether the shiny app should be opened in the browser or not
browserManager <- function(openinBrowser){
  if (openinBrowser){
    ##Use get to circumvent the global binding errors.
    options(shiny.launch.browser = get(".rs.invokeShinyWindowExternal", envir = as.environment("tools:rstudio")))
  }
  else{
    options(shiny.launch.browser = get(".rs.invokeShinyWindowViewer", envir = as.environment("tools:rstudio")))
  }
}

##' Parse inputs from strings.
##'
##'
##' Handle inputs from the time-varying transmission option, which are recieved as strings.
##'
##' @importFrom stringr str_length
##' @importFrom anytime anydate
##' @param valuesString the string we'd like to extract values from.
##' @param mode the type of data encoded in the string: either dates, symbols, or values
##' @export
justValues_f <- function(valuesString, mode){
  #Put a comma at the end to make life easier.
  valuesString <- paste(valuesString, ",", sep = "")
  justValueCounter <- 1
  justValueString <- ""
  justValues <- c()
  while (justValueCounter <= stringr::str_length(valuesString)){
    ##Grab the current letter
    currentLetter <- substr(valuesString, start = justValueCounter, stop = justValueCounter)
    ##If the current letter is not a comma or a space, add it to the justValueString string
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
      ##Clear the working string after we add the value we want to the list.
      justValueString = ""
    }
    justValueCounter <- justValueCounter + 1
  }
  return(justValues)
}

parameter.files <- c("CI_base.csv","CI_updApr1.csv","ICU1.csv", "ICU_diffs.csv")
default.parameter.file <- "ICU1.csv"
default.start.date <- "2020-01-01"
default.dropstates <- c("t","S","E","I","X")


##' Run the McMasterPandemic Shiny
##'
##' @import shiny
##' @importFrom shinythemes shinytheme
##' @importFrom anytime anytime
##' @importFrom ggplot2 scale_y_continuous theme_gray geom_step
##' @importFrom ggplot2 element_text
##' @importFrom directlabels direct.label
##' @importFrom scales log10_trans trans_breaks trans_format math_format
##' @importFrom shinyWidgets prettyCheckbox
##' @param useBrowser Open the shiny in the browser.
##' @return NULL
##' @export
run_shiny <- function(useBrowser = TRUE) {
    ## The ui (user interface) is what the user is shown when running
    ## the shiny.  The ui also gathers input that is the passed to the
    ## server (e.g., filling in boxes or sliders).
    ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
      #Set the title panel to be Heritage Maroon.
      h1(id = "heading", "McMasterPandemic Shiny"),
      tags$style(HTML("#heading{background-color: #7A003C; color: white;}")),
      #Set the colour of the sidebar panel to be Heritage Gold.
      tags$head(tags$style(HTML('#sidebar {background-color: #FDBF57;}'))),
      sidebarLayout(
        sidebarPanel(id = "sidebar",
          fluidRow(
            column(3,
              selectInput("fn",
                               label = "Default parameter file:",
                               choices = parameter.files, selected = default.parameter.file))),
          fluidRow(
            column(5,
                   textInput("sd", "Simulation Start Date (yyyy-mm-dd)",
                             value = default.start.date)),
            column(5,
                   textInput("ed", "Simulation End Date  (yyyy-mm-dd)",
                             value = "2020-08-01"))),
          ## Only show the selector to input parameters if that's selected.
          uiOutput("maintabPanel"),
          ##Colour the checkboxes to match the curves in the plot.
          ##Use the below object to call the HTML tags within the server function.
          htmlOutput("colourManager")),
        mainPanel(
          fluidRow(
            uiOutput("plotColumn"),
              column(2,
                   uiOutput("plotTogglePanel")
            )),
            fluidRow(
              column(2,
                   checkboxInput(inputId = "use_logYscale",
                                 label = "log-y scale",
                                 value = FALSE),
                   br(),
            fluidRow(
                column(6,
                        textOutput("summaryTitle"))),
            br(),
            fluidRow(
                column(6,
                        tableOutput("summary")))))
          )
        )
    )

#Everything else.
  server <- function(input, output, session){
    ## non-standard eval; circumvent 'no visible binding' check
    x <- Date <- Symbol <- Relative_value <- NULL
    ##Render the tab panel server-side to force tab changes the way we'd like, and give us the load-edit functionality we're after.
    output$maintabPanel <- renderUI({
      tabsetPanel(
      ##Use this to force the tab to change.
      id = "tabs",
      ##The parameters panel needs to be the default selected panel. This is so it Shiny can render the panel first. The input slots aren't created until the panel is created.
      ##So keep the default to be the parametersPanel to avoid ugly errors.
      selected = "plotaes",
      tabPanel(
        title = "Time changing transmission rates",
        value = "tcr",
        textOutput("trmsg"),
        br(),
        br(),
        textInput("timeParsDates", label = "Dates of changes, separated by commas", placeholder = "2020-02-20, 2020-05-20, 2020-07-02", value = "2020-02-20, 2020-05-20, 2020-07-02"),
        textInput("timeParsSymbols", label = "Parameter to change on each date", placeholder = "beta0, beta0, alpha", value = "beta0, beta0, alpha"),
        textInput("timeParsRelativeValues", label = "Relative change on each date", placeholder = "1, 1, 1", value = "1, 1, 1"),
        plotOutput("paramsPlot")
        ),
      tabPanel(
        title = "Process and Observation error",
        value = "procObsErr",
        textInput("procError", label = "Process error", value = 0),
        textInput("ObsError", label = "Observation error", value = 0)
      ),
      tabPanel(
        title = "Simulation Parameters",
               value = "parametersPanel",
        radioButtons(inputId = "showAll",
                     label = ("Show all params?"),
                     choices = list("No" = 1, "Yes" = 2),
                     selected = 1),
        conditionalPanel(condition = "input.showAll == 2",
               ##Using names to avoid factors getting passed as inputs to textInput_param.
                lapply(names(read_params("ICU1.csv"))[1:15],
                      FUN = textInput_param),
                lapply(names(read_params("ICU1.csv"))[16:length(names(read_params("ICU1.csv")))],
                      FUN = textInput_param))
               ),
      tabPanel(
        title = "Plot aesthetics",
        value = "plotaes",
        sliderInput("Globalsize", "Text size:",
                           min = 5, max = 45,
                           step = 0.25,
                           value = 25),
        sliderInput("lineThickness", "Line thickness:",
                           min = 0, max = 10,
                           step = 0.25,
                           value = 3),
          radioButtons(inputId = "automaticSize",
                            label = ("Change individual text elements size"),
                            choices = list("No" = 1, "Yes" = 2),
                            selected = 1),
        conditionalPanel(condition = "input.automaticSize == 2",
                            sliderInput("titleSize", "Title size:",
                                        min = 0, max = 25,
                                        value = 20),
                            sliderInput("XtextSize", "X axis title size:",
                                        min = 0, max = 25,
                                        value = 10),
                            sliderInput("YtextSize", "Y axis title size:",
                                        min = 0, max = 25,
                                        value = 10)
        ))
      )})
    ##Force the tab panel to be the parameters panel every time the default parameters drop down is changed.
    observeEvent(input$fn, {
      shiny::updateTabsetPanel(session, "tabs",
                        selected = "parametersPanel"
      )
    })
    get_factor_timePars <- reactive({
      relValues <- justValues_f(input$timeParsRelativeValues, mode = "values")
      dates <- justValues_f(input$timeParsDates, mode = "dates")
      symbols <- justValues_f(input$timeParsSymbols, mode = "symbols")
      currentPars <- data.frame("Date" = dates, "Symbol" = symbols, "Relative_value" = relValues, stringsAsFactors = FALSE)
      return(currentPars)
    })
    dp_1 <- describe_params(read_params("ICU1.csv"))
      output$trmsg <- renderText({"Transmission rate is constant by default but can be changed. You can have any number of parameters."})
    output$plot <- renderPlot({
        ##Detect changes from default values for time-changing transmission rates, and apply these changes in the simulation.
        time_pars <- get_factor_timePars()
        defaultTCParams <- data.frame("Date" = justValues_f(c("2020-02-20, 2020-05-20, 2020-07-02"), mode = "dates"), "Symbol"  = justValues_f(c("beta0, beta0, alpha"), mode = "symbols"), "Relative_value"= justValues_f(c(1, 1, 1), mode = "values"), stringsAsFactors = FALSE)
        ##If the length was changed from the default length, the parameters were definetly changed.
        if(nrow(time_pars) != nrow(defaultTCParams)){
          useTimeChanges <- TRUE
        }
        else{
          ##Are there any elements different from their default values?
          useTimeChanges <- sum(time_pars != defaultTCParams) != 0
        }
        ##Make the params.
        params <- makeParams()
        ##Throw in proc and obs error as zero by default.
        params <- update(params, c(proc_disp = justValues_f(input$procError, mode = "values"), obs_disp = justValues_f(input$ObsError, mode = "values")))
        if (useTimeChanges){
          sim = run_sim(params, start_date = justValues_f(input$sd, mode = "dates"), end_date = justValues_f(input$ed, mode = "dates"), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"), params_timevar = time_pars)
        }
        else{
          sim = run_sim(params, start_date = anytime::anydate(input$sd), end_date = anytime::anydate(input$ed), stoch = c(obs = input$ObsError != "0", proc = input$procError != "0"))
        }
        ##Allow for process and observation error, set to zero by default.
        p <- plot.pansim(sim, drop_states = getDropStates()) + labs(title = "Pandemic Simulation")
        if (input$automaticSize == 2){
        p <- p + ggplot2::theme(
          plot.title = element_text(color = "black", size = input$titleSize),
          axis.title.x = element_text(color = "black", size = input$XtextSize),
          axis.title.y = element_text(color = "black", size = input$YtextSize))
        }
        else{
          p <- p + theme_gray(base_size = input$Globalsize)
        }
        ##Line thickness is applied regardless.
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
        ##Use a positioning method to scale down the size for the direct labels, only label the last points for each graph, and add spacing beween the lines and the labels.
        p <- direct.label(p, list("last.points", cex = input$Globalsize/15, dl.trans(x = x + 0.05)))
        p
      })
      ##Take input and package it in a params_pansim object that run_sim and the like can read.
      makeParams <- function(){
        params <- c()
        for (param in describe_params(read_params("ICU1.csv"))$symbol){
          ##Grab the value of the input slot.
          paramValue <- eval(parse(text = paste0("input$", param)))
          ##If it's null beacause the panel hasn't been initialized yet, grab a default value to use a a placeholder.
          if (is.null(paramValue)){
            params <- c(params, read_params("ICU1.csv")[param])
          }
          else{
            ##Otherwise, package and make the params the way we'd like them. Changed for conciceness.
            params <- c(params, paramValue)
          }
        }
        paramNames <- describe_params(read_params("ICU1.csv"))$symbol
        ##I don't want numbers to be strings
        params <- vapply(params, function(z) eval(parse(text=z)), numeric(1))
        names(params) <- paramNames
        class(params) <- "params_pansim"
        ##Do this after because changing the numbers from strings removes their names.
        return(params)
      }
      loadParams <- function(param){
        Inputparams <- read_params(input$fn)
        numMissing <- sum(is.na(Inputparams))
        ##Also account for the fact that data might just be missing from the file entirely and not recorded as NA values.
        if (numMissing != 0 || length(Inputparams) < 26){
          #If the parameters file is missing info, fill in defaults for the missing values.
          DefaultParams <- read_params("ICU1.csv")
          NonMissingparams <- Inputparams[!is.na(Inputparams)]
          NonMissingparamNames <- names(NonMissingparams)
          DefaultParams[NonMissingparamNames] <- NonMissingparams
          params <- DefaultParams
        }
        else{
          params <- Inputparams
        }
        return(params[param])
      }
        output$summary <-renderTable({
            params <- makeParams()
            params <- update(params, c(proc_disp = justValues_f(input$procError, mode = "values"), obs_disp = justValues_f(input$ObsError, mode = "values")))
          data.frame("Symbol" = describe_params(summary(read_params("ICU1.csv")))$symbol, "Meaning" = describe_params(summary(read_params("ICU1.csv")))$meaning,"Value" = summary(params))
        })
        output$paramsPlot <- renderPlot({
          parameterChanges <- get_factor_timePars()
          ##Make the plot more descriptive by adding begining and ending points for every symbol.
          for (symbol in unique(parameterChanges$Symbol)){
            missingBegining <- sum(parameterChanges[parameterChanges$Symbol == symbol, "Date"] == input$sd) == 0
            missingEnd <- sum(parameterChanges[parameterChanges$Symbol == symbol, "Date"] ==  input$ed) == 0
            ##If we're missing the begining, populate it with the starting date and a default value of one.
            if (missingBegining){
              parameterChanges <- rbind(data.frame("Date" = justValues_f(input$sd, mode = "dates"), "Symbol" = symbol, "Relative_value" = 1, stringsAsFactors = FALSE), parameterChanges)
            }
            else {
            }
            ##Order by date to avoid confusion.
            parameterChanges <- parameterChanges[order(parameterChanges$Date),]
            ##If we're missing the end, populate the last date with the last relative value we put in.
            if (missingEnd){
              parameterChanges <- rbind(data.frame("Date" = justValues_f(input$ed, mode = "dates"), "Symbol" = symbol, "Relative_value" = parameterChanges[parameterChanges$Symbol == symbol,][nrow(parameterChanges[parameterChanges$Symbol == symbol,]), "Relative_value"], stringsAsFactors = FALSE), parameterChanges)
            }
            else {
            }
            ##Order by date again to avoid confusion.
            parameterChanges <- parameterChanges[order(parameterChanges$Date),]
          }
          p <- ggplot(parameterChanges,aes(anytime::anydate(Date), Relative_value, colour=Symbol)) + geom_step(size = 2) +
            theme_gray(base_size = input$Globalsize)
          p <- p + geom_vline(xintercept=parameterChanges$Date,lty=2) + labs(title = "Changes over time", x = "Date", y = "Relative value")
          p <- direct.label(p, list("last.points", cex = input$Globalsize/15))
          p
        })
        textInput_param <- function(param, dp = dp_1){
          maxVal <- 10
          return(sliderInput(param,
                           label = dp[dp$symbol == param,"meaning"],
                           value = loadParams(param),
                           min = 0,
                           max = maxVal,
                          step = 0.1))
          }
        ##Manage the states to drop.
        getDropStates <- function(){
          default.sim <- run_sim(read_params("ICU1.csv"))
          couldDropStates <- setdiff(colnames(default.sim)[2:length(default.sim)], default.dropstates)
          for (state in couldDropStates){
            stateVal <- eval(parse(text = paste0("input$", state)))
            ##Catch loading errors.
            if (is.null(stateVal)){
              return(default.dropstates)
            }
            else {
            }
            ##2 indicates we don't want to show the drop state.
            if (!stateVal){
              default.dropstates <- c(default.dropstates, state)
            }
            else{
            }
          }
          return(default.dropstates)
        }
        ##Create checkbuttons to display plots or not.
        checkButton_curve <- function(curve){
          #Don't show cum rep by default
          if (curve == "cumRep" || curve == "foi" || curve == "R"){
            showByDefault <- FALSE
          }
          else{
            showByDefault <- TRUE
          }
          theLabel <- paste0(curve)
          ##Change curve labels to be more intuitive.
          if (curve == "foi"){
            theLabel <- "force of \n infection"
          }
          if (curve == "cumRep"){
            theLabel <- "cumulative reports"
          }
          if (curve == "hosp"){
            theLabel <- "hospital \n admissions"
          }
          if (curve == "H"){
            theLabel <- "hospitalized"
          }
          return(shinyWidgets::prettyCheckbox(curve,
                           label = theLabel,
                           value = showByDefault))
        }
        create_togglePanel <- function(){
          ##Exclude the date as that's not a curve.
          defsim <- run_sim(read_params("ICU1.csv"))
          curves <- as.vector(colnames(defsim)[2:length(defsim)])
          ##Ignore curves that we're never going to show.
          curves <- setdiff(curves, c("t","S","E","I","X"))
            column(2,
                 lapply(curves,
                        FUN = checkButton_curve))
        }
        output$plotColumn <- renderUI({
          column(9,
                 plotOutput("plot"))
        })
        ##Panel to toggle curves showing
        output$plotTogglePanel <- renderUI({
          create_togglePanel()
        })
        output$colourManager <- renderUI({
          ##Grab all the curves.
          defsim <- run_sim(read_params("ICU1.csv"))
          curvesList <- as.vector(colnames(defsim)[2:length(defsim)])
          ##Remove the ones we're going to drop.
          curvesList <- setdiff(curvesList, getDropStates())
          ##Grab the corresponding tags
          theTags <- togglePanelColourManager(curves = curvesList, colourList = color_list(length(curvesList)))
          ##Render the tags.
          lapply(theTags, function(tag){
            return(tags$style(tag))
          })
          })
        output$summaryTitle <- renderText({"Summary characteristics"})
  }

  ##Set the viewing options first.
  browserManager(useBrowser)
  ##Run the shiny app. the default value of launch.browser looks for the option set by browserManager.
  shiny::runApp(appDir = shinyApp("ui" = ui, "server" = server), launch.browser = getOption("shiny.launch.browser", interactive()))
}
