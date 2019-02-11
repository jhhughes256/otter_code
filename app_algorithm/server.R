# server.R script for OTTER application
# Reactive objects (i.e., those dependent on widget input) are written here
# ------------------------------------------------------------------------------
server <- function(input, output, session) {

## Setup reactiveValues
# n - counter to determine number of boxes for reactive ui "time, conc" input
# io - on/off switch for processes to be used with data
# df - reactive data.frame that saves values form the "time, conc" input
  rv <- reactiveValues(
    n = 1, fit = FALSE, opt = FALSE,
    df = data.frame(
      time = double(1),
      conc = double(1)
    )  #data.frame
  )  #reactiveValues

## Observe Events
# When pressing add: check values in df are numbers
# If values are numbers then remove default renderUI values and increase values
# Also update dataframe for input - reference values for time, conc input
  observeEvent(input$addsamp, {
    rv$fit <- FALSE
    rv$opt <- FALSE
    rv$n <- rv$n + 1
    df <- ldply(seq_len(rv$n), function(i) {  #save current values
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
    rv$df <- arrange(df, time)
  })  #observeEvent

# When pressing remove: check to see if > 1 input box, if so decrease values
# Also update dataframe for input - reference values for time,dose,conc input
  observeEvent(input$remsamp, {
    rv$fit <- FALSE
    rv$opt <- FALSE
    if (rv$n > 1) {  
      rv$n <- rv$n - 1
    }  #if
    df <- ldply(seq_len(rv$n), function(i) {  #save current values 
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
    rv$df <- arrange(df, time)
  })  #observeEvent

# When pressing save: update dataframe for input
# Also print to console the state of the data (for debug purposes)
  observeEvent(input$savesamp, {
    df <- ldply(seq_len(rv$n), function(i){  #save current values 
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
    df$time[is.na(df$time)] <- 0
    df$conc[is.na(df$conc)] <- 0
    rv$df <- arrange(df, time)
    rv$fit <- TRUE
    print(rv$df)
    print(input$absorp)
  })  #observeEvent

## UI code for renderUI
# Shows "time, conc" ui as dictated by "n"
# Set up to allow saving of values as they are put in to prevent deletion of 
# progress when Add, Remove is pressed
  inputBoxes <- reactive({
    n <- rv$n
    df <- rv$df
    if(n > 0) {
      llply(seq_len(n), function(i) {
        div(class = "row-fluid",
        # Time
          div(class = "MyClass",  #TIME
            textInput(
              paste0("time", i),
              ifelse(i == 1, "Time", NA),
              df[i, 1]
            )  #textInput
          ),  #input$time-i
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"
          )),  #make inline
          tags$head(tags$style(type = "text/css", 
            paste0("#time", i, " {max-width: 110px}")
          )),  #change width
        # Conc
          div(class = "MyClass",
            textInput(
              paste0("conc", i),
              ifelse(i == 1, "Concentration", NA),
              df[i, 2]
            )  #textInput
          ),  #input$conc-i
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"
          )),  #make inline
          tags$head(tags$style(type = "text/css", 
            paste0("#conc", i, " {max-width: 110px}")
          ))  #change width
        )  #div
      })  #llply
    }  #if
  })  #textboxes

# Finally render the dataframe input UI so it can be used in ui body
  output$samptimeui <- renderUI({inputBoxes()})
  
# Reactive sumexp
  Rsumexp <- reactive({
    if (rv$fit) {
      df <- rv$df
      absorp <- as.logical(as.numeric(input$absorp))
      if (absorp) df[1, 1] <- 0
      optimResult <- optim.sumexp.new(df, oral = absorp)
      best.sumexp.aic(optimResult)
    } else {
      NULL
    }
  })  #Rsumexp
  
  RsumexpPred <- reactive({
    if (rv$fit) {
      rv$opt <- TRUE
      times <- seq(0, tail(rv$df[,1], 1), length.out = 97)
      data.frame(
        time = times,
        conc = pred.sumexp(Rsumexp()$sumexp, times)
      )
    } else {
      NULL
    }
  })  #RsumexpPred
  
  Rinterv <- reactive({
    tmaxStatus <- F
    tlastStatus <- F
    if (!is.null(input$adjunct)) {
      if (input$adjunct %in% "tmax") tmaxStatus <- T
      if (input$adjunct %in% "tlast") tlastStatus <- T
    }
    if (rv$opt) {
      nobs <- input$nobs
      tlast <- ifelse(tlastStatus, obs.tlast.lam(rv$df), tail(rv$df[,1], 1))
      times <- seq(0, tlast, length.out = nobs)
      optim.interv.dtmax(Rsumexp()$sumexp, times, tmax = tmaxStatus)
    } else {
      NULL
    }
  })  #Rinterv
  
# Reactive plot
# Shows data points when rv$io == TRUE
# Shows sumexp curve once available
  Rplot <- reactive({
    p <- NULL
    p <- ggplot()
    p <- p + geom_point(aes(x = time, y = conc), data = rv$df, na.rm = TRUE)
    if (rv$fit) p <- p + geom_line(aes(x = time, y = conc), data = RsumexpPred())
    if (rv$opt) p <- p + geom_vline(xintercept = Rinterv()$times, 
      colour = "green4", linetype = "dashed")
    p <- p + xlab("Time")
    p <- p + ylab("Concentration")
    p
  })  #Rplot
  
  output$plot <- renderPlot({Rplot()})

# Reactive text
  output$times <- renderUI({
    cssStyle <- paste0("color:green; text-align: center; font-size:", 26/input$nobs, "vw")
    if (rv$opt) {
      timeString <- paste(round(Rinterv()$times/0.5, 0)*0.5, collapse = ", ")
      div(timeString, style = cssStyle)
    } else {
      div("Waiting...", style = cssStyle)
    }
  })  #output$times

# Close the R session when Chrome closes
  session$onSessionEnded(function() {
    stopApp()
  })  #endsession
  
# Open debug console for R session
  # observe(label = "console", {
  #   if(input$console != 0) {
  #     options(browserNLdisabled = TRUE)
  #     # saved_console <- ".RDuetConsole"
  #     # if (file.exists(saved_console)) load(saved_console)
  #     isolate(browser())
  #     # save(file = saved_console, list = ls(environment()))
  #   }
  # })
  
}  #server