if (!require('devtools')){
  install.packages("devtools")
  library(devtools)
}

if (!require('shiny')){
  install.packages("shiny")
  library(shiny)
}

if (!require('MiRKATMC')){
  install_github("Zhiwen-Owen-Jiang/MiRKATMC")
  library(MiRKATMC)
}



ui <- fluidPage(
  titlePanel('MiRKAT-MC'),

  sidebarLayout(
    sidebarPanel(
      textInput(inputId = 'formula',
                label = 'Formula of Fix Effect',
                placeholder = 'Eg. outcome ~ age + sex'),
      textInput(inputId = 'random',
                label = 'Formula of Random Effect',
                value = '',
                placeholder = 'Eg. ~ 1 + time | ID or ~ 1 | ID'),
      selectInput(inputId = 'data.type',
                  label = 'Data Type',
                  choices = c('nominal', 'ordinal')),
      fileInput(inputId = 'Ks',
                label = 'RDA File of Kernel Matrix',
                accept = c('.rda', '.RData', '.RDA'),
                placeholder = 'Eg. test_kernel_list.rda'),
      fileInput(inputId = 'data',
                label = 'RDA File of Dataset',
                accept = c('.rda', '.RData', '.RDA'),
                placeholder = 'Eg. test_data.rda'),
      actionButton(inputId = 'go', label = 'Go')
    ),

      mainPanel(tableOutput('results'))
  )

)



server <- function(input, output){
  repeatInput <- eventReactive(input$go, {
    req(input$Ks)
    req(input$data)
    kernel <- get(load(input$Ks$datapath))
    data <- get(load(input$data$datapath))
    data.type <- input$data.type
    fix_formula <- as.formula(input$formula)
    if (input$random == ''){
      random_formula = NULL
    }else{
      random_formula = as.formula(input$random)
    }
    return(list(Ks = kernel, data = data, formula = fix_formula, random = random_formula,
           data.type = data.type))
  })

  output$results <- renderTable({
    input <- repeatInput()
    pvalues <- MiRKATMC(formula = input$formula, random = input$random,
                        Ks = input$Ks, data.type = input$data.type,
                        data = input$data)
    pvalues <- data.frame(as.character(pvalues), row.names = names(pvalues))
    return(t(pvalues))
  })

}

shinyApp(ui, server)




# generate some example

# set.seed(123)
# test.data <- data.frame(outcome = as.factor(sample(4, 100, replace = TRUE)),
#                         ID = gl(20, 5), time = rep(1:5, 20), age = rnorm(n = 100, mean = 30, sd = 5),
#                         sex = rbinom(100, 1, 1/2))
# D1 <- matrix(rbinom(10000, 2, 0.05), 100, 100)
# K1 <- crossprod(D1) # kernel matrix
# D2 <- matrix(rbinom(10000, 2, 0.1), 100, 100)
# K2 <- crossprod(D2) # kernel matrix
# K <- list(kernel1 = K1, kernel2 = K2)
# K_no_name <- list(K1, K2)
#
# save(test.data, file = 'test_data.rda')
# save(K, file = 'test_kernel_list.rda')
# save(K1, file = 'test_kernel_matrix.rda')
