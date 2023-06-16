library(shiny)
library(OncoSimulR)
library(readxl)
library(DT)
library(periscope)
allFitnessEffectsMatrix <- matrix()

# Function which receives the arguments to define the fitness of our simulation
# and the different parameters that characterize the cellular interactions.Then
# returns the result of the simulation performed by the oncoSimulIndiv function
oncoSimulIndivComplete <- function(allFitnessEffectsData, frequencyDependentFitness, mu, onlyCancer, initSize, finalTime, model, detectionSize, detectionDrivers,
                                   mutationPropGrowth, maxMemory, maxWallTime, maxNumTries,
                                   errorHitWallTime, errorHitMaxTries, AND_DrvProbExit) {
  print(allFitnessEffectsData)
  fitness <- allFitnessEffects(genotFitness = allFitnessEffectsData, frequencyDependentFitness = frequencyDependentFitness, frequencyType = "rel")
  result <- oncoSimulIndiv(fitness,
    model = model, mu = mu, initSize = initSize,
    finalTime = finalTime,
    onlyCancer = onlyCancer,
    detectionProb = "default",
    detectionSize = detectionSize,
    detectionDrivers = detectionDrivers,
    mutationPropGrowth = mutationPropGrowth,
    max.memory = maxMemory,
    max.wall.time = maxWallTime,
    max.num.tries = maxNumTries,
    errorHitWallTime = errorHitWallTime,
    errorHitMaxTries = errorHitMaxTries,
    AND_DrvProbExit = AND_DrvProbExit
  )

  result
}

# Function which given a file with table data, returns a table 
fileToTable <- function(file) {
  ext <- tools::file_ext(file$datapath)
  req(file)
  validate(need(ext == "txt", "Please upload a txt file"))
  read.table(file$datapath, header = TRUE, sep = "\t", fileEncoding="latin1")
}


detectBlankInput <- function(input){
  # When input$n is 3, filename is ./images/image3.jpeg
  filenameTick <- "./images/tick.jpeg"
  filenameCross <- "./images/x.jpeg"
  # Return a list containing the filename and alt text
  if (is.na(input) || input == "") {
    filename <- filenameCross
  } else {
    filename <- filenameTick
  }
  
  list(
    src = filename,
    height = "20px",
    width = "23px"
  )
}

# Server function that will manage how they interact the inputs and outputs(Placed in a UI file).
oncoSimulIndivLogic <- function(input, output) {
  # observer dedicated to cell edition of AllFitnessEffects table
  observeEvent(input$allFitnessEffectsDTEditable_cell_edit, {
    allFitnessEffectsMatrix <<- editData(allFitnessEffectsMatrix, input$allFitnessEffectsDTEditable_cell_edit, "allFitnessEffectsDTEditable", rownames = FALSE)
    
    })

  # Observer that saves input data for the user to decide when to run the simulation
  observeEvent(list(
    input$oncoSimulIndiv_frequencyDepedentFitnessCheckbox, input$oncoSimulIndiv_mutationRate, input$oncoSimulIndiv_onlyCancer,
    input$oncoSimulIndiv_initSize, input$oncoSimulIndiv_finalTime, input$oncoSimulIndiv_selectModel, input$oncoSimulIndiv_detectionSize,
    input$oncoSimulIndiv_detectionDrivers, input$oncoSimulIndiv_mutationPropGrowth, input$oncoSimulIndiv_maxMemory, input$oncoSimulIndiv_maxWallTime,
    input$oncoSimulIndiv_maxNumTries, input$oncoSimulIndiv_errorHitWallTime, input$oncoSimulIndiv_errorHitMaxTries, input$oncoSimulIndiv_AND_DrvProbExit,
    input$oncoSimulIndiv_plotType, input$oncoSimulIndiv_plotShow
  ), {
    frequencyDepedentFitnessCheckboxInstance <<- input$oncoSimulIndiv_frequencyDepedentFitnessCheckbox
    muInstance <<- input$oncoSimulIndiv_mutationRate
    onlyCancerInstance <<- input$oncoSimulIndiv_onlyCancer
    initSizeInstance <<- input$oncoSimulIndiv_initSize
    finalTimeInstance <<- input$oncoSimulIndiv_finalTime
    selectModelInstance <<- input$oncoSimulIndiv_selectModel
    detectionSizeInstance <<- input$oncoSimulIndiv_detectionSize
    detectionDriversInstance <<- input$oncoSimulIndiv_detectionDrivers
    mutationPropGrowthInstance <<- input$oncoSimulIndiv_mutationPropGrowth
    maxMemoryInstance <<- input$oncoSimulIndiv_maxMemory
    maxWallTimeInstance <<- input$oncoSimulIndiv_maxWallTime
    maxNumTriesInstance <<- input$oncoSimulIndiv_maxNumTries
    errorHitWallTimeInstance <<- input$oncoSimulIndiv_errorHitWallTime
    errorHitMaxTriesInstance <<- input$oncoSimulIndiv_errorHitMaxTries
    AND_DrvProbExitInstance <<- input$oncoSimulIndiv_AND_DrvProbExit
    oncosimulIndivPlotTypeInstance <<- input$oncoSimulIndiv_plotType
    oncosimulIndivPlotShow <<- input$oncoSimulIndiv_plotShow
  })

  # observer dedicated to execute the simulation and plot results
  observeEvent(input$oncoSimulIndiv_runExecutionButton, {
    output$oncoSimulIndiv_plot <- renderPlot(plot(
      (plotData <<- oncoSimulIndivComplete(
        allFitnessEffectsData = allFitnessEffectsDecider(input),
        frequencyDependentFitness = frequencyDepedentFitnessCheckboxInstance,
        mu = muInstance,
        onlyCancer = onlyCancerInstance,
        initSize = initSizeInstance,
        finalTime = finalTimeInstance,
        model = selectModelInstance,
        detectionSize = detectionSizeInstance,
        detectionDrivers = detectionDriversInstance,
        mutationPropGrowth = mutationPropGrowthInstance,
        maxMemory = maxMemoryInstance,
        maxWallTime = maxWallTimeInstance,
        maxNumTries = maxNumTriesInstance,
        errorHitWallTime = errorHitWallTimeInstance,
        errorHitMaxTries = errorHitMaxTriesInstance,
        AND_DrvProbExit = AND_DrvProbExitInstance
      )),
      type = oncosimulIndivPlotTypeInstance,
      show = oncosimulIndivPlotShow,
      addtot = TRUE
    ))
  })


  # reactive dedicated to decide which type of table is picked for AllFitnessEffects simulation and trasnform the data in case it is needed 
  allFitnessEffectsDecider <- function(input){
    if (input$oncoSimulIndiv_allFitnessEffectsType == 1) {
      print(fileToTable(input$oncoSimulIndiv_allFitnessEffectsFileInput))
      fileToTable(input$oncoSimulIndiv_allFitnessEffectsFileInput)
    } else {
      if (input$oncoSimulIndiv_allFitnessEffectsDataTableType == 1) {
        print(allFitnessEffectsMatrix)
        as.data.frame(allFitnessEffectsMatrix)
      } else {
        if (input$oncoSimulIndiv_frequencyDepedentFitnessCheckbox == FALSE) {
          d9 <- as.data.frame(allFitnessEffectsMatrix)
          D <- transform(d9, Birth = as.numeric(Birth))
        } else {
          colnames(allFitnessEffectsMatrix) <- c("Genotype", "Fitness")
          d9 <- as.data.frame(allFitnessEffectsMatrix)
          D <- transform(d9, Fitness = as.character.numeric_version(Fitness))
        }
      }
    }
  }
  
  ############Status REGION#############
  # Dedicated to show if a parameter has a value or not
  output$oncoSimulIndiv_modelStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_selectModel)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_detectionSizeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_detectionSize)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_mutationRateStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_mutationRate)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_detectionDriversStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_detectionDrivers)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_initSizeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_initSize)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_extraTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_extraTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_finalTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_finalTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_maxMemoryStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_maxMemory)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_maxWallTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_maxWallTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_maxTriesStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulIndiv_maxNumTries)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulIndiv_allFitnessEffectsStatus <- renderImage(
    {
      allFitnessData <- allFitnessEffectsDecider(input)
      # When input$n is 3, filename is ./images/image3.jpeg
      filenameTick <- "./images/tick.jpeg"
      filenameCross <- "./images/x.jpeg"
      # Return a list containing the filename and alt text
      if(input$oncoSimulIndiv_allFitnessEffectsType == 1 && is.data.frame(allFitnessData)){
        filename <- filenameTick
      }else if (input$oncoSimulIndiv_allFitnessEffectsType == 2 && input$oncoSimulIndiv_allFitnessEffectsDataTableType == 1 && input$oncoSimulIndiv_NumGenes <= 10) {
        filename <- filenameTick
      } else if (input$oncoSimulIndiv_allFitnessEffectsType == 2 && input$oncoSimulIndiv_allFitnessEffectsDataTableType == 2 && is.data.frame(allFitnessData)){
        filename <- filenameTick
      }else{
        filename <- filenameCross
      }
      
      list(
        src = filename,
        height = "20px",
        width = "23px"
      )
    },
    deleteFile = FALSE
  )
  
  ############END STATUS REGION#########
  
  # renders the data table that has been entered using a file
  output$oncoSimulIndiv_allfitnessResultTable <- renderDataTable(fileToTable(input$oncoSimulIndiv_allFitnessEffectsFileInput))

  # AllFitnessEffect editable table 
  output$allFitnessEffectsDTEditable <- renderDT({
    
    if (input$oncoSimulIndiv_allFitnessEffectsDataTableType == 1) {
      validate(need(input$oncoSimulIndiv_NumGenes <= 10, "Please specify less number of genes (Maximum 10)"))
      initialMatrix <- OncoSimulR:::generate_matrix_genotypes(input$oncoSimulIndiv_NumGenes)
      Birth <- c(0.0)
      finalMatrix <- cbind(initialMatrix, Birth)
      allFitnessEffectsMatrix <<- finalMatrix
      datatable(finalMatrix, editable = list(target = "column", disable = list(columns = c(0:(input$oncoSimulIndiv_NumGenes - 1)))))
    } else {
      initialMatrix <- matrix(data = "", nrow = input$oncoSimulIndiv_NumRowsDataTableOpt2, ncol = 1)
      Birth <- c(0.0)
      finalMatrix <- cbind(initialMatrix, Birth)
      colnames(finalMatrix) <- c("Genotype", "Birth")
      allFitnessEffectsMatrix <<- finalMatrix
      datatable(finalMatrix, editable = "all")
    }
  })

  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$oncoSimulIndiv_down <- downloadHandler(
    validate(need(is.null(plotData), "No plot Data")),
    filename = function() {
      paste("plotName", "png", sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      png(file) # open the png device
      plot(plotData,
        type = oncosimulIndivPlotTypeInstance, main = input$oncoSimulIndiv_plotTitle,
        addtot = TRUE
      )
      dev.off() # turn the device off
    }
  )
}
