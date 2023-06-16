library(shiny)
library(OncoSimulR)
library(readxl)
library(DT)
library(periscope)
library(vroom)
allFitnessEffectsMatrix <- matrix()

# Function which receives the arguments to define the fitness of our simulation
# and the different parameters that characterize the cellular interactions.Then
# returns the result of the simulation performed by the oncoSimulIndiv function
oncoSimulPopComplete <- function(allFitnessEffectsData, Nindiv, frequencyDependentFitness, mu, onlyCancer, initSize, finalTime, model, detectionSize, detectionDrivers,
                                   mutationPropGrowth, maxMemory, maxWallTime, maxNumTries,
                                   errorHitWallTime, errorHitMaxTries, AND_DrvProbExit) {
  print(allFitnessEffectsData)
  fitness <- allFitnessEffects(genotFitness = allFitnessEffectsData, frequencyDependentFitness = frequencyDependentFitness, frequencyType = "rel")
  result <- oncoSimulPop(fitness,Nindiv = Nindiv,
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
oncoSimulPopLogic <- function(input, output) {
  # observer dedicated to cell edition of AllFitnessEffects table
  observeEvent(input$allFitnessEffectsDTEditable_cell_edit, {
    allFitnessEffectsMatrix <<- editData(allFitnessEffectsMatrix, input$allFitnessEffectsDTEditable_cell_edit, "allFitnessEffectsDTEditable", rownames = FALSE)
    
  })
  
  # Observer that saves input data for the user to decide when to run the simulation
  observeEvent(list(
    input$oncoSimulPop_frequencyDepedentFitnessCheckbox, input$oncoSimulPop_mutationRate, input$oncoSimulPop_onlyCancer,
    input$oncoSimulPop_initSize, input$oncoSimulPop_finalTime, input$oncoSimulPop_selectModel, input$oncoSimulPop_detectionSize,
    input$oncoSimulPop_detectionDrivers, input$oncoSimulPop_mutationPropGrowth, input$oncoSimulPop_maxMemory, input$oncoSimulPop_maxWallTime,
    input$oncoSimulPop_maxNumTries, input$oncoSimulPop_errorHitWallTime, input$oncoSimulPop_errorHitMaxTries, input$oncoSimulPop_AND_DrvProbExit,
    input$oncoSimulPop_plotType, input$oncoSimulPop_plotShow, input$oncoSimulPop_nIndiv
  ), {
    nIndiv <<- input$oncoSimulPop_nIndiv
    frequencyDepedentFitnessCheckboxInstance <<- input$oncoSimulPop_frequencyDepedentFitnessCheckbox
    muInstance <<- input$oncoSimulPop_mutationRate
    onlyCancerInstance <<- input$oncoSimulPop_onlyCancer
    initSizeInstance <<- input$oncoSimulPop_initSize
    finalTimeInstance <<- input$oncoSimulPop_finalTime
    selectModelInstance <<- input$oncoSimulPop_selectModel
    detectionSizeInstance <<- input$oncoSimulPop_detectionSize
    detectionDriversInstance <<- input$oncoSimulPop_detectionDrivers
    mutationPropGrowthInstance <<- input$oncoSimulPop_mutationPropGrowth
    maxMemoryInstance <<- input$oncoSimulPop_maxMemory
    maxWallTimeInstance <<- input$oncoSimulPop_maxWallTime
    maxNumTriesInstance <<- input$oncoSimulPop_maxNumTries
    errorHitWallTimeInstance <<- input$oncoSimulPop_errorHitWallTime
    errorHitMaxTriesInstance <<- input$oncoSimulPop_errorHitMaxTries
    AND_DrvProbExitInstance <<- input$oncoSimulPop_AND_DrvProbExit
    oncosimulIndivPlotTypeInstance <<- input$oncoSimulPop_plotType
    oncosimulIndivPlotShow <<- input$oncoSimulPop_plotShow
  })
  
  # observer dedicated to execute the simulation and plot results
  observeEvent(input$oncoSimulPop_runExecutionButton, {
    output$oncoSimulPop_summary <- renderDataTable(summary(
      (oncoSimulPop_simulationData <<- oncoSimulPopComplete(
        allFitnessEffectsData = allFitnessEffectsDecider(input),
        Nindiv = nIndiv,
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
      ))
    ))
  })
  
  
  # reactive dedicated to decide which type of table is picked for AllFitnessEffects simulation and trasnform the data in case it is needed 
  allFitnessEffectsDecider <- function(input){
    if (input$oncoSimulPop_allFitnessEffectsType == 1) {
      print(fileToTable(input$oncoSimulPop_allFitnessEffectsFileInput))
      fileToTable(input$oncoSimulPop_allFitnessEffectsFileInput)
    } else {
      if (input$oncoSimulPop_allFitnessEffectsDataTableType == 1) {
        print(allFitnessEffectsMatrix)
        as.data.frame(allFitnessEffectsMatrix)
      } else {
        if (input$oncoSimulPop_frequencyDepedentFitnessCheckbox == FALSE) {
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
  output$oncoSimulPop_modelStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_selectModel)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_detectionSizeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_detectionSize)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_mutationRateStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_mutationRate)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_detectionDriversStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_detectionDrivers)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_initSizeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_initSize)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_extraTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_extraTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_finalTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_finalTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_maxMemoryStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_maxMemory)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_maxWallTimeStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_maxWallTime)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_maxTriesStatus <- renderImage(
    {
      detectBlankInput(input$oncoSimulPop_maxNumTries)
    },
    deleteFile = FALSE
  )
  
  output$oncoSimulPop_allFitnessEffectsStatus <- renderImage(
    {
      allFitnessData <- allFitnessEffectsDecider(input)
      # When input$n is 3, filename is ./images/image3.jpeg
      filenameTick <- "./images/tick.jpeg"
      filenameCross <- "./images/x.jpeg"
      # Return a list containing the filename and alt text
      if(input$oncoSimulPop_allFitnessEffectsType == 1 && is.data.frame(allFitnessData)){
        filename <- filenameTick
      }else if (input$oncoSimulPop_allFitnessEffectsType == 2 && input$oncoSimulPop_allFitnessEffectsDataTableType == 1 && input$oncoSimulPop_NumGenes <= 10) {
        filename <- filenameTick
      } else if (input$oncoSimulPop_allFitnessEffectsType == 2 && input$oncoSimulPop_allFitnessEffectsDataTableType == 2 && is.data.frame(allFitnessData)){
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
  output$oncoSimulPop_allfitnessResultTable <- renderDataTable(fileToTable(input$oncoSimulPop_allFitnessEffectsFileInput))
  
  # AllFitnessEffect editable table 
  output$allFitnessEffectsDTEditableOncoSimulPop <- renderDT({
    
    if (input$oncoSimulPop_allFitnessEffectsDataTableType == 1) {
      validate(need(input$oncoSimulPop_NumGenes <= 10, "Please specify less number of genes (Maximum 10)"))
      initialMatrix <- OncoSimulR:::generate_matrix_genotypes(input$oncoSimulPop_NumGenes)
      Birth <- c(0.0)
      finalMatrix <- cbind(initialMatrix, Birth)
      allFitnessEffectsMatrix <<- finalMatrix
      datatable(finalMatrix, editable = list(target = "column", disable = list(columns = c(0:(input$oncoSimulPop_NumGenes - 1)))))
    } else {
      initialMatrix <- matrix(data = "", nrow = input$oncoSimulPop_NumRowsDataTableOpt2, ncol = 1)
      Birth <- c(0.0)
      finalMatrix <- cbind(initialMatrix, Birth)
      colnames(finalMatrix) <- c("Genotype", "Birth")
      allFitnessEffectsMatrix <<- finalMatrix
      datatable(finalMatrix, editable = "all")
    }
  })
  
  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$oncoSimulPop_down <- downloadHandler(
    #validate(need(is.null(oncoSimulPop_simulationData), "No plot Data")),
    filename = function() {
      paste("plotName", "tsv", sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      vroom::vroom_write(summary(oncoSimulPop_simulationData), file)
      #write.csv(summary(oncoSimulPop_simulationData), file, row.names = FALSE)
      
    }
  )
}