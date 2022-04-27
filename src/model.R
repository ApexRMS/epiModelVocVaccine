# Script to wrap the voc_vaccines.R script and interface with SyncroSim

## Workspace setup ----
# Check for installed packages and install missing ones
packagesToLoad <- c("rsyncrosim", "purrr", "lubridate", "tidyr", "dplyr", "readr", "segmented")
packagesToInstall <- packagesToLoad[!(packagesToLoad %in% installed.packages()[,"Package"])]
if(length(packagesToInstall)) install.packages(packagesToInstall, repos = "https://cloud.r-project.org")
# Load packages
lapply(packagesToLoad, library, character.only = TRUE)

# Setup ----
transformerName <- "VOC + Vaccine Model: Run Model"

## Connect to SyncroSim ----
myScenario <- scenario()

packagePath <- ssimEnvironment()$PackageDirectory

inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F, optional = T) %>%
  mutate(TransformerID = as.character(TransformerID)) %>%
  replace_na(list(TransformerID = "Placeholder Transformer"))
runSettings <- datasheet(myScenario, "epiModelVocVaccine_RunSettings", lookupsAsFactors = F, optional = T)
jurisdictions <- datasheet(myScenario, "epiModelVocVaccine_RunJurisdictions", lookupsAsFactors = F, optional = T) %>% pull
vaccinationRates <- datasheet(myScenario, "epiModelVocVaccine_VaccinationRates", lookupsAsFactors = F, optional = T)

## Parse settings ----

sourceTransformer <- ifelse(is.na(runSettings$CaseSource), "Placeholder Transformer", runSettings$CaseSource)

if(length(jurisdictions) == 0)
  jurisdictions <- inputData %>%
    filter(TransformerID == sourceTransformer, Variable == "Cases - Daily") %>%
    pull(Jurisdiction) %>%
    unique

if(nrow(vaccinationRates) == 0)
  stop("Please provide some information on vaccination rates.")

# Use last day of data as first day of projection, if not provided
if(is.na(runSettings$MinimumTimestep))
  runSettings$MinimumTimestep <- inputData %>%
    filter(TransformerID == sourceTransformer) %>%
    pull(Timestep) %>%
    sort %>%
    tail(1)

# Use the past 180 days to fit the data if no fit start date is provided
if(is.na(runSettings$RegressionWindow))
  runSettings$RegressionWindow <- runSettings$MinimumTimestep %>%
    ymd %>%
    `-`(180) %>%
    format("%Y-%m-%d")

# Project 28 days forward if no projection end date is provided
if(is.na(runSettings$MaximumTimestep))
  runSettings$MaximumTimestep <- runSettings$MinimumTimestep %>%
    ymd %>%
    `+`(28) %>%
    format("%Y-%m-%d")

# Save run settings back to SyncroSim
saveDatasheet(myScenario, runSettings, "epiModelVocVaccine_RunSettings")

# Find the cumulative number of cases at the start of the regression window
initialCumulative <- 
  map_dbl(jurisdictions,
    ~{
      initialCumulativeJurisdiction <- inputData %>%
        filter(
          TransformerID == sourceTransformer,
          Jurisdiction == .x,
          Variable == "Cases - Cumulative",
          Timestep == runSettings$RegressionWindow) %>%
        pull(Value) %>%
        head(1)
      
      # Default to zero if cumulative cases data is not found
      if(length(initialCumulativeJurisdiction) == 0)
        initialCumulativeJurisdiction <- 0
      
      return(initialCumulativeJurisdiction)
    })

## Find and source the model ----
source(file.path(packagePath, "voc_vaccines.R"))
    
# Run model ----
outputData <- 
  map2_dfr(
    jurisdictions,
    initialCumulative,
    function(jurisdiction, initialCumulativeJurisdiction) {
      # Set up variables to loop over number of breakpoints in case of failure
      num_breakpoints <- runSettings$Breakpoints
      done_fit <- F
      output <- tibble()
      
      # Pull out the relevant vaccination rates for this jurisdiction
      jurisdictionVaccinationRates <- vaccinationRates %>%
          filter(is.na(Jurisdictions) | Jurisdictions == jurisdiction) %>%
          pull(Rate, name = Date)
      
      # Call the function
      while(!done_fit & num_breakpoints >= 0)
        tryCatch({
          # progress report
          print(paste0(jurisdiction, " - ", num_breakpoints))

          output <- voc_vaccines(
            transformerName = transformerName,
            jurisdiction = jurisdiction,
            cases = inputData %>% filter(Jurisdiction == jurisdiction, TransformerID == sourceTransformer, Variable == "Cases - Daily") %>% mutate(Value = Value + 1),
            vaccination_rates = jurisdictionVaccinationRates,
            vaccine_efficacy = runSettings$VaccineEfficacy,
            immunity_delay = runSettings$ImmunityDelay,
            startDate = runSettings$RegressionWindow %>% ymd,
            endDate = runSettings$MinimumTimestep %>% ymd,
            prediction_end_date = runSettings$MaximumTimestep %>% ymd,
            vocAdvantage = runSettings$VocAdvantage,
            voc_share_point = runSettings %>% pull(VocShare, name = VocShareDate),
            num_breakpoints = num_breakpoints,
            minimumIteration = runSettings$MinimumIteration,
            maximumIteration = runSettings$MaximumIteration,
            initialCumulative = initialCumulativeJurisdiction)
          done_fit <- TRUE},
          error = function(e) {
            if(grepl("Est.", e)) {
              warning("Could not estimate ", num_breakpoints, " breakpoints for ", jurisdiction, ". Decreasing the number of breakpoints.")
              num_breakpoints <<- num_breakpoints - 1
            } else if(grepl("L1 > L0", e)) {
              warning("No data found for ", jurisdiction, ". Excluding from projections.")
              done_fit <<- TRUE
              output <- NULL
            } else 
              stop("Unexpected error for ", jurisdiction, ".\n\n", e)
          })
        
        # If no appropriate number of breakpoints was found
        if(!done_fit)
          stop("Could not estimate breakpoints for one or more jurisdictions. Please change the number of breakpoints you would like to estimate.")
        
        return(output)
      }) %>%
  filter(
    TransformerID == transformerName,
    if(runSettings$HistoricalProjection == "No") Timestep >= runSettings$MinimumTimestep else TRUE) %>%
  dplyr::select(TransformerID, Iteration, Timestep, Variable, Jurisdiction, Value) %>%
  as.data.frame()

saveDatasheet(myScenario, outputData, "epi_DataSummary", append = T)
