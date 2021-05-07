# Script to wrap the voc_vaccines.R script and interface with SyncroSim

library(rsyncrosim)
library(dplyr)
library(purrr)
library(lubridate)

# Setup ----
transformerName <- "VOC + Vaccine Model: Run Model"

## Connect to SyncroSim ----
myScenario <- scenario()

packagePath <- ssimEnvironment()$PackageDirectory

pipeline <- datasheet(myScenario, "core_Pipeline", lookupsAsFactors = F)
inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F)
runSettings <- datasheet(myScenario, "epiModelVocVaccine_RunSettings", lookupsAsFactors = F, optional = T)
jurisdictions <- datasheet(myScenario, "epiModelVocVaccine_RunJurisdictions", lookupsAsFactors = F, optional = T) %>% pull
vaccinationRates <- datasheet(myScenario, "epiModelVocVaccine_VaccinationRates", lookupsAsFactors = F, optional = T)

## Decide which source transformer to use ----

# Find the position of the current transformer in the pipeline
currentRunOrder <- pipeline %>%
  filter(StageNameID == transformerName) %>%
  pull(RunOrder)

if(currentRunOrder > 1){
  sourceTransformer <- pipeline %>%
    filter(RunOrder == currentRunOrder - 1) %>%
    pull(StageNameID)
} else{
  sourceTransformer <- inputData %>%
    filter(Variable == "Cases - Daily") %>%
    pull(TransformerID) %>%
    tail(1)
}

## Parse settings ----
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

# Find the cumulative number of cases at the start of the regression window
initialCumulative <- inputData %>%
  filter(
    TransformerID == sourceTransformer,
    Variable == "Cases - Cumulative",
    Timestep == runSettings$RegressionWindow) %>%
  pull(Value)

# Default to zero if cumulative cases data is not found
if(length(initialCumulative) == 0)
  initialCumulative <- 0

## Find and source the model ----
source(file.path(packagePath, "voc_vaccines.R"))
    
# Run model ----
outputData <- 
  map_dfr(
    jurisdictions,
    function(jurisdiction) {
      # Pull out the relevant vaccination rates for this jurisdiction
      jurisdictionVaccinationRates <- vaccinationRates %>%
          filter(is.na(Jurisdictions) | Jurisdictions == jurisdiction) %>%
          pull(Rate, name = Date)
      
      # Call the function
      voc_vaccines(
        transformerName = transformerName,
        jurisdiction = jurisdiction,
        cases = inputData %>% filter(Jurisdiction == jurisdiction, Variable == "Cases - Daily"),
        vaccination_rates = jurisdictionVaccinationRates,
        vaccine_efficacy = runSettings$VaccineEfficacy,
        immunity_delay = runSettings$ImmunityDelay,
        startDate = runSettings$RegressionWindow %>% ymd,
        endDate = runSettings$MinimumTimestep %>% ymd,
        prediction_end_date = runSettings$MaximumTimestep %>% ymd,
        vocAdvantage = runSettings$VocAdvantage,
        voc_share_point = runSettings %>% pull(VocShare, name = VocShareDate),
        num_breakpoints = runSettings$Breakpoints,
        minimumIteration = runSettings$MinimumIteration,
        maximumIteration = runSettings$MaximumIteration,
        initialCumulative = initialCumulative)
      }) %>%
  filter(TransformerID == transformerName) %>%
  select(TransformerID, Iteration, Timestep, Variable, Jurisdiction, Value) %>%
  as.data.frame()

saveDatasheet(myScenario, outputData, "epi_DataSummary", append = T)
