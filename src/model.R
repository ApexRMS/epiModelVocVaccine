# Script to wrap the voc_vaccines.R script and interface with SyncroSim

library(rsyncrosim)
library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)

# Setup ----
transformerName <- "VOC + Vaccine Model: Run Model"

## Connect to SyncroSim ----
myScenario <- scenario()

packagePath <- ssimEnvironment()$PackageDirectory

inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F, optional = T) %>%
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
  try(map2_dfr(
    jurisdictions,
    initialCumulative,
    function(jurisdiction, initialCumulativeJurisdiction) {
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
        initialCumulative = initialCumulativeJurisdiction)
      }) %>%
  filter(TransformerID == transformerName) %>%
  select(TransformerID, Iteration, Timestep, Variable, Jurisdiction, Value) %>%
  as.data.frame())

if(class(outputData) == "try-error")
  case_when(
    grepl("Est.", outputData[1]) ~ stop("Could not estimate breakpoints for one or more jurisdictions. Please change the number of breakpoints you would like to estimate."),
    TRUE ~ stop("Unexpected error encountered while running the VOC and Vaccine model."))

saveDatasheet(myScenario, outputData, "epi_DataSummary", append = T)
