# voc_vaccines.R
# Runs a double exponential regression model with fixed growth rate advantage to forecast cases

library(tidyr)
library(dplyr)
library(lubridate)
library(readr)
library(segmented)

# Setup ----------
voc_vaccines <-
  function(
    transformerName,
    jurisdiction,
    cases,
    vaccination_rates,
    vaccine_efficacy,
    immunity_delay,
    startDate,
    endDate,
    prediction_end_date,
    vocAdvantage,
    voc_share_point,
    num_breakpoints,
    minimumIteration,
    maximumIteration,
    initialCumulative) {
    
  # The range of dates over which to project. The default is from the last day of
  # case data to 28 days later.
  # - You can specify any dates here in "yyyy-mm-dd" format (as a string)
  minimumTimestep <- startDate
  maximumTimestep <- prediction_end_date
  
  # Subset data by jurisdiction -------- 
  # This code would normally loop over jurisdictions, but I assume only one
  # jurisdiction to keep things simple
  
  # Setup jurisdiction-specific case data.frames
  jurisdictionCases <- cases %>% 
    subset(Jurisdiction == jurisdiction) 
  
  # Subset the case data to only days that should be fit by the regression
  jurisdictionRegressionCases <- jurisdictionCases %>%
    filter(Timestep >=startDate, Timestep<= endDate) %>%
    mutate(Day=difftime(Timestep,startDate,units = "day") %>% unclass) 
  
  # Set up vaccination and immunity function
  
  
  
  approx=tibble(Date=names(vaccination_rates),cvr=as.numeric(vaccination_rates)) %>%
    mutate(Day=difftime(Date,startDate,units = "day") %>% unclass)
  
  vcc.approx <- approxfun(approx %>% dplyr::select(Day,cvr))
  imm.approx <- approxfun(approx %>% dplyr::select(Day,cvr) %>%
                            mutate(Day=Day+immunity_delay,cvr=cvr*vaccine_efficacy))
  
  # VOC share
  voc_share <- function(date){
    ratio_to_share <- function(r)1/(1+1/r)
    share_to_ratio <- function(s)1/(1/s-1)
    r=share_to_ratio(as.numeric(voc_share_point)) %>% log
    t=as.numeric(difftime(date,as.Date(names(voc_share_point)),units = "days"))
    ratio_to_share(exp(r+vocAdvantage*t))
  }
  
  # Regression data
  regression_data <- jurisdictionRegressionCases %>%
    mutate(cir = imm.approx(Day)) %>%
    mutate(cir=coalesce(cir,0)) %>%
    mutate(S=1-cir) %>%
    mutate(voc=voc_share(Timestep)) %>%
    mutate(WT=Value*(1-voc)) %>%
    mutate(a=WT/S)
  
  model.lm <- lm(log(a) ~ Day,data=regression_data)
  model.s <- segmented::segmented(model.lm,npsi=num_breakpoints)
  
  break_points <- model.s$psi %>%
    as_tibble() %>%
    pull(Est.)
  
  c <- summary(model.s)$coefficients
  coeffs.s <- as_tibble(c) %>% 
    mutate(name=rownames(c)) %>%
    filter(grepl("Day",name),!grepl("psi",name)) %>%
    mutate(value=cumsum(Estimate))
  
  last_estimate <- coeffs.s %>% tail(1) %>% pull(Estimate)
  last_sd <- coeffs.s %>% tail(1) %>% pull(`Std. Error`)
  last_breakpoint <- break_points %>% tail(1)
  last_day <- regression_data %>% filter(Timestep==max(Timestep)) %>% pull(Day)
  #last_breakpoint_value <- predict(model.s,newdata=tibble(Day=last_breakpoint)) %>% as.numeric()
  
  project.s <- function(data,imm_lag=0,voc_lag=0,gr_var=0){
    data %>%
      mutate(a=(predict(model.s,newdata=.)+gr_var*(Day-last_day-7)) %>% exp) %>%
      mutate(cir = imm.approx(Day-imm_lag)) %>%
      mutate(cir=coalesce(cir,0)) %>%
      mutate(S=1-cir) %>%
      mutate(voc=voc_share(Timestep-voc_lag)) %>%
      mutate(WT=a*S) %>%
      mutate(Value=WT/(1-voc))
  }
  
  prediction_end_day=difftime(prediction_end_date,startDate,units="day") %>% unclass
  
  # projection <- tibble(Day=seq(0,difftime(prediction_end_date,startDate,units = "day") %>% unclass)) %>%
  #   mutate(Timestep=startDate+Day) %>%
  #   project.s()
  
  
  #td <- difftime(maximumTimestep,minimumTimestep,units="day") %>% unclass
  
  outputCases <- tibble(TransformerID=transformerName,
                                    Iteration = seq(1,maximumIteration),
                                    Variable = "Cases - Daily", 
                                    Jurisdiction = jurisdiction) %>%
    mutate(imm_lag=runif(maximumIteration,-7,7),
           voc_lag=runif(maximumIteration,-5,5),
           gr_var=rnorm(maximumIteration,0,last_sd),
           n=1) %>%
    full_join(tibble(n=1,Timestep=seq(minimumTimestep,maximumTimestep,by="day")) %>%
                mutate(Day=difftime(Timestep,startDate,units="day") %>% unclass),
              by="n") %>%
    mutate_at(c("imm_lag","voc_lag","gr_var"),function(d)ifelse(.$Day<last_day,0,d)) %>%
    mutate_at(c("imm_lag","voc_lag"),function(d)d*(.$Day-last_day)/(prediction_end_day-last_day)) %>%
    project.s(.,imm_lag,voc_lag,gr_var) %>%
    dplyr::select(-imm_lag,-voc_lag,-gr_var,-Day,-n) %>%
    bind_rows((.) %>%
                group_by(TransformerID,Jurisdiction,Iteration) %>%
                arrange(Timestep) %>%
                mutate(Variable="Cases - Cumulative",
                       Value=cumsum(Value),
                       Value=Value + initialCumulative))
  
  return(outputCases)
}
