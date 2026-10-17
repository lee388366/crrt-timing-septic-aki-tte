# =============================================================================
# Q6 grace-period sensitivity — build hourly time-varying covariates to hour 47
# (needed only for the 48-h grace arm; 12/24-h reuse the locked tv_hourly 0-23).
#
# Mirrors 03a_timevarying_covariates.R exactly, but:
#   - cohort source = v2 clean-t0 primary cohort_base.rds (N=2,463)
#   - GRACE = 48  (hour_index 0:47)
#   - output -> 数据溯源_grace_period/tv_hourly_0_47.rds
# =============================================================================
suppressPackageStartupMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr);library(lubridate);library(zoo)})

V2_DIR  <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data")
OUT_DIR <- file.path(Sys.getenv("CRRT_MANUSCRIPT_ROOT"), "11_母板重建_数据回补/数据溯源_grace_period")
GRACE <- 48L

con <- dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost", port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
on.exit(dbDisconnect(con), add=TRUE)

cohort <- readRDS(file.path(V2_DIR,"cohort_base.rds")) %>%
  group_by(subject_id) %>% slice_min(icu_intime, n=1, with_ties=FALSE) %>% ungroup()
cat("Cohort:", nrow(cohort), "\n")

time_grid <- cohort %>% select(stay_id,hadm_id,subject_id,t0) %>%
  crossing(hour_index = 0:(GRACE-1)) %>%
  mutate(window_end=t0+hours(hour_index+1), window_start=t0+hours(hour_index))

vitals <- dbGetQuery(con,"SELECT stay_id,charttime,mbp AS map FROM mimiciv_derived.vitalsign WHERE mbp IS NOT NULL AND mbp BETWEEN 20 AND 200")
map_hourly <- time_grid %>% select(stay_id,hour_index,window_start,window_end) %>%
  left_join(vitals,by="stay_id",relationship="many-to-many") %>%
  filter(charttime>=window_start,charttime<window_end) %>% group_by(stay_id,hour_index) %>%
  summarise(map_h=mean(map,na.rm=TRUE),.groups="drop")

vaso <- dbGetQuery(con,"SELECT stay_id,starttime,endtime, COALESCE(norepinephrine,0)+COALESCE(epinephrine,0)+COALESCE(dopamine,0)/100.0+COALESCE(phenylephrine,0)*0.1+(COALESCE(vasopressin,0)/60.0)*2.5 AS nee_dose FROM mimiciv_derived.vasoactive_agent")
nee_hourly <- time_grid %>% select(stay_id,hour_index,window_start,window_end) %>%
  left_join(vaso,by="stay_id",relationship="many-to-many") %>%
  filter(starttime<window_end,endtime>window_start) %>% group_by(stay_id,hour_index) %>%
  summarise(nee_h=mean(nee_dose,na.rm=TRUE),.groups="drop")

lac <- dbGetQuery(con,"SELECT hadm_id,charttime,lactate FROM mimiciv_derived.bg WHERE lactate IS NOT NULL AND lactate<30")
lactate_hourly <- time_grid %>% select(stay_id,hadm_id,hour_index,window_end) %>%
  left_join(lac,by="hadm_id",relationship="many-to-many") %>%
  filter(charttime>=window_end-hours(6),charttime<window_end) %>% group_by(stay_id,hour_index) %>%
  summarise(lactate_h=mean(lactate,na.rm=TRUE),.groups="drop")

sofa <- dbGetQuery(con,"SELECT stay_id,starttime,sofa_24hours AS sofa_total FROM mimiciv_derived.sofa WHERE sofa_24hours IS NOT NULL")
sofa_hourly <- time_grid %>% select(stay_id,hour_index,window_end) %>%
  left_join(sofa,by="stay_id",relationship="many-to-many") %>% filter(starttime<=window_end) %>%
  group_by(stay_id,hour_index) %>% slice_max(starttime,n=1,with_ties=FALSE) %>%
  summarise(sofa_h=first(sofa_total),.groups="drop")

kd <- dbGetQuery(con,"SELECT stay_id,charttime,aki_stage FROM mimiciv_derived.kdigo_stages WHERE aki_stage IS NOT NULL")
kdigo_hourly <- time_grid %>% select(stay_id,hour_index,window_end) %>%
  left_join(kd,by="stay_id",relationship="many-to-many") %>% filter(charttime<=window_end) %>%
  group_by(stay_id,hour_index) %>% slice_max(charttime,n=1,with_ties=FALSE) %>%
  summarise(kdigo_h=first(aki_stage),.groups="drop")

cr <- dbGetQuery(con,"SELECT hadm_id,charttime,creatinine FROM mimiciv_derived.chemistry WHERE creatinine IS NOT NULL AND creatinine<25")
creat_t0 <- cr %>% inner_join(cohort %>% select(stay_id,hadm_id,t0),by="hadm_id") %>%
  filter(charttime<=t0+hours(2),charttime>=t0-hours(6)) %>% group_by(stay_id) %>%
  slice_max(charttime,n=1,with_ties=FALSE) %>% ungroup() %>% select(stay_id,creatinine_t0=creatinine)
creat_hourly <- time_grid %>% select(stay_id,hadm_id,hour_index,window_end) %>%
  left_join(cr,by="hadm_id",relationship="many-to-many") %>%
  filter(charttime>=window_end-hours(6),charttime<window_end) %>% group_by(stay_id,hour_index) %>%
  summarise(creatinine_h=mean(creatinine,na.rm=TRUE),.groups="drop") %>%
  left_join(creat_t0,by="stay_id") %>% mutate(creatinine_delta=creatinine_h-coalesce(creatinine_t0,creatinine_h))

inp <- dbGetQuery(con,"SELECT stay_id,starttime,amount FROM mimiciv_icu.inputevents WHERE amount>0 AND amount<50000 AND amountuom IN ('mL','ml')")
outp <- dbGetQuery(con,"SELECT stay_id,charttime,value AS amount FROM mimiciv_icu.outputevents WHERE value>0 AND value<20000")
input_cumul <- time_grid %>% select(stay_id,hour_index,t0,window_end) %>%
  left_join(inp,by="stay_id",relationship="many-to-many") %>% filter(starttime>=t0,starttime<window_end) %>%
  group_by(stay_id,hour_index) %>% summarise(total_input=sum(amount,na.rm=TRUE),.groups="drop")
output_cumul <- time_grid %>% select(stay_id,hour_index,t0,window_end) %>%
  left_join(outp,by="stay_id",relationship="many-to-many") %>% filter(charttime>=t0,charttime<window_end) %>%
  group_by(stay_id,hour_index) %>% summarise(total_output=sum(amount,na.rm=TRUE),.groups="drop")
fluid_hourly <- full_join(input_cumul,output_cumul,by=c("stay_id","hour_index")) %>%
  mutate(fluid_balance_h=coalesce(total_input,0)-coalesce(total_output,0))

tv_hourly <- time_grid %>% select(stay_id,hour_index) %>%
  left_join(map_hourly,by=c("stay_id","hour_index")) %>%
  left_join(nee_hourly,by=c("stay_id","hour_index")) %>%
  left_join(lactate_hourly,by=c("stay_id","hour_index")) %>%
  left_join(sofa_hourly,by=c("stay_id","hour_index")) %>%
  left_join(kdigo_hourly,by=c("stay_id","hour_index")) %>%
  left_join(creat_hourly %>% select(stay_id,hour_index,creatinine_h,creatinine_delta),by=c("stay_id","hour_index")) %>%
  left_join(fluid_hourly %>% select(stay_id,hour_index,fluid_balance_h),by=c("stay_id","hour_index")) %>%
  group_by(stay_id) %>%
  mutate(across(c(map_h,nee_h,lactate_h,sofa_h,kdigo_h,creatinine_h,creatinine_delta,fluid_balance_h),
                ~zoo::na.locf(.,na.rm=FALSE))) %>%
  mutate(nee_h=coalesce(nee_h,0), fluid_balance_h=coalesce(fluid_balance_h,0)) %>% ungroup()

saveRDS(tv_hourly, file.path(OUT_DIR,"tv_hourly_0_47.rds"))
cat("Saved tv_hourly_0_47.rds | rows:",nrow(tv_hourly),"| hours:",paste(range(tv_hourly$hour_index),collapse="-"),"\n")
