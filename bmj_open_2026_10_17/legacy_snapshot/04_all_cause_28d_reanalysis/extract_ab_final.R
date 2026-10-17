suppressPackageStartupMessages({library(dplyr); library(DBI); library(RPostgres)})
con <- dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost", port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
base <- readRDS(file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data/baseline.rds")) %>% select(stay_id, subject_id, hadm_id, t0) %>% mutate(t0=as.POSIXct(t0,tz="UTC"))
g <- function(ids) dbGetQuery(con, sprintf("SELECT subject_id,hadm_id,charttime,valuenum FROM mimiciv_hosp.labevents WHERE itemid IN (%s) AND valuenum IS NOT NULL", paste(ids,collapse=",")))
near <- function(raw,nm,win=24) raw %>% inner_join(base, by=c("subject_id","hadm_id")) %>%
  mutate(dt=as.numeric(difftime(charttime,t0,units="hours"))) %>% filter(dt<=0, dt>=-win) %>%
  group_by(stay_id) %>% slice_max(charttime,n=1,with_ties=FALSE) %>% ungroup() %>% transmute(stay_id, !!nm := valuenum)
out <- base %>% select(stay_id) %>%
  left_join(near(g(50820),"ph_t0"), by="stay_id") %>%
  left_join(near(g(c(50971,50822)),"k_t0"), by="stay_id") %>%
  left_join(near(g(c(50882,50883,50884)),"hco3_t0"), by="stay_id")
for (v in c("ph_t0","k_t0","hco3_t0")) cat(sprintf("  %-8s 缺失率 %.1f%%\n", v, 100*mean(is.na(out[[v]]))))
saveRDS(out, file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/acidbase_t0_v2.rds")); dbDisconnect(con); cat("saved\n")
