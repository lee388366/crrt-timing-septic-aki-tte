# Read-only anchor feasibility audit for MIMIC-IV v3.1.
# Session-local temporary tables disappear on disconnect; no permanent database objects are created.
# No rescue-CRRT exclusion: the source pipeline contains no operational definition.
suppressPackageStartupMessages({library(DBI); library(RPostgres); library(dplyr); library(readr)})
out <- "analysis/anchor_feasibility_output"
dir.create(out, showWarnings=FALSE, recursive=TRUE)
con <- dbConnect(RPostgres::Postgres(),dbname="mimiciv31",host="localhost",port=5432,user="postgres",password = Sys.getenv("PGPASSWORD"))
on.exit(dbDisconnect(con),add=TRUE)
run <- function(q) dbExecute(con,q)

run("CREATE TEMP TABLE a AS WITH icu AS (SELECT i.stay_id,i.intime FROM mimiciv_icu.icustays i) SELECT 'KDIGO3_minus6h'::text anchor,k.stay_id,greatest(i.intime,k.t-interval '6 hours') t0 FROM (SELECT stay_id,min(charttime) t FROM mimiciv_derived.kdigo_stages WHERE aki_stage=3 GROUP BY stay_id) k JOIN icu i USING(stay_id) UNION ALL SELECT 'KDIGO3',stay_id,t FROM (SELECT stay_id,min(charttime) t FROM mimiciv_derived.kdigo_stages WHERE aki_stage=3 GROUP BY stay_id) k UNION ALL SELECT 'KDIGO2',stay_id,t FROM (SELECT stay_id,min(charttime) t FROM mimiciv_derived.kdigo_stages WHERE aki_stage=2 GROUP BY stay_id) k")
run("CREATE INDEX ON a(stay_id)")
run("CREATE INDEX ON a(anchor,stay_id)")
run("CREATE TEMP TABLE sepsis AS SELECT a.anchor,a.stay_id,bool_or(s.sepsis3 AND s.sofa_time<=a.t0) sepsis_ok FROM a LEFT JOIN mimiciv_derived.sepsis3 s USING(stay_id) GROUP BY a.anchor,a.stay_id")
run("CREATE INDEX ON sepsis(anchor,stay_id)")
run("CREATE TEMP TABLE vaso AS SELECT a.anchor,a.stay_id,count(v.*)>0 vaso_ok,avg(coalesce(v.norepinephrine,0)+coalesce(v.epinephrine,0)+coalesce(v.dopamine,0)/100.0+coalesce(v.phenylephrine,0)*0.1+coalesce(v.vasopressin,0)*2.5) nee_at_t0 FROM a LEFT JOIN mimiciv_derived.vasoactive_agent v ON v.stay_id=a.stay_id AND v.starttime<=a.t0 AND v.endtime>=a.t0-interval '6 hours' GROUP BY a.anchor,a.stay_id")
run("CREATE INDEX ON vaso(anchor,stay_id)")
run("CREATE TEMP TABLE crrt AS SELECT stay_id,min(charttime) crrt_first_time FROM mimiciv_derived.rrt WHERE dialysis_active=1 AND dialysis_type='CRRT' GROUP BY stay_id")
run("CREATE INDEX ON crrt(stay_id)")
run("CREATE TEMP TABLE f0 AS SELECT a.anchor,a.t0,i.subject_id,i.hadm_id,i.stay_id,i.intime,p.anchor_age age,ad.dischtime,ad.deathtime,cb.scr_baseline creatinine_baseline,cr.crrt_first_time,v.nee_at_t0 FROM a JOIN mimiciv_icu.icustays i USING(stay_id) JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions ad USING(hadm_id) LEFT JOIN mimiciv_derived.creatinine_baseline cb USING(hadm_id) JOIN sepsis s USING(anchor,stay_id) JOIN vaso v USING(anchor,stay_id) LEFT JOIN crrt cr USING(stay_id) WHERE p.anchor_age>=18 AND s.sepsis_ok AND v.vaso_ok AND (cr.crrt_first_time IS NULL OR cr.crrt_first_time>=a.t0) AND extract(epoch from(i.outtime-i.intime))/3600>=24 AND cb.scr_baseline IS NOT NULL AND greatest(0,least(28,coalesce(extract(epoch from(ad.deathtime-a.t0))/86400,extract(epoch from(ad.dischtime-a.t0))/86400,28)))>0 AND NOT EXISTS(SELECT 1 FROM mimiciv_hosp.diagnoses_icd d WHERE d.hadm_id=i.hadm_id AND ((d.icd_version=10 AND d.icd_code IN('N186','N185')) OR (d.icd_version=9 AND d.icd_code IN('5856','V4511','V4512')))) AND NOT(a.anchor='KDIGO2' AND EXISTS(SELECT 1 FROM mimiciv_derived.kdigo_stages k WHERE k.stay_id=a.stay_id AND k.charttime=a.t0 AND k.aki_stage=3))")
run("CREATE TEMP TABLE f AS SELECT DISTINCT ON(anchor,subject_id) * FROM f0 ORDER BY anchor,subject_id,intime")
run("CREATE INDEX ON f(stay_id)")
run("CREATE TEMP TABLE final AS SELECT f.*,sf.sofa_24hours sofa_t0,ap.apsiii,coalesce(f.crrt_first_time between f.t0 and f.t0+interval '24 hours',false) early_crrt FROM f LEFT JOIN LATERAL (SELECT sofa_24hours FROM mimiciv_derived.sofa x WHERE x.stay_id=f.stay_id AND x.starttime<=f.t0+interval '1 hour' ORDER BY x.starttime DESC LIMIT 1) sf ON true LEFT JOIN mimiciv_derived.apsiii ap USING(stay_id)")
d <- dbGetQuery(con,"SELECT *,sofa_t0 between 8 and 14 equipoise FROM final")
write_csv(d,file.path(out,"anchor_patient_level.csv"))
sm <- d %>% group_by(anchor) %>% summarise(final_n=n(),early_n=sum(early_crrt,na.rm=TRUE),early_pct=round(100*mean(early_crrt,na.rm=TRUE),1),equipoise_n=sum(equipoise,na.rm=TRUE),equipoise_early_n=sum(equipoise & early_crrt,na.rm=TRUE),equipoise_early_pct=round(100*mean(early_crrt[equipoise],na.rm=TRUE),1),.groups="drop")
write_csv(sm,file.path(out,"anchor_summary.csv"))
smd <- function(x,g){g[is.na(g)]<-FALSE;a<-x[g];b<-x[!g];if(sum(is.finite(a))<2||sum(is.finite(b))<2)return(NA_real_);abs(mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/sqrt((var(a,na.rm=TRUE)+var(b,na.rm=TRUE))/2)}
vv <- c("age","sofa_t0","apsiii","creatinine_baseline","nee_at_t0")
sp <- split(d[d$equipoise,],d$anchor[d$equipoise])
zz <- bind_rows(lapply(sp,function(x)data.frame(anchor=unique(x$anchor),variable=vv,smd=sapply(vv,function(v)smd(x[[v]],x$early_crrt)),early_n=sum(x$early_crrt),deferred_n=sum(!x$early_crrt))))
write_csv(zz,file.path(out,"unweighted_smd_equipoise.csv"))
ps_rows <- lapply(sp,function(x){
  z <- x[,c("early_crrt",vv)]
  z <- z[complete.cases(z),]
  if(nrow(z)<20 || length(unique(z$early_crrt))<2) return(data.frame(anchor=unique(x$anchor),n_complete=nrow(z),early_ps_min=NA,early_ps_max=NA,deferred_ps_min=NA,deferred_ps_max=NA,common_low=NA,common_high=NA))
  fit <- glm(early_crrt ~ age+sofa_t0+apsiii+creatinine_baseline+nee_at_t0,data=z,family=binomial())
  p <- predict(fit,type="response")
  data.frame(anchor=unique(x$anchor),n_complete=nrow(z),early_ps_min=min(p[z$early_crrt]),early_ps_max=max(p[z$early_crrt]),deferred_ps_min=min(p[!z$early_crrt]),deferred_ps_max=max(p[!z$early_crrt]),common_low=max(min(p[z$early_crrt]),min(p[!z$early_crrt])),common_high=min(max(p[z$early_crrt]),max(p[!z$early_crrt])))
})
ps <- bind_rows(ps_rows)
write_csv(ps,file.path(out,"propensity_overlap_equipoise.csv"))
lines <- c("# Three-anchor feasibility audit","","## Assumptions / ambiguities","","- Read-only session-local temporary tables only; no modeling.","- Rescue/urgent CRRT exclusion not applied: absent from source pipeline.","- KDIGO2 is first observed stage 2; simultaneous stage 3 excluded.","- First ICU is selected after source eligibility, following 03_ccw_ipcw.R.","","## Results","","| anchor | final n | early n | early % | SOFA 8-14 n | SOFA 8-14 early n | SOFA 8-14 early % |","|---|---:|---:|---:|---:|---:|---:|")
lines <- c(lines,apply(sm,1,function(x)sprintf("| %s | %s | %s | %s | %s | %s | %s |",x[[1]],x[[2]],x[[3]],x[[4]],x[[5]],x[[6]],x[[7]])),"","Unweighted SMDs: see unweighted_smd_equipoise.csv; >=0.10 flags imbalance.","","Logistic propensity-score range overlap: see propensity_overlap_equipoise.csv. This is diagnostic only, not used for weighting or inference.")
writeLines(lines,file.path(out,"anchor_feasibility_report.md"));print(sm);print(zz)
