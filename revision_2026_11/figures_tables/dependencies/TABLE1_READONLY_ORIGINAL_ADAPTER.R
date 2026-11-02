# CH70: SELECT-only adapter. It reuses the original Table1 SQL JOIN body verbatim.
# No new cohort, missing historical overall generator is NOT claimed reconstructed.
# Default --replay-old; --new PATH requires same key contract and prior replay PASS.
suppressPackageStartupMessages({library(DBI);library(RPostgres);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE)
arg<-function(k,default=NULL){i<-match(k,a);if(is.na(i))default else {stopifnot(i<length(a));a[i+1L]}}
mode<-if('--new'%in%a)'new'else'replay_old'
root<-normalizePath(arg('--root','.'),mustWork=TRUE)
out<-arg('--out',file.path(root,'private_data',paste0('rich_',mode)))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_sql<-arg('--source-sql',file.path(Sys.getenv('CRRT_REPO_DIR'),'upstream/table1/05_baseline_table1_v2.sql'))
source_helpers<-arg('--source-helpers',file.path(Sys.getenv('CRRT_REPO_DIR'),'upstream/table1/table1_rich_sofa8_14_v2.R'))
reference<-arg('--reference',Sys.getenv('CRRT_TABLE1_REFERENCE_CSV'))
keys_path<-if(mode=='new')arg('--new')else arg('--keys',file.path(Sys.getenv('CRRT_INPUT_DIR'),'ch24_q4_corrected_inputs/baseline.rds'))
sha<-function(p)unname(as.character(tools::md5sum(p))) # SHA256 recorded by root command; MD5 local drift aid only.
readkeys<-function(p){if(grepl('\\.rds$',p,ignore.case=TRUE))readRDS(p)else read.csv(p,check.names=FALSE)}
k<-readkeys(keys_path)
required<-c('stay_id','subject_id','hadm_id','icu_intime','icu_outtime','age','gender','icu_los_hours','kdigo3_time','t0')
stopifnot(all(required%in%names(k)),!anyDuplicated(k$stay_id),!anyNA(k$stay_id),!anyNA(k$t0))
if(mode=='replay_old')stopifnot(nrow(k)==2463L)
# JSON bound as a value. IDs and timestamps never interpolated into SQL or printed.
kbind<-k[,required]
for(v in c('icu_intime','icu_outtime','kdigo3_time','t0'))kbind[[v]]<-format(as.POSIXct(kbind[[v]],tz='UTC'),'%Y-%m-%d %H:%M:%OS6',tz='UTC')
kjson<-toJSON(kbind,dataframe='rows',auto_unbox=TRUE,na='null')
keys_cte<-"keys AS (SELECT * FROM jsonb_to_recordset($1::jsonb) AS k(stay_id bigint,subject_id bigint,hadm_id bigint,icu_intime timestamp,icu_outtime timestamp,age double precision,gender text,icu_los_hours double precision,kdigo3_time timestamp,t0 timestamp))"
# Original SQL has DDL before/after extraction. Extract only FROM cohort c SELECT.
s<-paste(readLines(source_sql,warn=FALSE),collapse='\n')
start<-regexpr('\nSELECT\n  -- ── identifiers',s,fixed=TRUE);stopifnot(start>0)
tail<-substring(s,start+1L);end<-regexpr('\n;\n',tail,fixed=TRUE);stopifnot(end>0)
body<-substring(tail,1,end-1L)
# Preserve original inner admissions join & all per-variable JOIN definitions.
# Supplied keys are already the authorized cohort: do not reselect first ICU here.
cohort_cte<-"cohort AS (SELECT cb.stay_id,cb.subject_id,cb.hadm_id,cb.icu_intime AS intime,cb.icu_outtime AS outtime,cb.age,cb.gender,adm.race,adm.insurance,cb.icu_los_hours/24.0 AS los_icu,cb.kdigo3_time,cb.t0 AS t0_time,cb.t0-INTERVAL '6 hour' AS baseline_start,adm.admittime FROM keys cb JOIN mimiciv_hosp.admissions adm ON adm.hadm_id=cb.hadm_id)"
query<-paste('WITH',keys_cte,',',cohort_cte,body)
# Separate *documented* rich-primary augmentation. Do not choose best-fitting source.
augment<-paste('WITH',keys_cte,"SELECT b.stay_id,fdw.weight AS first_day_weight, fds.sofa AS sofa_admission,fdg.gcs_min, n.nee_at_t0_corrected FROM keys b LEFT JOIN mimiciv_derived.first_day_weight fdw ON fdw.stay_id=b.stay_id LEFT JOIN mimiciv_derived.first_day_sofa fds ON fds.stay_id=b.stay_id LEFT JOIN mimiciv_derived.first_day_gcs fdg ON fdg.stay_id=b.stay_id LEFT JOIN LATERAL (SELECT AVG(coalesce(v.norepinephrine,0)+coalesce(v.epinephrine,0)+coalesce(v.dopamine,0)/100.0+coalesce(v.phenylephrine,0)*0.1+(coalesce(v.vasopressin,0)/60.0)*2.5) AS nee_at_t0_corrected FROM mimiciv_derived.vasoactive_agent v WHERE v.stay_id=b.stay_id AND v.starttime<=b.t0 AND v.endtime>=b.t0-interval '6 hours') n ON true")
# Label-driven arterial correction candidate: exact original JOIN7 with only specimen predicate added.
arterial_query<-paste('WITH',keys_cte,", cohort AS (SELECT *,t0 AS t0_time,t0-interval '6 hours' AS baseline_start FROM keys) SELECT c.stay_id,c.t0_time,b.ph,b.pao2,b.fio2,CASE WHEN b.pao2 IS NOT NULL AND b.fio2 IS NOT NULL AND b.fio2>0 THEN b.pao2/NULLIF(b.fio2/100.0,0) ELSE NULL END AS pf_ratio FROM cohort c LEFT JOIN LATERAL (SELECT b.ph,b.lactate,b.po2 AS pao2,b.fio2 FROM mimiciv_derived.bg b WHERE b.hadm_id=c.hadm_id AND b.charttime>c.baseline_start AND b.charttime<=c.t0_time AND b.specimen='ART.' ORDER BY b.charttime DESC LIMIT 1) b ON true")
insurance_query<-paste('WITH',keys_cte,"SELECT ad.insurance AS raw_category,COUNT(*)::int AS n,COUNT(ad.hadm_id)::int AS matched_admissions FROM keys b LEFT JOIN mimiciv_hosp.admissions ad ON ad.hadm_id=b.hadm_id GROUP BY ad.insurance ORDER BY ad.insurance NULLS FIRST")
schema_query<-"SELECT table_schema,table_name,column_name,data_type FROM information_schema.columns WHERE table_schema IN ('mimiciv_hosp','mimiciv_icu','mimiciv_derived') AND table_name IN ('admissions','icustays','weight_durations','charlson','vitalsign','chemistry','bg','enzyme','coagulation','complete_blood_count','gcs','apsiii','ventilation','vasoactive_agent','norepinephrine_equivalent_dose','sofa','diagnoses_icd','inputevents','prescriptions','urine_output','kdigo_stages','creatinine_baseline','first_day_weight','first_day_sofa','first_day_gcs') ORDER BY table_schema,table_name,ordinal_position"
stripcomments<-function(q)paste(sub('--.*$','',strsplit(q,'\n',fixed=TRUE)[[1]]),collapse='\n')
for(q in list(query,augment,schema_query,insurance_query,arterial_query)){qq<-stripcomments(q);stopifnot(!grepl(';',qq,fixed=TRUE),grepl('^\\s*(WITH|SELECT)\\b',qq),!grepl('\\b(INSERT|UPDATE|DELETE|CREATE|DROP|ALTER|TRUNCATE|COPY|CALL|DO|GRANT|REVOKE|ANALYZE|VACUUM)\\b',qq,ignore.case=TRUE))}
writeLines(arterial_query,file.path(out,'SQL_arterial_only_candidate.sql'));writeLines(insurance_query,file.path(out,'SQL_insurance_only.sql'));writeLines(query,file.path(out,'SQL_original_body_select.sql'));writeLines(augment,file.path(out,'SQL_documented_rich_primary_augmentation.sql'));writeLines(schema_query,file.path(out,'SQL_required_schema.sql'))
plan<-list(mode=mode,keys_path=normalizePath(keys_path),keys_n=nrow(k),source_sql=source_sql,source_helpers=source_helpers,reference=reference,source_MD5=list(sql=sha(source_sql),helpers=sha(source_helpers),keys=sha(keys_path),reference=sha(reference)),cohort_only_change='Original legacy cohort CTE replaced by bound authorized keys; all original main SELECT expressions/JOINS unchanged',augmentation_status='Separate explicitly documented source; never auto-selected to maximize match',not_recoverable_from_original_sql=c('Albumin','Total bilirubin','Exact original rich overall generator and historical augmentation receipt'),scientific_limit='APSIII first24h/CCI index-admission codes are descriptive, not pretreatment inputs; row-wise panel last observation is not analyte-wise last nonnull; MAP actual original generator NOT_REPRODUCED (raw6h missing801 vs displayed2); no blind64/65 replacement',execution_status='SQL_EMITTED_NOT_EXECUTED')
write_json(plan,file.path(out,'PLAN.json'),auto_unbox=TRUE,pretty=TRUE)
if('--emit-sql'%in%a){cat('SELECT-only SQL and plan written; no DB connection\n');quit(status=0)}
if(mode=='new'&&!any(c('--insurance-only','--arterial-only','--schema','--explain')%in%a)){
 gate<-arg('--old-gate');if(is.null(gate))stop('--new requires --old-gate from replay; never skip unresolved old reproduction')
 oldgate<-fromJSON(gate);if(!isTRUE(oldgate$data_select_ran_readonly)||!isTRUE(oldgate$keys_and_t0_exact))stop('Old SELECT replay/key-time validation must run first.');if(!isTRUE(oldgate$all_current_manuscript_rows_verified))warning('Old rich full-table unresolved: NEW may be extracted for diagnosis only, never promoted as verified complete Table1.')
}
# Never execute source file, never source() its connect/write statements.
con<-dbConnect(RPostgres::Postgres(),host=Sys.getenv('PGHOST','localhost'),port=as.integer(Sys.getenv('PGPORT','5432')),dbname=Sys.getenv('PGDATABASE','mimiciv31'),user=Sys.getenv('PGUSER','postgres'),password=Sys.getenv('PGPASSWORD'),options='-c default_transaction_read_only=on')
on.exit(dbDisconnect(con),add=TRUE)
dbExecute(con,'BEGIN TRANSACTION READ ONLY');on.exit(try(dbExecute(con,'ROLLBACK'),silent=TRUE),add=TRUE)
stopifnot(dbGetQuery(con,'SHOW transaction_read_only')[[1]][1]=='on')
dbExecute(con,"SET LOCAL statement_timeout = '10min'")
if('--schema'%in%a){sch<-dbGetQuery(con,schema_query);saveRDS(sch,file.path(out,'SCHEMA.rds'));write_json(sch,file.path(out,'SCHEMA.json'),dataframe='rows',pretty=TRUE);cat('Read-only schema saved; no patient rows queried\n');quit(status=0)}
if('--insurance-only'%in%a){
 tab<-dbGetQuery(con,insurance_query,params=list(as.character(kjson)));saveRDS(tab,file.path(out,'INSURANCE_RAW_CATEGORY_AGGREGATE.rds'));write.csv(tab,file.path(out,'INSURANCE_RAW_CATEGORY_AGGREGATE.csv'),row.names=FALSE)
 stopifnot(sum(tab$n)==nrow(k),all(tab$matched_admissions<=tab$n))
 catval<-as.character(tab$raw_category)
 null_n<-sum(tab$n[is.na(catval)]);empty_n<-sum(tab$n[!is.na(catval)&catval=='']);white_n<-sum(tab$n[!is.na(catval)&catval!=''&trimws(catval)==''])
 nonblank_n<-sum(tab$n[!is.na(catval)&trimws(catval)!=''])
 expected<-c(Medicare=1562L,Medicaid=301L,Private=512L,Other=55L,'No charge'=0L)
 got<-vapply(names(expected),function(n)sum(tab$n[!is.na(catval)&catval==n]),numeric(1))
 unmatched<-sum(tab$n-tab$matched_admissions)
 summary<-list(mode=mode,keys_n=nrow(k),joined_rows=sum(tab$n),matched_admissions=sum(tab$matched_admissions),unmatched_admissions=unmatched,sql_null=null_n,empty_string=empty_n,whitespace_only=white_n,recorded_nonblank=nonblank_n,absent_total=null_n+empty_n+white_n,categories=tab,old2463_reference_category_counts=as.list(expected),observed_category_counts=as.list(got),reference_categories_exact=isTRUE(all(got==expected)),old2463_raw_absence33_confirmed=mode=='replay_old'&&nrow(k)==2463L&&unmatched==0L&&(null_n+empty_n+white_n)==33L,transaction_read_only='on',heavy_rich_query_executed=FALSE,limit='Absent raw insurance is separately counted as NULL/empty/whitespace; omitted nonblank categories remain visible. No mapping or arithmetic-difference inference.')
 write_json(summary,file.path(out,'INSURANCE_RAW_AGGREGATE.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
 cat('Insurance-only READ ONLY aggregate saved; no heavy rich query and no patient rows printed.\n');quit(status=0)
}
if('--arterial-only'%in%a){
 if('--explain'%in%a){exp<-dbGetQuery(con,paste('EXPLAIN (FORMAT JSON)',arterial_query),params=list(as.character(kjson)));txt<-as.character(exp[[1]][1]);stopifnot(jsonlite::validate(txt));writeLines(txt,file.path(out,'EXPLAIN_arterial_only.json'));cat('Arterial candidate EXPLAIN saved; no extraction.\n');quit(status=0)}
 ar<-dbGetQuery(con,arterial_query,params=list(as.character(kjson)));stopifnot(nrow(ar)==nrow(k),!anyDuplicated(ar$stay_id),setequal(as.character(ar$stay_id),as.character(k$stay_id)))
 saveRDS(ar,file.path(out,'PRIVATE_arterial_only_candidate.rds'))
 agg<-lapply(c('ph','pao2','fio2','pf_ratio'),function(v){x<-ar[[v]];list(field=v,n=nrow(ar),missing=sum(is.na(x)),quartiles=as.list(quantile(x,c(.25,.5,.75),na.rm=TRUE,type=7)))})
 write_json(list(status='LABEL_DRIVEN_CORRECTION_CANDIDATE_NOT_OLD_REPRODUCTION',mode=mode,keys_n=nrow(k),transaction_read_only='on',change="Original JOIN7 plus b.specimen='ART.' only; original last row, prior6h and same-row PF retained; no new analyte filtering or tie breaker",original_scripts_have_ART_filter=FALSE,aggregates=agg,formal_files_modified=FALSE),file.path(out,'ARTERIAL_ONLY_CANDIDATE.json'),auto_unbox=TRUE,pretty=TRUE)
 cat('Arterial candidate aggregates saved; original reproduction remains a separate gate.\n');quit(status=0)
}
if('--explain'%in%a){
 exp<-dbGetQuery(con,paste('EXPLAIN (FORMAT JSON)',query),params=list(as.character(kjson)));saveRDS(exp,file.path(out,'EXPLAIN_original_select.rds'));plan_text<-as.character(exp[[1]][1]);stopifnot(jsonlite::validate(plan_text));writeLines(plan_text,file.path(out,'EXPLAIN_original_select.json'))
 cat('EXPLAIN only saved; ANALYZE not used and rich SELECT not executed.\n');quit(status=0)
}
if(!('--run-rich'%in%a))stop('Full original rich SELECT requires explicit --run-rich after --explain review; use --insurance-only for immediate independent insurance verification.')
reuse_original<-arg('--reuse-original-rds')
if(is.null(reuse_original)){d<-dbGetQuery(con,query,params=list(kjson))}else{stored_sql<-paste(readLines(file.path(dirname(reuse_original),'SQL_original_body_select.sql'),warn=FALSE),collapse='\n');stopifnot(identical(trimws(stored_sql),trimws(query)));d<-readRDS(reuse_original)}
saveRDS(d,file.path(out,'PRIVATE_original_select.rds'))
stopifnot(nrow(d)==nrow(k),!anyDuplicated(d$stay_id),setequal(as.character(d$stay_id),as.character(k$stay_id)))
j<-match(as.character(k$stay_id),as.character(d$stay_id));stopifnot(all(as.numeric(as.POSIXct(d$t0_time[j],tz='UTC'))==as.numeric(as.POSIXct(k$t0,tz='UTC'))))
u<-dbGetQuery(con,augment,params=list(kjson));saveRDS(u,file.path(out,'PRIVATE_documented_augmentation.rds'));stopifnot(nrow(u)==nrow(k),!anyDuplicated(u$stay_id))
# Reuse ONLY existing formatting helpers by extracting their definitions, never source().
e<-parse(source_helpers);want<-c('fmt_n','fmt_cont','fmt_cat','race_group');ee<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&is.symbol(x[[2]])&&as.character(x[[2]])%in%want,e);stopifnot(length(ee)==4L);for(ex in ee)eval(ex)
ref<-read.csv(reference,check.names=FALSE,stringsAsFactors=FALSE)
# Explicit label→original SQL field map; no alteration of field definitions.
cmap<-c('Age, years'='age','Body weight, kg'='weight_kg','Charlson Comorbidity Index'='cci','SOFA score at t0'='sofa_baseline','GCS score (minimum)'='gcs_baseline','APACHE III score'='apsiii_score','Heart rate, beats/min'='hr','Systolic blood pressure, mmHg'='sbp','Diastolic blood pressure, mmHg'='dbp','Mean arterial pressure, mmHg'='map','Respiratory rate, breaths/min'='rr','SpO2, %'='spo2','Body temperature, °C'='temp_c','White blood cell count, 10^9/L'='wbc','Hemoglobin, g/dL'='hb','Platelet count, 10^9/L'='plt','Sodium, mmol/L'='na','Chloride, mmol/L'='cl','Potassium, mmol/L'='k','Bicarbonate, mmol/L'='hco3','Serum creatinine (at t0), mg/dL'='scr','Blood urea nitrogen, mg/dL'='bun','AST, U/L'='ast','ALT, U/L'='alt','Arterial pH'='ph','Lactate, mmol/L'='lactate','PaO2, mmHg'='pao2','FiO2, %'='fio2','P/F ratio, mmHg'='pf_ratio','INR'='inr','Norepinephrine-equivalent dose, µg/kg/min'='ne_eq_baseline','Crystalloid volume, mL/kg'='crystalloid_ml_kg','Urine output (6 h before t0), mL/kg/h'='uo_mlkgph_t0_last6h','Baseline serum creatinine, mg/dL'='scr_baseline','Creatinine at t0, mg/dL'='scr')
bmap<-c('Diabetes mellitus, n (%)'='dm','Hypertension, n (%)'='htn','Congestive heart failure, n (%)'='chf','Chronic kidney disease (non-ESRD), n (%)'='ckd','Chronic liver disease, n (%)'='liver_disease','Active malignancy, n (%)'='malignancy','COPD, n (%)'='copd','Invasive mechanical ventilation, n (%)'='mv_baseline','Vasopressor use, n (%)'='vaso_baseline','Artificial colloid exposure, n (%)'='colloid_flag','Nephrotoxic drug exposure, n (%)'='nephrotoxic_flag')
miss<-function(x)sprintf('%s (%.1f)',fmt_n(sum(is.na(x))),100*mean(is.na(x)))
rows<-lapply(seq_len(nrow(ref)),function(i){lab<-ref$Variable[i];field<-unname(cmap[lab]);v<-NULL;val<-NA_character_;status<-'UNRESOLVED_ORIGINAL_OVERALL_GENERATOR_OR_MAPPING';missing<-NA_character_
 if(!is.na(field)){v<-d[[field]];dig<-if(grepl('^[0-9]+\\.',ref$Overall[i]))nchar(sub('^[0-9]+\\.([0-9]+).*','\\1',ref$Overall[i]))else 0L;val<-fmt_cont(v,dig);missing<-miss(v);status<-if(val==ref$Overall[i]&&missing==ref[['Missing n (%)']][i])'EXACT_ORIGINAL_SQL_REPLAY'else'MISMATCH_DO_NOT_PROMOTE'}
 bf<-unname(bmap[lab]);if(!is.na(bf)){field<-bf;v<-d[[bf]];val<-fmt_cat(v==1);missing<-miss(v);status<-if(val==ref$Overall[i]&&missing==ref[['Missing n (%)']][i])'EXACT_ORIGINAL_SQL_REPLAY'else'MISMATCH_DO_NOT_PROMOTE'}
 data.frame(variable=lab,source_field=ifelse(is.na(field),'',field),reference_value=ref$Overall[i],replay_value=val,reference_missing=ref[['Missing n (%)']][i],replay_missing=missing,status=status)
})
cmp<-do.call(rbind,rows);write.csv(cmp,file.path(out,'FIELD_BY_FIELD_ORIGINAL_SQL_REPLAY.csv'),row.names=FALSE)
# Raw insurance categories retain exact spelling, including blank/whitespace. No mapping exclusion inferred.
ins<-as.character(d$insurance);tab<-as.data.frame(table(ins,useNA='ifany'),stringsAsFactors=FALSE);names(tab)<-c('raw_category','n');write.csv(tab,file.path(out,'INSURANCE_RAW_CATEGORY_AGGREGATE.csv'),row.names=FALSE)
insurance<-list(n=nrow(d),sql_null=sum(is.na(ins)),empty_string=sum(!is.na(ins)&ins==''),whitespace_only=sum(!is.na(ins)&ins!=''&trimws(ins)==''),recorded_nonblank=sum(!is.na(ins)&trimws(ins)!=''),categories=tab,limit='Only raw categories; classification/mapping exclusions must be checked against missing original overall generator. NA+blank==33 can support raw absent33; unequal count means mapping mismatch, not automatic missing.')
write_json(insurance,file.path(out,'INSURANCE_RAW_AGGREGATE.json'),auto_unbox=TRUE,pretty=TRUE)
# Current Main display contract, separate from legacyCSV-only omitted analytes.
contract<-fromJSON(file.path(Sys.getenv('CRRT_REPO_DIR'),'revision_2026_11/figures_tables/dependencies/DISPLAYED_MAIN_TABLE1_CONTRACT.json'),simplifyVector=FALSE)
ui<-match(as.character(d$stay_id),as.character(u$stay_id));ki<-match(as.character(d$stay_id),as.character(k$stay_id))
get_source<-function(src){parts<-strsplit(src,'.',fixed=TRUE)[[1]];if(parts[1]=='raw')d[[parts[2]]]else if(parts[1]=='augmentation')u[[parts[2]]][ui]else if(parts[1]=='keys')k[[parts[2]]][ki]else NULL}
numtokens<-function(z)as.numeric(regmatches(gsub(',','',z,fixed=TRUE),gregexpr('[0-9]+(?:\\.[0-9]+)?',gsub(',','',z,fixed=TRUE),perl=TRUE))[[1]])
# Mechanical aliases of explicit ICU unit suffixes; save aggregate levels for audit.
care<-as.character(d$first_careunit);care_code<-ifelse(grepl('\\([^()]+\\)$',care),sub('.*\\(([^()]+)\\)$','\\1',care),care)
write.csv(as.data.frame(table(raw=care,display=care_code,useNA='ifany')),file.path(out,'ICU_RAW_DISPLAY_LEVEL_AGGREGATE.csv'),row.names=FALSE)
race<-race_group(d$race)
display_rows<-lapply(contract$rows,function(rr){
 kind<-rr$kind;src<-rr$source;calc<-NA_character_;passed<-FALSE
 if(kind=='section_or_category_header')return(data.frame(row=rr$row_index_zero_based,label=rr$label,expected=rr$current_display,calculated='',status='STRUCTURAL_NOT_A_DATA_GATE'))
 if(kind=='header_N'){calc<-as.character(nrow(k));passed<-nrow(k)==rr$expected_old_N}
 else if(kind=='admission_era'){
  if('anchor_year_group'%in%names(k)){yrs<-as.numeric(unlist(regmatches(k$anchor_year_group,gregexpr('[0-9]{4}',k$anchor_year_group))));calc<-paste(range(yrs,na.rm=TRUE),collapse='–');passed<-identical(calc,rr$current_display)}
 }else if(kind=='continuous'){
  x<-get_source(src);if(!is.null(x)){expected<-numtokens(rr$current_display);dg<-if(grepl('^[0-9]+\\.',rr$current_display))nchar(sub('^[0-9]+\\.([0-9]+).*','\\1',rr$current_display))else 0L;calc<-fmt_cont(x,dg);passed<-identical(numtokens(calc),expected)}
 }else if(kind=='missing'){
  x<-get_source(src);if(!is.null(x)){calc<-miss(x);passed<-identical(numtokens(calc),numtokens(rr$current_display))}
 }else if(kind=='binary'){
  x<-get_source(src);if(!is.null(x)){nn<-sum(x%in%c(TRUE,1));pct<-round(100*nn/nrow(k),1);calc<-paste(nn,pct);ex<-numtokens(rr$current_display);passed<-nn==ex[1]&&(if(grepl('<0.1%',rr$current_display,fixed=TRUE))nn>0&&100*nn/nrow(k)<.1 else pct==ex[2])}
 }else if(kind=='category'){
  sec<-rr$section;denom<-nrow(k);nn<-NA_real_
  if(sec=='ICU type, n (%)')nn<-if(rr$category=='Other')sum(!care_code%in%c('MICU','MICU/SICU','CCU','CVICU','SICU','TSICU'))else sum(care_code==rr$category,na.rm=TRUE)
  if(sec=='Sex, n (%)')nn<-sum(d$gender==if(rr$category=='Female')'F'else'M',na.rm=TRUE)
  if(sec=='Race/Ethnicity, n (%)'){lv<-c('Black or African American'='Black/African American','Hispanic or Latino'='Hispanic/Latino','Other / Unknown'='Other/Unknown');target<-if(rr$category%in%names(lv))unname(lv[rr$category])else rr$category;nn<-sum(race==target)}
  if(sec=='Insurance type, n (%)'){denom<-sum(!is.na(d$insurance));nn<-sum(d$insurance==rr$category,na.rm=TRUE)}
  if(sec=='Index-defining KDIGO stage 3 pathway, n (%)'){lv<-c('Urine output criterion only'='oliguria','Creatinine criterion only'='creatinine','Both criteria'='both','Indeterminate'='unknown');lv<-c('Urine output criterion only'='urine_output_only','Creatinine criterion only'='creatinine_only','Both criteria'='both_creatinine_and_uo','Indeterminate'='unknown');nn<-sum(k$stage3_path[ki]==lv[rr$category],na.rm=TRUE)}
  if(!is.na(nn)){pct<-round(100*nn/denom,1);calc<-paste(nn,pct);ex<-numtokens(rr$current_display);passed<-nn==ex[1]&&(if(grepl('<0.1%',rr$current_display,fixed=TRUE))nn>0&&100*nn/denom<.1 else pct==ex[2])}
 }
 data.frame(row=rr$row_index_zero_based,label=rr$label,expected=rr$current_display,calculated=calc,status=if(mode=='new')'NEW_VALUE_DIAGNOSTIC_NOT_OLD_COMPARISON'else if(isTRUE(passed))'EXACT_DISPLAYED_OLD'else if(is.na(calc))'UNRESOLVED_DISPLAYED_SOURCE'else'MISMATCH_DISPLAYED_DO_NOT_PROMOTE')
})
display_cmp<-do.call(rbind,display_rows);write.csv(display_cmp,file.path(out,'DISPLAYED_MAIN_FIELD_REPLAY.csv'),row.names=FALSE)
main_pass<-mode=='replay_old'&&all(display_cmp$status%in%c('EXACT_DISPLAYED_OLD','STRUCTURAL_NOT_A_DATA_GATE'))
write_json(list(mode=mode,current_main_display_rows=nrow(display_cmp),all_current_manuscript_rows_verified=main_pass,reference_scope='ActualMain117rows; omittedlegacyCSVfields do notfailMain',unresolved_or_mismatched=display_cmp[display_cmp$status%in%c('UNRESOLVED_DISPLAYED_SOURCE','MISMATCH_DISPLAYED_DO_NOT_PROMOTE'),],keys_n=nrow(k)),file.path(out,'DISPLAYED_MAIN_GATE.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
# Gate scopes displayedMain only; legacyCSV differences remain separately reported.
gate<-list(mode=mode,keys_n=nrow(k),data_select_ran_readonly=TRUE,keys_and_t0_exact=TRUE,all_current_manuscript_rows_verified=main_pass,original_sql_exact_rows=sum(cmp$status=='EXACT_ORIGINAL_SQL_REPLAY'),mismatch_rows=sum(cmp$status=='MISMATCH_DO_NOT_PROMOTE'),unresolved_rows=sum(grepl('UNRESOLVED',cmp$status)),reason='CurrentMaindisplaygate only;CSV-onlyalbumin/bilirubin/extra rows excluded. Displayed sourcebranches andcategoricalmappings must reproduceOLDexact beforepromotion; no nearestnumericfit.',formal_files_modified=FALSE)
write_json(gate,file.path(out,'REPLAY_GATE.json'),auto_unbox=TRUE,pretty=TRUE)
cat(if(main_pass)'Read-only replay saved; CURRENT MAIN DISPLAYED ROWS PASS. No scientific/model clearance claimed.\n'else'Read-only replay saved; displayedMain mismatches/unresolved rows HOLD. CSV-only omissions do notfailMain. No patient rows printed.\n')
