suppressPackageStartupMessages({library(dplyr);library(survival);library(readr);library(jsonlite)})
p<-'09_工作稿/out/ch24_q4_corrected_fig3';src<-'09_工作稿/out/run15/run15_NEW_anchors_L0_era.R';x<-parse(src,keep.source=TRUE)
asgn<-function(n){i<-which(vapply(x,function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&identical(as.character(e[[2]]),n),logical(1)));stopifnot(length(i)==1);x[[i]]}
for(n in c('mk','cox1'))eval(asgn(n))
selected<-which(vapply(x,function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&any(vapply(c('era_strat_eq','era_strat_all','l0_new_eq','l0_new_all'),function(k)identical(e[[2]],bquote(rows[[.(k)]])),logical(1))),logical(1)))
stopifnot(length(selected)==4)
writeLines(unlist(lapply(c(list(asgn('mk'),asgn('cox1')),as.list(x[selected])),deparse)),file.path(p,'BRIDGE_original_functions_calls_AST.R'))
dod<-readRDS('09_工作稿/out/ch24_q4_other_cohort_plan/dod_plan/UNION_subject_dod_FROM_DB.rds')
run<-function(canonfile,l0file){
 canon<-readRDS(canonfile);L<-readRDS(l0file);dc<-canon%>%distinct(clone_id,.keep_all=TRUE)%>%as.data.frame();dc$sf<-factor(dc$strategy,levels=c('B_deferred','A_early'),labels=c('Deferred CRRT','Early CRRT'));l0<-mk(L$ccw_l0);rows<-list();for(i in selected)eval(x[[i]]);bind_rows(rows)
}
a<-run('08_CH01/out/ccw_ch01_dod.rds',file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/results/strictK3_v2_l0_only_sensitivity.rds"))
# Keep exact original order of selected four rows and raw numeric strings in reference.
refpath<-'09_工作稿/out/run15/NEW/anchors_L0_era_rows.csv';raw<-readLines(refpath);use<-grepl('L0-only denominator (NEW CH-01)',raw,fixed=TRUE)|grepl('Calendar-era stratified Cox (NEW CH-01)',raw,fixed=TRUE)
writeLines(c(raw[1],raw[use]),file.path(p,'BRIDGE_OLD_REFERENCE_4.csv'));write_csv(a,file.path(p,'BRIDGE_OLD_4.csv'))
bytes<-function(f)readBin(f,'raw',n=file.info(f)$size);eq<-identical(bytes(file.path(p,'BRIDGE_OLD_REFERENCE_4.csv')),bytes(file.path(p,'BRIDGE_OLD_4.csv')))
write_json(list(status=if(eq)'PASS_FOUR_ROWS_BYTE_EXACT'else'FAIL_STOP',source=src,original_function_and_calls_AST_unchanged=TRUE,original_CH01_corrected_follow_columns_for_era=TRUE,original_mk_for_L0=TRUE,byte_exact=eq),file.path(p,'BRIDGE_OLD_GATE.json'),pretty=TRUE,auto_unbox=TRUE)
stopifnot(eq)
b<-run('09_工作稿/out/ch24_q4_corrected_primary/ccw_ch01_dod.rds','09_工作稿/out/ch24_q4_corrected_l0/NEW/results/strictK3_v2_l0_only_sensitivity.rds');stopifnot(identical(a$analysis,b$analysis),identical(a$population,b$population),identical(a$n_patients,b$n_patients))
write_csv(b,file.path(p,'BRIDGE_NEW_4.csv'))
print(b)
