#!/usr/bin/env Rscript
# I/O-only orchestration of original function/call AST; do not source make_variant.
suppressPackageStartupMessages({library(dplyr);library(readr);library(survival);library(tibble)})
code <- Sys.getenv('CRRT_CODE_ROOT')
root <- file.path(code,'current_sources/rescue')
ch01_dod <- readRDS(Sys.getenv('CRRT_DOD_RDS'))
stopifnot(inherits(ch01_dod$dod,'Date'),!anyDuplicated(ch01_dod$subject_id))
source(file.path(root,'run_cox_ORIGINAL_AST.R'))
source(file.path(root,'ch01_recode_ORIGINAL_AST.R'))
make_expr <- parse(file.path(root,'make_variant_ORIGINAL_AST.R'),keep.source=FALSE)[[1]]
original_call <- tail(as.list(make_expr[[3]][[3]]),1)[[1]]
stopifnot(identical(original_call[[1]],as.name('bind_rows')))
variant_expr <- parse(file.path(root,'variant_cohorts_ORIGINAL_AST.R'),keep.source=FALSE)[[1]]
branch_names <- names(as.list(variant_expr[[3]]))[-1]
branch_names <- branch_names[branch_names!='primary_R1_R2_R2b_R3']
stopifnot(length(branch_names)==5L)
results <- list()
for(name in branch_names){
 ccw <- readRDS(file.path(Sys.getenv('CRRT_RESCUE_WEIGHT_ROOT'),paste0('rescue_',name),'NEW/ccw_final.rds'))
 check <- ch01_recode(ccw)
 stopifnot(identical(ccw$died_28d,check$died_28d),identical(ccw$follow_time_28d,check$follow_time_28d))
 results[[name]] <- eval(original_call)
}
eval(parse(file.path(root,'out_ORIGINAL_AST.R'),keep.source=FALSE)[[1]])
dir.create(Sys.getenv('CRRT_RESCUE_OUT'),recursive=TRUE,showWarnings=FALSE)
write_csv(out,file.path(Sys.getenv('CRRT_RESCUE_OUT'),'strictK3_v2_rescue_sensitivity_model_results.csv'))
