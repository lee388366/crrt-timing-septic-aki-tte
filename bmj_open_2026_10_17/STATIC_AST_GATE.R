#!/usr/bin/env Rscript
# Parse/structural check only. Does not evaluate any analysis/model/database expression.
suppressPackageStartupMessages(library(jsonlite))
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
map<-fromJSON(args[1]);code<-args[2];output<-args[3]
config<-c('OUT','B','V2_DIR','OUT_DIR','base_dir','strict_dir','data_dir','results_dir','fig_dir','ccw_path','secondary_csv_path','dod','a_raw','k_raw','con')
lhs<-function(x)if(is.call(x)&&paste(deparse(x[[1]]),collapse='')%in%c('<-','='))paste(deparse(x[[2]]),collapse='')else''
normalize<-function(x){
 if(is.character(x)&&length(x)==1L&&grepl('/',x,fixed=TRUE))return('__IO_PATH__')
 if(!is.call(x))return(x)
 if(paste(deparse(x[[1]]),collapse='')%in%c('file.path','Sys.getenv'))return('__IO_PATH__')
 as.call(lapply(as.list(x),normalize))
}
iskeep<-function(x){
 if(lhs(x)%in%config)return(FALSE)
 if(is.call(x)&&paste(deparse(x[[1]]),collapse='')%in%c('on.exit','dir.create','readr::write_csv'))return(FALSE)
 TRUE
}
rows<-list()
for(i in seq_len(nrow(map))){
 if(grepl('current_sources/analysis/',map$destination[i],fixed=TRUE))next
 a<-as.list(parse(map$source[i],keep.source=FALSE));b<-as.list(parse(file.path(code,map$destination[i]),keep.source=FALSE))
 funcs<-function(p)Filter(function(x)is.call(x)&&paste(deparse(x[[1]]),collapse='')=='<-'&&is.call(x[[3]])&&paste(deparse(x[[3]][[1]]),collapse='')=='function',p)
 fa<-funcs(a);fb<-funcs(b)
 function_exact<-identical(fa,fb)
 # Paths are normalized only; all other eligible expressions are compared in original sequence.
 aa<-lapply(Filter(iskeep,a),normalize);bb<-lapply(Filter(iskeep,b),normalize)
 whole<-identical(aa,bb)
 rows[[i]]<-list(destination=map$destination[i],function_definitions=length(fa),function_AST_exact=function_exact,non_IO_expression_sequence_exact=whole,source_expression_count=length(a),candidate_expression_count=length(b))
}
rows<-Filter(Negate(is.null),rows)
passed<-all(vapply(rows,function(x)x$function_AST_exact&&x$non_IO_expression_sequence_exact,logical(1)))
write_json(list(status=if(passed)'PASS_STATIC_AST_IO_ONLY'else'FAIL_REVIEW',passed=passed,rows=rows,analysis_evaluated=FALSE),output,auto_unbox=TRUE,pretty=TRUE)
if(!passed)quit(status=1L)
