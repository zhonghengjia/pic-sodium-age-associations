suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(jsonlite);library(digest)})
entry <- normalizePath(sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]),winslash='/',mustWork=TRUE)
R <- normalizePath(file.path(dirname(entry),'../..'),winslash='/',mustWork=TRUE)
P <- dirname(R); C <- dirname(entry); Q <- file.path(R,'reviews/figure_worker'); F <- file.path(R,'figures')
dir.create(F,showWarnings=FALSE)
sha <- function(p)digest(file=p,algo='sha256',serialize=FALSE)
wj <- function(x,n)write_json(x,file.path(Q,n),auto_unbox=TRUE,pretty=TRUE,digits=17)
protected <- list.files(file.path(P,'history_clinical_language_20261011_v16'),recursive=TRUE,full.names=TRUE)
protected <- protected[!dir.exists(protected)]; ph <- vapply(protected,sha,'')
b <- fromJSON(file.path(P,'analysis_daily_mean_20261009_v1/figure_association_roundticks_20261009_v10/control/SOURCE_BINDINGS.json'),simplifyVector=FALSE)
paths <- vapply(b$sources,function(x)x$path,'')
all_inputs <- c(paths[c('full_curves','rug','joint_tests','support','accepted_review','caption')],
 file.path(P,'source/Figure_1_flow_nodes.csv'),
 file.path(P,'KM_main_annotation_removed_20261009_v7/source_data',c('logrank_tests.csv','km_estimates_full.csv','alive_numbers.csv','group_counts_overall_age.csv')),
 file.path(P,'side_day1_day3_combined_ROC_age_20261010_v1/results',c('ROC_coordinates.csv','Three_model_AUC_age.csv')),
 file.path(P,'sodium_day3_CIF_color_revision_20261009_v2/source_data',c('CIF_full_followup.csv','risk_table.csv')),
 file.path(P,'sodium_day3_CIF_color_revision_20261009_v2/results/display_counts.csv'))
ih <- vapply(all_inputs,sha,'')
wj(as.list(setNames(ih,all_inputs)),'input_sha256.json')
newenv <- function(name){
 e<-new.env(parent=globalenv());e$paths<-paths;e$sha<-sha;e$before<-vapply(paths,sha,'');e$base<-P
 if(name=='Figure4')e$P<-file.path(P,'KM_main_annotation_removed_20261009_v7/source_data')
 if(name=='Supplementary_Figure1')e$O<-file.path(P,'sodium_day3_CIF_color_revision_20261009_v2')
 e
}
compare <- function(a,b,figure,panel,max_layers=NULL){
 x<-ggplot_build(a);y<-ggplot_build(b);n<-if(is.null(max_layers))length(x$data) else max_layers
 stopifnot(length(y$data)>=n)
 for(i in seq_len(n))stopifnot(identical(x$data[[i]],y$data[[i]]))
 stopifnot(identical(x$layout$panel_params[[1]]$x.range,y$layout$panel_params[[1]]$x.range),identical(x$layout$panel_params[[1]]$y.range,y$layout$panel_params[[1]]$y.range))
 saveRDS(x$data[seq_len(n)],file.path(Q,paste0(figure,'_',panel,'_baseline_geoms.rds')))
 saveRDS(y$data[seq_len(n)],file.path(Q,paste0(figure,'_',panel,'_revised_geoms.rds')))
 data.frame(figure=figure,panel=panel,layers=n,geom_rows=sum(vapply(x$data[seq_len(n)],nrow,0L)),identical_geoms=TRUE,identical_numeric_coordinate_ranges=TRUE)
}
saveplot <- function(fig,name,width,height){
 ggsave(file.path(F,paste0(name,'.png')),fig,width=width,height=height,units='in',dpi=300,device=ragg::agg_png,bg='white',limitsize=FALSE)
 ggsave(file.path(F,paste0(name,'.pdf')),fig,width=width,height=height,units='in',device=cairo_pdf,family='Arial',bg='white',limitsize=FALSE)
}
checks<-list();warnings<-character()
withCallingHandlers({
 for(name in c('Figure1','Figure3','Figure4','Figure5','Supplementary_Figure1')){
  a<-newenv(name);z<-newenv(name)
  sys.source(file.path(C,paste0(name,'_baseline.R')),envir=a);sys.source(file.path(C,paste0(name,'_revised.R')),envir=z)
  if(name=='Figure1'){
   flow<-read.csv(file.path(P,'source/Figure_1_flow_nodes.csv'))
   old<-a$flow_plot(flow);fig<-z$flow_plot(flow)
   x<-ggplot_build(old);y<-ggplot_build(fig)
   stopifnot(nrow(x$data[[1]])==10,nrow(y$data[[1]])==9)
   for(i in 1:2){
    retained<-x$data[[i]][1:9,]; revised<-y$data[[i]]
    # ggplot2 attaches a subset-related 'n' attribute to group; values and every geom property must be exact.
    retained$group<-as.integer(retained$group);revised$group<-as.integer(revised$group)
    stopifnot(identical(retained,revised))
   }
   stopifnot(identical(x$data[[3]],y$data[[3]]),nrow(y$data[[3]])==8)
   saveRDS(list(baseline=x$data,revised=y$data),file.path(Q,'Figure1_flow_geometry.rds'))
   checks[[length(checks)+1]]<-data.frame(figure=name,panel='flow',layers=3,geom_rows=26,identical_geoms=TRUE,identical_numeric_coordinate_ranges=FALSE)
   # Height reduced only to remove the age-box region; all retained node and edge coordinates are exact.
   saveplot(fig,name,10,6.8)
  } else if(name=='Figure3'){
   for(i in 1:18)checks[[length(checks)+1]]<-compare(a$plots[[i]],z$plots[[i]],name,paste0('panel',i))
   stopifnot(all(vapply(z$plots,function(p)p$labels$x=='Sodium (mmol/L)',TRUE)))
   saveplot(z$p,name,210/25.4,285/25.4)
  } else if(name=='Figure4'){
   for(i in 1:6){
    checks[[length(checks)+1]]<-compare(a$panels[[i]]$patches$plots[[1]],z$panels[[i]]$patches$plots[[1]],name,paste0('curve',i))
    checks[[length(checks)+1]]<-compare(a$panels[[i]],z$panels[[i]],name,paste0('alive_table',i))
   }
   stopifnot(identical(a$rrall,z$rrall),identical(a$ccall,z$ccall))
   write.csv(z$rrall,file.path(Q,'Figure4_alive_numbers_verification.csv'),row.names=FALSE)
   saveplot(z$fig,name,15.5,10.7)
  } else if(name=='Figure5'){
   for(i in 1:6)checks[[length(checks)+1]]<-compare(a$plots[[i]]$patches$plots[[1]],z$plots[[i]],name,paste0('ROC',i),2)
   stopifnot(identical(a$curves,z$curves),identical(a$tab,z$tab),nrow(z$tab)==12)
   write.csv(z$tab,file.path(Q,'Figure5_AUC_source_values_including_preserved_CI.csv'),row.names=FALSE)
   saveplot(z$fig,name,12.5,8.65)
  } else {
   for(i in 1:6){
    checks[[length(checks)+1]]<-compare(a$panels[[i]]$patches$plots[[1]],z$panels[[i]]$patches$plots[[1]],name,paste0('CIF',i))
    checks[[length(checks)+1]]<-compare(a$panels[[i]],z$panels[[i]],name,paste0('risk_table',i))
   }
   stopifnot(identical(a$cv,z$cv),identical(a$rt,z$rt))
   write.csv(z$rt,file.path(Q,'Supplementary_Figure1_risk_table_verification.csv'),row.names=FALSE)
   saveplot(z$fig,name,7.2,7.9)
  }
  cat(name,'R export complete; fixed-data geometry validated.\n')
 }
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
stopifnot(identical(ph,vapply(protected,sha,'')),identical(ih,vapply(all_inputs,sha,'')))
write.csv(do.call(rbind,checks),file.path(Q,'geom_invariance_checks.csv'),row.names=FALSE)
outputs<-list.files(F,pattern='^(Figure[1345]|Supplementary_Figure1)\\.(png|pdf)$',full.names=TRUE)
wj(list(status='R_RENDER_GEOMETRY_VERIFIED_PENDING_VISUAL_REVIEW',mode='IMPLEMENT_APPROVED',source_based_R_render=TRUE,
 no_model_fit=TRUE,no_new_CI=TRUE,no_resampling=TRUE,no_patient_records_read=TRUE,png_dpi=300,vector_pdf_device='cairo_pdf',
 source_inputs_unchanged=length(all_inputs),V16_files_hash_unchanged=length(protected),geom_comparison_count=length(checks),
 warnings=warnings,output_hashes=as.list(setNames(vapply(outputs,sha,''),basename(outputs))),sessionInfo=capture.output(sessionInfo())), 'render_receipt.json')
cat('All source and V16 hashes unchanged. No fitting or uncertainty estimation executed.\n')
