# Display-only consumer of independently reviewed aggregate results.
# CLI: helper render_figures.R template|preflight|render
# No patient records, MI objects, model fitting, ROC estimation or recalibration.
suppressPackageStartupMessages({
  library(jsonlite); library(digest); library(ggplot2); library(patchwork); library(ragg)
})
sha <- function(path) digest(file=path, algo='sha256', serialize=FALSE)
writej <- function(x,path) write_json(x,path,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null')
check <- function(x,message) if(!isTRUE(x)) stop(message,call.=FALSE)
scalar_equal <- function(x,y) is.atomic(x)&&length(x)==1L&&!is.na(x)&&
  is.atomic(y)&&length(y)==1L&&!is.na(y)&&isTRUE(x==y)
same_map <- function(x,y) is.list(x)&&!is.null(names(x))&&!anyDuplicated(names(x))&&
  length(x)==length(y)&&setequal(names(x),names(y))&&
  all(vapply(names(y),function(n)is.character(x[[n]])&&scalar_equal(x[[n]],y[[n]]),logical(1)))
check_release <- function(gate,expected) {
  for(n in setdiff(names(expected),'csv_hashes'))
    check(scalar_equal(gate[[n]],expected[[n]]),paste('Display release mismatch:',n))
  check(same_map(gate$csv_hashes,expected$csv_hashes),'Display release CSV map mismatch')
  invisible(TRUE)
}
ptext <- function(x) ifelse(x<.001,'P<0.001',sprintf('P=%.3f',x))
models <- c('D1','D3','D1D3')
model_names <- c(D1='Clinical + day 1',D3='Clinical + day 3',D1D3='Clinical + days 1 and 3')
model_colours <- c(D1='#9DD1FA',D3='#FF6666',D1D3='#459F91')
model_types <- c(D1='dashed',D3='solid',D1D3='dotdash')
layers <- paste0('L0',1:5)
age_names <- c('Neonates (\u226428 d)','Infants (>28 d-<1 y)','1-4 y','5-9 y','10-14 y')
target_names <- c(L_vs_N='Day-3 low vs normal',H_vs_N='Day-3 high vs normal')
theme_plot <- function() theme_classic(base_size=10,base_family='Arial')+
  theme(axis.text=element_text(size=9,colour='black'),axis.title=element_text(size=10),
        plot.title=element_text(size=11,face='bold'),plot.subtitle=element_text(size=9),
        legend.text=element_text(size=9),legend.title=element_blank(),
        plot.margin=margin(4,4,4,4,'mm'))
probability_axes <- function() list(
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=c(0,0)),
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=c(0,0)),
  coord_equal(expand=FALSE))
validate_probability <- function(x,name) check(all(is.finite(x))&&all(x>=0 & x<=1),paste(name,'outside [0,1]'))
validate_columns <- function(x,required,name) check(all(required %in% names(x)),paste(name,'missing required columns'))
axis_from_ci <- function(z) {
  check(all(is.finite(as.matrix(z[c('OR','ci_low','ci_high')]))),'Non-finite OR/CI')
  check(all(z$ci_low>0 & z$ci_low<=z$OR & z$OR<=z$ci_high),'Invalid positive OR interval')
  lo <- min(c(1,z$ci_low))*.94; hi <- max(c(1,z$ci_high))*1.06
  nice <- sort(unique(as.vector(outer(c(1,2,5),10^seq(floor(log10(lo))-1,ceiling(log10(hi))+1)))))
  limits <- c(max(nice[nice<=lo]),min(nice[nice>=hi]))
  check(all(z$ci_low>=limits[1] & z$ci_high<=limits[2])&&limits[1]<=1&&limits[2]>=1,'OR axis clips source intervals/null')
  list(limits=limits,breaks=nice[nice>=limits[1]&nice<=limits[2]])
}
roc_plot <- function(roc,performance) {
  labels <- setNames(vapply(models,function(m){z<-performance[performance$model_id==m,]
    sprintf('%s: AUC %.3f (%.3f-%.3f)',model_names[[m]],z$AUC,z$auc_low,z$auc_high)},''),models)
  ggplot(roc,aes(false_positive_rate,sensitivity,colour=model_id,linetype=model_id,group=model_id))+
    geom_abline(intercept=0,slope=1,colour='#459F91',linetype='dotted',linewidth=.45)+
    geom_path(linewidth=.85)+probability_axes()+
    scale_colour_manual(values=model_colours,breaks=models,labels=labels)+
    scale_linetype_manual(values=model_types,breaks=models,labels=labels)+
    labs(x='False-positive rate (1 - specificity)',y='Sensitivity',title='Exploratory OOF discrimination',
         subtitle='Same 72-hour landmark cohort; repeated cross-validation risk scores')+
    theme_plot()+theme(legend.position='bottom',legend.direction='vertical')
}
calibration_panel <- function(bins,performance,model) {
  z <- bins[bins$model_id==model,];z <- z[order(z$bin),]
  a <- performance[performance$model_id==model,]
  p <- ggplot(z,aes(mean_pred,observed_rate))+
    geom_abline(intercept=0,slope=1,colour='#459F91',linetype='dotted',linewidth=.45)+
    geom_ribbon(aes(ymin=ci_low,ymax=ci_high),fill=model_colours[[model]],alpha=.18)+
    geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.015,colour=model_colours[[model]],linewidth=.4)+
    geom_point(size=1.9,colour=model_colours[[model]])+probability_axes()+
    labs(x='Mean predicted risk',y='Observed mortality',title=model_names[[model]],
         subtitle=sprintf('Calibration: intercept %.3f; slope %.3f',a$calibration_intercept,a$calibration_slope))+
    theme_plot()+theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=8),
                       axis.title=element_text(size=9),axis.text=element_text(size=8))
  # Explicit bin counts are below the axes, not annotations over nearby risk points.
  z$x <- (z$bin-1)%%5+1;z$y <- ifelse(z$bin<=5,2,1)
  strip <- ggplot(z,aes(x,y,label=paste0('Bin ',bin,'\nn=',N)))+
    geom_text(size=2.7,family='Arial',colour='#3C7AAB')+
    scale_x_continuous(limits=c(.4,5.6),expand=c(0,0))+
    scale_y_continuous(limits=c(.4,2.6),expand=c(0,0))+theme_void()+
    theme(plot.margin=margin(0,3,0,3,'mm'))
  wrap_plots(list(wrap_elements(full=p),wrap_elements(full=strip)),ncol=1,heights=c(5,1.25))
}
association_panel <- function(z,title,joint,axis) {
  # gtable construction also occurs before ggsave; never create an implicit PDF device.
  if(grDevices::dev.cur()==1L){
    ragg::agg_capture(width=230,height=230,units='mm',res=96)
    assembly_device<-grDevices::dev.cur();on.exit(grDevices::dev.off(assembly_device),add=TRUE)
  }
  z <- z[match(c('ALL',layers),z$age_layer),];z$y <- 6:1
  z$label <- c('Overall common association',age_names)
  z$shape <- ifelse(z$age_layer=='ALL','Common','Representative age')
  yl <- c(.4,6.6)
  left <- ggplot(z,aes(.02,y,label=label))+geom_text(hjust=0,size=3.05,family='Arial')+
    scale_x_continuous(limits=c(0,1),expand=c(0,0))+
    scale_y_continuous(limits=yl,expand=c(0,0))+theme_void()+theme(plot.margin=margin(3,0,3,3,'mm'))
  p <- ggplot(z,aes(OR,y))+
    geom_vline(xintercept=1,colour='#459F91',linetype='dashed',linewidth=.45)+
    geom_segment(aes(x=ci_low,xend=ci_high,yend=y),colour='#F5B6B6',linewidth=1.1)+
    geom_point(aes(shape=shape),colour='#FF6666',size=2.5)+
    scale_shape_manual(values=c(Common=15,'Representative age'=16),guide='none')+
    scale_x_log10(limits=axis$limits,breaks=axis$breaks,labels=function(x)format(x,trim=TRUE,scientific=FALSE))+
    scale_y_continuous(limits=yl,expand=c(0,0))+
    labs(x='Adjusted OR (95% CI)',y=NULL)+theme_plot()+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          plot.margin=margin(3,2,3,2,'mm'))
  z$support <- paste0(z$target_N,'/',z$target_deaths,'  |  ',z$reference_N,'/',z$reference_deaths)
  right <- ggplot(z,aes(.5,y,label=support))+geom_text(size=3,family='Arial',colour='#3C7AAB')+
    scale_x_continuous(limits=c(0,1),expand=c(0,0))+
    scale_y_continuous(limits=yl,expand=c(0,0))+theme_void()+theme(plot.margin=margin(3,1,3,1,'mm'))
  # Shared gtable row heights include the middle x-axis space in both text columns.
  # Wrapping independent text plots without it shifts later labels away from their CI.
  aligned <- lapply(list(left,p,right),ggplotGrob)
  panel_rows <- lapply(aligned,function(g)g$layout[g$layout$name=='panel',c('t','b')])
  check(all(vapply(aligned,function(g)length(g$heights),0L)==length(aligned[[2]]$heights))&&
    identical(panel_rows[[1]],panel_rows[[2]])&&identical(panel_rows[[2]],panel_rows[[3]]),
    'Association column panel structures do not match')
  shared_heights <- grid::unit.pmax(aligned[[1]]$heights,aligned[[2]]$heights,aligned[[3]]$heights)
  aligned <- lapply(aligned,function(g){g$heights<-shared_heights;g})
  body <- wrap_plots(lapply(aligned,function(g)wrap_elements(full=g)),
                     nrow=1,widths=c(2.8,3.4,2.4))
  header <- ggplot()+annotate('text',x=0,y=.78,label=title,hjust=0,size=3.6,fontface='bold',family='Arial')+
    annotate('text',x=0,y=.27,label=paste('Continuous-age interaction (2 df; Holm):',ptext(joint$p_holm)),hjust=0,size=3,family='Arial')+
    annotate('text',x=.91,y=.27,label='Target | Reference\nN/deaths',size=2.8,family='Arial',colour='#3C7AAB')+
    xlim(0,1.1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(3,3,0,3,'mm'))
  wrap_plots(list(wrap_elements(full=header),wrap_elements(full=body)),ncol=1,heights=c(1,5))
}
reference_receipt <- function() list(
  scope='Form references only; prior full-text reading reused, with current targeted re-read and actual Figure 6a view. No new 20-paper reading claim.',
  local_paper=list(id='D144',title='Predicting Risk of Mortality in Pediatric ICU Based on Ensemble Step-Wise Feature Selection',
    DOI='10.34133/2021/9365125',source=Sys.getenv("PIC_REFERENCE_XML"),
    source_sha256='32a97daa3740f094a61000b258ed9575f67761d538a56da4966fbee6bbb82f34',
    native_pdf_sha256='cf82f695646ac38b4eb2cde627310ac58055457f4310e17f3750523434ddc36b',
    figure='Figure 6a, native PDF page 6',read_units='Current re-read: performance-evaluation Methods, Figure 6 caption/Results, PDF page 6 actual pixels; prior full body card D144 reused.',
    applicable='Multiple model risk-score ROC curves on common false-positive-rate/sensitivity axes, model-specific AUC legend.',
    not_transferred='Holdout design, feature ranking, median filling, PRISM comparison, ROC-AUC SD and selection logic are not adopted.'),
  official_docs=list(
    coordinate=list(URL='https://ggplot2.tidyverse.org/reference/coord_fixed.html',read='ratio=1 and expand=FALSE arguments',applied='Equal physical unit scales and exact [0,1] boundaries.'),
    scale=list(URL='https://ggplot2.tidyverse.org/reference/scale_continuous.html',read='limits and expand arguments',applied='Exact shared [0,1] limits; no automatic expansion. Out-of-limits data must fail before plotting.')))
main <- function() {
  args <- commandArgs(trailingOnly=TRUE);phase <- if(length(args))args[1]else 'template'
  check(phase %in% c('template','preflight','render','preflight-figure3','render-figure3'),'Unknown display phase')
  figure3_only <- phase %in% c('preflight-figure3','render-figure3')
  entry <- normalizePath(sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]),winslash='/',mustWork=TRUE)
  root <- normalizePath(file.path(dirname(entry),'..'),winslash='/',mustWork=TRUE)
  notes <- file.path(root,'figures/notes');delivery <- file.path(root,'交付审核图形')
  dir.create(notes,recursive=TRUE,showWarnings=FALSE)
  controls <- c(plan=file.path(root,'control/ANALYSIS_PLAN.md'),config=file.path(root,'control/ANALYSIS_CONFIG.json'),
                bindings=file.path(root,'control/SOURCE_BINDINGS.json'),protection=file.path(root,'control/PROTECTED_MANIFEST.json'))
  control_hashes <- vapply(controls,sha,'');cfg <- fromJSON(controls[['config']])
  check(cfg$cohort_N==3459&&cfg$deaths==306&&cfg$folds==5&&cfg$repeats==5&&cfg$prediction_M==10&&cfg$association_M==50,'Plan population/validation configuration drift')
  pm <- fromJSON(controls[['protection']],simplifyVector=FALSE)$protected_files
  check(is.list(pm)&&length(pm)>0L&&!anyDuplicated(names(pm)),'Invalid protected map')
  if(phase=='template') {
    writej(reference_receipt(),file.path(notes,'form_reference_receipt.json'))
    writej(list(status='TEMPLATE_READY_NO_RESULT_READ',code_sha256=sha(entry),control_hashes=as.list(control_hashes),
      protected_files=length(pm),source_schemas=c('roc_curves.csv','performance.csv','calibration_bins.csv','association_contrasts.csv','joint_tests.csv'),
      packages=as.list(vapply(c('ggplot2','patchwork','ragg','jsonlite','digest'),function(p)as.character(packageVersion(p)),'')),
      no_patient_read=TRUE,no_result_read=TRUE,no_render=TRUE,no_model_fit=TRUE,no_MI=TRUE),file.path(notes,'template_receipt.json'))
    cat('Template ready; no result CSV read or figure rendered.\n');return(invisible(NULL))
  }
  protbefore <- vapply(names(pm),sha,'');check(identical(unname(protbefore),unname(unlist(pm))),'Protected sources changed before display')
  receipt_path <- file.path(root,'results/analysis_receipt.json');review_path <- file.path(root,'reviews/result_review.json')
  rr <- fromJSON(receipt_path,simplifyVector=FALSE);review <- fromJSON(review_path,simplifyVector=FALSE)
  check(scalar_equal(rr$status,'COMPLETED_PENDING_INDEPENDENT_RESULT_REVIEW'),'Results not complete')
  check(scalar_equal(rr$fits_attempted,850)&&scalar_equal(rr$fits_passed,850),'Expected 850 successful fits not established')
  check(scalar_equal(rr$prediction_folds,25)&&scalar_equal(rr$association_models_completed,2)&&
    scalar_equal(rr$expected_attempts,850),'Completed prediction/association task counts mismatch')
  check(scalar_equal(review$status,'ACCEPT')&&scalar_equal(review$analysis_receipt_sha256,sha(receipt_path)),'Independent actual result review missing/mismatched')
  check(is.list(rr$source_hashes)&&is.list(rr$source_paths)&&length(rr$source_hashes)>0&&
    length(rr$source_hashes)==length(rr$source_paths)&&!anyDuplicated(names(rr$source_hashes))&&
    setequal(names(rr$source_hashes),names(rr$source_paths)),'Invalid statistical source path/hash map')
  for(n in names(rr$source_hashes))check(scalar_equal(rr$source_hashes[[n]],sha(rr$source_paths[[n]])),paste('Statistical source changed:',n))
  rel <- c(roc='roc_curves.csv',performance='performance.csv',calibration='calibration_bins.csv',association='association_contrasts.csv',joint='joint_tests.csv')
  paths <- setNames(file.path(root,'results',rel),names(rel));csv_hashes <- setNames(vapply(paths,sha,''),rel)
  for(n in rel)check(scalar_equal(rr$output_hashes[[n]],csv_hashes[[n]]),paste('Result CSV hash mismatch:',n))
  roc <- read.csv(paths[['roc']],stringsAsFactors=FALSE);performance <- read.csv(paths[['performance']],stringsAsFactors=FALSE)
  bins <- read.csv(paths[['calibration']],stringsAsFactors=FALSE);cs <- read.csv(paths[['association']],stringsAsFactors=FALSE)
  jt <- read.csv(paths[['joint']],stringsAsFactors=FALSE)
  validate_columns(roc,c('model_id','threshold','false_positive_rate','sensitivity','specificity'),'ROC')
  validate_columns(performance,c('model_id','N','events','AUC','auc_low','auc_high','Brier','calibration_intercept','calibration_slope','metric_basis','CI_type'),'Performance')
  validate_columns(bins,c('model_id','bin','mean_pred','N','deaths','observed_rate','ci_low','ci_high'),'Calibration')
  validate_columns(cs,c('model_id','contrast','age_layer','age_days','OR','ci_low','ci_high','target_N','target_deaths','reference_N','reference_deaths'),'Association')
  validate_columns(jt,c('test_id','role','p','p_holm'),'Joint tests')
  check(nrow(performance)==3&&setequal(performance$model_id,models)&&!anyDuplicated(performance$model_id),'Performance model identity mismatch')
  check(all(performance$N==3459 & performance$events==306)&&all(performance$metric_basis=='mean_5_repeat_OOF_probability')&&
    all(performance$CI_type=='conditional_fixed_OOF_stratified_patient_bootstrap'),'Unexpected population or performance estimand')
  validate_probability(unlist(performance[c('AUC','auc_low','auc_high','Brier')],use.names=FALSE),'Performance probability')
  check(all(performance$auc_low<=performance$AUC & performance$AUC<=performance$auc_high),'Invalid AUC CI')
  check(all(is.finite(performance$calibration_intercept)&is.finite(performance$calibration_slope)),'Calibration scalar missing')
  check(setequal(roc$model_id,models),'ROC model mismatch')
  for(n in c('false_positive_rate','sensitivity','specificity'))validate_probability(roc[[n]],paste('ROC',n))
  check(all(abs(roc$false_positive_rate+roc$specificity-1)<1e-12),'ROC x is not 1-specificity')
  for(m in models){z<-roc[roc$model_id==m,];check(any(z$false_positive_rate==0&z$sensitivity==0)&&any(z$false_positive_rate==1&z$sensitivity==1),'ROC lacks exact (0,0)/(1,1)');check(all(diff(z$false_positive_rate)>=0 & diff(z$sensitivity)>=0),'ROC source row order is not increasing risk-score curve')}
  check(nrow(bins)==30&&!anyDuplicated(paste(bins$model_id,bins$bin)),'Calibration must contain 30 unique model/bin aggregates')
  for(m in models){z<-bins[bins$model_id==m,];check(setequal(z$bin,1:10)&&sum(z$N)==3459&&sum(z$deaths)==306,'Calibration bins/population mismatch')}
  for(n in c('mean_pred','observed_rate','ci_low','ci_high'))validate_probability(bins[[n]],paste('Calibration',n))
  check(all(bins$N>0 & bins$deaths>=0 & bins$deaths<=bins$N)&&all(abs(bins$observed_rate-bins$deaths/bins$N)<1e-12)&&
    all(bins$ci_low<=bins$observed_rate & bins$ci_high>=bins$observed_rate),'Invalid calibration group interval/count')
  display <- cs[(cs$model_id=='day3_common'&cs$age_layer=='ALL')|(cs$model_id=='day3_age'&cs$age_layer %in% layers),]
  check(nrow(display)==12&&setequal(display$contrast,names(target_names))&&!anyDuplicated(paste(display$contrast,display$age_layer)),'Association display requires ALL + five representative ages for both contrasts')
  for(k in names(target_names))check(setequal(display$age_layer[display$contrast==k],c('ALL',layers)),'Missing association age estimate')
  check(all(display$target_N>=0 & display$reference_N>=0 & display$target_deaths<=display$target_N & display$reference_deaths<=display$reference_N),'Invalid support counts')
  medians <- setNames(vapply(layers,function(l){a<-unique(display$age_days[display$age_layer==l]);check(length(a)==1&&is.finite(a)&&a>=0,'Representative age missing/inconsistent');a},0),layers)
  joint <- jt[jt$test_id %in% paste0(names(target_names),'_age2df')&jt$role=='target_age_Holm2',]
  check(nrow(joint)==2&&!anyDuplicated(joint$test_id),'Missing Holm2 age-interaction tests')
  validate_probability(joint$p_holm,'Holm interaction P');axis <- axis_from_ci(display)
  identity <- list(status='ACCEPT',code_sha256=sha(entry),plan_sha256=control_hashes[['plan']],config_sha256=control_hashes[['config']],
    source_bindings_sha256=control_hashes[['bindings']],protected_manifest_sha256=control_hashes[['protection']],
    result_review_sha256=sha(review_path),analysis_receipt_sha256=sha(receipt_path),csv_hashes=as.list(csv_hashes))
  keep_paths <- character();keep_hashes <- character()
  if(figure3_only){
    history <- file.path(root,'figures/history/row_alignment_attempt01')
    history_manifest <- file.path(history,'copy_manifest.json')
    hm <- fromJSON(history_manifest,simplifyVector=FALSE)
    check(scalar_equal(hm$status,'BYTE_IDENTICAL_HISTORY_COPY_COMPLETE'),'First render archive not verified')
    for(f in hm$files)check(scalar_equal(f$sha256,sha(file.path(history,f$relative_path))),'Historical file changed')
    keep_rel <- c('Prediction_ROC.png','Prediction_ROC.pdf','Prediction_ROC_caption.txt',
                  'Prediction_calibration.png','Prediction_calibration.pdf','Prediction_calibration_caption.txt')
    keep_paths <- file.path(delivery,keep_rel);keep_hashes <- vapply(keep_paths,sha,'')
    for(i in seq_along(keep_rel))check(scalar_equal(keep_hashes[i],sha(file.path(history,'交付审核图形',keep_rel[i]))),'ROC/calibration changed since first render')
    identity$render_scope <- 'Figure3_only_row_alignment'
    identity$history_manifest_sha256 <- sha(history_manifest)
  }
  unchanged <- function(){check(identical(control_hashes,vapply(controls,sha,'')),'Controls changed during display');
    check(identical(protbefore,vapply(names(pm),sha,'')),'Protected artifacts changed during display');
    check(identical(csv_hashes,setNames(vapply(paths,sha,''),rel)),'Results changed during display');
    if(figure3_only)check(identical(keep_hashes,vapply(keep_paths,sha,'')),'ROC/calibration changed during Figure3-only operation')}
  captions <- list(
    Prediction_ROC=paste('Exploratory internal-validation ROC curves for three clinical-plus-sodium models in the same 72-hour landmark cohort. Day-1 and day-3 sodium are arithmetic means of valid distinct measurement times in [0,24) and [48,72) hours after first ICU admission. The curves use mortality-risk probabilities, not unidirectional sodium values or dichotomized sodium. Each score averages ten fold-local imputation predictions and then five out-of-fold repeats. Legend AUC intervals are patient-paired, death-stratified bootstrap intervals conditional on the saved OOF scores; they do not include the full uncertainty of training, fold assignment or imputation parameters. These are exploratory results, not external validation or confirmatory evidence of prediction superiority.'),
    Prediction_calibration=paste('Exploratory out-of-fold calibration of the same three models. Each panel has the same 0-1 axes. The dotted teal line is ideal calibration. Points show observed mortality against mean predicted risk within ten equal-rank groups; vertical intervals and light bands are pointwise Wilson 95% intervals, not simultaneous bands or smoothed/recalibrated curves. Bin sample sizes appear below each panel. Calibration intercept and slope are copied from the reviewed score-evaluation results and were not used to recalibrate predictions. OOF scores average ten fold-local imputation predictions across five repeats; the cohort comprises children alive and in hospital through the 72-hour landmark.'),
    Prediction_D3_age=paste('Adjusted associations of day-3 mean sodium state with subsequent in-hospital mortality. Low is <135 mmol/L, normal is 135-145 inclusive and high is >145, using unrounded means. The overall square is the common-association estimate from the model without day-3-state-by-age interactions. The five circles are evaluated from one whole-cohort continuous-age interaction model at representative ages, not five independent stratum regressions or stratum-average effects. Groups are equally spaced for display and points are not connected. Models adjust for age layer, continuous age, sex, first ICU type, first-24-hour albumin, creatinine, derived BUN, lactate, glucose and pH, plus day-1 sodium and its continuous-age interaction. Intervals are pointwise Rubin-pooled 95% CIs from the existing fifty clinical-variable imputations. Target and reference support entries are N/deaths; normal is the reference. Joint continuous-age interaction P values have two degrees of freedom and Holm adjustment across the two contrasts. These are association estimates, not prediction performance or treatment effects. The two author-approved deterministic sodium reconstructions remain; their uncertainty is not propagated. Representative median ages, in group order, are',paste(sprintf('%.1f',medians),collapse=', '),'days.'))
  pf_path <- file.path(notes,if(figure3_only)'display_preflight_row_alignment.json'else 'display_preflight.json')
  if(phase %in% c('preflight','preflight-figure3')) {
    check(!file.exists(pf_path),'Preflight already exists; preserve history before an authorized replacement')
    if(figure3_only){
      # Preflight needs font metrics, but must not open the default Rplots.pdf device.
      ragg::agg_capture(width=230,height=230,units='mm',res=96)
      memory_device <- grDevices::dev.cur();gtable_warnings <- character()
      tryCatch(withCallingHandlers({for(k in names(target_names))association_panel(display[display$contrast==k,],
        target_names[[k]],joint[joint$test_id==paste0(k,'_age2df'),],axis)},
        warning=function(w){gtable_warnings<<-c(gtable_warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
        finally=grDevices::dev.off(memory_device))
      check(length(gtable_warnings)==0L,'Association gtable preflight emitted warnings')
    }
    for(n in if(figure3_only)'Prediction_D3_age'else names(captions))writeLines(captions[[n]],file.path(notes,paste0(n,'_caption.txt')),useBytes=TRUE)
    unchanged();writej(list(status='DISPLAY_PREFLIGHT_COMPLETE_NO_RENDER',bindings=identity,OR_axis=axis,
      calibration_bins=30,association_points=12,representative_age_days=as.list(medians),
      source_row_counts=as.list(vapply(list(roc=roc,performance=performance,calibration=bins,association=display,joint=joint),nrow,0L)),
      protected_files_unchanged=length(pm),no_patient_read=TRUE,no_model_fit=TRUE,no_MI=TRUE,no_render=TRUE),pf_path)
    cat('Display preflight completed; independent display release required.\n');return(invisible(NULL))
  }
  pf <- fromJSON(pf_path,simplifyVector=FALSE)
  check(scalar_equal(pf$status,'DISPLAY_PREFLIGHT_COMPLETE_NO_RENDER'),'Invalid display preflight')
  check_release(pf$bindings,identity)
  expected <- identity;expected$preflight_sha256 <- sha(pf_path)
  release_path <- file.path(root,'reviews/display_code_release.json');check_release(fromJSON(release_path,simplifyVector=FALSE),expected)
  if(dir.exists(delivery)&&!figure3_only)check(length(list.files(delivery,all.files=TRUE,no..=TRUE))==0L,'Delivery folder not empty; do not overwrite')
  dir.create(delivery,recursive=TRUE,showWarnings=FALSE)
  if(figure3_only)plots <- list(Prediction_D3_age=wrap_plots(lapply(names(target_names),function(k)wrap_elements(full=association_panel(
      display[display$contrast==k,],target_names[[k]],joint[joint$test_id==paste0(k,'_age2df'),],axis))),ncol=1))
  else plots <- list(Prediction_ROC=roc_plot(roc,performance),
    Prediction_calibration=wrap_plots(lapply(models,function(m)wrap_elements(full=calibration_panel(bins,performance,m))),nrow=1),
    Prediction_D3_age=wrap_plots(lapply(names(target_names),function(k)wrap_elements(full=association_panel(
      display[display$contrast==k,],target_names[[k]],joint[joint$test_id==paste0(k,'_age2df'),],axis))),ncol=1))
  sizes <- list(Prediction_ROC=c(180,180),Prediction_calibration=c(210,115),Prediction_D3_age=c(230,230))
  outputs <- list()
  for(n in names(plots)) {
    for(ext in c('png','pdf')){
      path<-file.path(delivery,paste0(n,'.',ext))
      if(figure3_only)check(scalar_equal(sha(path),sha(file.path(history,'交付审核图形',basename(path)))),'Figure3 original differs from preserved attempt')
      else check(!file.exists(path),'Output already exists')
      if(ext=='png')ggsave(path,plots[[n]],device=ragg::agg_png,width=sizes[[n]][1],height=sizes[[n]][2],units='mm',dpi=300,bg='white')
      else ggsave(path,plots[[n]],device=grDevices::cairo_pdf,width=sizes[[n]][1],height=sizes[[n]][2],units='mm',family='Arial',bg='white')
      outputs[[basename(path)]]<-sha(path)
    }
    cp<-file.path(delivery,paste0(n,'_caption.txt'));writeLines(captions[[n]],cp,useBytes=TRUE);outputs[[basename(cp)]]<-sha(cp)
  }
  if(!figure3_only)dir.create(file.path(delivery,'source_data'),showWarnings=FALSE)
  if(!figure3_only)for(n in names(rel)){dest<-file.path(delivery,'source_data',rel[[n]]);check(file.copy(paths[[n]],dest,overwrite=FALSE),'Source copy failed');
    check(scalar_equal(sha(dest),csv_hashes[[rel[[n]]]]),'Copied source hash mismatch');outputs[[paste0('source_data/',rel[[n]])]]<-sha(dest)}
  unchanged();writej(list(status='FIGURES_RENDERED_PENDING_ACTUAL_VISUAL_REVIEW',bindings=expected,
    display_release_sha256=sha(release_path),output_hashes=outputs,protected_files_unchanged=length(pm),
    no_model_fit=TRUE,no_MI=TRUE,no_patient_read=TRUE,producer_exit_status=0,
    untouched_ROC_calibration_hashes=as.list(setNames(keep_hashes,basename(keep_paths))),
    sessionInfo=capture.output(sessionInfo())),file.path(notes,if(figure3_only)'render_row_alignment_receipt.json'else 'render_receipt.json'))
  cat(if(figure3_only)'Figure3 only rendered; actual row-alignment review pending.\n'else 'Three real figures rendered; actual PNG/PDF visual review pending.\n')
}
if(sys.nframe()==0L)main()
