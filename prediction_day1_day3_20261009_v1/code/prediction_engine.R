# Frozen-plan consumer. source() defines functions only; CLI ROOT preflight|models.
p_labs <- c('albumin','creatinine','derived_BUN_mg_dL','lactate','glucose','ph')
p_preds <- c('D1','D3','D1D3')
p_assocs <- c('day3_common','day3_age')
p_check <- function(x,msg) if(!isTRUE(x))stop(msg,call.=FALSE)
p_context <- function(root) {
  root<-normalizePath(root,winslash='/',mustWork=TRUE)
  .libPaths(c(file.path(dirname(dirname(root)),'landmark72_reanalysis_20261006/runtime/R-library'),.libPaths()))
  cfg<-jsonlite::read_json(file.path(root,'control/ANALYSIS_CONFIG.json'),simplifyVector=TRUE)
  cfg$base <- Sys.getenv("PIC_ASSOCIATION_ROOT")
  p_check(nzchar(cfg$base),"Set PIC_ASSOCIATION_ROOT to the authorized association workspace")
  source(file.path(cfg$base,'code/outcome_engine.R'))
  base<-context(cfg$base)
  paths<-c(code=file.path(root,'code/prediction_engine.R'),plan=file.path(root,'control/ANALYSIS_PLAN.md'),
    config=file.path(root,'control/ANALYSIS_CONFIG.json'),bindings=file.path(root,'control/SOURCE_BINDINGS.json'),
    protected=file.path(root,'control/PROTECTED_MANIFEST.json'),plan_review=file.path(root,'reviews/plan_review.json'),
    setNames(base$paths,paste0('BASE_',names(base$paths))))
  p_check(all(file.exists(paths)),'Missing new source or plan release')
  bindings<-jsonlite::read_json(paths[['bindings']],simplifyVector=TRUE)
  p_check(identical(vapply(bindings$inputs,function(v)v$sha256,character(1)),base$input_hashes),
    'New binding differs from 56 accepted inputs')
  protect<-unlist(jsonlite::read_json(paths[['protected']],simplifyVector=TRUE)$protected_files,use.names=TRUE)
  p_check(is.character(protect)&&length(protect)>0&&!anyNA(protect)&&
    identical(vapply(names(protect),sha,character(1)),protect),'Protected source mismatch')
  review<-jsonlite::read_json(paths[['plan_review']],simplifyVector=TRUE)
  p_check(isTRUE(review$status=='ACCEPT')&&isTRUE(review$plan_sha256==sha(paths[['plan']]))&&
    isTRUE(review$config_sha256==sha(paths[['config']])),'Actual plan approval/hash absent')
  p_check(cfg$cohort_N==3459&&cfg$deaths==306&&cfg$association_M==50&&cfg$folds==5&&cfg$repeats==5&&
    cfg$prediction_M==10&&cfg$prediction_maxit==10&&cfg$pmm_donors==5&&cfg$bootstrap_B==2000&&
    cfg$seed==20261009&&cfg$bootstrap_seed==20261029&&
    identical(cfg$primary_prediction,'mean_patient_OOF_probability_over_repeats')&&
    identical(cfg$prediction_covariates,p_labs)&&identical(cfg$prediction_models,p_preds)&&
    identical(cfg$association_models,p_assocs)&&identical(cfg$original_config,base$cfg),'Unexpected frozen configuration')
  list(root=root,cfg=cfg,base=base,paths=paths,hashes=vapply(paths,sha,character(1)),
    protected=protect,input_hashes=base$input_hashes,N=base$N,events=base$events)
}
p_assert <- function(c) {
  assert_unchanged(c$base)
  p_check(identical(vapply(c$paths,sha,character(1)),c$hashes),'New source changed')
  p_check(identical(vapply(names(c$protected),sha,character(1)),c$protected),'Protected source changed')
}
p_setup <- function(d,cfg) {
  p_check(all(complete.cases(d[p_labs])),'Model clinical values incomplete')
  d$age_layer<-factor(d$age_layer,levels=paste0('L0',1:5))
  d$sex<-factor(d$sex,levels=c('F','M'))
  d$first_careunit<-factor(d$first_careunit,levels=c('PICU','CICU','General ICU','NICU','SICU'))
  a<-rcs3(log1p(d$age_days),cfg$age_log_knots,log1p(cfg$reference_age_days));d$a1<-a[,1];d$a2<-a[,2]
  d$z_albumin<-d$albumin/5
  for(v in c('creatinine','lactate','derived_BUN_mg_dL'))p_check(all(is.finite(d[[v]]))&&all(d[[v]]>=0),'Invalid clinical log input')
  d$z_creatinine<-log1p(d$creatinine);d$z_BUN<-log1p(d$derived_BUN_mg_dL);d$z_lactate<-log1p(d$lactate)
  for(v in c('glucose','ph')) {k<-cfg[[paste0(v,'_knots')]];z<-rcs3(d[[v]],k,k[2],if(v=='glucose')5 else .1)
    for(j in 1:2)d[[paste0('z_',v,j)]]<-z[,j]}
  for(day in c(1,3)) {z<-rcs3(d[[paste0('D',day)]],cfg$sodium_knots,cfg$reference_sodium,cfg$sodium_scale)
    for(j in 1:2){nm<-paste0('n',day,j);d[[nm]]<-z[,j];for(a in 1:2)d[[paste0(nm,'_a',a)]]<-z[,j]*d[[paste0('a',a)]]}}
  d$D3status<-factor(ifelse(d$D3<135,'L',ifelse(d$D3>145,'H','N')),levels=c('N','L','H'))
  contrasts(d$D3status)<-contr.treatment(c('N','L','H'),base=1)
  for(level in c('L','H'))for(a in 1:2)d[[paste0('D3status',level,'_a',a)]]<-(d$D3status==level)*d[[paste0('a',a)]]
  d
}
p_na_terms <- function(day)c(paste0('n',day,1:2),as.vector(t(outer(paste0('n',day,1:2),paste0('_a',1:2),paste0))))
p_formula <- function(id) {
  terms<-clinical_terms
  if(id%in%p_preds){for(day in if(id=='D1')1 else if(id=='D3')3 else c(1,3))terms<-c(terms,p_na_terms(day))}
  else {p_check(id%in%p_assocs,'Unknown model');terms<-c(terms,p_na_terms(1),'D3status')
    if(id=='day3_age')terms<-c(terms,as.vector(t(outer(paste0('D3status',c('L','H')),paste0('_a',1:2),paste0))))}
  reformulate(terms,response='death')
}
p_expected <- function(id)switch(id,D1=26L,D3=26L,D1D3=32L,day3_common=28L,day3_age=32L)
p_folds <- function(raw,cfg) {
  rows<-list();strata<-interaction(raw$age_layer,raw$death,drop=TRUE)
  for(r in 1:cfg$repeats){set.seed(cfg$seed+r);f<-integer(nrow(raw))
    for(s in levels(strata)){ii<-which(strata==s);f[ii]<-sample(rep(1:cfg$folds,length.out=length(ii)))}
    rows[[r]]<-data.frame(row=seq_len(nrow(raw)),id=raw$id,repeat_id=r,fold=f)}
  do.call(rbind,rows)
}
p_imputer_data <- function(d) {
  z<-d[p_labs];z$log_age<-log1p(d$age_days)
  z$age_layer<-factor(d$age_layer,levels=paste0('L0',1:5));z$sex<-factor(d$sex,levels=c('F','M'))
  z$first_careunit<-factor(d$first_careunit,levels=c('PICU','CICU','General ICU','NICU','SICU'))
  z
}
p_impute <- function(raw,train,cfg,seed,audit=FALSE) {
  x<-p_imputer_data(raw);ignore<-!train
  p_check(is.logical(train)&&length(train)==nrow(x)&&all(colSums(!is.na(x[train,p_labs]))>=cfg$pmm_donors),
    'Insufficient original training donors')
  p_check(!any(c('death','D1','D2','D3','urea','hb','wbc','platelet')%in%names(x)),'Forbidden imputer feature')
  set.seed(seed);init<-x
  for(v in p_labs){ii<-is.na(x[[v]]);donors<-which(train&!ii);init[[v]][ii]<-sample(x[[v]][donors],sum(ii),replace=TRUE)}
  pred<-matrix(0L,ncol(x),ncol(x),dimnames=list(names(x),names(x)))
  pred[p_labs,]<-1L;diag(pred)<-0L
  method<-setNames(rep('',ncol(x)),names(x));method[p_labs]<-if(audit)'audit_prediction_pmm' else 'pmm'
  if(audit){.p_audit_env<<-new.env(parent=emptyenv());.p_audit_env$ignore<-ignore;.p_audit_env$calls<-0L
    .p_audit_env$train_only<-TRUE;.p_audit_env$training_donors_only<-TRUE}
  elapsed<-system.time(imp<-mice::mice(x,m=cfg$prediction_M,maxit=cfg$prediction_maxit,donors=cfg$pmm_donors,
    predictorMatrix=pred,method=method,visitSequence=p_labs,ignore=ignore,where=is.na(x),data.init=init,
    remove.constant=FALSE,remove.collinear=FALSE,seed=seed+1L,printFlag=FALSE))['elapsed']
  p_check(identical(imp$ignore,ignore),'Imputer ignore rows changed')
  p_check(isTRUE(all.equal(imp$predictorMatrix,pred,check.attributes=TRUE)),
    'Fixed predictor matrix values or column identities changed')
  ds<-lapply(seq_len(cfg$prediction_M),function(j){z<-mice::complete(imp,j);d<-raw
    p_check(all(complete.cases(z[p_labs])),'Clinical imputation incomplete')
    for(v in p_labs){obs<-!is.na(raw[[v]]);p_check(identical(z[[v]][obs],raw[[v]][obs]),'Observed clinical value changed');d[[v]]<-z[[v]]}
    p_setup(d,cfg$original_config)})
  list(datasets=ds,imp=imp,elapsed=unname(elapsed),seed=seed,initialization_seed=seed,
    sampler_seed=seed+1L,training_rows=which(train),validation_rows=which(ignore),initialization=init,
    audit=if(audit)as.list(.p_audit_env) else NULL)
}
mice.impute.audit_prediction_pmm <- function(y,ry,x,wy=NULL,...) {
  .p_audit_env$calls<-.p_audit_env$calls+1L
  .p_audit_env$train_only<-.p_audit_env$train_only&&!any(ry[.p_audit_env$ignore])
  out<-mice::mice.impute.pmm(y,ry,x,wy=wy,...)
  .p_audit_env$training_donors_only<-.p_audit_env$training_donors_only&&all(out%in%y[ry]);out
}
p_synthetic <- function(c) {
  set.seed(123818);n<-3459L
  x<-data.frame(id=sprintf('SYN%05d',seq_len(n)),age_days=sample(1:5000,n,TRUE),
    age_layer=sample(paste0('L0',1:5),n,TRUE),sex=sample(c('F','M'),n,TRUE),
    first_careunit=sample(c('PICU','CICU','General ICU','NICU','SICU'),n,TRUE),
    albumin=rnorm(n,38,5),creatinine=exp(rnorm(n,3.4,.4)),derived_BUN_mg_dL=exp(rnorm(n,2.6,.4)),
    lactate=exp(rnorm(n,.5,.3)),glucose=exp(rnorm(n,1.8,.3)),ph=rnorm(n,7.37,.06),
    D1=rnorm(n,138,5),D3=rnorm(n,139,5),death=rbinom(n,1,.09))
  for(i in seq_along(p_labs))x[[p_labs[i]]][sample(n,c(781,534,534,3,32,1)[i])]<-NA_real_
  train<-seq_len(n)%%5!=0
  z<-p_impute(x,train,c$cfg,20261200L,audit=TRUE)
  p_check(isTRUE(z$audit$train_only)&&isTRUE(z$audit$training_donors_only)&&z$audit$calls>0,
    'Synthetic training/validation isolation failed')
  q<-setNames(c(.1,.2,.3),c('(Intercept)','x','z'));u<-diag(c(.2,.1,.1));dimnames(u)<-list(names(q),names(q))
  fit<-structure(list(q=q,u=u,formula='death ~ x + z'),class='prediction_numeric_fit')
  pp<-pool_model(rep(list(fit),50),400)
  p_check(identical(unname(pp$Q[,1]),unname(q))&&max(abs(pp$W-u))<1e-14&&max(abs(pp$B))<1e-14,
    'Numeric fit accessors do not reproduce full covariance pooling')
  sc<-scalar(pp,setNames(c(0,1,0),names(q)))
  p_check(all(is.finite(as.matrix(sc))),'Finite-df synthetic scalar failed')
  au<-p_auc(c(0,1,0,1),c(.1,.9,.3,.8));p_check(au==1,'AUC orientation/ties implementation failed')
  list(status='PASS',N=n,M=c$cfg$prediction_M,maxit=c$cfg$prediction_maxit,elapsed_seconds=z$elapsed,
    estimated_25_fold_imputation_seconds=25*z$elapsed,train_only=z$audit$train_only,
    training_donors_only=z$audit$training_donors_only,PMM_calls=z$audit$calls,
    predictor_names=names(z$imp$data),loggedEvents=z$imp$loggedEvents,seed=z$seed,
    no_real_outcome_fit=TRUE,no_environment_in_numeric_fit=TRUE,synthetic_scalar=sc)
}
coef.prediction_numeric_fit <- function(object,...)object$q
vcov.prediction_numeric_fit <- function(object,...)object$u
formula.prediction_numeric_fit <- function(x,...)as.formula(paste(x$formula,collapse=' '),env=baseenv())
p_write_design <- function(X,y,folder,j) {
  bin<-sprintf('design_%03d.bin',j);con<-file(file.path(folder,bin),'wb')
  tryCatch(writeBin(as.double(c(y,as.vector(t(X)))),con,endian='little'),finally=close(con))
  list(imputation=j,path=bin,rows=nrow(X),columns=ncol(X),sha256=sha(file.path(folder,bin)))
}
p_lp <- function(c,folder,meta) {
  writej(list(model=basename(folder),designs=meta),file.path(folder,'design_metadata.json'))
  exit<-system2('pwsh',c('-NoProfile','-File',shQuote(c$base$paths[['python_helper']]),
    shQuote(c$base$paths[['separation']]),shQuote(folder)),stdout=file.path(folder,'separation_stdout.log'),stderr=file.path(folder,'separation_stderr.log'))
  receipt<-if(file.exists(file.path(folder,'separation_validation.json')))jsonlite::read_json(file.path(folder,'separation_validation.json'),simplifyVector=TRUE) else NULL
  list(exit=exit,all_pass=identical(exit,0L)&&isTRUE(receipt$all_pass),receipt=receipt)
}
p_age_C <- function(nms,level,ag=c(0,0),joint=FALSE) {
  if(joint){out<-matrix(0,2,length(nms),dimnames=list(NULL,nms));for(a in 1:2)out[a,paste0('D3status',level,'_a',a)]<-1;return(out)}
  out<-setNames(rep(0,length(nms)),nms);out[paste0('D3status',level)]<-1
  for(a in 1:2){nm<-paste0('D3status',level,'_a',a);if(nm%in%nms)out[nm]<-ag[a]};out
}
p_preflight <- function(c) {
  out<-file.path(c$root,'results/preflight');private<-file.path(out,'restricted')
  p_check(!dir.exists(private)&&!file.exists(file.path(out,'preflight.json')),'Preflight previously attempted; preserve evidence')
  dir.create(private,recursive=TRUE);writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  synth<-p_synthetic(c);writej(synth,file.path(out,'synthetic_validation.json'))
  folds<-p_folds(c$base$raw,c$cfg);saveRDS(folds,file.path(private,'fold_indices.rds'))
  fc<-merge(folds,transform(c$base$raw[c('id','death','age_layer')],row=seq_len(c$N)),by=c('id','row'))
  write.csv(aggregate(rep(1,nrow(fc)),fc[c('repeat_id','fold','age_layer','death')],sum),file.path(out,'fold_counts.csv'),row.names=FALSE)
  rows<-list();meta<-list();for(key in c('association','prediction_probe')){meta[[key]]<-list();dir.create(file.path(private,key))}
  for(j in 1:50){d<-load_completed(c$base,j);z<-p_setup(d,c$cfg$original_config)
    legacy<-setup(d,c$base,1L);cols<-setdiff(clinical_terms,c('age_layer','sex','first_careunit'))
    p_check(identical(z[cols],legacy[cols]),'New six-clinical transformations differ from BASE')
    for(id in p_assocs){X<-model.matrix(p_formula(id),z);k<-length(meta$association)+1L
      rows[[length(rows)+1L]]<-data.frame(scope='association',model_id=id,task=j,N=nrow(X),parameters=ncol(X),rank=qr(X)$rank)
      meta$association[[k]]<-p_write_design(X,z$death,file.path(private,'association'),k)}}
  for(r in 1:5)for(f in 1:5){train<-folds$fold[folds$repeat_id==r]!=f;d<-c$base$raw
    # This is a rank/LP placeholder, never the production imputer or a prediction.
    for(v in p_labs)d[[v]][is.na(d[[v]])]<-median(d[[v]][train],na.rm=TRUE)
    z<-p_setup(d,c$cfg$original_config)
    for(id in p_preds){X<-model.matrix(p_formula(id),z[train,]);k<-length(meta$prediction_probe)+1L
      rows[[length(rows)+1L]]<-data.frame(scope='prediction_train_median_probe',model_id=id,task=(r-1)*5+f,N=nrow(X),parameters=ncol(X),rank=qr(X)$rank)
      meta$prediction_probe[[k]]<-p_write_design(X,z$death[train],file.path(private,'prediction_probe'),k)}}
  checks<-do.call(rbind,rows);expected<-vapply(checks$model_id,p_expected,integer(1))
  fullrank<-all(checks$rank==checks$parameters)&&all(checks$parameters==expected)
  write.csv(checks,file.path(out,'design_checks.csv'),row.names=FALSE)
  lp<-lapply(names(meta),function(key)p_lp(c,file.path(private,key),meta[[key]]));names(lp)<-names(meta)
  p_assert(c)
  writej(list(status='NO_OUTCOME_FIT_PREFLIGHT_COMPLETE_PENDING_RELEASE',no_real_outcome_fit=TRUE,
    N=c$N,events=c$events,association_designs=100L,prediction_probe_designs=75L,all_designs_full_rank=fullrank,
    all_LP_pass=all(vapply(lp,function(x)x$all_pass,logical(1))),LP=lp,synthetic_validation=synth,
    prediction_probe_qualification='Training-only median placeholders establish design diagnostics, not actual fold imputation/predictions. Runtime checks all750 completed designs.',
    source_paths=as.list(c$paths),source_hashes=as.list(c$hashes),input_hashes=as.list(c$input_hashes),
    protected_files=length(c$protected),protected_unchanged=TRUE,fold_indices_sha256=sha(file.path(private,'fold_indices.rds')),
    models_command=p_models_command(c),sessionInfo_sha256=sha(file.path(out,'sessionInfo.txt')),
    package_versions=as.list(c$base$package_versions),utc=utc()),file.path(out,'preflight.json'))
  cat('NO_OUTCOME_FIT_PREFLIGHT_COMPLETE_PENDING_RELEASE\n')
}
p_models_command <- function(c)paste0("& '",c$base$paths[['R_helper']],"' '",c$paths[['code']],"' '",c$root,"' models")
p_gate <- function(c) {
  release_path<-file.path(c$root,'reviews/code_release.json');pre_path<-file.path(c$root,'results/preflight/preflight.json')
  p_check(all(file.exists(c(release_path,pre_path))),'Code release/preflight absent')
  r<-jsonlite::read_json(release_path,simplifyVector=TRUE);pre<-jsonlite::read_json(pre_path,simplifyVector=TRUE)
  p_check(isTRUE(r$status=='ACCEPT')&&isTRUE(r$model_execution_allowed)&&same_hashes(r$source_hashes,c$hashes)&&
    same_hashes(r$input_hashes,c$input_hashes)&&isTRUE(r$preflight_sha256==sha(pre_path))&&
    isTRUE(r$plan_review_sha256==sha(c$paths[['plan_review']])), 'Actual code release/hash mismatch')
  p_check(isTRUE(pre$no_real_outcome_fit)&&isTRUE(pre$all_designs_full_rank)&&isTRUE(pre$all_LP_pass)&&
    same_hashes(pre$source_hashes,c$hashes)&&same_hashes(pre$input_hashes,c$input_hashes)&&
    same_hashes(pre$package_versions,c$base$package_versions),'No-fit preflight not fully passed or bound')
  reg_path<-file.path(c$root,'control/HARNESS_NATIVE_REGISTRATION.json')
  p_check(file.exists(reg_path)&&isTRUE(r$harness_registration_sha256==sha(reg_path)),'Native registration not release-bound')
  reg<-jsonlite::read_json(reg_path,simplifyVector=FALSE);p_check(is.list(reg$runs)&&length(reg$runs)==1L,'Exactly one registered run required')
  entry<-reg$runs[[1]];resolve<-function(p){p_check(is.character(p)&&length(p)==1L&&nzchar(p),'Native path missing')
    if(!grepl('^[A-Za-z]:|^/',p))p<-file.path(dirname(c$root),p);normalizePath(p,winslash='/',mustWork=TRUE)}
  run_path<-resolve(entry$path);native<-jsonlite::read_json(run_path,simplifyVector=FALSE)
  p_check(isTRUE(native$run_id==entry$run_id)&&isTRUE(native$input_snapshot_sha256==entry$input_snapshot_sha256)&&
    is.character(native$input_snapshot_sha256)&&nchar(native$input_snapshot_sha256)==64L&&length(native$input_snapshot_errors)==0L,
    'Native input snapshot invalid')
  p_check(isTRUE(native$command==p_models_command(c))&&isTRUE(entry$command==native$command)&&
    isTRUE(pre$models_command==native$command),'Native command differs from frozen command')
  for(key in c('code','BASE_code'))p_check(any(vapply(native$code_inputs,function(v)isTRUE(resolve(v$path)==c$paths[[key]])&&
    isTRUE(v$sha256==c$hashes[[key]]),logical(1))),'Native code dependency hash absent')
  p_check(isTRUE(resolve(native$configuration_path)==c$paths[['config']])&&isTRUE(native$configuration_sha256==c$hashes[['config']])&&
    isTRUE(resolve(native$environment_path)==normalizePath(file.path(c$root,'results/preflight/sessionInfo.txt'),winslash='/'))&&
    isTRUE(native$environment_sha256==sha(resolve(native$environment_path))),'Native configuration/environment mismatch')
  old<-jsonlite::read_json(file.path(c$cfg$base,'control/HARNESS_NATIVE_REGISTRATION.json'),simplifyVector=FALSE)
  p_check(setequal(unlist(native$data_ids),vapply(old$data,function(v)v$data_id,character(1))),'Registered source identities mismatch')
  list(code_release_sha256=sha(release_path),preflight_sha256=sha(pre_path),native_run_id=entry$run_id,
    native_run_path=run_path,harness_registration_sha256=sha(reg_path),native_run_sha256=sha(run_path))
}
p_fit <- function(d,c,id,j,folder) {
  ap<-file.path(folder,sprintf('attempt_%03d.json',j));aa<-list(model_id=id,imputation=j,status='ATTEMPTED',utc=utc(),source_hashes=as.list(c$hashes));writej(aa,ap)
  warnings<-character();tryCatch({X<-model.matrix(p_formula(id),d);N<-nrow(X);p<-ncol(X)
    p_check(p==p_expected(id)&&qr(X)$rank==p&&all(is.finite(X)),'Actual completed design rank/finite failed')
    meta<-p_write_design(X,d$death,folder,j)
    f<-withCallingHandlers(glm(p_formula(id),data=d,family=binomial(),x=TRUE,y=TRUE,control=glm.control(maxit=100)),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    u<-vcov(f);mineig<-min(eigen(u,symmetric=TRUE,only.values=TRUE)$values)
    p_check(length(warnings)==0&&isTRUE(f$converged)&&!isTRUE(f$boundary)&&f$rank==p&&f$df.residual==N-p&&
      all(is.finite(coef(f)))&&all(is.finite(u))&&mineig>0,'GLM diagnostic failure')
    cook<-cooks.distance(f);hat<-hatvalues(f);dfb<-dfbetas(f);p_check(all(is.finite(c(cook,hat,dfb))),'Nonfinite influence')
    saveRDS(data.frame(id=d$id,cook=cook,leverage=hat,max_abs_DFBETAS=apply(abs(dfb),1,max)),file.path(folder,sprintf('influence_%03d.rds',j)))
    diag<-data.frame(model_id=id,imputation=j,N=N,events=sum(d$death),parameters=p,rank=f$rank,converged=f$converged,
      boundary=f$boundary,iterations=f$iter,warnings='',cov_min_eigen=mineig,condition_number=kappa(X,exact=TRUE),
      probability_min=min(fitted(f)),probability_max=max(fitted(f)),max_cook=max(cook),max_leverage=max(hat),
      max_abs_DFBETAS=max(abs(dfb)),LP_status='PENDING',fit_pass=FALSE)
    compact<-structure(list(q=coef(f),u=u,formula=deparse(formula(f)),design_columns=colnames(X),diagnostic=diag,
      N=N,events=sum(d$death),rank=f$rank,complete_df=N-p,iterations=f$iter,converged=f$converged,boundary=f$boundary,
      deviance=f$deviance,AIC=AIC(f)),class='prediction_numeric_fit')
    saveRDS(compact,file.path(folder,sprintf('numeric_fit_%03d.rds',j)));aa$status<-'GLM_PASSED_LP_PENDING';writej(aa,ap)
    list(fit=compact,diagnostic=diag,meta=meta)
  },error=function(e){aa$status<-'FAILED';aa$error<-conditionMessage(e);aa$warnings<-warnings;writej(aa,ap)
    list(error=conditionMessage(e),diagnostic=data.frame(model_id=id,imputation=j,N=nrow(d),events=sum(d$death),parameters=p_expected(id),rank=NA,
      converged=FALSE,boundary=NA,iterations=NA,warnings=paste(warnings,collapse=' | '),cov_min_eigen=NA,condition_number=NA,
      probability_min=NA,probability_max=NA,max_cook=NA,max_leverage=NA,max_abs_DFBETAS=NA,LP_status='NOT_RUN_GLM_FAILED',fit_pass=FALSE))})
}
p_finish_group <- function(attempts,c,folder) {
  ok<-vapply(attempts,function(x)is.null(x$error),logical(1));lp<-p_lp(c,folder,lapply(attempts[ok],`[[`,'meta'))
  for(j in seq_along(attempts)){row<-attempts[[j]]$diagnostic
    if(ok[j]&&!is.null(lp$receipt)){k<-match(j,lp$receipt$results$imputation)
      if(!is.na(k)){row$LP_status<-lp$receipt$results$status[k];row$fit_pass<-isTRUE(lp$receipt$results$pass[k])&&lp$exit%in%c(0L,2L)}}
    attempts[[j]]$diagnostic<-row
    if(ok[j]){attempts[[j]]$fit$diagnostic<-row;saveRDS(attempts[[j]]$fit,file.path(folder,sprintf('numeric_fit_%03d.rds',j)))}
    ap<-file.path(folder,sprintf('attempt_%03d.json',j));aa<-jsonlite::read_json(ap);aa$LP_exit_code<-lp$exit;aa$LP_status<-row$LP_status
    aa$fit_pass<-row$fit_pass;aa$status<-if(isTRUE(row$fit_pass))'PASSED_ALL_DIAGNOSTICS' else 'FAILED';writej(aa,ap)}
  list(attempts=attempts,all_pass=all(ok)&&isTRUE(lp$all_pass))
}
p_auc <- function(y,p) {
  n1<-sum(y==1);n0<-sum(y==0);p_check(n1>0&&n0>0&&all(is.finite(p)),'AUC score/outcome invalid')
  (sum(rank(p,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)
}
p_score <- function(y,p,calibration=TRUE) {
  p_check(all(is.finite(p))&&all(p>0&p<1),'Risk scores outside open probability interval')
  out<-data.frame(AUC=p_auc(y,p),Brier=mean((y-p)^2))
  if(calibration){lp<-qlogis(p);a<-glm(y~1+offset(lp),family=binomial());b<-glm(y~lp,family=binomial())
    p_check(isTRUE(a$converged)&&isTRUE(b$converged)&&all(is.finite(c(coef(a),coef(b)))),'Calibration evaluation not estimable')
    out$calibration_intercept<-unname(coef(a)[1]);out$calibration_slope<-unname(coef(b)[2])};out
}
p_curves <- function(y,p,id) {
  o<-order(p,decreasing=TRUE);pp<-p[o];yy<-y[o];ends<-c(which(diff(pp)!=0),length(pp));tp<-cumsum(yy)[ends];fp<-cumsum(1-yy)[ends]
  data.frame(model_id=id,threshold=c(Inf,pp[ends]),false_positive_rate=c(0,fp/sum(y==0)),
    sensitivity=c(0,tp/sum(y==1)),specificity=c(1,1-fp/sum(y==0)))
}
p_bins <- function(y,p,id) {
  bins<-ceiling(rank(p,ties.method='first')*10/length(p));rows<-list()
  for(b in 1:10){ii<-bins==b;n<-sum(ii);k<-sum(y[ii]);phat<-k/n;z<-qnorm(.975);den<-1+z*z/n
    cen<-(phat+z*z/(2*n))/den;half<-z*sqrt(phat*(1-phat)/n+z*z/(4*n*n))/den
    rows[[b]]<-data.frame(model_id=id,bin=b,mean_pred=mean(p[ii]),N=n,deaths=k,observed_rate=phat,ci_low=cen-half,ci_high=cen+half)}
  do.call(rbind,rows)
}
p_joint <- function(a,C,id,test_id,role) {
  k<-nrow(C);p_check(ncol(C)==a$p&&qr(t(C))$rank==k,'Association constraints rank mismatch')
  q<-C%*%a$Q;u<-array(unlist(lapply(a$u,function(v)C%*%v%*%t(C))),dim=c(k,k,a$M))
  rownames(q)<-paste0('b',1:k);dimnames(u)<-list(rownames(q),rownames(q),NULL)
  z<-mitml::testConstraints(qhat=q,uhat=u,constraints=rownames(q),method='D1',df.com=a$N-a$p)$test
  p_check(all(is.finite(z[1,c('F.value','df1','df2','P(>F)')])),'D1 joint test nonfinite')
  data.frame(model_id=id,test_id=test_id,role=role,method='Rubin_full_covariance_D1_F',statistic=z[1,'F.value'],
    df1=z[1,'df1'],df2=z[1,'df2'],p=z[1,'P(>F)'],p_holm=NA_real_,N=a$N,events=a$events,parameters=a$p)
}
p_models <- function(c) {
  gate<-p_gate(c);out<-file.path(c$root,'results');private<-file.path(out,'restricted_models')
  p_check(!dir.exists(private),'Actual models already attempted; no restart')
  dir.create(private,recursive=TRUE);writej(c(list(status='ATTEMPTING_FIXED_850_MODELS',utc=utc(),source_hashes=as.list(c$hashes)),gate),file.path(out,'models_started.json'))
  pre<-jsonlite::read_json(file.path(out,'preflight/preflight.json'),simplifyVector=TRUE)
  foldpath<-file.path(out,'preflight/restricted/fold_indices.rds');p_check(isTRUE(sha(foldpath)==pre$fold_indices_sha256),'Fold indices changed')
  folds<-readRDS(foldpath);raw<-c$base$raw;OOF<-array(NA_real_,dim=c(c$N,5,3),dimnames=list(raw$id,paste0('R',1:5),p_preds))
  diagnostics<-list();failures<-character();fit_lists<-list();imputation_rows<-list()
  for(r in 1:5)for(f in 1:5){train<-folds$fold[folds$repeat_id==r]!=f;key<-sprintf('R%d_F%d',r,f);folder<-file.path(private,key);dir.create(folder)
    seed<-c$cfg$seed+10000L*r+100L*f
    mi<-tryCatch(p_impute(raw,train,c$cfg,seed),error=function(e)e)
    if(inherits(mi,'error')){writej(list(status='FAILED',error=conditionMessage(mi),seed=seed),file.path(folder,'imputation_receipt.json'))
      failures<-c(failures,key);next}
    saveRDS(mi$imp,file.path(folder,'prediction_mids.rds'));saveRDS(mi$initialization,file.path(folder,'training_initialization.rds'))
    writej(list(status='COMPLETE',seed=seed,initialization_seed=mi$initialization_seed,sampler_seed=mi$sampler_seed,
      m=10,maxit=10,donors=5,training_N=sum(train),validation_N=sum(!train),elapsed=mi$elapsed,
      predictorMatrix=mi$imp$predictorMatrix,methods=mi$imp$method,loggedEvents=mi$imp$loggedEvents,
      qualification='Whole-data constant/collinear selection disabled; iterative training-only lindep retained'),file.path(folder,'imputation_receipt.json'))
    imputation_rows[[key]]<-data.frame(repeat_id=r,fold=f,seed=seed,elapsed=mi$elapsed,logged_events=if(is.null(mi$imp$loggedEvents))0 else nrow(mi$imp$loggedEvents))
    for(id in p_preds){ff<-file.path(folder,id);dir.create(ff)
      ats<-lapply(1:10,function(j)p_fit(mi$datasets[[j]][train,],c,id,j,ff));done<-p_finish_group(ats,c,ff)
      for(j in 1:10){row<-done$attempts[[j]]$diagnostic;row$scope<-'prediction';row$repeat_id<-r;row$fold<-f;diagnostics[[length(diagnostics)+1L]]<-row}
      fit_lists[[paste(key,id,sep='_')]]<-lapply(done$attempts,`[[`,'fit')
      if(done$all_pass){pred<-sapply(1:10,function(j){X<-model.matrix(p_formula(id),mi$datasets[[j]][!train,]);a<-done$attempts[[j]]$fit
          p_check(identical(colnames(X),names(a$q)),'Validation matrix coefficient order mismatch');plogis(drop(X%*%a$q))})
        OOF[!train,r,id]<-rowMeans(pred)} else failures<-c(failures,paste(key,id,sep='_'))
      writej(list(status=if(done$all_pass)'PASS' else 'FAILED',model_id=id,attempted=10),file.path(ff,'model_receipt.json'))}
    saveRDS(OOF,file.path(private,'OOF.rds'));write.csv(do.call(rbind,diagnostics),file.path(out,'fit_diagnostics.csv'),row.names=FALSE)
    cat('Prediction fold actual attempts',key,'completed\n')}
  pooled<-list();for(id in p_assocs){folder<-file.path(private,id);dir.create(folder)
    ats<-lapply(1:50,function(j)p_fit(p_setup(load_completed(c$base,j),c$cfg$original_config),c,id,j,folder))
    done<-p_finish_group(ats,c,folder)
    for(j in 1:50){row<-done$attempts[[j]]$diagnostic;row$scope<-'association';row$repeat_id<-NA_integer_;row$fold<-NA_integer_;diagnostics[[length(diagnostics)+1L]]<-row}
    fit_lists[[id]]<-lapply(done$attempts,`[[`,'fit')
    if(done$all_pass){a<-pool_model(fit_lists[[id]],c$N);a$events<-c$events;pooled[[id]]<-a} else failures<-c(failures,id)
    cat('Association actual attempts',id,'all_pass=',done$all_pass,'\n')}
  saveRDS(fit_lists,file.path(private,'fits.rds'));saveRDS(pooled,file.path(private,'pooled_models.rds'))
  write.csv(do.call(rbind,diagnostics),file.path(out,'fit_diagnostics.csv'),row.names=FALSE)
  if(length(imputation_rows))write.csv(do.call(rbind,imputation_rows),file.path(out,'imputation_diagnostics.csv'),row.names=FALSE)
  ac<-list();jt<-list();for(id in names(pooled)){a<-pooled[[id]]
    for(level in c('L','H'))for(layer in if(id=='day3_common')'ALL' else paste0('L0',1:5)){
      age<-if(layer=='ALL')NA_real_ else median(raw$age_days[raw$age_layer==layer]);ag<-if(layer=='ALL')c(0,0) else as.numeric(rcs3(log1p(age),c$cfg$original_config$age_log_knots,log1p(c$cfg$original_config$reference_age_days)))
      z<-scalar(a,p_age_C(rownames(a$Q),level,ag));z$p<-2*pt(-abs(z$logOR/z$SE),z$df)
      state<-ifelse(raw$D3<135,'L',ifelse(raw$D3>145,'H','N'));eligible<-if(layer=='ALL')rep(TRUE,c$N) else raw$age_layer==layer
      ac[[length(ac)+1L]]<-cbind(data.frame(model_id=id,contrast=paste0(level,'_vs_N'),age_layer=layer,age_days=age,
        target_N=sum(eligible&state==level),target_deaths=sum(raw$death[eligible&state==level]),reference_N=sum(eligible&state=='N'),
        reference_deaths=sum(raw$death[eligible&state=='N']),N=c$N,events=c$events,parameters=a$p),z)}
    if(id=='day3_age'){for(level in c('L','H'))jt[[length(jt)+1L]]<-p_joint(a,p_age_C(rownames(a$Q),level,joint=TRUE),id,paste0(level,'_vs_N_age2df'),'target_age_Holm2')
      jt[[length(jt)+1L]]<-p_joint(a,rbind(p_age_C(rownames(a$Q),'L',joint=TRUE),p_age_C(rownames(a$Q),'H',joint=TRUE)),id,'all_D3_state_age4df','support')}}
  if(length(ac))write.csv(do.call(rbind,ac),file.path(out,'association_contrasts.csv'),row.names=FALSE,na='')
  if(length(jt)){j<-do.call(rbind,jt);ii<-j$role=='target_age_Holm2';j$p_holm[ii]<-p.adjust(j$p[ii],'holm');write.csv(j,file.path(out,'joint_tests.csv'),row.names=FALSE,na='')}
  if(all(is.finite(OOF))){scores<-apply(OOF,c(1,3),mean);performance<-lapply(p_preds,function(id)cbind(data.frame(model_id=id,N=c$N,events=c$events),p_score(raw$death,scores[,id])))
    performance<-do.call(rbind,performance);rp<-list();for(r in 1:5)for(id in p_preds)rp[[length(rp)+1L]]<-cbind(data.frame(repeat_id=r,model_id=id),p_score(raw$death,OOF[,r,id]))
    write.csv(do.call(rbind,rp),file.path(out,'repeat_performance.csv'),row.names=FALSE)
    set.seed(c$cfg$bootstrap_seed);bi<-array(NA_real_,c(c$cfg$bootstrap_B,2,3),dimnames=list(NULL,c('AUC','Brier'),p_preds))
    for(b in 1:c$cfg$bootstrap_B){idx<-c(sample(which(raw$death==0),sum(raw$death==0),TRUE),sample(which(raw$death==1),sum(raw$death==1),TRUE))
      for(id in p_preds)bi[b,,id]<-as.numeric(p_score(raw$death[idx],scores[idx,id],calibration=FALSE))}
    saveRDS(bi,file.path(private,'conditional_bootstrap.rds'));performance$auc_low<-sapply(p_preds,function(id)quantile(bi[,1,id],.025));performance$auc_high<-sapply(p_preds,function(id)quantile(bi[,1,id],.975))
    performance$brier_low<-sapply(p_preds,function(id)quantile(bi[,2,id],.025));performance$brier_high<-sapply(p_preds,function(id)quantile(bi[,2,id],.975))
    performance$metric_basis<-'mean_5_repeat_OOF_probability';performance$CI_type<-'conditional_fixed_OOF_stratified_patient_bootstrap'
    write.csv(performance,file.path(out,'performance.csv'),row.names=FALSE)
    comparisons<-list();for(id in c('D3','D1D3'))for(metric in c('AUC','Brier')){delta<-bi[,metric,id]-bi[,metric,'D1'];point<-performance[[metric]][match(id,performance$model_id)]-performance[[metric]][match('D1',performance$model_id)]
      comparisons[[length(comparisons)+1L]]<-data.frame(contrast=paste0(id,'_vs_D1'),metric=metric,difference=point,ci_low=quantile(delta,.025),ci_high=quantile(delta,.975),bootstrap_B=2000,CI_type='conditional_fixed_OOF',direction=if(metric=='AUC')'positive_favors_added_or_D3' else 'negative_favors_added_or_D3')}
    write.csv(do.call(rbind,comparisons),file.path(out,'paired_comparisons.csv'),row.names=FALSE)
    write.csv(do.call(rbind,lapply(p_preds,function(id)p_curves(raw$death,scores[,id],id))),file.path(out,'roc_curves.csv'),row.names=FALSE)
    write.csv(do.call(rbind,lapply(p_preds,function(id)p_bins(raw$death,scores[,id],id))),file.path(out,'calibration_bins.csv'),row.names=FALSE)}
  else failures<-c(failures,'INCOMPLETE_OOF_NO_PERFORMANCE_ESTIMATION')
  p_assert(c);files<-c('association_contrasts.csv','joint_tests.csv','performance.csv','paired_comparisons.csv','repeat_performance.csv',
    'roc_curves.csv','calibration_bins.csv','fit_diagnostics.csv','imputation_diagnostics.csv','restricted_models/fits.rds','restricted_models/pooled_models.rds','restricted_models/OOF.rds','restricted_models/conditional_bootstrap.rds')
  files<-files[file.exists(file.path(out,files))]
  writej(c(list(status=if(length(failures))'PARTIAL_FAILED_PENDING_REVIEW' else 'COMPLETED_PENDING_INDEPENDENT_RESULT_REVIEW',
    failures=failures,prediction_folds=length(imputation_rows),association_models_completed=length(pooled),
    fits_attempted=length(diagnostics),fits_passed=sum(vapply(diagnostics,function(v)isTRUE(v$fit_pass),logical(1))),
    expected_attempts=850,N=c$N,events=c$events,source_paths=as.list(c$paths),source_hashes=as.list(c$hashes),input_hashes=as.list(c$input_hashes),
    protected_unchanged=TRUE,output_hashes=as.list(setNames(vapply(file.path(out,files),sha,character(1)),files)),utc=utc(),
    qualification='Exploratory 72h-landmark comparison. Conditional fixed-OOF bootstrap omits training/fold/imputation uncertainty; no confirmatory superiority claim.'),gate),file.path(out,'analysis_receipt.json'))
  p_check(length(failures)==0&&length(diagnostics)==850,'Actual failure/attempt count mismatch; no fallback')
  cat('COMPLETED_PENDING_INDEPENDENT_RESULT_REVIEW\n')
}
p_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  p_check(length(args)==2&&args[2]%in%c('preflight','models'),'Usage prediction_engine.R ROOT preflight|models')
  c<-p_context(args[1]);if(args[2]=='preflight')p_preflight(c) else p_models(c)
}
if(sys.nframe()==0L)tryCatch(p_main(),error=function(e){a<-commandArgs(trailingOnly=TRUE)
  if(length(a)>0&&dir.exists(a[1])&&requireNamespace('jsonlite',quietly=TRUE)){
    out<-file.path(a[1],'results');dir.create(out,showWarnings=FALSE)
    jsonlite::write_json(list(status='FAILED_NOT_PASS',error=conditionMessage(e),utc=format(Sys.time(),tz='UTC',format='%Y-%m-%dT%H:%M:%SZ')),
      file.path(out,paste0('execution_failure_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.json')),auto_unbox=TRUE,pretty=TRUE)}
  message(conditionMessage(e));quit(status=1L)})
