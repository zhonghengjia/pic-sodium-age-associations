# Accepted daily-mean plan consumer. CLI: outcome_engine.R ROOT preflight|models.
# source() defines functions only. No new MI, cleaning, case deletion or fallback.
sha <- function(p) digest::digest(file=p,algo='sha256',serialize=FALSE)
writej <- function(x,p) jsonlite::write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=17,na='null')
utc <- function() format(Sys.time(),tz='UTC',format='%Y-%m-%dT%H:%M:%SZ')
check <- function(ok,message) if(!isTRUE(ok)) stop(message,call.=FALSE)
same_hashes <- function(x,y) {
  x <- unlist(x,use.names=TRUE)
  is.character(x)&&length(x)==length(y)&&!anyNA(x)&&
    !anyDuplicated(names(x))&&setequal(names(x),names(y))&&
    identical(unname(x[names(y)]),unname(y))
}
rcs3 <- function(x,k,ref,scale=1) {
  check(length(k)==3L&&all(is.finite(k))&&all(diff(k)>0)&&
          all(is.finite(x))&&is.finite(ref)&&scale>0,'Invalid fixed RCS input')
  h <- function(v) (pmax(v-k[1],0)^3-(k[3]-k[1])/(k[3]-k[2])*pmax(v-k[2],0)^3+
                     (k[2]-k[1])/(k[3]-k[2])*pmax(v-k[3],0)^3)/(k[3]-k[1])^2
  cbind((x-ref)/scale,(h(x)-h(ref))/scale)
}
clinical_terms <- c('age_layer','a1','a2','sex','first_careunit',
                    'z_albumin','z_creatinine','z_BUN','z_lactate','z_glucose1','z_glucose2','z_ph1','z_ph2')
initial_terms <- c('s1_a1','s1_a2','s2_a1','s2_a2')
dynamic_terms <- c('delta12_a1','delta12_a2','delta23_a1','delta23_a2')
cbc_terms <- c('z_hb','z_wbc','z_platelet1','z_platelet2')
target_labs <- c('albumin','creatinine','urea','lactate','glucose','ph','hb','wbc','platelet')
model_formula <- function(id) {
  terms <- clinical_terms
  if(startsWith(id,'dynamic_')) {
    terms <- c(terms,'s1','s2','delta12','delta23',initial_terms,dynamic_terms)
    if(id=='dynamic_CBC') terms <- c(terms,cbc_terms)
    if(id=='dynamic_frequency') terms <- c(terms,paste0('obs',1:3))
  } else {
    check(grepl('^daily[123]_(overall|age)$',id),'Unexpected model identifier')
    terms <- c(terms,'s1','s2')
    if(endsWith(id,'_age')) terms <- c(terms,initial_terms)
  }
  reformulate(terms,response='death')
}
read_patient <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,
                                    colClasses=c(id='character'))
setup <- function(d,c,day=1L) {
  d$age_layer <- factor(d$age_layer,levels=paste0('L0',1:5))
  d$sex <- factor(d$sex,levels=c('F','M'))
  d$first_careunit <- factor(d$first_careunit,levels=c('PICU','CICU','General ICU','NICU','SICU'))
  a <- rcs3(log1p(d$age_days),c$cfg$age_log_knots,log1p(c$cfg$reference_age_days))
  d$a1 <- a[,1];d$a2 <- a[,2]
  s <- rcs3(d[[paste0('D',day)]],c$cfg$sodium_knots,c$cfg$reference_sodium,c$cfg$sodium_scale)
  d$s1 <- s[,1];d$s2 <- s[,2]
  d$delta12 <- (d$D2-d$D1)/c$cfg$sodium_scale
  d$delta23 <- (d$D3-d$D2)/c$cfg$sodium_scale
  for(v in c('s1','s2','delta12','delta23')) for(j in 1:2)
    d[[paste0(v,'_a',j)]] <- d[[v]]*d[[paste0('a',j)]]
  for(j in 1:3) {
    n <- d[[paste0('nD',j)]]-d[[paste0('D',j,'_reconstructed_n')]]
    check(all(is.finite(n))&&all(n>=0)&&all(n==floor(n)),'Invalid original measurement counts')
    d[[paste0('obs',j)]] <- log1p(n)
  }
  d$z_albumin <- d$albumin/5;d$z_hb <- d$hb/10
  for(v in c('creatinine','lactate','wbc')) {
    check(all(is.finite(d[[v]]))&&all(d[[v]]>=0),'Invalid log-transform input')
    d[[paste0('z_',v)]] <- log1p(d[[v]])
  }
  check(all(is.finite(d$derived_BUN_mg_dL))&&all(d$derived_BUN_mg_dL>=0),'Invalid BUN input')
  d$z_BUN <- log1p(d$derived_BUN_mg_dL)
  for(v in c('glucose','ph','platelet')) {
    k <- if(v=='platelet') c$platelet_knots else c$cfg[[paste0(v,'_knots')]]
    z <- rcs3(d[[v]],k,k[2],if(v=='glucose') 5 else if(v=='ph') .1 else 1)
    d[[paste0('z_',v,'1')]] <- z[,1];d[[paste0('z_',v,'2')]] <- z[,2]
  }
  d
}
hash_inputs <- function(bindings) vapply(bindings$inputs,function(v)sha(v$path),character(1))
context <- function(root) {
  root <- normalizePath(root,winslash='/',mustWork=TRUE)
  study <- dirname(dirname(root))
  .libPaths(c(file.path(study,'landmark72_reanalysis_20261006/runtime/R-library'),.libPaths()))
  packages <- c('digest','jsonlite','mice','mitml')
  check(all(vapply(packages,requireNamespace,logical(1),quietly=TRUE)),'Missing required R package; no installation allowed')
  paths <- c(code=file.path(root,'code/outcome_engine.R'),plan=file.path(root,'control/ANALYSIS_PLAN.md'),
             config=file.path(root,'control/ANALYSIS_CONFIG.json'),bindings=file.path(root,'control/SOURCE_BINDINGS.json'),
             plan_review=file.path(root,'reviews/plan_review_round01.json'),
             separation=file.path(study,'daily_min_analysis_20261008/code/check_separation.py'),
             python_helper=Sys.getenv("PIC_PYTHON_HELPER"),
             R_helper=Sys.getenv("PIC_R_HELPER"))
  check(all(file.exists(paths)),'Missing source, release or dependency')
  hashes <- vapply(paths,sha,character(1))
  review <- jsonlite::read_json(paths['plan_review'],simplifyVector=TRUE)
  check(isTRUE(review$status=='ACCEPT')&&isTRUE(review$plan_sha256==hashes[['plan']])&&
          isTRUE(review$config_sha256==hashes[['config']]),'Independent plan acceptance/hash gate failed')
  for(nm in names(review$source_hashes))
    check(isTRUE(sha(file.path(root,nm))==review$source_hashes[[nm]]),'Plan-review source binding mismatch')
  cfg <- jsonlite::read_json(paths['config'],simplifyVector=TRUE)
  check(cfg$M==50&&cfg$expected_N==3459&&cfg$expected_events==306&&
          isTRUE(cfg$no_new_MI)&&cfg$daily_interaction=='sodium_RCS2_x_continuous_age_RCS2_4df',
        'Unexpected accepted model configuration')
  ids <- c('dynamic_main','dynamic_CBC','dynamic_frequency',paste0('daily',1:3,'_overall'),paste0('daily',1:3,'_age'))
  check(identical(cfg$model_ids,ids),'Unexpected model set/order')
  bindings <- jsonlite::read_json(paths['bindings'],simplifyVector=TRUE)
  check(length(bindings$inputs)==56L,'Expected 56 single-authority input bindings')
  check(all(vapply(bindings$inputs,function(v)file.exists(v$path),logical(1))),'Missing bound input')
  input_hashes <- hash_inputs(bindings)
  expected <- vapply(bindings$inputs,function(v)v$sha256,character(1))
  check(identical(input_hashes,expected),'Actual source input hash mismatch')
  d <- read_patient(bindings$inputs$observed$path)
  check(nrow(d)==cfg$expected_N&&!anyDuplicated(d$id)&&!anyNA(d$id)&&sum(d$death)==cfg$expected_events&&
          all(d$death%in%c(0,1)),'Accepted cohort/count mismatch')
  mi <- readRDS(bindings$inputs$MI_object$path)
  check(mi$imp$m==cfg$M&&mi$imp$iteration==30&&identical(mi$targets,target_labs)&&
          nrow(mi$input)==nrow(d)&&identical(as.character(mi$master$id),d$id),'MI object correspondence mismatch')
  check(isTRUE(unname(mi$source_hashes['input'])==input_hashes[['observed']]),'MI original input hash mismatch')
  release <- jsonlite::read_json(bindings$inputs$MI_release$path,simplifyVector=TRUE)
  check(isTRUE(release$status=='ACCEPT'),'MI data-stage acceptance missing')
  bun <- jsonlite::read_json(bindings$inputs$BUN_binding$path,simplifyVector=TRUE)
  check(isTRUE(bun$decision=='ACCEPT')&&bun$binding$source_column=='derived_BUN_mg_dL', 'BUN binding mismatch')
  for(v in c('sodium','age_log')) {
    key <- paste0(v,'_knots')
    micfg <- if(v=='sodium') 'sodium_knots' else 'age_log_knots'
    check(isTRUE(all.equal(as.numeric(cfg[[key]]),as.numeric(jsonlite::read_json(
      file.path(dirname(bindings$inputs$data_authority$path),'MI_config.json'),simplifyVector=TRUE)[[micfg]]))),
      'Consumer Na/age knots differ from accepted MI')
  }
  pk <- as.numeric(quantile(d$platelet,c(.1,.5,.9),na.rm=TRUE,type=7))
  check(all(is.finite(pk))&&all(diff(pk)>0),'Observed platelet knots not distinct')
  list(root=root,paths=paths,hashes=hashes,cfg=cfg,bindings=bindings,input_hashes=input_hashes,
       raw=d,mi=mi,N=nrow(d),events=sum(d$death),platelet_knots=pk,
       package_versions=setNames(vapply(packages,function(v)as.character(packageVersion(v)),character(1)),packages))
}
load_completed <- function(c,j) {
  key <- sprintf('completed_%02d',j);p <- c$bindings$inputs[[key]]$path
  check(isTRUE(sha(p)==c$input_hashes[[key]]),'Completed input changed during run')
  d <- read_patient(p);raw <- c$raw
  check(identical(names(d),names(raw))&&identical(d$id,raw$id)&&nrow(d)==c$N,'Completed row/column correspondence mismatch')
  allowed <- c(target_labs,'derived_BUN_mg_dL')
  check(identical(d[setdiff(names(raw),allowed)],raw[setdiff(names(raw),allowed)]),'Non-target input changed')
  z <- mice::complete(c$mi$imp,j)
  for(v in target_labs) {
    check(all(is.finite(d[[v]])),'Completed target has nonfinite/missing value')
    obs <- !is.na(raw[[v]])
    check(identical(d[[v]][obs],raw[[v]][obs]),'Observed target changed')
    check(isTRUE(all.equal(d[[v]],z[[v]],tolerance=1e-12)),'Export/object target mismatch')
  }
  check(all(is.finite(d$derived_BUN_mg_dL))&&
          max(abs(d$derived_BUN_mg_dL-d$urea*2.8))<1e-9,'Per-dataset BUN formula mismatch')
  d
}
support_ranges <- function(c) {
  rows <- list();d <- c$raw
  for(j in 1:3)for(layer in c('ALL',paste0('L0',1:5))) {
    z <- if(layer=='ALL') d else d[d$age_layer==layer,]
    x <- z[[paste0('D',j)]];q <- quantile(x,c(.025,.975),type=7)
    rows[[length(rows)+1L]] <- data.frame(day=j,age_layer=layer,N=nrow(z),events=sum(z$death),
      representative_age_days=if(layer=='ALL') NA_real_ else median(z$age_days),
      min=min(x),p025=q[1],p975=q[2],max=max(x),
      local_reference_n=sum(x>=c$cfg$local_reference_support_window[1]&x<=c$cfg$local_reference_support_window[2]),
      reference_in_display_range=c$cfg$reference_sodium>=q[1]&&c$cfg$reference_sodium<=q[2],
      reconstructed_contributions=sum(z[[paste0('D',j,'_reconstructed_n')]]))
  }
  do.call(rbind,rows)
}
assert_unchanged <- function(c) {
  check(identical(vapply(c$paths,sha,character(1)),c$hashes),'Code/control/dependency changed during execution')
  check(identical(hash_inputs(c$bindings),c$input_hashes),'Protected bound inputs changed during execution')
}
models_command <- function(c) paste0("& '",c$paths[['R_helper']],"' '",c$paths[['code']],"' '",c$root,"' models")
preflight <- function(c) {
  out <- file.path(c$root,'results/preflight');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  rows <- list()
  for(j in seq_len(c$cfg$M)) {
    d <- load_completed(c,j)
    for(id in c$cfg$model_ids) {
      day <- if(startsWith(id,'daily')) as.integer(substr(id,6,6)) else 1L
      z <- setup(d,c,day);f <- model_formula(id);X <- model.matrix(f,z)
      check(nrow(X)==c$N&&all(is.finite(X))&&qr(X)$rank==ncol(X),'Preflight design not finite/full rank')
      rows[[length(rows)+1L]] <- data.frame(model_id=id,imputation=j,N=nrow(X),parameters=ncol(X),
        rank=qr(X)$rank,complete_df=c$N-ncol(X),condition_number=kappa(X,exact=TRUE))
    }
    if(j%%10==0) cat('Preflight completed input/design',j,'of',c$cfg$M,'\n')
  }
  designs <- do.call(rbind,rows)
  write.csv(designs,file.path(out,'design_checks.csv'),row.names=FALSE)
  ranges <- support_ranges(c);write.csv(ranges,file.path(out,'support_ranges.csv'),row.names=FALSE)
  probe <- file.path(out,'dependency_probe.py')
  writeLines(c('import sys,json,hashlib,numpy,scipy',
    'from scipy.optimize import linprog',
    'from pathlib import Path',
    'files=[Path(numpy.__file__),Path(scipy.__file__),Path(sys.executable)]',
    'print(json.dumps({"executable":sys.executable,"python":sys.version,"numpy":numpy.__version__,"scipy":scipy.__version__,"files":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}))'),probe)
  status <- system2('pwsh',c('-NoProfile','-File',shQuote(c$paths[['python_helper']]),shQuote(probe)),
                    stdout=file.path(out,'python_probe_stdout.txt'),stderr=file.path(out,'python_probe_stderr.txt'))
  check(identical(status,0L),'LP Python dependency preflight failed')
  py <- jsonlite::fromJSON(paste(readLines(file.path(out,'python_probe_stdout.txt'),warn=FALSE),collapse='\n'))
  check(is.character(py$executable)&&length(py$executable)==1L,'Python dependency receipt invalid')
  assert_unchanged(c)
  p <- list(status='PREFLIGHT_COMPLETE_PENDING_CODE_RELEASE',created_utc=utc(),no_outcome_fit=TRUE,no_MI_run=TRUE,
    N=c$N,events=c$events,M=c$cfg$M,models=c$cfg$model_ids,designs_checked=nrow(designs),
    source_hashes=as.list(c$hashes),input_hashes=as.list(c$input_hashes),R_version=R.version.string,
    models_command=models_command(c),
    R_package_versions=as.list(c$package_versions),python_dependencies=py,python_probe_exit_code=status,
    platelet_knots=c$platelet_knots,consumer_basis=list(sodium=c$cfg$sodium_knots,age=c$cfg$age_log_knots,
      glucose=c$cfg$glucose_knots,ph=c$cfg$ph_knots),design_summary=do.call(rbind,lapply(split(designs,designs$model_id),function(z)
        data.frame(model_id=z$model_id[1],parameters=unique(z$parameters),min_rank=min(z$rank),
                   max_condition=max(z$condition_number)))),
    support_values_include_author_reconstructions=TRUE,
    preflight_files_sha256=as.list(vapply(file.path(out,c('sessionInfo.txt','design_checks.csv','support_ranges.csv',
      'dependency_probe.py','python_probe_stdout.txt','python_probe_stderr.txt')),sha,character(1))),
    limitations=c('Design coverage does not prove SMC or MAR','Raw PMM then clinical log/RCS remains approximate FCS',
      'Accepted tail donors and deterministic reconstruction uncertainty remain'))
  writej(p,file.path(out,'preflight.json'))
  cat('PREFLIGHT_COMPLETE_PENDING_CODE_RELEASE; no outcome model fitted.\n')
}
code_gate <- function(c) {
  pf <- file.path(c$root,'results/preflight/preflight.json')
  rf <- file.path(c$root,'reviews/code_release.json')
  check(file.exists(pf)&&file.exists(rf),'Independent code release/preflight absent; no outcome fit')
  pre <- jsonlite::read_json(pf,simplifyVector=TRUE);r <- jsonlite::read_json(rf,simplifyVector=TRUE)
  check(isTRUE(r$status=='ACCEPT')&&isTRUE(r$preflight_sha256==sha(pf))&&
          same_hashes(r$source_hashes,c$hashes)&&same_hashes(r$input_hashes,c$input_hashes),
        'Independent code release/hash gate failed; no outcome fit')
  check(isTRUE(pre$status=='PREFLIGHT_COMPLETE_PENDING_CODE_RELEASE')&&isTRUE(pre$no_outcome_fit)&&
          same_hashes(pre$source_hashes,c$hashes)&&same_hashes(pre$input_hashes,c$input_hashes),
        'Bound preflight source/input mismatch; no outcome fit')
  check(same_hashes(pre$R_package_versions,c$package_versions),'R package version drift since preflight')
  for(p in names(pre$python_dependencies$files)) check(isTRUE(sha(p)==pre$python_dependencies$files[[p]]),
                                                     'Python dependency drift since preflight')
  regfile <- file.path(c$root,'control/HARNESS_NATIVE_REGISTRATION.json')
  check(file.exists(regfile)&&isTRUE(r$harness_registration_sha256==sha(regfile)),
        'Native harness registration not bound to independent code release')
  reg <- jsonlite::read_json(regfile,simplifyVector=FALSE)
  check(is.list(reg$runs)&&length(reg$runs)==1L,'Exactly one frozen native model run registration required')
  entry <- reg$runs[[1]]
  check(is.character(entry$run_id)&&length(entry$run_id)==1L&&nzchar(entry$run_id), 'Native run_id missing')
  resolve <- function(p) {
    check(is.character(p)&&length(p)==1L&&nzchar(p),'Native run path missing')
    if(!grepl('^[A-Za-z]:|^/',p)) p <- file.path(dirname(c$root),p)
    normalizePath(p,winslash='/',mustWork=TRUE)
  }
  runpath <- resolve(entry$path)
  check(startsWith(runpath,paste0(normalizePath(file.path(dirname(c$root),'.research-harness/runs'),
    winslash='/',mustWork=TRUE),'/')),'Native run path outside parent harness')
  native <- jsonlite::read_json(runpath,simplifyVector=FALSE)
  check(isTRUE(native$run_id==entry$run_id)&&is.character(native$input_snapshot_sha256)&&
          length(native$input_snapshot_sha256)==1L&&nchar(native$input_snapshot_sha256)==64L&&
          isTRUE(native$input_snapshot_sha256==entry$input_snapshot_sha256)&&length(native$input_snapshot_errors)==0L,
        'Native run snapshot not frozen successfully')
  check(isTRUE(native$command==models_command(c))&&isTRUE(entry$command==native$command)&&
          isTRUE(pre$models_command==native$command),'Native models command does not match the frozen invocation')
  check(is.list(native$code_inputs)&&length(native$code_inputs)>=1L,'Native code inputs missing')
  codeok <- vapply(native$code_inputs,function(v)
    isTRUE(resolve(v$path)==c$paths[['code']])&&isTRUE(v$sha256==c$hashes[['code']]),logical(1))
  check(any(codeok),'Native run does not bind current outcome engine')
  check(isTRUE(resolve(native$configuration_path)==c$paths[['config']])&&
          isTRUE(native$configuration_sha256==c$hashes[['config']]),'Native configuration binding mismatch')
  session <- normalizePath(file.path(c$root,'results/preflight/sessionInfo.txt'),winslash='/',mustWork=TRUE)
  check(isTRUE(resolve(native$environment_path)==session)&&isTRUE(native$environment_sha256==sha(session)),
        'Native environment/sessionInfo binding mismatch')
  check(is.list(reg$data)&&length(reg$data)==56L&&
          setequal(vapply(reg$data,function(v)v$sha256,character(1)),unname(c$input_hashes))&&
          setequal(unlist(native$data_ids),vapply(reg$data,function(v)v$data_id,character(1))),
        'Native registered data inputs do not match the 56 frozen sources')
  list(preflight_sha256=sha(pf),code_release_sha256=sha(rf),harness_registration_sha256=sha(regfile),
       native_run_id=entry$run_id,native_run_path=runpath,native_run_sha256=sha(runpath),models_command=native$command)
}
pool_model <- function(fits,N) {
  q <- lapply(fits,coef);u <- lapply(fits,vcov);nms <- names(q[[1]])
  check(all(vapply(q,function(v)identical(names(v),nms),logical(1))),'Coefficient names differ across imputations')
  Q <- do.call(cbind,q);M <- ncol(Q);W <- Reduce('+',u)/M;B <- cov(t(Q));T <- W+(1+1/M)*B
  check(all(is.finite(T))&&min(eigen(T,symmetric=TRUE,only.values=TRUE)$values)>0,'Pooled covariance not positive definite')
  list(q=q,u=u,Q=Q,W=W,B=B,T=T,N=N,p=length(nms),M=M,formula=deparse(formula(fits[[1]])))
}
scalar <- function(a,C) {
  check(identical(names(C),rownames(a$Q)),'Scalar contrast name/order mismatch')
  if(all(C==0)) return(data.frame(logOR=0,OR=1,ci_low=1,ci_high=1,SE=0,df=Inf))
  q <- vapply(a$q,function(z)sum(z*C),numeric(1));u <- vapply(a$u,function(z)drop(t(C)%*%z%*%C),numeric(1))
  v <- mice::pool.scalar(q,u,n=a$N,k=a$p);s <- sqrt(v$t);crit <- qt(.975,v$df)
  ans <- data.frame(logOR=v$qbar,OR=exp(v$qbar),ci_low=exp(v$qbar-crit*s),ci_high=exp(v$qbar+crit*s),SE=s,df=v$df)
  check(all(is.finite(as.matrix(ans)))&&ans$df>0,'Nonfinite scalar contrast; no automatic suppression')
  ans
}
joint <- function(a,terms,id,test_id,role,holm_family='none') {
  ii <- match(terms,rownames(a$Q));check(!anyNA(ii),'Joint constraint name missing')
  q <- a$Q[ii,,drop=FALSE]
  u <- array(unlist(lapply(a$u,function(z)z[ii,ii,drop=FALSE])),dim=c(length(ii),length(ii),a$M))
  rownames(q) <- paste0('b',seq_along(ii));dimnames(u) <- list(rownames(q),rownames(q),NULL)
  tt <- mitml::testConstraints(qhat=q,uhat=u,constraints=rownames(q),method='D1',df.com=a$N-a$p)$test
  check(all(is.finite(tt[1,c('F.value','df1','df2','P(>F)')])),'D1 returned nonfinite results')
  data.frame(model_id=id,test_id=test_id,role=role,df1=tt[1,'df1'],df2=tt[1,'df2'],
    statistic=tt[1,'F.value'],p=tt[1,'P(>F)'],p_holm=NA_real_,holm_family=holm_family,
    N=a$N,events=306L,parameters=a$p,complete_df=a$N-a$p)
}
fit_one <- function(d,c,id,j,folder) {
  day <- if(startsWith(id,'daily')) as.integer(substr(id,6,6)) else 1L
  z <- setup(d,c,day);f <- model_formula(id);X <- model.matrix(f,z)
  warnings <- character();attempt <- list(model_id=id,imputation=j,started_utc=utc(),
    input_sha256=c$input_hashes[[sprintf('completed_%02d',j)]],source_hashes=as.list(c$hashes),status='ATTEMPTED')
  ap <- file.path(folder,sprintf('attempt_%02d.json',j));writej(attempt,ap)
  result <- tryCatch({
    check(nrow(X)==c$N&&qr(X)$rank==ncol(X)&&all(is.finite(X)),'Fit design gate failed')
    fit <- withCallingHandlers(glm(f,data=z,family=binomial(),x=TRUE,y=TRUE,
                                  control=glm.control(maxit=100)),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
    })
    V <- vcov(fit);eig <- min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)
    check(length(warnings)==0&&isTRUE(fit$converged)&&fit$rank==ncol(X)&&
            all(is.finite(coef(fit)))&&all(is.finite(V))&&eig>0&&fit$df.residual==c$N-ncol(X),
          'GLM rank/convergence/warnings/finite/covariance gate failed')
    cook <- cooks.distance(fit);hat <- hatvalues(fit);dfb <- dfbetas(fit)
    check(all(is.finite(cook))&&all(is.finite(hat))&&all(is.finite(dfb)),'Nonfinite influence diagnostic')
    infl <- data.frame(id=d$id,imputation=j,cook=cook,leverage=hat,max_abs_DFBETAS=apply(abs(dfb),1,max))
    saveRDS(infl,file.path(folder,sprintf('influence_%02d.rds',j)))
    bin <- sprintf('design_%02d.bin',j);con <- file(file.path(folder,bin),'wb')
    tryCatch(writeBin(as.double(c(d$death,as.vector(t(X)))),con,endian='little'),finally=close(con))
    diag <- data.frame(model_id=id,imputation=j,N=c$N,events=c$events,parameters=ncol(X),rank=fit$rank,
      complete_df=fit$df.residual,converged=fit$converged,iterations=fit$iter,warnings='',
      cov_min_eigen=eig,condition_number=kappa(X,exact=TRUE),max_abs_coefficient=max(abs(coef(fit))),
      probability_min=min(fitted(fit)),probability_max=max(fitted(fit)),
      probability_below_1e8=sum(fitted(fit)<1e-8),probability_above_1minus1e8=sum(fitted(fit)>1-1e-8),
      max_cook=max(cook),max_leverage=max(hat),max_abs_DFBETAS=max(abs(dfb)),
      cook_screen_count=sum(cook>4/fit$df.residual),leverage_screen_count=sum(hat>2*ncol(X)/c$N),
      dfbetas_screen_count=sum(apply(abs(dfb),1,max)>2/sqrt(c$N)),LP_status='PENDING',fit_pass=FALSE)
    attempt$status <- 'GLM_DIAGNOSTICS_PASSED_LP_PENDING';attempt$ended_utc <- utc()
    writej(attempt,ap)
    list(fit=fit,diagnostic=diag,design=list(imputation=j,path=bin,rows=c$N,columns=ncol(X),
      sha256=sha(file.path(folder,bin))))
  },error=function(e) {
    attempt$status <- 'FAILED';attempt$error <- conditionMessage(e);attempt$warnings <- warnings;attempt$ended_utc <- utc()
    writej(attempt,ap)
    list(error=conditionMessage(e),diagnostic=data.frame(model_id=id,imputation=j,N=c$N,events=c$events,
      parameters=ncol(X),rank=qr(X)$rank,complete_df=c$N-ncol(X),converged=FALSE,iterations=NA_integer_,
      warnings=paste(warnings,collapse=' | '),cov_min_eigen=NA_real_,condition_number=kappa(X,exact=TRUE),
      max_abs_coefficient=NA_real_,probability_min=NA_real_,probability_max=NA_real_,
      probability_below_1e8=NA_integer_,probability_above_1minus1e8=NA_integer_,max_cook=NA_real_,
      max_leverage=NA_real_,max_abs_DFBETAS=NA_real_,cook_screen_count=NA_integer_,leverage_screen_count=NA_integer_,
      dfbetas_screen_count=NA_integer_,LP_status='NOT_RUN_GLM_FAILED',fit_pass=FALSE))
  })
  result
}
curve_outputs <- function(c,pooled,ranges) {
  dayrows <- list();dynrows <- list()
  for(id in names(pooled)) {
    a <- pooled[[id]]
    if(startsWith(id,'daily')) {
      j <- as.integer(substr(id,6,6));layers <- if(endsWith(id,'_overall')) 'ALL' else paste0('L0',1:5)
      for(layer in layers) {
        rr <- ranges[ranges$day==j&ranges$age_layer==layer,]
        age <- rr$representative_age_days
        ag <- if(layer=='ALL') c(0,0) else as.numeric(rcs3(log1p(age),c$cfg$age_log_knots,log1p(c$cfg$reference_age_days)))
        grid <- sort(unique(c(seq(rr$p025,rr$p975,length.out=161),c$cfg$reference_sodium)))
        grid <- grid[grid>=rr$p025&grid<=rr$p975]
        for(x in grid) {
          s <- as.numeric(rcs3(x,c$cfg$sodium_knots,c$cfg$reference_sodium,c$cfg$sodium_scale))
          C <- setNames(rep(0,nrow(a$Q)),rownames(a$Q));C[c('s1','s2')] <- s
          if(layer!='ALL') for(k in 1:2)for(t in 1:2)C[paste0('s',k,'_a',t)] <- s[k]*ag[t]
          dayrows[[length(dayrows)+1L]] <- cbind(data.frame(model_id=id,day=j,
            scope=if(layer=='ALL') 'overall' else 'representative_age',age_layer=layer,age_days=age,sodium=x),scalar(a,C))
        }
      }
    } else {
      ar <- data.frame(point_role='grid',age_layer='ALL',age_days=expm1(seq(log1p(min(c$raw$age_days)),log1p(max(c$raw$age_days)),length.out=201)))
      reps <- ranges[ranges$day==1&ranges$age_layer!='ALL',]
      ar <- rbind(ar,data.frame(point_role='representative_age',age_layer=reps$age_layer,age_days=reps$representative_age_days))
      for(i in seq_len(nrow(ar)))for(segment in c('delta12','delta23')) {
        ag <- as.numeric(rcs3(log1p(ar$age_days[i]),c$cfg$age_log_knots,log1p(c$cfg$reference_age_days)))
        C <- setNames(rep(0,nrow(a$Q)),rownames(a$Q));C[segment] <- 1
        C[paste0(segment,'_a',1:2)] <- ag
        dynrows[[length(dynrows)+1L]] <- cbind(data.frame(model_id=id,segment=segment,
          point_role=ar$point_role[i],age_layer=ar$age_layer[i],age_days=ar$age_days[i]),scalar(a,C))
      }
    }
  }
  list(day=do.call(rbind,dayrows),dynamic=do.call(rbind,dynrows))
}
models <- function(c) {
  gate <- code_gate(c);out <- file.path(c$root,'results');private <- file.path(out,'restricted')
  start <- file.path(out,'models_started.json')
  check(!file.exists(start)&&!dir.exists(private),'Models already attempted; no overwrite/restart without new authorization')
  dir.create(private,recursive=TRUE);writej(c(list(started_utc=utc(),status='ATTEMPTING_FROZEN_MODELS',
    source_hashes=as.list(c$hashes),input_hashes=as.list(c$input_hashes)),gate),start)
  datasets <- lapply(seq_len(c$cfg$M),function(j)load_completed(c,j))
  allfits <- list();pooled <- list();diags <- list();failed <- character();tests <- list();coeffs <- list()
  for(id in c$cfg$model_ids) {
    folder <- file.path(private,id);dir.create(folder)
    attempts <- lapply(seq_len(c$cfg$M),function(j)fit_one(datasets[[j]],c,id,j,folder))
    ok <- vapply(attempts,function(a)is.null(a$error),logical(1));designs <- lapply(attempts[ok],`[[`,'design')
    writej(list(model=id,designs=designs),file.path(folder,'design_metadata.json'))
    lpstatus <- if(length(designs)) system2('pwsh',c('-NoProfile','-File',shQuote(c$paths[['python_helper']]),
      shQuote(c$paths[['separation']]),shQuote(folder)),stdout=file.path(folder,'separation_stdout.txt'),
      stderr=file.path(folder,'separation_stderr.txt')) else NA_integer_
    lp <- if(file.exists(file.path(folder,'separation_validation.json')))
      jsonlite::read_json(file.path(folder,'separation_validation.json'),simplifyVector=TRUE) else NULL
    model_pass <- all(ok)&&identical(lpstatus,0L)&&isTRUE(lp$all_pass)&&length(designs)==c$cfg$M
    for(j in seq_len(c$cfg$M)) {
      row <- attempts[[j]]$diagnostic
      if(ok[j]&&!is.null(lp)) {
        ix <- match(j,lp$results$imputation)
        if(!is.na(ix)) {row$LP_status <- lp$results$status[ix];row$fit_pass <- isTRUE(lp$results$pass[ix])&&identical(lpstatus,0L)}
        else row$LP_status <- 'LP_RECEIPT_MISSING'
      } else if(ok[j]) row$LP_status <- 'LP_EXECUTION_FAILED'
      ap <- file.path(folder,sprintf('attempt_%02d.json',j));aa <- jsonlite::read_json(ap)
      aa$LP_exit_code <- lpstatus;aa$LP_status <- row$LP_status;aa$fit_pass <- row$fit_pass
      aa$status <- if(isTRUE(row$fit_pass)) 'PASSED_ALL_DIAGNOSTICS' else 'FAILED';writej(aa,ap)
      diags[[length(diags)+1L]] <- row
    }
    allfits[[id]] <- lapply(attempts,`[[`,'fit')
    saveRDS(allfits,file.path(private,'fits.rds'))
    write.csv(do.call(rbind,diags),file.path(out,'fit_diagnostics.csv'),row.names=FALSE)
    if(model_pass) {
      a <- pool_model(allfits[[id]],c$N);pooled[[id]] <- a
      for(term in rownames(a$Q)) {
        C <- setNames(rep(0,nrow(a$Q)),rownames(a$Q));C[term] <- 1
        z <- scalar(a,C);z$term <- term;z$model_id <- id;z$MCSE <- sqrt(a$B[term,term]/a$M)
        z$relative_MCSE <- z$MCSE/z$SE;coeffs[[length(coeffs)+1L]] <- z
      }
      if(startsWith(id,'dynamic')) {
        tests[[length(tests)+1L]] <- joint(a,dynamic_terms,id,'dynamic_age_4df',if(id=='dynamic_main') 'primary' else 'support')
        tests[[length(tests)+1L]] <- joint(a,c('delta12','delta23',dynamic_terms),id,'all_dynamic_6df','support')
        tests[[length(tests)+1L]] <- joint(a,initial_terms,id,'initial_age_4df','secondary')
      } else if(endsWith(id,'_overall')) tests[[length(tests)+1L]] <- joint(a,c('s1','s2'),id,'daily_Na_2df','secondary','daily_overall_3')
      else tests[[length(tests)+1L]] <- joint(a,initial_terms,id,'daily_age_4df','secondary','daily_age_3')
    } else failed <- c(failed,id)
    writej(list(status=if(model_pass) 'MODEL_ALL50_DIAGNOSTICS_PASSED' else 'FAILED_MODEL_NO_POOL_RELEASE',
      model_id=id,attempted=c$cfg$M,GLM_passed=sum(ok),LP_exit_code=lpstatus,ended_utc=utc()),file.path(folder,'model_receipt.json'))
    cat('Actual attempts finished for',id,'; all50 diagnostic pass:',model_pass,'\n')
  }
  assert_unchanged(c)
  if(length(failed)) {
    writej(c(list(status='FAILED_MODELS_NO_COMPLETE_RELEASE',failed_models=failed,fits_attempted=nrow(do.call(rbind,diags)),
      fits_passed=sum(vapply(diags,function(z)isTRUE(z$fit_pass),logical(1))),source_hashes=as.list(c$hashes),
      input_hashes=as.list(c$input_hashes),ended_utc=utc(),no_fallback=TRUE),gate),file.path(out,'analysis_receipt.json'))
    stop('One or more frozen models failed; all attempted fits retained; no complete result release',call.=FALSE)
  }
  jt <- do.call(rbind,tests)
  for(family in c('daily_overall_3','daily_age_3')) {
    ii <- jt$holm_family==family;check(sum(ii)==3L,'Holm family incomplete');jt$p_holm[ii] <- p.adjust(jt$p[ii],method='holm')
  }
  ranges <- support_ranges(c);curves <- curve_outputs(c,pooled,ranges)
  write.csv(jt,file.path(out,'joint_tests.csv'),row.names=FALSE,na='')
  write.csv(do.call(rbind,coeffs),file.path(out,'pooled_coefficients.csv'),row.names=FALSE)
  write.csv(ranges,file.path(out,'support_ranges.csv'),row.names=FALSE,na='')
  write.csv(curves$day,file.path(out,'day_curves.csv'),row.names=FALSE,na='')
  write.csv(curves$dynamic,file.path(out,'dynamic_curves.csv'),row.names=FALSE,na='')
  saveRDS(pooled,file.path(private,'pooled_models.rds'))
  outputs <- c('joint_tests.csv','pooled_coefficients.csv','fit_diagnostics.csv','support_ranges.csv',
               'day_curves.csv','dynamic_curves.csv','restricted/pooled_models.rds','restricted/fits.rds')
  assert_unchanged(c)
  writej(c(list(status='COMPLETED_PENDING_INDEPENDENT_RESULT_REVIEW',ended_utc=utc(),N=c$N,events=c$events,M=c$cfg$M,
    models_completed=length(pooled),fits_attempted=length(diags),fits_passed=sum(vapply(diags,function(z)isTRUE(z$fit_pass),logical(1))),
    model_formulas=lapply(pooled,`[[`,'formula'),source_hashes=as.list(c$hashes),input_hashes=as.list(c$input_hashes),
    output_hashes=as.list(setNames(vapply(file.path(out,outputs),sha,character(1)),outputs)),
    platelet_knots=c$platelet_knots,interpretation='Conditional odds at accepted72h landmark; no causal treatment interpretation',
    uncertainty_limits=c('Raw PMM then nonlinear clinical transforms is approximate FCS','Tail donor limitations persist',
      'Deterministic reconstruction/model-selection uncertainty not included'),protected_bound_inputs_unchanged=TRUE),gate),
    file.path(out,'analysis_receipt.json'))
  cat('COMPLETED_PENDING_INDEPENDENT_RESULT_REVIEW\n')
}
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  check(length(args)==2L&&args[2]%in%c('preflight','models'),'Usage: outcome_engine.R ROOT preflight|models')
  c <- context(args[1]);set.seed(c$cfg$seed)
  if(args[2]=='preflight') preflight(c) else models(c)
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  tryCatch(main(args),error=function(e) {
    if(length(args)>=1L&&dir.exists(args[1])&&requireNamespace('jsonlite',quietly=TRUE)) {
      folder <- file.path(args[1],'results');dir.create(folder,showWarnings=FALSE)
      writej(list(status='FAILED_NOT_PASS',phase=if(length(args)>1L)args[2] else 'invalid_cli',
        error=conditionMessage(e),utc=utc()),file.path(folder,paste0('execution_failure_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.json')))
    }
    message(conditionMessage(e));quit(status=1L)
  })
}
