# Data-stage MI only; all numerical inputs must have current independent cleaning acceptance.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2,args[2]%in%c('preflight','mi'))
root<-normalizePath(args[1],winslash='/',mustWork=TRUE)
if(nzchar(Sys.getenv("PIC_R_LIBRARY"))) .libPaths(c(Sys.getenv("PIC_R_LIBRARY"),.libPaths()))
sha<-function(p)digest::digest(file=p,algo='sha256',serialize=FALSE)
js<-function(x,p)jsonlite::write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=17,na='null')
same_hashes<-function(candidate,current){
 z<-unlist(candidate,use.names=TRUE)
 is.character(z)&&length(z)==length(current)&&!anyNA(z)&&setequal(names(z),names(current))&&identical(unname(z[names(current)]),unname(current))
}
p<-c(input=file.path(root,'data/restricted/clean_observed_patient.csv'),ready=file.path(root,'data/clean_data_ready.json'),
 policy=file.path(root,'control/AUXILIARY_QC_POLICY.json'),plan=file.path(root,'control/MI_PLAN.md'),config=file.path(root,'control/MI_config.json'),code=file.path(root,'code/run_clean_MI.R'),
 export_code=file.path(root,'code/preserve_export_literals.py'),clean_release=file.path(root,'supervision/cleaned_data_release.json'))
stopifnot(all(file.exists(p)));hashes<-vapply(p,sha,character(1))
cleanrelease<-jsonlite::read_json(file.path(root,'supervision/cleaned_data_release.json'),simplifyVector=TRUE)
stopifnot(isTRUE(cleanrelease$status=='ACCEPT'),isTRUE(cleanrelease$input_sha256==unname(hashes['input'])))
ready<-jsonlite::read_json(p['ready'],simplifyVector=TRUE)
stopifnot(isTRUE(ready$input_sha256==unname(hashes['input'])),isTRUE(ready$N==3459),isTRUE(ready$Na_exact_preserved))
d<-read.csv(p['input'],check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(d)==3459,!anyDuplicated(d$id),sum(d$death)==306)
cfg<-jsonlite::read_json(p['config'],simplifyVector=TRUE);labs<-cfg$clinical_labs
stopifnot(cfg$M==50,cfg$iterations==30,cfg$donors==5,cfg$seed==20261009)
stopifnot(all(labs%in%names(d)),all(is.finite(as.matrix(d[c('D1','D2','D3','age_days')]))))
complete_labs<-labs[!vapply(d[labs],anyNA,logical(1))];targets<-setdiff(labs,complete_labs)
knots<-cfg
rcs3<-function(x,k,ref,scale=1){
 h<-function(v)(pmax(v-k[1],0)^3-(k[3]-k[1])/(k[3]-k[2])*pmax(v-k[2],0)^3+(k[2]-k[1])/(k[3]-k[2])*pmax(v-k[3],0)^3)/(k[3]-k[1])^2
 cbind((x-ref)/scale,(h(x)-h(ref))/scale)
}
a<-rcs3(log1p(d$age_days),knots$age_log_knots,log1p(365.25));d$a1<-a[,1];d$a2<-a[,2]
d$age_layer<-factor(d$age_layer,levels=paste0('L0',1:5));d$sex<-factor(d$sex,levels=c('F','M'))
d$first_careunit<-factor(d$first_careunit,levels=c('PICU','CICU','General ICU','NICU','SICU'))
aux<-c('age_layer','sex','first_careunit','death','a1','a2',complete_labs)
for(j in 1:3){
 observed_times<-d[[paste0('nD',j)]]-d[[paste0('D',j,'_reconstructed_n')]]
 stopifnot(all(is.finite(observed_times)),all(observed_times>=0),all(observed_times==floor(observed_times)))
 nm<-paste0('obs',j);d[[nm]]<-log1p(observed_times);aux<-c(aux,nm)
 s<-rcs3(d[[paste0('D',j)]],knots$sodium_knots,140,5)
 for(k in 1:2){
  ns<-paste0('na',j,'_',k);d[[ns]]<-s[,k];aux<-c(aux,ns)
  for(t in 1:2){ni<-paste0(ns,'_a',t);d[[ni]]<-s[,k]*a[,t];aux<-c(aux,ni)}
 }
 for(l in paste0('L0',2:5)){ni<-paste0('na',j,'_1_',l);d[[ni]]<-s[,1]*as.numeric(d$age_layer==l);aux<-c(aux,ni)}
}
stopifnot(!anyNA(d[aux]));mm<-model.matrix(~.,d[aux]);mq<-qr(mm);stopifnot(mq$rank==ncol(mm))
qq<-qr.Q(mq)[,-1,drop=FALSE]*sqrt(nrow(mm)-1);colnames(qq)<-sprintf('auxQR_%02d',seq_len(ncol(qq)))
newmm<-cbind('(Intercept)'=1,qq);back<-qr.coef(qr(newmm),mm);forward<-qr.coef(mq,newmm)
err<-max(abs(newmm%*%back-mm));relerr<-sqrt(sum((newmm%*%back-mm)^2)/sum(mm^2))
stopifnot(err<1e-7,relerr<1e-10,max(abs(mm%*%forward-newmm))<1e-7)
x<-cbind(d[targets],as.data.frame(qq));methods<-setNames(rep('',ncol(x)),names(x));methods[targets]<-'pmm'
pm<-matrix(0L,ncol(x),ncol(x),dimnames=list(names(x),names(x)));pm[targets,]<-1L;diag(pm)<-0L
target_rank<-sapply(targets,function(v)qr(newmm[!is.na(x[[v]]),,drop=FALSE])$rank)
stopifnot(all(target_rank==ncol(mm)),all(vapply(x,function(v)length(unique(v[!is.na(v)]))>1,logical(1))))
pre<-list(status='PREFLIGHT_COMPLETE_PENDING_CODE_RELEASE',N=nrow(d),deaths=sum(d$death),labs=labs,targets=targets,
 complete_labs=complete_labs,missing_cells=sum(is.na(d[labs])),M=50L,iterations=30L,donors=5L,seed=20261009L,
 original_auxiliary=aux,auxiliary_design_rank=ncol(mm),target_aux_rank=as.list(target_rank),QR_space_error=err,QR_relative_space_error=relerr,
 predictor_names=names(x),methods=as.list(methods),predictor_matrix=unname(pm),source_hashes=as.list(hashes),
 excluded=c('INR','LOS','K_source5226','ALT','bilirubin','pCO2','perioperative'),no_outcome_models=TRUE)
if(args[2]=='preflight'){dir.create(file.path(root,'results'),showWarnings=FALSE);js(pre,file.path(root,'results/MI_preflight.json'));cat('Preflight only; MI not started.\n');quit(status=0)}
release<-jsonlite::read_json(file.path(root,'supervision/MI_code_release.json'),simplifyVector=TRUE)
stopifnot(isTRUE(release$status=='ACCEPT'),same_hashes(release$source_hashes,hashes))
prefile<-file.path(root,'results/MI_preflight.json');stopifnot(file.exists(prefile),isTRUE(release$preflight_sha256==sha(prefile)))
pref<-jsonlite::read_json(prefile,simplifyVector=TRUE)
stopifnot(isTRUE(pref$status=='PREFLIGHT_COMPLETE_PENDING_CODE_RELEASE'),same_hashes(pref$source_hashes,hashes))
out<-file.path(root,'results/MI');private<-file.path(root,'results/restricted_MI');stopifnot(!dir.exists(out),!dir.exists(private))
dir.create(out,recursive=TRUE);dir.create(private,recursive=TRUE);js(pre,file.path(out,'started.json'))
warns<-character();imp<-withCallingHandlers(mice::mice(x,m=50,maxit=30,seed=20261009,method=methods,predictorMatrix=pm,
 donors=5,printFlag=FALSE,remove.constant=FALSE,remove.collinear=FALSE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(list(imp=imp,input=x,master=d,predictor=pm,labs=labs,targets=targets,complete_labs=complete_labs,source_hashes=hashes,
 original_auxiliary=d[aux],original_design=mm,QR_design=qq,back=back,forward=forward,auxiliary_space_error=err),file.path(private,'mi50_clean.rds'))
checks<-list();dist<-list();chain<-list()
for(j in 1:50){
 z<-mice::complete(imp,j);stopifnot(!anyNA(z))
 for(v in names(x)){obs<-!is.na(x[[v]]);stopifnot(identical(z[[v]][obs],x[[v]][obs]))}
 for(v in targets){mis<-is.na(x[[v]]);stopifnot(all(z[[v]][mis]%in%x[[v]][!mis]))}
 master<-read.csv(p['input'],check.names=FALSE,stringsAsFactors=FALSE)
 for(v in targets)master[[v]]<-z[[v]]
 master$derived_BUN_mg_dL<-master$urea*2.8
 write.csv(master,file.path(private,sprintf('completed_%02d.csv',j)),row.names=FALSE,na='')
 checks[[j]]<-list(imputation=j,complete=TRUE,observed_unchanged=TRUE,donor_member=TRUE)
 for(v in targets){
  obs<-!is.na(x[[v]]);u<-x[[v]][obs];b<-z[[v]][!obs]
  dist[[length(dist)+1]]<-data.frame(imputation=j,variable=v,group='overall',level='ALL',observed_n=length(u),imputed_n=length(b),
   observed_median=median(u),imputed_median=median(b),observed_min=min(u),observed_max=max(u),imputed_min=min(b),imputed_max=max(b))
  for(g in c('age_layer','first_careunit'))for(level in levels(d[[g]])){
   ix<-d[[g]]==level;uo<-x[[v]][ix&obs];ui<-z[[v]][ix&!obs]
   dist[[length(dist)+1]]<-data.frame(imputation=j,variable=v,group=g,level=level,observed_n=length(uo),imputed_n=length(ui),
    observed_median=if(length(uo))median(uo)else NA_real_,imputed_median=if(length(ui))median(ui)else NA_real_,
    observed_min=if(length(uo))min(uo)else NA_real_,observed_max=if(length(uo))max(uo)else NA_real_,
    imputed_min=if(length(ui))min(ui)else NA_real_,imputed_max=if(length(ui))max(ui)else NA_real_)
  }
 }
}
for(v in targets){
 a<-imp$chainMean[v,,];b<-imp$chainVar[v,,];taila<-a[16:30,,drop=FALSE]
 within<-mean(apply(taila,2,var));between<-15*var(colMeans(taila));rh<-if(is.finite(within)&&within>0)sqrt((14/15*within+between/15)/within)else NA_real_
 chain[[v]]<-data.frame(variable=v,missing_n=sum(is.na(x[[v]])),rhat_style_last15=rh,variance_chain_estimable=all(is.finite(b)),
  first_half_mean=mean(a[1:15,]),last_half_mean=mean(a[16:30,]))
 grDevices::png(file.path(out,paste0('trace_',v,'.png')),width=1500,height=650,res=140)
 par(mfrow=c(1,2),mar=c(4,4,3,1));matplot(a,type='l',lty=1,col=gray(seq(.1,.8,length.out=50)),xlab='Iteration',ylab='Imputed mean',main=v)
 if(any(is.finite(b)))matplot(b,type='l',lty=1,col=gray(seq(.1,.8,length.out=50)),xlab='Iteration',ylab='Imputed variance',main='Variance chain')else{plot.new();text(.5,.5,'Variance not estimable: one missing value')}
 dev.off()
}
write.csv(do.call(rbind,chain),file.path(out,'chain_diagnostics.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,dist),file.path(out,'observed_imputed_distributions.csv'),row.names=FALSE,na='')
saveRDS(imp$chainMean,file.path(private,'chain_means.rds'));saveRDS(imp$chainVar,file.path(private,'chain_variances.rds'))
js(list(status='MI_EXECUTED_PENDING_INDEPENDENT_REVIEW',M=imp$m,iterations=imp$iteration,seed=imp$seed,
 warnings=warns,loggedEvents=imp$loggedEvents,checks=checks,targets=targets,complete_labs=complete_labs,
 missing_cells=sum(is.na(x[targets])),source_hashes=as.list(hashes),object_sha256=sha(file.path(private,'mi50_clean.rds')),
 imputation_methods=as.list(imp$method),no_outcome_models=TRUE,reconstruction_uncertainty_propagated=FALSE),file.path(out,'execution_receipt.json'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));cat('MI executed; awaiting actual independent review.\n')
