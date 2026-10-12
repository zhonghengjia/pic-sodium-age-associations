args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
base <- args[1]; out <- args[2]
projlib <- file.path(dirname(base),'landmark72_reanalysis_20261006/runtime/R-library')
if(dir.exists(projlib)) .libPaths(c(projlib,.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(jsonlite);library(digest)})
sha <- function(p) digest(file=p,algo='sha256',serialize=FALSE)
check <- function(x,m) if(!isTRUE(x))stop(m,call.=FALSE)
run <- file.path(base,'analysis_daily_mean_20261009_v1')
paths <- c(observed=file.path(run,'results/display/Figure2A_outcome_points.csv'),
           adjusted=file.path(run,'图表审核_20261009_v5/source/equal_age_group_points.csv'),
           original=file.path(run,'results/dynamic_curves.csv'),
           joint=file.path(run,'results/joint_tests.csv'),
           observed_receipt=file.path(run,'results/display/Figure2A_outcome_preflight_receipt.json'),
           figure_bindings=file.path(base,'figures_tables_mainline_20261009_v1/figures/SOURCE_BINDINGS.json'),
           old_main_png=file.path(base,'图表重排_20261009_v3/主图/Figure3.png'),
           old_main_pdf=file.path(base,'图表重排_20261009_v3/主图/Figure3.pdf'),
           old_supp_png=file.path(base,'图表重排_20261009_v3/补充图/Supplementary_Figure1.png'),
           old_supp_pdf=file.path(base,'图表重排_20261009_v3/补充图/Supplementary_Figure1.pdf'))
check(all(file.exists(paths)),'Missing source')
before <- vapply(paths,sha,'')
bd <- fromJSON(paths[['figure_bindings']],simplifyVector=FALSE)$sources
check(sha(paths[['adjusted']])==bd$dynamic_points$sha256,'Accepted dynamic point source changed')
check(sha(paths[['original']])==bd$dynamic_original$sha256,'Original dynamic source changed')
pf <- fromJSON(paths[['observed_receipt']],simplifyVector=FALSE)
check(sha(paths[['observed']])==pf$points_sha256,'Accepted observed point source changed')
obs <- read.csv(paths[['observed']],stringsAsFactors=FALSE)
adj <- read.csv(paths[['adjusted']],stringsAsFactors=FALSE)
orig <- read.csv(paths[['original']],stringsAsFactors=FALSE)
jt <- read.csv(paths[['joint']],stringsAsFactors=FALSE)
layers <- c('ALL',paste0('L0',1:5))
Ns <- c(3459,822,1223,809,371,234); deaths <- c(306,94,94,67,32,19)
age_names <- c('Overall','Neonates (\u226428 d)','Infants (>28 d to <1 y)','Children 1\u20134 y','Children 5\u20139 y','Children 10\u201314 y')
check(nrow(obs)==36L&&!anyDuplicated(paste(obs$age_layer,obs$outcome,obs$day)),'Observed36 mismatch')
check(setequal(obs$age_layer,layers)&&setequal(obs$day,1:3)&&setequal(obs$outcome,0:1),'Observed design mismatch')
for(i in seq_along(layers)){
 z<-obs[obs$age_layer==layers[i],]
 check(all(z$n[z$outcome==0]==Ns[i]-deaths[i])&&all(z$n[z$outcome==1]==deaths[i]),'Outcome/group counts changed')
}
check(all(is.finite(as.matrix(obs[c('mean','ci_low','ci_high')]))),'Observed nonfinite')
check(all(obs$ci_low<=obs$mean & obs$ci_high>=obs$mean),'Observed interval mismatch')
z <- orig[orig$model_id=='dynamic_main' & orig$point_role=='representative_age',]
check(nrow(adj)==10L&&nrow(z)==10L,'Adjusted10 mismatch')
key <- function(d)paste(d$segment,d$age_layer)
z<-z[match(key(adj),key(z)),]
for(k in c('OR','ci_low','ci_high','age_days'))check(all(abs(adj[[k]]-z[[k]])<1e-12),paste('Adjusted original mismatch',k))
check(all(adj$ci_low>0 & adj$ci_low<=adj$OR & adj$ci_high>=adj$OR),'Adjusted intervals invalid')
primary <- jt[jt$model_id=='dynamic_main' & jt$test_id=='dynamic_age_4df',]
check(nrow(primary)==1L&&primary$df1==4&&primary$N==3459&&primary$events==306,'Joint test mismatch')
p_label<-if(primary$p<.001)'P < 0.001' else sprintf('P = %.3f',primary$p)
write.csv(obs,file.path(out,'results/Observed_daily_points.csv'),row.names=FALSE)
write.csv(adj,file.path(out,'results/Adjusted_dynamic_points.csv'),row.names=FALSE)
write.csv(primary,file.path(out,'results/Primary_joint_interaction.csv'),row.names=FALSE)
for(k in c('mean','ci_low','ci_high'))check(all(abs(read.csv(file.path(out,'results/Observed_daily_points.csv'))[[k]]-obs[[k]])<1e-12),'Observed display source changed')
for(k in c('OR','ci_low','ci_high','age_days'))check(all(abs(read.csv(file.path(out,'results/Adjusted_dynamic_points.csv'))[[k]]-adj[[k]])<1e-12),'Adjusted display source changed')
facet_labels<-setNames(paste0(age_names,'\nSurvived n=',Ns-deaths,'; died n=',deaths),layers)
obs$facet <- factor(facet_labels[obs$age_layer],levels=facet_labels)
obs$group <- factor(ifelse(obs$outcome==0,'Survived hospitalization','Subsequent in-hospital death'),
                    levels=c('Survived hospitalization','Subsequent in-hospital death'))
pal<-c('Survived hospitalization'='#3C8DCC','Subsequent in-hospital death'='#FF6268')
yr<-c(floor(min(obs$ci_low)/5)*5,ceiling(max(obs$ci_high)/5)*5)
top<-ggplot(obs,aes(day,mean,group=group,colour=group,shape=group,linetype=group))+
 geom_hline(yintercept=c(135,145),colour='#C9CDD3',linetype='dotted',linewidth=.35)+
 geom_line(linewidth=.8)+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.09,linewidth=.55)+
 geom_point(size=2.3)+facet_wrap(~facet,ncol=3,axes='all',axis.labels='all')+
 scale_colour_manual(values=pal)+scale_linetype_manual(values=c('solid','dashed'))+
 scale_shape_manual(values=c(16,17))+
 scale_x_continuous(breaks=1:3,labels=paste0('Day ',1:3),limits=c(.9,3.1),expand=c(0,0))+
 scale_y_continuous(limits=yr,breaks=seq(yr[1],yr[2],5),expand=expansion(mult=.025))+
 labs(x=NULL,y='Mean sodium (mmol/L)',colour=NULL,shape=NULL,linetype=NULL)+
 theme_classic(base_size=11,base_family='Arial')+
 theme(strip.background=element_rect(fill='white',colour='#AAB0B7',linewidth=.4),strip.text=element_text(size=10.5),
       axis.text=element_text(colour='black'),axis.title.y=element_text(size=11.5),
       legend.position='bottom',legend.text=element_text(size=10.5),legend.key.width=grid::unit(15,'mm'),
       panel.spacing=grid::unit(6,'mm'),plot.margin=margin(4,12,3,8))
age_labels<-c('Neonates\n\u226428 d','Infants\n>28 d to <1 y','Children\n1\u20134 y','Children\n5\u20139 y','Children\n10\u201314 y')
or_limits<-c(.3,2.5)
check(all(adj$ci_low>=or_limits[1] & adj$ci_high<=or_limits[2]),'Adjusted CI clipped')
make_bottom<-function(seg,title){
 d<-adj[adj$segment==seg,];d$x<-match(d$age_layer,paste0('L0',1:5));check(nrow(d)==5&&!anyNA(d$x),'Segment age mismatch')
 ggplot(d,aes(x,OR))+
 geom_hline(yintercept=1,colour='#459F91',linetype='dashed',linewidth=.55)+
 geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.12,colour='#FF6268',linewidth=.75)+
 geom_point(colour='#FF6268',size=3)+
 scale_x_continuous(breaks=1:5,labels=age_labels,limits=c(.65,5.35),expand=c(0,0))+
 scale_y_log10(limits=or_limits,breaks=c(.3,.5,1,1.5,2,2.5),labels=c('0.3','0.5','1','1.5','2','2.5'),expand=expansion(mult=.025))+
 labs(title=title,x='Age group',y='Adjusted OR (95% CI)')+
 theme_classic(base_size=11,base_family='Arial')+
 theme(plot.title=element_text(face='bold',size=12),axis.text=element_text(colour='black'),
       axis.text.x=element_text(size=9.5,lineheight=1.05),axis.title.y=element_text(size=11.5),
       axis.title.x=element_text(margin=margin(t=8)),plot.margin=margin(8,13,8,9))
}
head<-function(title,subtitle){ggplot()+annotate('text',x=0,y=.78,label=title,hjust=0,fontface='bold',size=4.8,family='Arial')+
 annotate('text',x=0,y=.22,label=subtitle,hjust=0,size=3.6,family='Arial')+
 scale_x_continuous(limits=c(0,1),expand=c(0,0))+scale_y_continuous(limits=c(0,1),expand=c(0,0))+
 theme_void()+theme(plot.margin=margin(0,8,0,8))}
ha<-head('A  Observed three-day sodium changes','Unadjusted group means of patient daily means; error bars: 95% CI')
hb<-head('B  Adjusted associations of sodium changes',paste0('OR per 5 mmol/L greater change; joint age-by-change interaction: ',p_label))
bottom<-make_bottom('delta12','Day 1 to Day 2 change')+make_bottom('delta23','Day 2 to Day 3 change')
fig<-wrap_plots(list(ha,top,hb,bottom),ncol=1,heights=c(.6,5.9,.65,3.7))
stem<-file.path(out,'图形/Figure3_dynamic_composite')
ggsave(paste0(stem,'.png'),fig,width=12,height=11.6,units='in',dpi=300,device=ragg::agg_png,bg='white')
ggsave(paste0(stem,'.pdf'),fig,width=12,height=11.6,units='in',device=grDevices::cairo_pdf,bg='white')
caption<-paste('Figure 3. Observed three-day sodium changes and adjusted associations with subsequent in-hospital mortality.',
 'A, Group averages of patient daily mean sodium by subsequent hospitalization outcome, overall and in five descriptive age groups. Days 1, 2 and 3 represent 0-<24, 24-<48 and 48-<72 hours after first ICU admission. Each child contributes one daily mean per window without weighting by measurement count. Error bars are the original pointwise 95% t confidence intervals; lines connect group summaries. Dotted horizontal lines indicate 135 and 145 mmol/L.',
 'B, Adjusted odds ratios per 5 mmol/L greater change for the Day 1-to-2 and Day 2-to-3 intervals. The five equally spaced points in each panel are estimates from one common continuous-age model evaluated at median ages of 1.25, 121.00, 769.24, 2517.57 and 4464.29 days; they are not independently fitted subgroup effects. Error bars are the original pointwise 95% confidence intervals pooled across fifty clinical-variable imputations. The OR axis uses a standard logarithmic scale; age points are not connected.',
 'The Day 1-to-2 contrast holds Day 1 sodium and the Day 2-to-3 change fixed, so Day 2 and Day 3 sodium both increase by 5 mmol/L. The Day 2-to-3 contrast holds Day 1 sodium and the Day 1-to-2 change fixed. Adjustment includes age layer and continuous age, sex, first ICU type, first-24-hour albumin, creatinine, BUN, lactate, glucose and pH, with Day 1 sodium and the defined age interactions. The displayed P value is the primary joint age-by-change interaction test (4 df).',
 'The cohort includes 3,459 children alive and hospitalized at 72 hours, with 306 subsequent in-hospital deaths. Daily means retain the two author-approved deterministic sodium reconstructions. These displays describe observed changes and adjusted associations, not treatment effects.')
writeLines(caption,paste0(stem,'_caption.txt'),useBytes=TRUE)
check(identical(vapply(paths,sha,''),before),'Protected inputs/old figures changed')
write_json(list(source_hashes=as.list(before),source_paths=as.list(paths),source_unchanged=TRUE,
 observed_points=36,adjusted_points=10,N=3459,deaths=306,joint_p=primary$p,
 observed_y_limits=yr,adjusted_OR_limits=or_limits,all_original_intervals_displayed=TRUE,
 age_points_connected=FALSE,models_refitted=FALSE,imputation_rerun=FALSE,
 independent_agent_review=FALSE,status='AUTHOR_PREVIEW_PENDING_ACTUAL_VISUAL_CHECK'),
 file.path(out,'control/render_receipt.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(out,'control/R_session.txt'))
cat('FIGURE3_PREVIEW_RENDERED; OBSERVED=36; ADJUSTED=10; ALL_SOURCE_VALUES_UNCHANGED\n')
