args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grep('^--file=', args)])
O <- dirname(dirname(normalizePath(script, winslash='/')))
pkgs <- c('survival','ggplot2','patchwork','ragg','svglite','jsonlite')
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)))
library(survival); library(ggplot2); library(patchwork)
d <- read.csv(file.path(O,'data/restricted/time_input.csv'))
stopifnot(nrow(d)==3459,sum(d$event=='Death')==306,all(d$time_days>0))
d$event <- factor(d$event, levels=c('Censored','Death','Discharged alive'))
d$sodium_status <- factor(d$sodium_status,levels=c('Low','Normal','High'))
layers <- c('Overall','L01','L02','L03','L04','L05')
titles <- c('Overall','Neonates (\u226428 d)','Infants (>28 d\u2013<1 y)',
            'Children 1\u20134 y','Children 5\u20139 y','Children 10\u201314 y')
colours <- c(Low='#F16B6B',Normal='#009D9A',High='#4B8EC9')
labels <- c('Low (<135)','Normal (135\u2013145)','High (>145)')
ticks <- c(0,7,14,21,30)
curves <- list(); risk <- list(); counts <- list(); checks <- list(); fits <- list()
for (j in seq_along(layers)) {
 z <- if (j==1) d else d[d$age_layer==layers[j],]
 stopifnot(all(table(z$sodium_status)>0))
 fits[[layers[j]]] <- list()
 for(g in levels(d$sodium_status)) {
  q <- z[z$sodium_status==g,]
  fit <- survfit(Surv(time_days,event)~1,data=q,conf.int=.95,conf.type='log-log')
  stopifnot(inherits(fit,'survfitms'),all(c('Death','Discharged alive')%in%fit$states))
  k <- match('Death',fit$states); cidx <- match('Discharged alive',fit$states)
  pp <- data.frame(time=fit$time,estimate=fit$pstate[,k],lower=fit$lower[,k],upper=fit$upper[,k],
                   discharge=fit$pstate[,cidx],in_hospital=fit$pstate[,1])
  stopifnot(all(is.finite(as.matrix(pp))),all(pp$lower<=pp$estimate+1e-12),all(pp$upper>=pp$estimate-1e-12),
            all(abs(pp$estimate+pp$discharge+pp$in_hospital-1)<1e-10),all(diff(pp$estimate)>=-1e-12))
  # Independent event-count AJ recursion (not an independent-discharge KM).
  ts <- sort(unique(q$time_days)); p0 <- 1; pdth <- 0; pdis <- 0; alt <- numeric(length(ts))
  for(i in seq_along(ts)) {
   nr <- sum(q$time_days>=ts[i]); nd <- sum(q$time_days==ts[i]&q$event=='Death')
   nl <- sum(q$time_days==ts[i]&q$event=='Discharged alive')
   pdth <- pdth+p0*nd/nr; pdis <- pdis+p0*nl/nr; p0 <- p0*(1-(nd+nl)/nr); alt[i] <- pdth
  }
  err <- max(abs(alt-pp$estimate)); stopifnot(err<1e-10)
  # With complete terminal outcomes and no administrative censoring, final CIF equals observed death proportion.
  stopifnot(abs(tail(pp$estimate,1)-mean(q$event=='Death'))<1e-10)
  checks[[length(checks)+1]] <- data.frame(panel=layers[j],group=g,n=nrow(q),deaths=sum(q$event=='Death'),
    aj_recursion_max_error=err,final_CIF=tail(pp$estimate,1),observed_death_fraction=mean(q$event=='Death'))
  pp <- rbind(data.frame(time=0,estimate=0,lower=0,upper=0,discharge=0,in_hospital=1),pp)
  pp$panel <- layers[j]; pp$group <- g
  curves[[length(curves)+1]] <- pp
  risk[[length(risk)+1]] <- data.frame(panel=layers[j],group=g,time=ticks,n_risk=vapply(ticks,function(t)sum(q$time_days>=t),integer(1)))
  counts[[length(counts)+1]] <- data.frame(panel=layers[j],group=g,n=nrow(q),deaths=sum(q$event=='Death'))
  fits[[layers[j]]][[g]] <- fit
 }
}
cv <- do.call(rbind,curves); rt <- do.call(rbind,risk); cc <- do.call(rbind,counts); qa <- do.call(rbind,checks)
write.csv(cv,file.path(O,'results/CIF_full_followup.csv'),row.names=FALSE)
write.csv(rt,file.path(O,'results/risk_table.csv'),row.names=FALSE)
write.csv(cc,file.path(O,'results/panel_group_counts.csv'),row.names=FALSE)
write.csv(qa,file.path(O,'results/AJ_numerical_checks.csv'),row.names=FALSE)
saveRDS(fits,file.path(O,'results/AJ_fits.rds'))
# Preserve stepwise confidence bands: duplicate the prior ordinate at each new event time.
step_band <- function(x) {
 x <- x[x$time<=30,]; last <- x[nrow(x),]; last$time <- 30; x <- rbind(x,last)
 stopifnot(nrow(x)>=2)
 a <- x[-nrow(x),]; a$time <- x$time[-1]
 all <- rbind(x,a); all$order <- c(seq(2,2*nrow(x),2),seq(3,2*nrow(x)-1,2))
 all[order(all$order),]
}
steps <- do.call(rbind,lapply(split(cv,interaction(cv$panel,cv$group,drop=TRUE)),step_band))
max_upper <- max(steps$upper)
ymax <- if(max_upper<=.5) .5 else min(1,ceiling(max_upper*10)/10)
stopifnot(ymax>=max_upper)
yb <- if(ymax<=.5) seq(0,ymax,.1) else pretty(c(0,ymax),n=5)
panels <- list()
for (j in seq_along(layers)) {
 x <- steps[steps$panel==layers[j],]; r <- rt[rt$panel==layers[j],]
 z <- cc[cc$panel==layers[j],]
 p <- ggplot(x,aes(time,estimate,colour=group,fill=group,group=group))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA)+geom_line(linewidth=.65)+
  scale_colour_manual(values=colours,breaks=names(colours),labels=labels,drop=FALSE)+
  scale_fill_manual(values=colours,breaks=names(colours),labels=labels,drop=FALSE)+
  scale_x_continuous(breaks=ticks,limits=c(0,30),expand=expansion(mult=c(.015,.015)))+
  scale_y_continuous(breaks=yb,labels=function(x)paste0(round(100*x),'%'),limits=c(0,ymax),expand=expansion(mult=c(0,.035)))+
  labs(title=paste0(LETTERS[j],'  ',titles[j]),subtitle=sprintf('n = %s; subsequent in-hospital deaths = %s',format(sum(z$n),big.mark=','),sum(z$deaths)),
       x='Days after the 72-h landmark',y='Cumulative in-hospital mortality',colour='Day 3 mean Na (mmol/L)',fill='Day 3 mean Na (mmol/L)')+
  theme_classic(base_family='Arial',base_size=11)+
  theme(plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),legend.position='bottom',
        axis.title.y=element_text(size=10.5),axis.title.x=element_text(size=10.5),plot.margin=margin(8,12,0,5))
 r$group <- factor(r$group,levels=c('High','Normal','Low'))
 tab <- ggplot(r,aes(time,group,label=n_risk,colour=group))+
  geom_text(aes(hjust=ifelse(time==0,0,ifelse(time==30,1,.5))),size=3.1,family='Arial',show.legend=FALSE)+
  scale_colour_manual(values=colours,drop=FALSE)+
  scale_x_continuous(breaks=ticks,limits=c(0,30),expand=expansion(mult=c(.015,.015)))+
  scale_y_discrete(labels=c(High='High',Normal='Normal',Low='Low'),expand=expansion(add=.6))+
  labs(title='Number at risk',x=NULL,y=NULL)+
  theme_classic(base_family='Arial',base_size=10)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_text(size=9,margin=margin(r=10)),plot.title=element_text(size=10),plot.margin=margin(3,12,8,5))
 panels[[j]] <- (p/tab)+plot_layout(heights=c(4,1.05))
}
caption <- 'Unadjusted Aalen\u2013Johansen estimates; shaded bands: pointwise 95% CI. Discharge alive is a competing event.\nDay 3 uses the mean sodium during 48\u2013<72 h. Full follow-up is used; only the first 30 days after the landmark are shown.'
fig <- wrap_plots(panels,ncol=2,guides='collect')+
 plot_annotation(caption=caption,theme=theme(plot.caption=element_text(family='Arial',size=10,hjust=0,margin=margin(t=10)))) &
 theme(legend.position='bottom',legend.title=element_text(size=11),legend.text=element_text(size=11))
out <- file.path(O,'output'); dir.create(out,showWarnings=FALSE)
ggsave(file.path(out,'Sodium_day3_CIF_age_preview.png'),fig,width=13,height=15.5,dpi=300,device=ragg::agg_png,bg='white',limitsize=FALSE)
ggsave(file.path(out,'Sodium_day3_CIF_age_preview.pdf'),fig,width=13,height=15.5,device=cairo_pdf,bg='white',limitsize=FALSE)
ggsave(file.path(out,'Sodium_day3_CIF_age_preview.svg'),fig,width=13,height=15.5,device=svglite::svglite,bg='white',limitsize=FALSE)
ggsave(file.path(out,'Sodium_day3_CIF_age_preview.tiff'),fig,width=13,height=15.5,dpi=300,device=ragg::agg_tiff,compression='lzw',bg='white',limitsize=FALSE)
writeLines(capture.output(sessionInfo()),file.path(O,'control/R_session.txt'))
report <- list(status='NUMERICAL_AND_RENDER_EXECUTED',N=3459,deaths=306,method='Aalen-Johansen with death and live discharge',
 confidence_method='survival::survfit multistate, log-log pointwise 95%',display_days=c(0,30),risk_times=ticks,
 risk_convention='T >= t: in hospital event-free immediately before events at t',ymax=ymax,
 full_followup_max_days=max(d$time_days),groups=as.list(table(d$sodium_status)),panels=length(panels),
 independent_AJ_max_error=max(qa$aj_recursion_max_error),no_new_tests=TRUE,author_review_only=TRUE)
jsonlite::write_json(report,file.path(O,'control/render_QA.json'),pretty=TRUE,auto_unbox=TRUE)
cat(jsonlite::toJSON(report,auto_unbox=TRUE))
