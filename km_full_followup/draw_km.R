args<-commandArgs(trailingOnly=TRUE)
O<-normalizePath(args[1],winslash="/",mustWork=TRUE)
library(survival);library(ggplot2);library(patchwork)
d<-read.csv(file.path(O,"restricted/km_input.csv"))
d$sodium_group<-factor(d$sodium_group,levels=c("Low","Normal","High"))
ids<-c("ALL","L01","L02","L03","L04","L05")
titles<-c("Overall","Neonates (<=28 days)","Infants (>28 days to <1 year)","Children 1-4 years","Children 5-9 years","Children 10-14 years")
fits<-list();test<-list();curves<-list();risks<-list();counts<-list()
for(i in seq_along(ids)){
  k<-ids[i];z<-if(k=="ALL")d else subset(d,age_layer==k)
  fit<-survfit(Surv(fup_days,death)~sodium_group,data=z,conf.type="log-log")
  lr<-survdiff(Surv(fup_days,death)~sodium_group,data=z,rho=0)
  fits[[k]]<-fit
  test[[i]]<-data.frame(panel=k,n=nrow(z),deaths=sum(z$death),chisq=unname(lr$chisq),df=length(lr$n)-1,p_raw=pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE))
  sm<-summary(fit,censored=TRUE)
  cc<-data.frame(panel=k,group=sub("sodium_group=","",as.character(sm$strata)),
                 time=sm$time,surv=sm$surv,lower=sm$lower,upper=sm$upper,
                 n_risk=sm$n.risk,n_event=sm$n.event,n_censor=sm$n.censor)
  for(g in levels(d$sodium_group)){
    cc<-rbind(data.frame(panel=k,group=g,time=0,surv=1,lower=1,upper=1,n_risk=sum(z$sodium_group==g),n_event=0,n_censor=0),cc)
  }
  cc<-cc[order(cc$group,cc$time),];curves[[i]]<-cc
  rr<-expand.grid(group=levels(d$sodium_group),time=c(0,7,14,21,30),stringsAsFactors=FALSE)
  rr$n_risk<-mapply(function(g,t)sum(z$sodium_group==g & z$fup_days>=t),rr$group,rr$time)
  rr$panel<-k;risks[[i]]<-rr
  counts[[i]]<-data.frame(panel=k,group=levels(d$sodium_group),n=as.numeric(table(z$sodium_group)),
                         deaths=sapply(levels(d$sodium_group),function(g)sum(z$death[z$sodium_group==g])))
}
tests<-do.call(rbind,test);tests$p_holm_age<-NA_real_
tests$p_holm_age[-1]<-p.adjust(tests$p_raw[-1],method="holm")
write.csv(tests,file.path(O,"logrank_tests.csv"),row.names=FALSE)
write.csv(do.call(rbind,curves),file.path(O,"km_estimates_full.csv"),row.names=FALSE)
write.csv(do.call(rbind,risks),file.path(O,"risk_numbers.csv"),row.names=FALSE)
write.csv(do.call(rbind,counts),file.path(O,"group_counts_overall_age.csv"),row.names=FALSE)
pformat<-function(p)if(p<.001)"<0.001" else paste0("= ",sprintf("%.3f",p))
panels<-list()
for(i in seq_along(ids)){
  k<-ids[i];cc<-curves[[i]];rr<-risks[[i]];nn<-counts[[i]];tr<-tests[i,]
  cc$group<-factor(cc$group,levels=levels(d$sodium_group))
  # Step confidence intervals change only at observed follow-up times.
  stepped<-do.call(rbind,lapply(split(cc,cc$group),function(a){
    a<-a[order(a$time),]
    s<-a[rep(seq_len(nrow(a)),each=2),]
    s$time<-as.vector(rbind(a$time,c(a$time[-1],tail(a$time,1))))
    s
  }))
  ptext<-paste0("Three-group log-rank P ",pformat(tr$p_raw))
  if(k!="ALL")ptext<-paste0(ptext,"; Holm P ",pformat(tr$p_holm_age))
  p<-ggplot(cc,aes(time,surv,linetype=group,colour=group))+
    geom_ribbon(data=subset(stepped,is.finite(lower)&is.finite(upper)),mapping=aes(x=time,ymin=lower,ymax=upper,group=group),fill="grey45",alpha=.11,colour=NA,inherit.aes=FALSE)+
    geom_step(linewidth=.75,direction="hv")+
    geom_point(data=subset(cc,n_censor>0 & time<=30),shape=3,size=1.1,stroke=.3,show.legend=FALSE)+
    scale_linetype_manual(values=c(Low="dashed",Normal="solid",High="dotdash"),labels=c("Low (<135)","Normal (135-145)","High (>145)"),name="Day-3 mean sodium (mmol/L)")+
    scale_colour_manual(values=c(Low="black",Normal="black",High="grey40"),labels=c("Low (<135)","Normal (135-145)","High (>145)"),name="Day-3 mean sodium (mmol/L)")+
    scale_x_continuous(breaks=c(0,7,14,21,30),expand=c(0,0))+
    scale_y_continuous(breaks=c(0,.25,.5,.75,1),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0))+
    coord_cartesian(xlim=c(0,30),ylim=c(0,1),clip="on")+
    labs(title=paste0(LETTERS[i],"  ",titles[i]),subtitle=paste0("N = ",format(tr$n,big.mark=","),"; in-hospital deaths = ",tr$deaths,"\n",ptext),
         x=NULL,y="KM survival estimate")+
    theme_classic(base_size=10.5,base_family="Arial")+
    theme(plot.title=element_text(size=12,face="bold"),plot.subtitle=element_text(size=9,lineheight=1.1),
          legend.position="bottom",plot.margin=margin(5,10,0,6),axis.text=element_text(colour="black"))
  rr$group<-factor(rr$group,levels=rev(levels(d$sodium_group)))
  rowlabs<-setNames(paste0(nn$group," (d=",nn$deaths,")"),nn$group)
  rt<-ggplot(rr,aes(time,group,label=n_risk))+geom_text(aes(hjust=ifelse(time==0,0,ifelse(time==30,1,.5))),size=3.1,family="Arial")+
    scale_x_continuous(breaks=c(0,7,14,21,30),limits=c(0,30),expand=c(0,0))+
    scale_y_discrete(labels=rowlabs,drop=FALSE)+
    labs(title="Number at risk (d = total deaths)",x="Days after the 72-hour landmark",y=NULL)+
    coord_cartesian(clip="off")+
    theme_classic(base_size=10.5,base_family="Arial")+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(colour="black"),
          plot.title=element_text(size=9),plot.margin=margin(2,10,2,6),axis.text.y=element_text(size=9,margin=margin(r=10)))
  panels[[i]]<-p/rt+plot_layout(heights=c(3.2,1.1))
}
fig<-wrap_plots(panels,ncol=3,guides="collect")+
  plot_annotation(title="Day-3 sodium and subsequent in-hospital survival across ages",
    subtitle="Unadjusted Kaplan-Meier estimates in children alive and hospitalized at 72 hours",
    caption=paste0("Low: <135; normal: 135-145; high: >145 mmol/L. Categories use the accepted day-3 arithmetic daily mean.\n",
      "Time zero: first ICU admission + 72 hours. Recorded in-hospital death is the event; alive discharge is censored.\n",
      "Shading: pointwise 95% CIs; +: censoring. Curves show days 0-30; log-rank tests use the complete recorded hospital follow-up.\n",
      "Holm correction covers the five age-stratum log-rank tests. These unadjusted comparisons do not test sodium-by-age interaction."),
    theme=theme(text=element_text(family="Arial"),plot.title=element_text(size=17,face="bold"),
                plot.subtitle=element_text(size=11),plot.caption=element_text(size=9.2,hjust=0)))&
  theme(legend.position="bottom")
ragg::agg_png(file.path(O,"KM_day3_sodium_age_review.png"),width=15.5,height=11.8,units="in",res=260)
print(fig);dev.off()
cairo_pdf(file.path(O,"KM_day3_sodium_age_review.pdf"),width=15.5,height=11.8,family="Arial")
print(fig);dev.off()
writeLines(capture.output(sessionInfo()),file.path(O,"R_session_info.txt"))
print(tests,row.names=FALSE)
