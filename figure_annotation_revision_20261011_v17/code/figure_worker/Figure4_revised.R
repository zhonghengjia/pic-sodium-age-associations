ids<-c("ALL","L01","L02","L03","L04","L05")
titles<-c("Overall","Neonates (<=28 days)","Infants (>28 days to <1 year)","Children 1-4 years","Children 5-9 years","Children 10-14 years")
d<-data.frame(sodium_group=factor(c("Low","Normal","High"),levels=c("Low","Normal","High")))
tests<-read.csv(file.path(P,"logrank_tests.csv"))
ccall<-read.csv(file.path(P,"km_estimates_full.csv"))
rrall<-read.csv(file.path(P,"alive_numbers.csv"))
nnall<-read.csv(file.path(P,"group_counts_overall_age.csv"))
curves<-lapply(ids,function(k)subset(ccall,panel==k))
risks<-lapply(ids,function(k)subset(rrall,panel==k))
counts<-lapply(ids,function(k)subset(nnall,panel==k))
counts<-lapply(counts,function(x)x[match(c("Low","Normal","High"),x$group),])
pformat<-function(p)if(p<.001)"<0.001" else paste0("= ",sprintf("%.3f",p))
panels<-list()
render_checks<-list()
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
    geom_step(linewidth=.75,direction="hv")+
    geom_point(data=subset(cc,n_censor>0 & time<=30),shape=3,size=1.1,stroke=.3,show.legend=FALSE)+
    scale_linetype_manual(values=c(Low="dashed",Normal="solid",High="dotdash"),labels=c("Low (<135)","Normal (135-145)","High (>145)"),name="Day-3 sodium (mmol/L)")+
    scale_colour_manual(values=c(Low="#0072B2",Normal="#009E73",High="#D55E00"),labels=c("Low (<135)","Normal (135-145)","High (>145)"),name="Day-3 sodium (mmol/L)")+
    scale_fill_manual(values=c(Low="#0072B2",Normal="#009E73",High="#D55E00"),guide="none")+
    scale_x_continuous(breaks=c(0,7,14,21,30),expand=c(0,0))+
    scale_y_continuous(breaks=c(0,.25,.5,.75,1),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0))+
    coord_cartesian(xlim=c(0,30),ylim=c(0,1),clip="on")+
    labs(title=paste0(LETTERS[i],"  ",titles[i]),subtitle=" \n ",
         x="Days after the 72-hour landmark",y="KM survival estimate")+
    theme_classic(base_size=10.5,base_family="Arial")+
    theme(plot.title=element_text(size=12,face="bold"),plot.subtitle=element_text(size=9,lineheight=1.1),
          legend.position="bottom",plot.margin=margin(5,10,0,6),axis.text=element_text(colour="black"))
  stopifnot(identical(p$data$time,cc$time),identical(p$data$surv,cc$surv),identical(p$data$lower,cc$lower),identical(p$data$upper,cc$upper))
  render_checks[[i]]<-data.frame(panel=k,curve_estimates_unchanged=TRUE,ribbon_layers=sum(vapply(p$layers,function(l)inherits(l$geom,"GeomRibbon"),logical(1))),step_layers=sum(vapply(p$layers,function(l)inherits(l$geom,"GeomStep"),logical(1))))
  rr$group<-factor(rr$group,levels=rev(levels(d$sodium_group)))
  rowlabs<-setNames(as.character(nn$group),nn$group)
  rowlabel_data<-data.frame(time=-1.2,group=factor(nn$group,levels=rev(levels(d$sodium_group))),label=unname(rowlabs[nn$group]))

  rt<-ggplot(rr,aes(time,group,label=n_alive))+geom_text(aes(hjust=ifelse(time==0,0,ifelse(time==30,1,.5))),size=3.1,family="Arial")+
    scale_x_continuous(breaks=c(0,7,14,21,30),expand=c(0,0))+
    geom_text(data=rowlabel_data,mapping=aes(x=time,y=group,label=label,colour=group),inherit.aes=FALSE,hjust=1,size=3.1,family="Arial",show.legend=FALSE)+
    scale_colour_manual(values=c(Low="#0072B2",Normal="#009E73",High="#D55E00"),guide="none")+
    scale_y_discrete(labels=rowlabs,drop=FALSE)+
    labs(title="Number alive",x=NULL,y=NULL)+
    coord_cartesian(xlim=c(0,30),clip="off")+
    theme_classic(base_size=10.5,base_family="Arial")+
    theme(axis.line.y=element_blank(),axis.line.x=element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text=element_text(colour="black"),
          plot.title=element_text(size=9),plot.margin=margin(2,10,2,6),axis.text.y=element_text(size=9,colour="transparent",margin=margin(r=10)))
  label_check<-ggplot_build(rt)$data[[2]]
  stopifnot(nrow(label_check)==3,!anyNA(label_check$x),identical(label_check$colour,c("#0072B2","#009E73","#D55E00")))
  render_checks[[i]]$risk_labels_correctly_colored<-TRUE
  panels[[i]]<-p/rt+plot_layout(heights=c(3.2,1.1))
}
fig<-wrap_plots(panels,ncol=3,guides="collect")+
  plot_annotation(title=NULL,
    subtitle=NULL,
    caption=NULL,
    theme=theme(text=element_text(family="Arial"),plot.title=element_text(size=17,face="bold"),
                plot.subtitle=element_text(size=11),plot.caption=element_text(size=9.2,hjust=0)))&
  theme(legend.position="bottom")
