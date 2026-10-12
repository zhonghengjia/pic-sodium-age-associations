cv <- read.csv(file.path(O,'source_data/CIF_full_followup.csv'))
rt <- read.csv(file.path(O,'source_data/risk_table.csv'))
counts <- read.csv(file.path(O,'results/display_counts.csv'))
layers <- c('Overall','L01','L02','L03','L04','L05')
titles <- c('Overall','Neonates (\u226428 d)','Infants (>28 d\u2013<1 y)',
            'Children 1\u20134 y','Children 5\u20139 y','Children 10\u201314 y')
groups <- c('Low','Normal','High')
colours <- c(Low='#0072B2',Normal='#009E73',High='#D55E00')
linetypes <- c(Low='dashed',Normal='solid',High='dotdash')
labels <- c('Low (<135)','Normal (135\u2013145)','High (>145)')
ticks <- c(0,7,14,21,30)
stopifnot(nrow(rt)==90,nrow(counts)==18,all(cv$lower<=cv$estimate+1e-12),all(cv$upper>=cv$estimate-1e-12))
step_band <- function(x) {
 x <- x[x$time<=30,];last <- x[nrow(x),];last$time <- 30;x <- rbind(x,last)
 a <- x[-nrow(x),];a$time <- x$time[-1]
 both <- rbind(x,a);both$order <- c(seq(2,2*nrow(x),2),seq(3,2*nrow(x)-1,2))
 both[order(both$order),]
}
steps <- do.call(rbind,lapply(split(cv,interaction(cv$panel,cv$group,drop=TRUE)),step_band))
stopifnot(max(steps$upper)<.5)
panels <- list()
for (j in seq_along(layers)) {
 x <- steps[steps$panel==layers[j],];r <- rt[rt$panel==layers[j],]
 z <- counts[counts$panel==layers[j],]
 x$group <- factor(x$group,levels=groups)
 p <- ggplot(x,aes(time,estimate,colour=group,fill=group,linetype=group,group=group))+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA,linetype=0)+geom_line(linewidth=.55)+
  scale_colour_manual(values=colours,breaks=groups,labels=labels,drop=FALSE)+
  scale_fill_manual(values=colours,breaks=groups,labels=labels,drop=FALSE)+
  scale_linetype_manual(values=linetypes,breaks=groups,labels=labels,drop=FALSE)+
  scale_x_continuous(breaks=ticks,limits=c(0,30),expand=expansion(mult=c(.015,.015)))+
  scale_y_continuous(breaks=seq(0,.5,.1),labels=function(x)round(100*x),limits=c(0,.5),expand=expansion(mult=c(0,.035)))+
  labs(title=paste0(LETTERS[j],'  ',titles[j]),subtitle=' ',
       x='Days after 72-h landmark',y='Cumulative mortality (%)',colour='Day 3 mean Na (mmol/L)',fill='Day 3 mean Na (mmol/L)',linetype='Day 3 mean Na (mmol/L)')+
  theme_classic(base_family='Arial',base_size=8)+
  theme(plot.title=element_text(face='bold',size=9),plot.subtitle=element_text(size=7.4),legend.position='bottom',
        axis.title.y=element_text(size=7.5),axis.title.x=element_text(size=7.5),axis.text=element_text(size=7.5),
        plot.margin=margin(5,3,0,3))
 r$group <- factor(r$group,levels=rev(groups))
 tab <- ggplot(r,aes(time,group,label=n_risk,colour=group))+
  geom_text(aes(hjust=ifelse(time==0,0,ifelse(time==30,1,.5))),size=2.6,family='Arial',show.legend=FALSE)+
  scale_colour_manual(values=colours,drop=FALSE)+
  scale_x_continuous(breaks=ticks,limits=c(0,30),expand=expansion(mult=c(.015,.015)))+
  scale_y_discrete(labels=c(High='High',Normal='Normal',Low='Low'),expand=expansion(add=.6))+
  labs(title='Number at risk',x=NULL,y=NULL)+theme_classic(base_family='Arial',base_size=7.5)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_text(size=7.5,margin=margin(r=4)),plot.title=element_text(size=7.5),plot.margin=margin(4,3,5,3))
 panels[[j]] <- (p/tab)+plot_layout(heights=c(4,1.15))
}
caption <- paste0('Unadjusted Aalen\u2013Johansen estimates with pointwise 95% CI; discharge alive is a competing event.\n',
 'Day 30 is counted from the 72-h landmark. Full follow-up is used; only days 0\u201330 are displayed.\n',
 'Before the first death, zero-width intervals do not establish zero underlying risk. Late risk sets can be small.')
fig <- wrap_plots(panels,ncol=3,guides='collect')+
 plot_annotation(title=NULL,caption=caption,
 theme=theme(plot.title=element_text(family='Arial',size=11,face='bold'),
             plot.caption=element_text(family='Arial',size=7.3,hjust=0,margin=margin(t=5)))) &
 theme(legend.position='bottom',legend.title=element_text(size=7.5),legend.text=element_text(size=7.5),
       legend.key.width=grid::unit(10,'pt'),legend.margin=margin(0,0,0,0))
