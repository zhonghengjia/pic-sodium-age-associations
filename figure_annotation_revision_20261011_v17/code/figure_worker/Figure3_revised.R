dc<-read.csv(paths[['full_curves']],stringsAsFactors=FALSE);rug<-read.csv(paths[['rug']],stringsAsFactors=FALSE)
jt<-read.csv(paths[['joint_tests']],stringsAsFactors=FALSE);sp<-read.csv(paths[['support']],stringsAsFactors=FALSE)
review<-fromJSON(paths[['accepted_review']],simplifyVector=FALSE)
stopifnot(review$status=='ACCEPT',review$full_curves_sha256==before[['full_curves']],review$source_hashes$rug==before[['rug']],
 nrow(dc)==7267,nrow(rug)==4334,all(dc$ci_low>=.01),all(dc$ci_high<=1000))
layers<-c('ALL',paste0('L0',1:5));titles<-c('Overall','Neonates (\u226428 d)','Infants (>28 d-<1 y)','Children 1-4 y','Children 5-9 y','Children 10-14 y')
for(i in seq_len(nrow(sp))){z<-dc[dc$day==sp$day[i]&dc$age_layer==sp$age_layer[i],];r<-rug[rug$day==sp$day[i]&rug$age_layer==sp$age_layer[i],];stopifnot(min(z$sodium)==sp$min[i],max(z$sodium)==sp$max[i],sum(r$n)==sp$N[i])}
hist_source<-rug
hist_source$bin_left<-2*floor(hist_source$sodium/2)
# Include the exact upper display boundary in the final visible interval [158,160].
hist_source$bin_left[hist_source$sodium==160]<-158
histogram<-aggregate(n~day+age_layer+bin_left,hist_source,sum)
histogram$bin_right<-histogram$bin_left+2
stopifnot(all(histogram$n==floor(histogram$n)),all(histogram$n>0))
for(i in seq_len(nrow(sp))){h<-histogram[histogram$day==sp$day[i]&histogram$age_layer==sp$age_layer[i],];stopifnot(sum(h$n)==sp$N[i])}
ref<-dc[dc$sodium==140,];stopifnot(nrow(ref)==18,all(ref$OR==1))
pp<-vapply(1:3,function(j){z<-jt[jt$model_id==paste0('daily',j,'_age')&jt$test_id=='daily_age_4df',];stopifnot(nrow(z)==1);z$p_holm},0)
stopifnot(identical(sprintf('%.3f',pp),c('0.022','0.005','0.002')))
caption<-paste(readLines(paths[['caption']],warn=FALSE,encoding='UTF-8'),collapse='\n')

cfg<-list(x_view=c(120,160),x_breaks=seq(120,160,10),y_limits=c(.1,10),y_scale='standard_log10',coordinate_zoom_only=TRUE,data_filter=FALSE,
 width_mm=210,height_mm=285,dpi=600,font='Arial',minimum_text_pt=7.5,coral='#FF6262',CI_fill='#FADBD5',teal='#009E9C',histogram_blue='#65AEEF',bin_width=2,histogram_log10_height=1,right_axis_title='Count',right_axis_minimum_tick_count=3,right_axis_ticks='round integers',right_axis_tick_angle=90)
clipping<-do.call(rbind,lapply(layers,function(l)do.call(rbind,lapply(1:3,function(j){z<-dc[dc$age_layer==l&dc$day==j,];z<-z[order(z$sodium),];inside<-z$sodium>=120&z$sodium<=160;bound<-c(120,160);bound<-bound[bound>=min(z$sodium)&bound<=max(z$sodium)];upper<-c(z$ci_high[inside],approx(z$sodium,z$ci_high,xout=bound)$y);lower<-c(z$ci_low[inside],approx(z$sodium,z$ci_low,xout=bound)$y);estimate<-c(z$OR[inside],approx(z$sodium,z$OR,xout=bound)$y);data.frame(age_layer=l,day=j,curve_above10=any(estimate>10),CI_above10=any(upper>10),curve_below0_1=any(estimate<.1),CI_below0_1=any(lower<.1))}))))
stopifnot(all(dc$OR>0),all(dc$ci_low>0),all(dc$ci_high>0))
round_count_ticks<-function(peak){top<-2*peak;steps<-c(1,2,5,10)*10^floor(log10(top/2));step<-max(steps[steps<=top/2]);z<-seq(0,floor(top/step)*step,by=step);stopifnot(length(z)>=3,all(z==floor(z)));z}
stopifnot(identical(round_count_ticks(591),c(0,500,1000)))
panel<-function(l,j){z<-dc[dc$age_layer==l&dc$day==j,];r<-rug[rug$age_layer==l&rug$day==j,];q<-ref[ref$age_layer==l&ref$day==j,]
 h<-histogram[histogram$age_layer==l&histogram$day==j,];visible_h<-h[h$bin_left>=120&h$bin_right<=160,]
 count_max<-max(visible_h$n);stopifnot(count_max>=max(visible_h$n),count_max==floor(count_max))
 h$y0<-.1;h$y1<-.1*10^(h$n/count_max)
 stopifnot(abs(max(h$y1[h$bin_left>=120&h$bin_right<=160])-1)<1e-12)
 g<-ggplot(z,aes(sodium,OR))+geom_ribbon(aes(ymin=ci_low,ymax=ci_high),fill=cfg$CI_fill)+
 geom_rect(data=h,aes(xmin=bin_left,xmax=bin_right,ymin=y0,ymax=y1),inherit.aes=FALSE,fill=cfg$histogram_blue,colour='white',linewidth=.1,alpha=.65)+geom_line(colour=cfg$coral,linewidth=.45,linetype='solid')+
 geom_hline(yintercept=1,linetype='dashed',colour=cfg$teal,linewidth=.3)+geom_vline(xintercept=135,linetype='dotted',colour='#717171',linewidth=.3)+
 geom_point(data=q,shape=21,fill='white',colour=cfg$teal,size=1.4,stroke=.4)+
 scale_x_continuous(breaks=cfg$x_breaks,expand=expansion(mult=0))+
 scale_y_log10(breaks=c(.1,.25,.5,1,2,4,10),labels=c('0.1','0.25','0.5','1','2','4','10'),expand=expansion(mult=0),sec.axis=sec_axis(~log10(./.1)*count_max,name='Count',breaks=round_count_ticks(count_max),labels=format(round_count_ticks(count_max),trim=TRUE,scientific=FALSE)))+
 coord_cartesian(xlim=cfg$x_view,ylim=cfg$y_limits,expand=FALSE)+
 labs(title=titles[match(l,layers)],x='Sodium (mmol/L)',y='Adjusted OR (log scale)')+theme_classic(base_size=8,base_family='Arial')+
 theme(plot.title=element_text(size=7.5,hjust=.5,face='plain'),axis.text=element_text(size=7.5,colour='black'),axis.title=element_text(size=8),
 axis.title.y=element_text(angle=90),axis.title.y.right=element_text(angle=90,colour='#3976A7',size=7.5),axis.text.y.right=element_text(colour='#3976A7',size=7.5,angle=90,hjust=.5,vjust=.5),axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),axis.ticks.length=grid::unit(.8,'mm'),plot.margin=margin(1.5,1.5,1.5,1.5,'mm'))
 flag<-clipping[clipping$age_layer==l&clipping$day==j,];if(flag$curve_above10||flag$CI_above10) g<-g+annotate('text',x=159,y=9.8,label='\u2191',colour=cfg$coral,size=3.5,hjust=1,vjust=1)
 if(flag$curve_below0_1||flag$CI_below0_1) g<-g+annotate('text',x=159,y=.105,label='\u2193',colour=cfg$coral,size=3.5,hjust=1,vjust=0)
 g
}
plots<-unlist(lapply(layers,function(l)lapply(1:3,function(j)panel(l,j))),recursive=FALSE)
built<-lapply(plots,ggplot_build)
stopifnot(length(built)==18,all(vapply(built,function(z)isTRUE(all.equal(z$layout$panel_params[[1]]$x.range,c(120,160))),logical(1))))
stopifnot(sum(vapply(built,function(z)nrow(z$data[[1]]),integer(1)))==7267,all(vapply(built,function(z)isTRUE(all.equal(z$layout$panel_params[[1]]$y.range,c(-1,1))),logical(1))))
stopifnot(all(vapply(built,function(z)length(z$layout$panel_params[[1]]$y.sec$get_labels())>=3,logical(1))),all(vapply(plots,function(z)z$theme$axis.text.y.right$angle==90,logical(1))))
headers<-lapply(1:3,function(j)ggplot()+annotate('text',x=.5,y=.5,label=paste0(LETTERS[j],'  Day ',j,'\n '),family='Arial',fontface='bold',size=7.5/ggplot2::.pt,lineheight=1.05)+xlim(0,1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(0,0,0,0)))
p<-wrap_plots(headers,ncol=3)/wrap_plots(plots,ncol=3)+plot_layout(heights=c(.38,6))
