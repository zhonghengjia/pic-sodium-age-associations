# Publication refinement of existing figures only. No model fitting.
# Layout references: Lex 2014 Fig4; BaSICS 2021 Fig4; Peltan 2016 Fig2.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root <- normalizePath(args[1],winslash='/',mustWork=TRUE)
options(device=function(...) grDevices::cairo_pdf(tempfile(fileext='.pdf'),...))
project <- dirname(root)
lib <- file.path(dirname(project),'landmark72_reanalysis_20261006/runtime/R-library')
if(dir.exists(lib)) .libPaths(c(lib,.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(svglite);library(jsonlite);library(digest)})
# Measure on the same Cairo font backend as the submitted PDF. The global helper is read-only.
alignment_source <- readLines(Sys.getenv("PIC_PANEL_ALIGNMENT_R"),warn=FALSE)
alignment_source <- sub('grDevices::pdf(probe_path, width = width_in, height = height_in, useDingbats = FALSE)',
 'grDevices::cairo_pdf(probe_path, width = width_in, height = height_in)',alignment_source,fixed=TRUE)
eval(parse(text=alignment_source))
# Patchwork 1.3.2 retains unresolved null units and fixed-aspect panel names.
# Measure actual rendered viewports, rather than treating their null widths as zero.
write_patchwork_panel_layout <- function(plot,manifest_path,width_in,height_in,panel_ids,row_groups=NULL,column_groups=NULL){
 g <- patchwork::patchworkGrob(plot)
 cells <- g$layout[startsWith(g$layout$name,'panel') & g$layout$name!='panel-area',,drop=FALSE]
 cells <- cells[order(cells$t,cells$l),,drop=FALSE]
 stopifnot(nrow(cells)==length(panel_ids))
 probe <- tempfile(fileext='.pdf');grDevices::cairo_pdf(probe,width=width_in,height=height_in)
 grid::grid.newpage();grid::grid.draw(g);grid::grid.force()
 vps <- grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name
 panels <- lapply(seq_len(nrow(cells)),function(i){
  r <- cells[i,];vp <- sprintf('%s.%d-%d-%d-%d',r$name,r$t,r$l,r$b,r$r)
  stopifnot(vp %in% vps);grid::seekViewport(vp)
  lo <- grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)
  hi <- grid::deviceLoc(grid::unit(1,'npc'),grid::unit(1,'npc'),valueOnly=TRUE)
  box <- c(lo$x,lo$y,hi$x,hi$y)*72;stopifnot(box[3]>box[1],box[4]>box[2])
  grid::upViewport(0)
  list(id=panel_ids[i],bbox_pt=unname(box),grid_id='patchwork-grid-1',row_start=as.integer(r$t-1),row_stop=as.integer(r$b),col_start=as.integer(r$l-1),col_stop=as.integer(r$r))
 })
 grDevices::dev.off();unlink(probe)
 exemptions <- list()
 if(length(panel_ids)>2){
  exemptions <- list(list(panels=panel_ids,checks='panel-width',reason='Intentional unequal semantic columns: status matrix/count/deaths or labels/counts/numeric estimates/forest. Corresponding columns between A and B have identical width; row/column alignment is not exempted.'))
 }
 manifest <- list(schema_version=1L,backend='r-patchwork',measurement='Cairo rendered panel viewport deviceLoc',figure=list(width_pt=width_in*72,height_pt=height_in*72),panels=panels,exemptions=exemptions)
 if(!is.null(row_groups))manifest$row_groups <- row_groups
 if(!is.null(column_groups))manifest$column_groups <- column_groups
 jsonlite::write_json(manifest,manifest_path,auto_unbox=TRUE,pretty=TRUE,digits=NA)
}
src <- file.path(root,'source_data');out <- file.path(root,'figures');qa <- file.path(root,'reviews')
paths <- list.files(src,full.names=TRUE);before <- vapply(paths,function(p)digest(file=p,algo='sha256'),'')
histories <- read.csv(file.path(src,'histories.csv'),colClasses=c(history='character'))
assoc <- read.csv(file.path(src,'associations.csv'))
tests <- read.csv(file.path(src,'joint_tests.csv'))
bins <- read.csv(file.path(src,'calibration_bins.csv'));perf <- read.csv(file.path(src,'performance.csv'))
blue <- '#377EB8';coral <- '#E96752';ink <- '#252525';gray <- '#E7E7E7'
font <- 'Arial';pt <- 8.5;txt <- pt/ggplot2::.pt
base <- theme_classic(base_size=pt,base_family=font)+theme(axis.text=element_text(colour=ink,size=pt),axis.title=element_text(colour=ink,size=pt),plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=8),plot.margin=margin(4,4,4,4))
tabtheme <- theme_void(base_size=pt,base_family=font)+theme(plot.margin=margin(4,2,4,2))
writefig <- function(g,i,w,h,ids=NULL,rg=NULL,cg=NULL){
 stem <- file.path(out,paste0('Supplementary_Figure',i))
 ggsave(paste0(stem,'.png'),g,width=w,height=h,units='in',dpi=600,device=ragg::agg_png,bg='white')
 ggsave(paste0(stem,'.pdf'),g,width=w,height=h,units='in',device=grDevices::cairo_pdf,bg='white')
 svglite::svglite(paste0(stem,'.svg'),width=w,height=h,bg='white');print(g);grDevices::dev.off()
 ggsave(paste0(stem,'.tiff'),g,width=w,height=h,units='in',dpi=600,device=ragg::agg_tiff,compression='lzw',bg='white')
 write_patchwork_panel_layout(g,file.path(qa,paste0('S',i,'_panel_layout.json')),w,h,panel_ids=ids,row_groups=rg,column_groups=cg)
}

# S1: adapt the existing three-column status matrix and horizontal cardinality bars.
stopifnot(nrow(histories)==16)
historyplots <- list()
for(k in seq_along(c('Low sodium','High sodium'))){
 typ <- c('Low sodium','High sodium')[k];z <- histories[histories$type==typ,];z <- z[order(z$history),]
 stopifnot(identical(z$history,sprintf('%03d',c(0,1,10,11,100,101,110,111))),sum(z$n)==3459,sum(z$events)==306)
 stopifnot(max(abs(z$crude_mortality_pct-100*z$events/z$n))<1e-10)
 z$y <- 8:1;col <- c(blue,coral)[k]
 griddata <- do.call(rbind,lapply(1:8,function(r)data.frame(y=z$y[r],day=1:3,state=strsplit(z$history[r],'')[[1]]=='1')))
 matrixplot <- ggplot(griddata,aes(day,y))+geom_tile(data=data.frame(day=2,y=c(2,4,6,8)),width=3.18,height=.9,fill='#F3F3F3')+
  geom_point(aes(fill=state),shape=21,size=3.4,stroke=.35,colour='#A0A0A0')+
  scale_fill_manual(values=c('FALSE'='white','TRUE'=col),guide='none')+
  scale_x_continuous(limits=c(.4,3.6),breaks=1:3,labels=c('1','2','3'),expand=c(0,0))+
  scale_y_continuous(limits=c(.4,9.5),breaks=z$y,labels=z$history,expand=c(0,0))+
  annotate('text',x=2,y=9.15,label='Day',size=txt,family=font,fontface='bold')+
  labs(x=NULL,y=NULL)+base+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=pt),plot.margin=margin(4,2,4,2))
 bars <- ggplot(z,aes(n,y))+geom_col(orientation='y',width=.52,fill=col,alpha=.78)+
  geom_text(aes(label=format(n,big.mark=',',trim=TRUE)),hjust=-.15,size=txt,family=font)+
  annotate('text',x=1850,y=9.15,label='Children, n',size=txt,family=font,fontface='bold')+
  scale_x_continuous(limits=c(0,3700),breaks=c(0,1000,2000,3000),labels=c('0','1,000','2,000','3,000'),expand=c(0,0))+
  scale_y_continuous(limits=c(.4,9.5),expand=c(0,0))+labs(x=NULL,y=NULL)+base+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(4,2,4,2))
 death <- ggplot(z,aes(.02,y))+geom_text(aes(label=sprintf('%d/%s (%.1f%%)',events,format(n,big.mark=',',trim=TRUE),crude_mortality_pct)),hjust=0,size=txt,family=font)+
  annotate('text',x=.02,y=9.15,label='Deaths/N (%)',hjust=0,size=txt,fontface='bold',family=font)+
  scale_x_continuous(limits=c(0,1),expand=c(0,0))+scale_y_continuous(limits=c(.4,9.5),expand=c(0,0))+tabtheme
 historyplots <- c(historyplots,list(matrixplot,bars,death))
}
g1 <- wrap_plots(historyplots,ncol=3,widths=c(1.6,4.4,2.1))+plot_annotation(title='Observed three-day sodium histories',subtitle='A  Low sodium (<135 mmol/L)                         B  High sodium (>145 mmol/L)',caption='Filled circles: threshold met; open circles: threshold not met. Digits follow days 1, 2 and 3.',theme=theme(plot.title=element_text(size=10,face='bold',family=font),plot.subtitle=element_blank(),plot.caption=element_text(size=8,hjust=0,family=font)))
# Panel labels belong to each block; no connecting line across independent daily states.
historyplots[[1]] <- historyplots[[1]]+labs(title='A  Low sodium (<135 mmol/L)')+theme(plot.title=element_text(size=9,face='bold',hjust=0))
historyplots[[4]] <- historyplots[[4]]+labs(title='B  High sodium (>145 mmol/L)')+theme(plot.title=element_text(size=9,face='bold',hjust=0))
g1 <- wrap_plots(historyplots,ncol=3,widths=c(1.6,4.4,2.1))+plot_annotation(caption='Filled circles: threshold met; open circles: threshold not met. Digits follow days 1, 2 and 3.',theme=theme(plot.caption=element_text(size=8,hjust=0,family=font)))
writefig(g1,1,183/25.4,164/25.4,ids=c('A_states','A_counts','A_deaths','B_states','B_counts','B_deaths'),rg=list(c('A_states','A_counts','A_deaths'),c('B_states','B_counts','B_deaths')),cg=list(c('A_states','B_states'),c('A_counts','B_counts'),c('A_deaths','B_deaths')))
write.csv(histories,file.path(qa,'S1_displayed_values.csv'),row.names=FALSE)

# S2: the existing estimates arranged as a clinical results table plus forest plot.
rows <- c('ALL',paste0('L0',1:5)); labels <- c('Overall','Neonates (<=28 d)','Infants (>28 d-<1 y)','Children 1-4 y','Children 5-9 y','Children 10-14 y')
forests <- list();displayed <- list()
for(k in 1:2){
 con <- c('L_vs_N','H_vs_N')[k];z <- assoc[assoc$contrast==con & ((assoc$model_id=='day3_common' & assoc$age_layer=='ALL')|assoc$model_id=='day3_age'),]
 z <- z[match(rows,z$age_layer),];stopifnot(nrow(z)==6,all(z$N==3459),all(z$events==306),all(z$ci_low>0),all(z$ci_high<20),all(z$OR>=z$ci_low & z$OR<=z$ci_high))
 z$y <- c(7,5,4,3,2,1);z$label <- labels;z$shape <- ifelse(z$age_layer=='ALL','common','age');displayed[[k]] <- z
 col <- c(blue,coral)[k];p <- tests$p_holm[tests$test_id==paste0(c('L','H')[k],'_vs_N_age2df')];stopifnot(length(p)==1)
 blank <- ggplot()+scale_x_continuous(limits=c(0,1),expand=c(0,0))+scale_y_continuous(limits=c(.4,8.8),expand=c(0,0))+tabtheme
 column <- function(vals,heading,hjust=0){column_data <- data.frame(y=z$y,label=vals);blank+geom_text(data=column_data,aes(x=.02,y=y,label=label),hjust=hjust,size=txt,family=font)+annotate('text',x=.02,y=8.2,label=heading,hjust=0,fontface='bold',size=txt,family=font)+geom_hline(yintercept=7.5,linewidth=.25,colour='#787878')}
 lab <- column(z$label,'Population')+annotate('text',x=.02,y=6.1,label='Representative age',hjust=0,fontface='italic',size=7.8/ggplot2::.pt,family=font)
 targets <- column(sprintf('%d/%s',z$target_deaths,format(z$target_N,big.mark=',',trim=TRUE)),paste0(c('Low','High')[k],'\nDeaths/N'))
 refs <- column(sprintf('%d/%s',z$reference_deaths,format(z$reference_N,big.mark=',',trim=TRUE)),'Normal\nDeaths/N')
 nums <- column(sprintf('%.2f (%.2f-%.2f)',z$OR,z$ci_low,z$ci_high),'Adjusted OR\n(95% CI)')
 forest <- ggplot(z,aes(OR,y))+annotate('segment',x=1,xend=1,y=.5,yend=7.45,linetype='dashed',linewidth=.35,colour='#7A7A7A')+
  geom_segment(aes(x=ci_low,xend=ci_high,yend=y),linewidth=.55,colour=col)+geom_point(aes(shape=shape),size=2.4,colour=col)+
  scale_shape_manual(values=c(common=15,age=16),guide='none')+scale_x_log10(limits=c(.2,20),breaks=c(.2,.5,1,2,5,10,20),labels=c('0.2','0.5','1','2','5','10','20'),expand=expansion(mult=.035))+
  scale_y_continuous(limits=c(.4,8.8),expand=c(0,0))+annotate('text',x=2,y=8.2,label='Adjusted OR',fontface='bold',size=txt,family=font)+geom_hline(yintercept=7.5,linewidth=.25,colour='#787878')+
  labs(x=NULL,y=NULL)+base+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.text.x=element_text(size=8),plot.margin=margin(4,2,4,2))
 lab <- lab+labs(title=paste0(c('A  Low','B  High')[k],' vs normal sodium'),subtitle=sprintf('Age interaction: Holm P = %.3f',p))+theme(plot.title=element_text(size=10,face='bold',family=font),plot.subtitle=element_text(size=8,family=font))
 forests <- c(forests,list(lab,targets,refs,nums,forest))
}
g2 <- wrap_plots(forests,ncol=5,widths=c(2.55,1.4,1.5,2.65,3.3))+plot_annotation(caption='Square: common association. Circles: median-age estimates from one continuous-age interaction model.',theme=theme(plot.caption=element_text(size=8,family=font,hjust=0)))
writefig(g2,2,183/25.4,162/25.4,ids=c('A_labels','A_target','A_normal','A_numbers','A_forest','B_labels','B_target','B_normal','B_numbers','B_forest'),rg=list(c('A_labels','A_target','A_normal','A_numbers','A_forest'),c('B_labels','B_target','B_normal','B_numbers','B_forest')),cg=list(c('A_labels','B_labels'),c('A_target','B_target'),c('A_normal','B_normal'),c('A_numbers','B_numbers'),c('A_forest','B_forest')))
write.csv(do.call(rbind,displayed),file.path(qa,'S2_displayed_estimates.csv'),row.names=FALSE)

# S3: retain existing ten-group observations, Wilson intervals and calibration summaries.
stopifnot(nrow(bins)==30)
calplots <- list();shown <- bins[bins$model_id %in% c('D1','D3'),]
for(k in 1:2){
 m <- c('D1','D3')[k];z <- shown[shown$model_id==m,];z <- z[order(z$mean_pred),];q <- perf[perf$model_id==m,]
 stopifnot(nrow(z)==10,sum(z$N)==3459,sum(z$deaths)==306,nrow(q)==1,max(abs(z$observed_rate-z$deaths/z$N))<1e-12,all(z$ci_high<=.4),all(z$ci_low>=0))
 col <- c(blue,coral)[k]
 g <- ggplot(z,aes(mean_pred,observed_rate))+geom_abline(aes(linetype='Ideal',intercept=0,slope=1),linewidth=.4,colour='#505050')+
  geom_line(colour=col,linewidth=.38)+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.006,colour=col,linewidth=.4)+
  geom_point(shape=21,fill='white',colour=col,stroke=.65,size=2.1)+
  scale_linetype_manual(name=NULL,values=c(Ideal='dashed'))+
  scale_x_continuous(limits=c(0,.4),breaks=seq(0,.4,.1),labels=function(x)paste0(round(100*x),'%'),expand=c(0,0))+
  scale_y_continuous(limits=c(0,.4),breaks=seq(0,.4,.1),labels=function(x)paste0(round(100*x),'%'),expand=c(0,0))+
  coord_fixed(1)+labs(title=paste0(c('A','B')[k],'  Clinical + day ',c(1,3)[k]),subtitle=sprintf('Intercept %.3f; slope %.3f',q$calibration_intercept,q$calibration_slope),x='Mean predicted mortality',y='Observed mortality')+
  base+theme(axis.text.x=element_text(margin=margin(t=5)),axis.text.y=element_text(margin=margin(r=5)),legend.text=element_text(size=pt),legend.position='inside',legend.position.inside=c(.77,.18),legend.background=element_rect(fill='white',colour=NA),legend.key.width=grid::unit(9,'mm'),plot.margin=margin(6,8,6,6))
 calplots[[k]] <- g
}
g3 <- wrap_plots(calplots,nrow=1)
writefig(g3,3,183/25.4,93/25.4,ids=c('A_calibration','B_calibration'),rg=list(c('A_calibration','B_calibration')),cg=NULL)
write.csv(shown,file.path(qa,'S3_displayed_calibration.csv'),row.names=FALSE)
stopifnot(identical(before,vapply(paths,function(p)digest(file=p,algo='sha256'),'')))
write_json(list(source_aggregates_unchanged=TRUE,models_refitted=FALSE,main_figures_modified=FALSE,S1=list(rows=16,N_each=3459,events_each=306),S2=list(estimates=12,intervals=12,age_p_Holm=tests$p_holm[1:2]),S3=list(points=20,intervals=20,axis=c(0,.4),interval_method='Existing Wilson'),figure_width_mm=183),file.path(qa,'figure_numeric_receipt.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(qa,'R_session.txt'))
for(i in 1:3){
 status <- system2('pwsh',c('-NoProfile','-File',shQuote(Sys.getenv("PIC_PYTHON_HELPER")),shQuote(Sys.getenv("PIC_PANEL_ALIGNMENT_AUDIT")),shQuote(file.path(qa,paste0('S',i,'_panel_layout.json'))),'--json-out',shQuote(file.path(qa,paste0('S',i,'_alignment_audit.json'))),'--strict'))
 if(status!=0)stop('Rendered panel-alignment gate failed: figure ',i)
}
cat('Three supplementary figures rendered; original aggregate values and models unchanged.\n')
