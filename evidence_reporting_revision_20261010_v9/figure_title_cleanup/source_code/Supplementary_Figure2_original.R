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
font <- 'Arial';pt <- 10;txt <- pt/ggplot2::.pt
base <- theme_classic(base_size=pt,base_family=font)+theme(axis.text=element_text(colour=ink,size=pt),axis.title=element_text(colour=ink,size=pt),plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=9.5),plot.margin=margin(4,4,4,4))
tabtheme <- theme_void(base_size=pt,base_family=font)+theme(plot.margin=margin(4,2,4,2))
writefig <- function(g,i,w,h,ids=NULL,rg=NULL,cg=NULL){
 stem <- file.path(out,paste0('Supplementary_Figure',i))
 ggsave(paste0(stem,'.png'),g,width=w,height=h,units='in',dpi=600,device=ragg::agg_png,bg='white')
 ggsave(paste0(stem,'.pdf'),g,width=w,height=h,units='in',device=grDevices::cairo_pdf,bg='white')
 svglite::svglite(paste0(stem,'.svg'),width=w,height=h,bg='white');print(g);grDevices::dev.off()
 ggsave(paste0(stem,'.tiff'),g,width=w,height=h,units='in',dpi=600,device=ragg::agg_tiff,compression='lzw',bg='white')
 write_patchwork_panel_layout(g,file.path(qa,paste0('S',i,'_panel_layout.json')),w,h,panel_ids=ids,row_groups=rg,column_groups=cg)
}

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
 lab <- column(z$label,'Population')+annotate('text',x=.02,y=6.1,label='Representative age',hjust=0,fontface='italic',size=9.5/ggplot2::.pt,family=font)
 targets <- column(sprintf('%d/%s',z$target_deaths,format(z$target_N,big.mark=',',trim=TRUE)),paste0(c('Low','High')[k],'\nDeaths/N'))
 refs <- column(sprintf('%d/%s',z$reference_deaths,format(z$reference_N,big.mark=',',trim=TRUE)),'Normal\nDeaths/N')
 nums <- column(sprintf('%.2f (%.2f-%.2f)',z$OR,z$ci_low,z$ci_high),'Adjusted OR\n(95% CI)')
 forest <- ggplot(z,aes(OR,y))+annotate('segment',x=1,xend=1,y=.5,yend=7.45,linetype='dashed',linewidth=.35,colour='#7A7A7A')+
  geom_segment(aes(x=ci_low,xend=ci_high,yend=y),linewidth=.55,colour=col)+geom_point(aes(shape=shape),size=2.8,colour=col)+
  scale_shape_manual(values=c(common=15,age=16),guide='none')+scale_x_log10(limits=c(.2,20),breaks=c(.2,.5,1,2,5,10,20),labels=c('0.2','0.5','1','2','5','10','20'),expand=expansion(mult=.035))+
  scale_y_continuous(limits=c(.4,8.8),expand=c(0,0))+annotate('text',x=2,y=8.2,label='Adjusted OR',fontface='bold',size=txt,family=font)+geom_hline(yintercept=7.5,linewidth=.25,colour='#787878')+
  labs(x=NULL,y=NULL)+base+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.text.x=element_text(size=9.5),plot.margin=margin(4,2,4,2))
 lab <- lab+labs(title=paste0(c('A  Low','B  High')[k],' vs normal sodium'))+theme(plot.title=element_text(size=12,face='bold',family=font),plot.subtitle=element_text(size=9.5,family=font))
 forests <- c(forests,list(lab,targets,refs,nums,forest))
}
g2 <- wrap_plots(forests,ncol=5,widths=c(2.75,1.25,1.25,2.65,5.6))
writefig(g2,2,183/25.4,170/25.4,ids=c('A_labels','A_target','A_normal','A_numbers','A_forest','B_labels','B_target','B_normal','B_numbers','B_forest'),rg=list(c('A_labels','A_target','A_normal','A_numbers','A_forest'),c('B_labels','B_target','B_normal','B_numbers','B_forest')),cg=list(c('A_labels','B_labels'),c('A_target','B_target'),c('A_normal','B_normal'),c('A_numbers','B_numbers'),c('A_forest','B_forest')))
write.csv(do.call(rbind,displayed),file.path(qa,'S2_displayed_estimates.csv'),row.names=FALSE)


stopifnot(identical(before,vapply(paths,function(p)digest(file=p,algo='sha256'),'')))
write_json(list(source_aggregates_unchanged=TRUE,models_refitted=FALSE,other_figures_modified=FALSE,estimates=12,intervals=12,column_widths=c(2.75,1.25,1.25,2.65,5.6),body_font_pt=10,smallest_requested_font_pt=9.5,figure_mm=c(183,170)),file.path(qa,'figure_numeric_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(qa,'R_session.txt'))
status <- system2('pwsh',c('-NoProfile','-File',shQuote(Sys.getenv("PIC_PYTHON_HELPER")),shQuote(Sys.getenv("PIC_PANEL_ALIGNMENT_AUDIT")),shQuote(file.path(qa,'S2_panel_layout.json')),'--json-out',shQuote(file.path(qa,'S2_alignment_audit.json')),'--strict'))
if(status!=0)stop('Rendered S2 alignment failed')
cat('Supplementary Figure 2 only: wider forest and larger font, existing estimates unchanged.\n')
