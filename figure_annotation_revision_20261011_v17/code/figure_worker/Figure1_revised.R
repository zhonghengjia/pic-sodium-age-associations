flow_plot <- function(flow) {
  # Same accepted nodes/geometry/style; only the internal approval phrase is replaced.
  p <- ggplot(flow[seq_len(9),])+geom_rect(aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill='white',colour='black',linewidth=.4)+
    geom_text(aes(x,y,label=label),family='Arial',size=3.1,lineheight=.95)+theme_void(base_family='Arial')+
    coord_cartesian(xlim=c(.7,11.3),ylim=c(2.05,8.85),clip='off')+
    labs(title=NULL,caption=NULL)+
    theme(plot.title=element_text(face='bold',size=13),plot.caption=element_text(size=9,hjust=0),plot.margin=margin(12,12,12,12))
  edge <- data.frame(x=c(rep(4,4),rep(4,4)),xend=c(rep(4,4),rep(7.85,4)),y=c(7.9,6.5,5.1,3.7,7.6,6.2,4.8,3.4),yend=c(7.3,5.9,4.5,3.1,7.6,6.2,4.8,3.4))
  p+geom_segment(data=edge,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,linewidth=.4,arrow=grid::arrow(length=grid::unit(1.5,'mm')))
}
