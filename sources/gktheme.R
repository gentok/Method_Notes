## Default Theme for ggplot2 04/25/2017##

gktheme <-
  theme(axis.text=element_text(size=10, colour="black"),
        axis.title.x=element_text(size=12,face="bold", vjust=-1.5),
        axis.title.y=element_text(size=12,face="bold", vjust=1.5),
        plot.title=element_text(size=12,face="bold", vjust=2,hjust=0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.background = element_rect(fill=NA, colour="black", size=0.5, linetype=1),
        legend.background = element_rect(fill=NA,colour=NA),
        legend.position = c(0.75,0.68),
        legend.key.width = unit(1.5, "cm"))
