library(ggplot2)


# Polygone coordinates
las_poly <- data.frame(x=c(0.75,1,1),y=c(0.75,0.5,1))
ada_poly <- data.frame(x=c(0.875,1,1),y=c(0.875,0.75,1))


ppoly <- ggplot() +
	theme_minimal() +
	theme(panel.grid.major = element_line(colour = "grey50"), legend.box = "horizontal",
		  legend.position="bottom", plot.title = element_text(vjust=2, size=12)) +
	geom_abline(intercept=3/2 ,slope=-1,  linetype=2) +
	geom_abline(intercept=0,   slope=1, linetype=2) +
	#geom_abline(intercept=1 ,slope=-1/2,  linetype=3)+
	geom_abline(intercept=7/4,   slope=-1, linetype=3)+
	xlim(0,1)+ylim(0,1)  + xlab(label = 'a') + ylab(label = 'd')+
	geom_polygon(data=las_poly,aes(x=x,y=y),fill='gray25',alpha=0.5) + 
	geom_polygon(data=ada_poly,aes(x=x,y=y),fill='gray25',alpha=0.5)
print(ppoly)

ggsave(ppoly,filename = 'admissible.pdf',width=16,height=10,units='cm')

