#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


library(ggplot2)

print(paste("Working in : ", getwd()))

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

data <- read.csv(file=args[1], sep='\t', header=T, comment.char = "#")
data["y1"]<-1
data["y2"]<-2


##########################################################################
# all queries in 1 pdf
#p <- ggplot() + 
#  scale_x_continuous(name="x", minor_breaks=seq(0,30000,by=1000), guide = guide_axis(minor.ticks = TRUE)) + 
#  scale_y_continuous(name=data$query) +
#  #uncomment and adapt to select types colours
#  #scale_fill_manual(values = c(
#  #    "typeA" = "#00ff00",
#  #    "typeB" = "#ff0000",
#  #    "typeC" = "#00ffff",
#  #    "typeD" = "#0000ff",
#  #    "typeF" = "#000000"
#  #  )) +         
#  geom_rect(data=data, mapping=aes(xmin=position_start, xmax=position_end, ymin=y1, ymax=y2, fill=predicted_type)) +
#  facet_grid(query ~ . , switch="y") +
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())
#h=2+length(unique(data$query)) 
#ggsave("recomb_pattern.pdf", units="in", dpi=300, width = 10, height=h)

##########################################################################
# or 1 image per query
for (i in unique(data$query)) {
  p <- ggplot() + 
    scale_x_continuous(name="x", limits = c(0,max(data$position_end)), minor_breaks=seq(0,max(data$position_end),by=1000), guide = guide_axis(minor.ticks = TRUE)) + 
    scale_y_continuous(name=data[data$query==i,]$query) +
    #uncomment and adapt to select types colours
    #scale_fill_manual(values = c(
    #    "typeA" = "#00ff00",
    #    "typeB" = "#ff0000",
    #    "typeC" = "#00ffff",
    #    "typeD" = "#0000ff",
    #    "typeF" = "#000000"
    #  )) +         
    geom_rect(data=data[data$query==i,], mapping=aes(xmin=position_start, xmax=position_end, ymin=y1, ymax=y2, fill=predicted_type)) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  ggsave(paste("recomb_pattern__",i,".pdf",sep=""), units="in", dpi=300, width = 10, height=2)
}

