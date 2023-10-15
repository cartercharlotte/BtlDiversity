#Just graphing the ANI values so they are prettier
library(ggplot2)
library(reshape2)
library(ggdendro)
library(viridis)
setwd("")

ANI_data <- read.csv("ANI_matrix.csv", header=T)
colnames(ANI_data)[1]<- "Strain"
ANI_matrix <- as.matrix(ANI_data[,-1])
rownames(ANI_matrix) <- ANI_data$Strain
ANI_data2 <- melt(ANI_data, na.rm = F)
colnames(ANI_data2) <- c("Strain1", "Strain2", "ANI")
ANI_data2$ANI<-(1-ANI_data2$ANI)*100


ANI_dendro <- as.dendrogram(hclust(d = dist(x = ANI_matrix)))
ANI_order <- order.dendrogram(ANI_dendro)
ANI_data2$Strain1 <- factor(x = ANI_data2$Strain1,
                               levels = ANI_data$Strain[ANI_order], 
                               ordered = TRUE)
ANI_data2$Strain2 <- factor(x = ANI_data2$Strain2,
                            levels = ANI_data$Strain[ANI_order], 
                            ordered = TRUE)

ANI_plot<- ggplot(ANI_data2, aes(x=Strain1, y=Strain2, fill= ANI))+
  geom_tile()+
  scale_fill_gradient(low="#ffffff", high="#5671d4", name='% ANI')+
  theme(
    panel.background=element_blank(),
    axis.text.x=element_text(angle = 45, size=12, hjust=0.95),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=12),
    axis.title.y=element_blank(),
    #legend.key=element_blank(),     #removes the border
    legend.text=element_text(size=10),
    legend.title=element_text(size=14, face="bold"),
    legend.title.align = 0)+
  annotate("rect", xmin = 16.5, xmax = 28.5, ymin =16.5, ymax = 28.5,
     fill = "transparent", color="#5a396a", lwd=1, linetype = "dotted")+
  annotate("rect", xmin = 0.5, xmax = 16.5, ymin = 0.5, ymax = 16.5,
     fill = "transparent", color="#547772", linetype = "dotted", lwd=1)+
  annotate("text", x = 20.25, y = 27.75, label = "M. endofungorum", color="white", size=7, fontface = "italic")+
  annotate("text", x = 13.5, y = 1.5, label = "M. rhizoxinica", color="white", size=7, fontface = "italic")
ANI_plot
  
