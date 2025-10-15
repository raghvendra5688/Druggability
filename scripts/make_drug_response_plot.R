library(data.table)
library(ggplot2)
library(ggrepel)
library(tibble)

setwd("/home/rmall/TII/Projects/Raghav/Druggability/")

nlrp3_df <- fread("Results/Raj_Analysis/NLRP3_Drug_Predictions_Subset.csv",header=TRUE)
nlrp3_df <- as.data.frame(nlrp3_df)
nlrp3_df$Size <- 1
nlrp3_df$Color <- "blue"
nlrp3_df[nlrp3_df$Drug_Name%in%c("ZAPi","LOC14"),]$Color <- "red"
nlrp3_df[nlrp3_df$Drug_Name%in%c("ZAPi","LOC14"),]$Size <- 3
subset_nlrp3_df <- nlrp3_df[nlrp3_df$Drug_Name%in%c("ZAPi","LOC14"),]


g_nlrp3 <- ggplot(nlrp3_df, aes(x=Morgan_AAC_KIBA,y=MPNN_CNN_BindingDB))+
  geom_point(aes(color=Color,size=Size)) + xlab("DeepPurpose Method 1") + ylab("DeepPurpose Method 2") + 
  ggtitle("Drugs targeting NLRP3") + scale_color_manual(values=c("red"="red","blue"="blue"))+
  xlim(c(10,13))+ylim(c(0,20))+
  # geom_label_repel(data = subset_nlrp3_df, aes(label = Drug_Name, x = Morgan_AAC_KIBA, y = MPNN_CNN_BindingDB),
  #                  box.padding = 0.25,
  #                  point.padding = 0.05,
  #                  nudge_x = .25,
  #                  nudge_y = .25,
  #                  segment.linetype = 1,
  #                  max.overlaps = 20,
  #                  segment.curvature = -1e-20,
  #                  arrow = arrow(length = unit(0.01, "npc")),
  #                  segment.color = 'grey50') +
  theme_light()+
  theme(legend.position="none")+
  #theme(legend.text=element_text(size=20),
  #      legend.title=element_text(size=20))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title=element_text(color="grey20", size=20, face="plain"))

ggsave(filename="Results/Raj_Analysis/NLRP3_LOC14_Hit_without_annotation.pdf",plot = g_nlrp3, device = pdf(), height = 4, width = 6, units="in")
dev.off()