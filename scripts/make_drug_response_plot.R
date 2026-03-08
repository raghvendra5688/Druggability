library(data.table)
library(ggplot2)
library(ggrepel)
library(tibble)

setwd("/home/rmall/TII/Projects/Raghav/Druggability/")

ZBP1_df <- fread("Results/Raj_Analysis/ZBP1_Drug_Predictions_Subset.csv",header=TRUE)
ZBP1_df <- as.data.frame(ZBP1_df)
ZBP1_df$Size <- 1
ZBP1_df$Color <- "blue"
ZBP1_df[ZBP1_df$Drug_Name%in%c("ZAPi","LOC14"),]$Color <- "blue"
ZBP1_df[ZBP1_df$Drug_Name%in%c("ZAPi","LOC14"),]$Size <- 3
subset_ZBP1_df <- ZBP1_df[ZBP1_df$Drug_Name%in%c("ZAPi","LOC14"),]


g_ZBP1 <- ggplot(ZBP1_df, aes(x=Morgan_AAC_BindingDB_IC50,y=Morgan_CNN_BindingDB_IC50))+
  geom_point(aes(color=Color,size=Size)) + xlab(expr("DeepPurpose Method 1 \n Binding Affinity Score"))+ ylab(expr("DeepPurpose Method 2 \n Binding Affinity Score")) + 
  ggtitle("Drugs targeting ZBP1") + scale_color_manual(values=c("red"="red","blue"="blue"))+
  #xlim(c(10,13))+ylim(c(0,20))+
  # geom_label_repel(data = subset_ZBP1_df, aes(label = Drug_Name, x = Morgan_AAC_KIBA, y = MPNN_CNN_BindingDB),
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

ggsave(filename="Results/Raj_Analysis/ZBP1_ZAPi_Hit_without_annotation.pdf",plot = g_ZBP1, device = pdf(), height = 4, width = 6, units="in")
dev.off()