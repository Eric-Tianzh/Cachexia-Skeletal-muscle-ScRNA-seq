####################################################
####------------------figure 6------------------####
####################################################

###code in R###
library(pcutils)
library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)
library(colorspace)
library(corrplot)
library(paletteer)

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratObject)

library(ROGUE)
library(reshape2)


###--Proportions of subclusters in LLC and PF muscle groups--###
color_dict <- c(
  "Activated-1"= "#E4A031",
  "Activated-2" = "#F6E2C1",
  "Activated-3" = "#D75B4E",
  "Activated-4" = "#DBE4FB",
  "Transitional" = "#A1D4A2",
  "Myoblast-1" = "#7C4D77",
  "Myoblast-2" = "#B6B3D6",
  "Myocyte" = "#26445E",
  "Quiescent" = "#1f77b4"
)

data <- read.table('paired_group_name_perct_t.xls', header = TRUE, sep = "\t",row.names=NULL)
rownames(data)<-data$labels
data<-data[,-1]
pdf("figure6c.pdf",width = 7,height=7)
stackplot(data,flow = T,bar_params = list(width = 0.8, position = "stack"),topN = 10,number = T)+
  scale_fill_manual(values = color_dict)+ 
  theme_bw()+ 
  ggtitle("")+ xlab("")+ 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_blank())
dev.off()   


###--MuSCs trajectory analysis--###
##Monocle3
library(monocle3)
library(ggplot2)
library(viridis)
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
library(ggridges)

p1<-plot_cells(cds1, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_roots = FALSE, label_branch_points = FALSE,
           cell_size = 1.5,trajectory_graph_color = "white",trajectory_graph_segment_size = 1)+
  theme(axis.text.x = element_text(angle=0,vjust = 1,colour = "black",size = 15),
        axis.text.y = element_text(angle=0,hjust = 1,colour = "black",size = 15),
        axis.title.x = element_text(colour = "black",size = 15),
        axis.title.y = element_text(colour = "black",size = 15))+
  scale_color_viridis(discrete=FALSE, name="pseudotime")
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 14, name = "Spectral")))
ggsave(p1, file='figure6e-left.pdf', width=10, height=10)

#pseudotime ridge plot
ps_tim <- data.frame(pseudotime=pseudotime(cds1),cell=names(pseudotime(cds1)))
meta<-data.frame(pData(cds))
metacluster<-data.frame(cluster=meta$labels,cell=rownames(meta))
psdata<-merge(metacluster,ps_tim)
psdata<-psdata[!is.infinite(psdata$pseudotime),]
psdata$cluster<-factor(psdata$cluster,levels = c("Quiescent","Activated-1","Activated-2","Activated-3","Activated-4",
                                               "Transitional","Myoblast-1","Myoblast-2","Myocyte"))

p2<-ggplot(psdata, aes(x = pseudotime, y = cluster, fill = cluster)) +
  geom_density_ridges(scale = 2) +
  theme_ridges(font_size = 16, grid = TRUE,center_axis_labels=TRUE) +
  #scale_fill_brewer(palette = 4)
  scale_fill_viridis(discrete=TRUE,name="cluster")+
  theme(legend.position = "None")
ggsave(p2,file='figure6e-right.pdf', width=10, height=10)


##The activity of 12 signaling pathways in LLC and PF groups, predicted by PROGENy
result<-read.csv("summarized_progeny_scores_matrix_llcpf.csv")
result<-melt(result)
colnames(result)<-c("celltype","pathway","Activity")
split_cols <- do.call(rbind, strsplit(as.character(result$celltype), "_"))
colnames(split_cols) <- c("subtype", "group")
result <- cbind(result, split_cols)
resultplot <- result[-which(result$pathway %in% c("Androgen","Estrogen")),]
resultplot$pathway<-gsub("JAK.STAT","JAK-STAT",resultplot$pathway)
resultplot$celltype<-factor(resultplot$celltype,levels = rev(c("Quiescent_LLC","Quiescent_PF",
                                                           "Activated-1_LLC","Activated-1_PF","Activated-2_LLC","Activated-2_PF",
                                                           "Activated-3_LLC","Activated-3_PF","Activated-4_LLC","Activated-4_PF","Intermediate-state_LLC","Intermediate-state_PF",
                                                           "Myoblast-1_LLC","Myoblast-1_PF","Myoblast-2_LLC","Myoblast-2_PF",
                                                           "Myocyte_LLC","Myocyte_PF")))
resultplot$subtype<-factor(resultplot$subtype,levels=c("Quiescent","Activated-1","Activated-2","Activated-3","Activated-4",
                                                       "Intermediate-state","Myoblast-1","Myoblast-2","Myocyte"))

groupMean <- resultplot %>%
  dplyr::group_by(pathway) %>%
  dplyr::mutate(MeanV = mean(Activity)) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::select(pathway, MeanV)

pdf("figure6f.pdf",width = 13,height = 6)
ggplot(resultplot, aes(x = celltype, y = Activity,color = subtype,fill=subtype)) +
  #geom_bar(aes(fill=celltype),stat = "identity",width=0.1)+
  geom_vline(aes(xintercept = celltype,color = celltype),linewidth=0.6)+
  geom_hline(data = groupMean, aes(yintercept = MeanV), linetype = "dashed",linewidth=0.3) +
  geom_point(aes(shape=group),size = 3) +#,color="black"
  coord_flip() +
  theme_bw()+
  scale_color_manual(values = color_dict)+
  scale_fill_manual(values = color_dict)+
  scale_shape_manual(values = c(15,19))+
  #theme_gray(base_size = 12) +
  theme(axis.text.x = element_text(color = "black",size = 7),
        axis.text.y = element_text(color = "black"),
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #panel.grid.major=element_line(colour=NA),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        strip.background = element_rect(color = "grey") 
  ) +
  facet_grid(~ pathway,scales = "free")#scales = "free"
dev.off()


group="Hypoxia"
result<-read.csv("progeny_scores_df.csv")
result$cellcluster<-paste(result$CellType,result$CellGroup,sep="_")
result$CellType<-factor(result$CellType,levels = rev(c("Quiescent","Activated-1","Activated-2","Activated-3","Activated-4",
                                                       "Intermediate-state","Myoblast-1","Myoblast-2","Myocyte")))

resultsub<-result[grep(group,result$Pathway),]

p<-ggplot(resultsub, aes(
  x = CellType,           
  y = Activity,
  #group=CellGroup,
  # fill = Activity
  fill = CellGroup)) +
  geom_boxplot(aes(fill=CellGroup),staplewidth=0.5,outliers = TRUE,outlier.colour = "gray",outlier.alpha = 0.6)+
  geom_hline(yintercept = 0,color="#333333",linetype="dashed")+
  #scale_fill_manual(values=c("#D75B4E","#4292c6"))+,
  scale_fill_manual(values=c("#8D2F25","#4A5F7E"))+
  stat_compare_means(aes(group = CellGroup),method = "t.test",label = "p.signif")+
  #scale_fill_gradient2(low="#4292c6",mid="white", high="#D75B4E")+
  #scale_fill_gradient2(low="#238b45", high="#E4A031",midpoint = 0)+
  #scale_color_gradient2(low="#238b45", high="#E4A031",midpoint = 0)+
  #coord_flip()+# x轴与y轴互换位置
  labs(x="",y=group,title = '')+
  #coord_flip()+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(color="black",angle=45,vjust =0.6,hjust=0.5),
        axis.text.y = element_text(color='black'),
        panel.grid.major=element_line(colour=NA),
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_rect(colour="#4292c6", fill="#4292c6"),
        strip.text.y=element_text(colour="white",face = "bold"))+                                          # 标签大小
  # geom_text(aes(label=dup,hjust=ifelse(depth=="2x",-2.5,3)),size=2)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) 
pdf("figure6g.pdf",width = 6,height = 8)
print(p)
dev.off()