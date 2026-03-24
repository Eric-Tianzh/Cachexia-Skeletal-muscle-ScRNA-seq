####################################################
####------------------figure 1------------------####
####################################################


###--Proportions of subclusters in LLC and PF muscle groups--###
###code in R###
library(pcutils)
library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)
library(colorspace)
color_dict <- c(
  'B' = "#1f77b4",                     
  'EC-1' = "#38917E",        
  'EC-2' = "#66c2a4",                 
  'EC-3' = "#56b4b0",                 
  'EC-4' = "#238b45",                                 
  'FAPs-1' = "#dd7c4f",
  'FAPs-2' = "#CC5200",                
  'FAPs-3' = "#FF9933",                
  'FAPs-4' = "#FFD699",                
  'Glial_cell' = "#af8dc3",            
  'Mono/Macro' = "#c51b7d",            
  'MuSCs' = "#d0dd97",
  'Myonuclei' = "#4D261A",             
  'Neutrophil' = "#17becf",            
  'Schwann cell' ="#963B79",  
  'Smooth M-1' = "#f1a0c4",            
  'Smooth M-2' = "#fbb4b9",            
  'Smooth M-3' = "#F08080",            
  'T' = "#FFD700",                     
  'Tenocytes-1' = "#7895C1",
  'Tenocytes-2' = "#2171b5",           
  'Tenocytes-3' = "#4292c6"               
)

qmobs<-read.csv('qM_labels0.5_obs.csv', header = TRUE,row.names=NULL)

cluster_proportions <- qmobs %>%
  group_by(labels, paired_group_name) %>%
  summarise(count = n()) %>%
  group_by(labels) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
llc_proportions <- qmobs %>%
  filter(paired_group_name == "LLC") %>%
  group_by(labels) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

pf_proportions <- qmobs %>%
  filter(paired_group_name == "PF") %>%
  group_by(labels) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

llc_proportions$group<-"LLC"
pf_proportions$group<-"PF"
qmprop<-rbind(llc_proportions,pf_proportions)

qmprop$tissue<-"quadriceps"
prop<-qmprop
prop$labels<-factor(prop$labels,levels=c("T","B","Mono/Macro","Neutrophil","Schwann cell","Glial_cell","EC-1","EC-2","EC-3","EC-4","MuSCs","Myonuclei",
                                         "FAPs-1","FAPs-2","FAPs-3","FAPs-4","Smooth M-1" ,"Smooth M-2","Smooth M-3","Tenocytes-1","Tenocytes-2","Tenocytes-3"))

pdf("figure1g.pdf",width = 6,height=9)
ggplot(prop,aes(group,proportion,fill=as.factor(labels)))+ 
  geom_bar(stat="identity",
           position="fill")+
  scale_fill_manual(values = color_dict)+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes = list(shape = 21,size=3)))+
  ggtitle("")+ xlab("")+ guides(fill=guide_legend(title=NULL))+ 
  theme(axis.text.x = element_text(angle=0,vjust = 1),
        panel.grid.major=element_line(colour=NA),
        legend.key = element_rect(fill = NA),
        panel.grid.minor = element_blank())
dev.off()