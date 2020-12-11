library(tidyverse)
library(readxl)
library(writexl)
library(plyr)
library(zoo)
library(lubridate)
library(ggsci)
library(RColorBrewer)

############ 1 DATA ADDRESS###################
############1-1 Merge data in 2013-2018 to be a dataset called "NoV1318"################### 

setwd("$Datadir")

a = list.files("$Datadir/final.seq1/")   #####Set 98.9% for OTU clustering

dir=paste("$Datadir/final.seq1/",a,sep="")

NoV1318=read.csv(file=dir[1],header=T,sep=",")

n=length(dir)

for(i in 2:n){
  
  new.data=read.csv(file=dir[i],header=T,sep=",")
  
  NoV1318=rbind(NoV1318,new.data)
  
}

str(NoV1318)

colnames(NoV1318)[1] <- "date"

str(NoV1318)

GI1617<-NoV1318[grep("-GI",NoV1318$date),]

GI1617<-GI1617[grep("GI.[1-9]",GI1617$type),]

GI1617$date <- gsub("-GI", "", GI1617$date, fixed = TRUE)

NoV1318<-rbind(NoV1318,GI1617)

NoV1318$date<-as.yearmon(NoV1318$date,"%Y-%m")%>% as.Date()

class(NoV1318$date)

write.csv(NoV1318,file="$Datadir/NoV1318.csv",row.names=F)

NoV1318<-read.csv("$Datadir/NoV1318.csv",header=T)



######### 1-2 Export 2013-2018 GI and GII data respectively ########## 

GI1318<-NoV1318[grep("GI.[1-9]",NoV1318$type),]

GII1318<-NoV1318[grep("GII",NoV1318$type),]

write.csv(GII1318,file="$Datadir/GII1318.csv",row.names=F)

write.csv(GI1318,file="$Datadir/GI1318.csv",row.names=F)



######################### 2 Analyse and Plot###################
####################2-1 GI genetic distribution #########################
#########################2-1-1 GI genotype distribution in monthly sample from 2013-2018 ####################

GI1318 <-GI1318 %>%select(-Seq) %>% filter(reads>1000) ######Based on plasmid assay, the genotype with reads < 1000 is at a low ratio and may have PCR bias when counting final no. of reads;

str(GI1318)

GI1318_1<-aggregate(reads~date+type,data= GI1318,FUN="sum", na.rm=T)

GI1318_1<-ddply(GI1318_1,"date",transform,per=reads/sum(reads) * 100,options(digits = 2))%>%select(-digits)%>%arrange(date) 

GI1318_1$year <- year(GI1318_1$date)

GI1318_1$month<- month(GI1318_1$date)

GI1318_1$month<-as.numeric(GI1318_1$month)

str(GI1318_1)

p1<-ggplot(GI1318_1,aes(x=month,y=per,fill=type))+
  
  scale_x_continuous(breaks = GI1318_1$month)+ facet_wrap(~year)+
  
  geom_bar(stat="identity",width=0.9,alpha=0.8)+
  
  scale_fill_npg()+labs(x="",y="Percentage(%)")+
  
  guides(fill = guide_legend(title = "Genotype"))+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.title.y=element_text(family="Times",size=8),legend.title=element_text(family="Times",size=6),
        
        legend.text=element_text(family="Times",size=6),
        
        legend.key= element_rect(fill = "transparent",color="transparent"),
        
        legend.key.size=unit(10,"pt"))


p1 #####Figure 3B

ggsave("$Datadir/PDF/GI1318.pdf",p1,width=5,height=2.5)


#########################2-1-2 GI genotype distributed annually from 2013-2018 ####################

GI1318_1$year<-factor(GI1318_1$year,levels=c("2013","2014","2015","2016","2017","2018",ordered=TRUE))

GI1318_2<-aggregate(reads~year+type,data= GI1318_1,FUN="sum", na.rm=T)

GI1318_2 <- ddply(GI1318_2,"year",transform,per=reads/sum(reads) * 100,options(digits = 2))%>%select(-digits)%>%arrange(year)

p2 <- ggplot(GI1318_2,aes(x=type,y=per,color=year))+geom_jitter(width = 0.1,height=0.1,alpha = 0.8)+coord_flip()+labs(x="Genotype",y="Percentage(%)")+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.title=element_text(family="Times",size=8),legend.title=element_text(family="Times",size=8),
        
        legend.text=element_text(family="Times",size=8),
        
        legend.key= element_rect(fill = "transparent",color="transparent"),
        
        legend.key.size=unit(14,"pt"))

p2

####################2-2 GII genetic distribution #########################

#########################2-2-1 GII genotype distribution in monthly sample from 2013-2018 ####################

GII1318<-GII1318%>% filter(reads>1000) 

GII1318_1<-ddply(GII1318,"date",transform,per=reads/sum(reads) * 100,options(digits = 2)) %>% select(-digits)%>%arrange(date)

GII.4_1318<-GII1318_1[grep("GII.4",GII1318_1$type),] ###### Select GII.4 sequences to construct phylogenetic tree

GII.17_1318<-GII1318_1[grep("GII.17",GII1318_1$type),]  ###### Select GII.17 sequences to construct phylogenetic tree

GII.2_1318<-GII1318_1%>% filter(type=="GII.2")  ###### Select GII.2 sequences to construct phylogenetic tree

write.csv(GII.4_1318,file="$Datadir/GII.4_1318.csv",row.names=F)

write.csv(GII.17_1318,file="$Datadir/GII.17_1318.csv",row.names=F)

write.csv(GII.2_1318,file="$Datadir/GII.2_1318.csv",row.names=F)

str(GII1318)

GII1318_2<-aggregate(reads~date+type,data= GII1318,FUN="sum", na.rm=T)

GII1318_2<-ddply(GII1318_2,"date",transform,per=reads/sum(reads) * 100,options(digits = 2)) %>%select(-digits)%>%arrange(date) ####add variable:percentage of each genotype

GII1318_2$year <- year(GII1318_2$date)

GII1318_2$month<- month(GII1318_2$date)

GII1318_2$month<-as.numeric(GII1318_2$month)

unique(GII1318_2$type)

GII1318_2$type<-factor(GII1318_2$type,levels=c("GII.1","GII.2","GII.3","GII.4","GII.5","GII.6","GII.7","GII.9","GII.10","GII.12","GII.13","GII.14","GII.17","GII.21"),ordered = TRUE)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

GIIcolourCount = length(unique(GII1318_2$type))
p3<-ggplot(GII1318_2,aes(x=month,y=per,fill=type))+
  
  scale_x_continuous(breaks = GII1318_2$month)+ facet_wrap(~year)+
  
  geom_bar(stat="identity",width=0.9,alpha=0.8)+
  
  scale_fill_manual(values = getPalette(GIIcolourCount))+labs(x="",y="Percentage(%)")+
  
  guides(fill = guide_legend(title = "Genotype"))+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.title.y=element_text(family="Times",size=8),legend.title=element_text(family="Times",size=6),
        
        legend.text=element_text(family="Times",size=6),
        
        legend.key= element_rect(fill = "transparent",color="transparent"),
        
        legend.key.size=unit(10,"pt"))


p3 ###Figure 3C

ggsave("$Datadir/PDF/GII1318.pdf",p2,width=6,height=3)

#########################2-2-2 GII genotype distributed annually from 2013-2018 ####################

GII1318_3<-aggregate(reads~year+type,data= GII1318_2,FUN="sum", na.rm=T)

GII1318_3<-ddply(GII1318_3,"year",transform,per=reads/sum(reads) * 100,options(digits = 2)) %>%select(-digits)

GII1318_3$year<-factor(GII1318_3$year,levels=c("2013","2014","2015","2016","2017","2018",ordered=TRUE))

p4 <- ggplot(GII1318_3,aes(x=type,y=per,color=year))+geom_jitter(width = 0.1,height=0.1,alpha = 0.8)+coord_flip()+
  
  labs(x="Genotype",y="Percentage(%)")+guides(color="none")+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.title=element_text(family="Times",size=8),legend.title=element_text(family="Times",size=8),
        
        legend.text=element_text(family="Times",size=8),
        
        legend.key= element_rect(fill = "transparent",color="transparent"))


ggsave("$Datadir/PDF/NoV type-year1.pdf",grid.arrange(p2,p4,ncol = 2,widths=c(5/9,4/9)),width=6,height=3)######Figure 3A


###############3 Correlation between viral load and the number of outbreak ##################
library(ggplot2)
library(gtable)
library(grid) 
library(ggsci)
library(reshape2)
library(zoo)

setwd("$Datadir\outbreakdata")

copout<-read.csv("copies-outbreak.csv",header=T)

copout$totcop<-ifelse(copout$totcop==0,1,copout$totcop)

str(copout)

copout$date<-as.Date(copout$date,origin="%y-%m-%d")

datebreak<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

copout$genotype<-factor(copout$genotype,levels=c("GI","GII"))

copout<-subset(copout,totcop!="1")

copout$totcop<-log10(copout$totcop)

GII<-copout[grep("GII",copout$genotype),]

out<- read.csv("GD_NoV_Outbreaks_13-18.csv",header = T)

names(out)

str(out)

colnames(out)[1] <- "date"

names(out)<-c("date","GII.4","GII.17","GII.2","GII.3","Other genotype")

out<-melt(out,id="date")

names(out)<-c("date","Genotype","outbreaks")

out$date<-as.Date(out$date,origin="%y-%m-%d")

datebreak1<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

out$Genotype<- factor(out$Genotype,levels=c("GII.2","GII.3","GII.4" ,"GII.17","Other genotype"),ordered=TRUE)



g1 <- ggplot(copout, aes(date, totcop,color=genotype)) +geom_point(shape=1,position=position_jitter(width = 0.1,height=0.1),alpha=0.8)+
  
  stat_smooth(method ="loess",span=0.2,size=0.3,se=FALSE,alpha=0.8)+scale_color_manual(values=c("#0B3A42","#FF404A")) +
  
  scale_x_date(breaks = datebreak)+labs(x="",y = 'Log 10(copies/L)') +  theme(panel.grid.major = element_line(color="ghostwhite"), 
                                                                              
                                                                              panel.grid.minor = element_line(linetype="blank"), 
                                                                              
                                                                              axis.text = element_text(family="Times",size = 8), 
                                                                              
                                                                              panel.background = element_rect(fill = "transparent",color="gray"),
                                                                              
                                                                              panel.border = element_rect(fill = "transparent",color="gray"),
                                                                              
                                                                              axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(size=8))
g1

g2<- ggplot(out, aes(date, outbreaks,fill=Genotype)) + geom_bar(stat="identity",position="stack",alpha=0.8)+
  
  scale_x_date(breaks = datebreak1)+scale_fill_npg()+
  
  labs(x="",y = 'No. of Outbreak') + theme(panel.grid.major = element_line(linetype="blank"), 
                                           
                                           panel.grid.minor = element_line(linetype="blank"), 
                                           
                                           axis.text = element_text(family="Times",size = 8), 
                                           
                                           panel.background = element_rect(fill = "transparent",color="gray"),
                                           
                                           panel.border = element_rect(fill = "transparent",color="gray"),
                                           
                                           axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(size=8))
g2


ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

ggplot2.two_y_axis(g1, g2)########Figure 1A




###################4 Phylogenetic tree visualization#################

library(ggtree)
library(treeio)
library(ape)
library(phangorn)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(tidytree)
library(viridis)
library(lubridate)
library(zoo)
library(gridExtra)
library(gtable)
library(grid) 

##################4-1-1 GII.4 tree ################

setwd("$Datadir/GII.4")

G4ML<- read.newick("$Datadir/GII.4/GII.4-rooted.newick")

ggtree(G4ML)+geom_text2(aes(subset=!isTip,label=node),hjust=-.3,color="red",size=2)

G41318<-read.csv("GII.4_1318.csv",header=T,stringsAsFactors = F,sep = ",")%>%select(-Seq,-OTU,-type)

str(G41318)

G41318$date<-as.Date(G41318$date)

G4data <- as_tibble(G4ML)

str(G4data)

G4data <- left_join(G4data,G41318,by=c("label"="tip.lable"))

G4data$year <- year(G4data$date)

head(G4data,10)

G4data$per[is.na(G4data$per)]<- 0

G4data$year[is.na(G4data$year)]<- "2013"

str(G4data)

G4tree <- as.treedata(G4data) ####switch to treedata with traits

str(G4tree)

str(G4tree@data$year)

p3 <- ggtree(G4tree,color="black",size=0.5,linetype=1,ladderize=T)+
  
  geom_tippoint(aes(color=year,size=per))+geom_tiplab(size=3,hjust=0)+
  
  theme(legend.position =c(0.1,0.75),legend.box="horizontal",
        
        legend.box.background =element_rect(colour = "#3C4155", size =0.1 ),
        
        legend.box.margin = margin(6, 6, 6, 6))+xlim(NA,0.15)+
  
  scale_size(breaks=c(10,30,50,70,90),limits=c(0.12,100))+
  
  geom_treescale(width = 0.01,fontsize=3, linesize=0.1)+
  
  geom_cladelabel(120,label="GII.4_Syd12_Var1",barsize = 1.5,
                  
                  color="#00BEC3",fontsize = 6,offset = 0.005,align=TRUE,offset.text = 0.001)+
  
  geom_cladelabel(148,label="GII.4_Syd12_Var2",barsize = 1.5,
                  
                  color="#CA83FF",fontsize =6,offset = 0.005,align=TRUE,offset.text = 0.001)+
  
  geom_strip(147,184,label="GII.4_Syd12",barsize = 1.5,
             
             color="#7CAE00",fontsize = 6,offset = 0.005,align=TRUE,offset.text = 0.001)

p3

ggsave("$Datadir/GII.4/GII.4 tree.pdf",p3,width=10,height=5)



###############################4-1-2 GII.4 outbreak+GII.4 variant distribution################

setwd("C:/Lu Boss/Sewage NoV/data")

ob <- read.csv("GD_NoV_Outbreaks_13-18.csv",header = T)

str(ob)

colnames(ob)[1] <- "date"

ob$date<-as.Date(ob$date,origin="%y-%m-%d")

datebreak1<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

ob <- ob%>% gather("GII.4","GII.17","GII.2" ,"GII.3","other", key="type",value="outbreaks")

g1 <- ggplot(filter(ob,type=="GII.4"),aes(x=date,y=outbreaks))+geom_bar(stat = "identity",fill="gray",alpha=0.7)+
  
  scale_x_date(breaks = datebreak1,date_labels = "%m")+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1),
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8),
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))+
  
  labs(x="",y="No. of GII.4-outbreak")

g1


GII.41318<-read.csv("$Datadir/GII.4/GII.4_1318.csv",header=T) %>% select(-Seq,-OTU,-type,-tip.label)

GII.41318$date<-as.Date(GII.41318$date)

datebreak<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

GII.41318_1<-aggregate(per~date+clade,data= GII.41318,FUN="sum", na.rm=T)

g2<-ggplot(GII.41318_1,aes(x=date,y=per,color=clade))+
  
  geom_point(shape=1,size=2,position=position_jitter(width = 0.1,height=0.1))+
  
  scale_color_manual(values=c("#00BEC3","#CA83FF","#7CAE00"),guide=FALSE)+
  
  scale_x_date(breaks = datebreak, date_labels = "%m")+
  
  labs(x='',y = 'Percentage in Sewage') +
  
  theme(panel.grid.major = element_line(linetype="blank"), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray"),
        
        panel.border = element_rect(fill = "transparent",color="gray"),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))

g2

ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

ggsave("$Datadir/GII.4/GII.4 outbreak-clade.pdf",ggplot2.two_y_axis(g1, g2),width=6,height=3)


###################4-2-1 GII.17 tree#############

setwd("$Datadir/GII.17")

G17ML<- read.newick("$Datadir/GII.17/GII17-rooted.newick")

ggtree(G17ML)+geom_text2(aes(subset=!isTip,label=node),hjust=-.3,color="red",size=2)

G171318<-read.csv("GII.17_1318.csv",header=T,stringsAsFactors = F,sep = ",")%>%select(-Seq,-OTU,-type)

str(G171318)

G171318$date<-as.Date(G171318$date)

G17data <- as_tibble(G17ML)

str(G17data)

G17data <- left_join(G17data,G171318,by=c("label"="tip.label"))

G17data$year <- year(G17data$date)

head(G17data,10)

G17data$per[is.na(G17data$per)]<- 0

G17data$year[is.na(G17data$year)]<- "2013"

str(G17data)

G17tree <- as.treedata(G17data) ####switch to treedata with traits

str(G17tree)

str(G17tree@data$year)

p4 <- ggtree(G17tree,color="black",size=0.1,linetype=1,ladderize=T)+
  
  geom_tippoint(aes(color=year,size=per))+###geom_tiplab(size=3,hjust=0)+
  
  theme(legend.position =c(0.15,0.75),legend.box="horizontal",
        
        legend.box.background =element_rect(colour = "#3C4155", size =0.1 ),
        
        legend.box.margin = margin(4, 4, 4, 4))+xlim(NA,0.3)+
  
  scale_size(breaks=c(10,30,50,70,90),limits=c(0.57,100))+
  
  geom_treescale(width = 0.02,fontsize=3, linesize=0.1)+
  
  geom_cladelabel(108,label="GII.17_preEpi(Kawaski323)",barsize = 1.5,
                  
                  color="#F8756C",fontsize = 4,offset = 0.008,align=TRUE,offset.text = 0.001)+
  
  geom_strip(80,117,label="GII.17_Epi(Kawaski308)",barsize = 1.5,
             
             color="#00BA38",fontsize = 4,offset = 0.008,align=TRUE,offset.text = 0.001)+
  
  geom_strip(83,107,label="GII.17_Epi(Kawaski308)",barsize = 1.5,
             
             color="#00BA38",fontsize = 4,offset = 0.008,align=TRUE,offset.text = 0.001)+
  
  geom_cladelabel(133,label="GII.17_variant 1",barsize = 1.5,
                  
                  color="#619CFF",fontsize = 4,offset = 0.008,align=TRUE,offset.text = 0.001)

p4

ggsave("$Datadir/GII.17/GII.17 tree.pdf",p4,width=12,height=10)


###############################4-2-2 GII.17 outbreak+GII.17 variant distribution################

g1 <- ggplot(filter(ob,type=="GII.17"),aes(x=date,y=outbreaks))+geom_bar(stat = "identity",fill="gray",alpha=0.7)+
  
  scale_x_date(breaks = datebreak1,date_labels = "%m")+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1),
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8),
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))+
  
  labs(x="",y="No. of GII.17-outbreak")

g1


GII.171318<-read.csv("$Datadir/GII.17/GII.17_1318.csv",header=T) %>% select(-Seq,-OTU,-type,-tip.label)

unique(GII.171318$clade)

GII.171318$clade<-factor(GII.171318$clade,levels=c("GII.17_Epi(Kawaski308)","GII.17_preEpi(Kawaski323)","GII.17_variant 1"),order=T)

GII.171318$date<-as.Date(GII.171318$date)

datebreak<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

GII.171318_1<-aggregate(per~date+clade,data= GII.171318,FUN="sum", na.rm=T)

g2<-ggplot(GII.171318_1,aes(x=date,y=per,color=clade))+
  
  geom_point(shape=1,size=2,position=position_jitter(width = 0.1,height=0.1))+
  
  scale_color_manual(values=c("#00BA38","#F8756C","#619CFF"),guide=FALSE)+
  
  scale_x_date(breaks = datebreak, date_labels = "%m")+
  
  labs(x='',y = 'Percentage in Sewage') +
  
  theme(panel.grid.major = element_line(linetype="blank"), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray"),
        
        panel.border = element_rect(fill = "transparent",color="gray"),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))

g2

ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

ggsave("$Datadir/GII.17/GII.17 outbreak-clade.pdf",ggplot2.two_y_axis(g1, g2),width=6,height=3)



###############4-3-1 GII.2 tree ##############

setwd("$Datadir/GII.2/")

G2ML<- read.tree("$Datadir/GII.2/GII.2-rooted.newick")

ggtree(G2ML)+geom_tiplab(size=2)+geom_text2(aes(subset=!isTip,label=node, hjust=-0.5),size=2)+geom_rootedge()

G21318<-read.csv("GII.2_1318.csv",header=T,stringsAsFactors = F,sep = ",")%>%select(-Seq,-OTU,-type)

str(G21318)

G21318$date<-as.Date(G21318$date)


G2data <- as_tibble(G2ML)

str(G2data)

G2data <- left_join(G2data,G21318,by=c("label"="tip.label"))

G2data$year <- year(G2data$date)

head(G2data,10)

G2data$per[is.na(G2data$per)]<- 0

G2data$year[is.na(G2data$year)]<- "2013"

str(G2data)

G2tree <- as.treedata(G2data) ####switch to treedata with traits

str(G2tree)

str(G2tree@data$year)

p1 <- ggtree(G2tree,color="black",size=0.1,linetype=1,ladderize=T)+
  
  geom_tippoint(aes(color=year,size=per))+####geom_tiplab(size=1,hjust=0)+
  
  theme(legend.position =c(0.15,0.75),legend.box="horizontal",
        
        legend.box.background =element_rect(colour = "#3C4155", size =0.1 ),
        
        legend.box.margin = margin(4, 4, 4, 4))+xlim(NA,0.06)+
  
  scale_size(breaks=c(10,30,50,70,90),limits=c(0.41,100))+
  
  geom_treescale(width = 0.01,fontsize=3, linesize=0.1)+
  
  geom_cladelabel(62,label="Var1(GII.P16-GII.2_Epi)",family="Times",barsize = 1.5,
                  
                  color="#7CAE00",fontsize = 3,align=TRUE,offset.text = 0.0002)+
  
  geom_cladelabel(93,label="Var2(GII.P16-GII.2_12-15)",family="Times",barsize = 1.5,
                  
                  color="#CA83FF",fontsize = 3,align=TRUE,offset.text = 0.0002)+
  
  geom_strip(115,116,label="Var3(GII.P2-GII.2)",family="Times",barsize = 1.5,
             
             color="#00BEC3",fontsize = 8,align=TRUE,offset.text = 0.0002)

p1


ggsave("$Datadir/GII.2/GII.2 tree-new.pdf",p1,width=12,height=6)


#####################4-3-2 GII.2 outbreak+GII.2 variant distribution#########

setwd("C:/Lu Boss/Sewage NoV/data")

ob <- read.csv("GD_NoV_Outbreaks_13-18.csv",header = T)

str(ob)

colnames(ob)[1] <- "date"

ob$date<-as.Date(ob$date,origin="%y-%m-%d")

datebreak1<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

ob <- ob%>% gather("GII.4","GII.17","GII.2" ,"GII.3","other", key="type",value="outbreaks")

g1 <- ggplot(filter(ob,type=="GII.2"),aes(x=date,y=outbreaks))+geom_bar(stat = "identity",fill="gray",alpha=0.7)+
  
  scale_x_date(breaks = datebreak1,date_labels = "%m")+
  
  theme(panel.grid.major = element_line(color="ghostwhite",size=0.1),
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8),
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray",size=0.1),
        
        panel.border = element_rect(fill = "transparent",color="gray",size=0.1),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))+
  
  labs(x="",y="No. of Novs outbreak")

g1


GII.21318<-read.csv("$Datadir/GII.2/GII.2_1318.csv",header=T) %>% select(-Seq,-OTU,-type,-tip.label)

GII.21318$date<-as.Date(GII.21318$date)

datebreak<-seq(as.Date("2013-01-01"),as.Date("2018-12-01"),by="1 month")

GII.21318_1<-aggregate(per~date+clade,data= GII.21318,FUN="sum", na.rm=T)

g2<-ggplot(GII.21318_1,aes(x=date,y=per,color=clade))+
  
  geom_point(shape=1,size=2,position=position_jitter(width = 0.1,height=0.1))+
  
  scale_color_manual(values=c("#7CAE00","#CA83FF","#00BEC3"),guide=FALSE)+###scale_x_date(breaks = datebreak)+
  
  scale_x_date(breaks = datebreak, date_labels = "%m")+
  
  labs(x='',y = 'Percentage in Sewage') +
  
  theme(panel.grid.major = element_line(linetype="blank"), 
        
        panel.grid.minor = element_line(linetype="blank"), 
        
        axis.text = element_text(family="Times",size = 8), 
        
        axis.ticks = element_line(size=0.1),
        
        panel.background = element_rect(fill = "transparent",color="gray"),
        
        panel.border = element_rect(fill = "transparent",color="gray"),
        
        axis.text.x = element_text(angle=45, hjust=1),axis.title.y=element_text(family="Times",size=8))

g2

ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

ggsave("$Datadir/GII.2/GII.2 outbreak-clade.pdf",p2,width=6,height=3)






