#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("ggplot2")
library("cowplot")
library("dplyr")

specimen<-args[2]
P1<-paste0(args[3],"X",args[3])
P2<-paste0(args[4],"X",args[4])
HE<-paste0(args[3],"X",args[4])

aims<-read.table(args[1], sep=",", header= T)

colour_scale<-c("darkblue","orange","darkred")
colour_scale_hist<-c("darkblue","orange","darkred")

names(colour_scale)<-c(P1,HE,P2)
names(colour_scale_hist)<-c(P1,HE,P2)

LABLES<-c(P1,HE, P2,"mix")
LABLES_hist<-c(P1,HE, P2)

png_aims<-paste0("./", specimen, "_AIMs.png")
  
aims_auto<-aims[!aims$chromo %in% c("chrZ", "chrX", "chrW", "chrY"), ]
aims_sex<-aims[aims$chromo %in%  c("chrZ", "chrX", "chrW", "chrY"), ]

t<-levels(reorder(aims_auto$chromo, aims_auto$position))
order<-c("chrZ",t)

aimsplot<-ggplot()+
  geom_point(data=aims, aes(x=position,y=chromo, color=geno), size=1, shape=18)+
  theme_bw() +
  ggtitle(specimen) +  xlab("Position") +  ylab("") +
  theme(axis.text.y = element_text(angle=15, vjust=-0.1))+
  scale_colour_manual(labels=LABLES,values=colour_scale,name="Ancestry")+
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=6)) +guides(col=guide_legend(override.aes = list(size = 5)))+ theme(legend.position =  c(0.85, 0.35)) + 
  scale_y_discrete(limits=order)
  
aims_hist<-aims[!aims$chromo %in% c("chrZ", "chrX", "chrW", "chrY"), ]

histplot<-ggplot(aims_hist, aes(geno, fill=geno))+geom_bar()+theme_bw()+theme(axis.title.x =element_blank(), axis.text.x = element_text(angle = 25, vjust = 0.5, hjust=0.7))  + geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.1)  + ggtitle("")+scale_fill_manual(labels=LABLES_hist,values=colour_scale_hist,name="Ancestry")+scale_x_discrete(labels=c("","",""))

comb<-plot_grid(aimsplot, histplot, ncol = 1,nrow = 2, rel_heights =c(1,1) )
ggsave(png_aims, plot=comb, dpi=600, device="png", width=6, height = 8, units = "in")

pdf_aims<-paste0("./", specimen, "_AIMs.pdf")
ggsave(pdf_aims, plot=comb, dpi=600, device="pdf", width=6, height = 8, units = "in")

########################################################################################################
#hybrid indices and heterozygosity

AIMs_auto<-aims[!aims$chromo %in% c("chrZ", "chrX", "chrW", "chrY"),]

P1_count<-length(AIMs_auto[AIMs_auto$geno == P1,]$geno)
P2_count<-length(AIMs_auto[AIMs_auto$geno == P2,]$geno)
HE_count<-length(AIMs_auto[AIMs_auto$geno == HE,]$geno)
tot_count<-length(AIMs_auto$geno)

#interspecific heterozygosity
inteHE<-HE_count/tot_count

#Hybrid Index parental allele frequencies
P1_freq<-(P1_count*2 + HE_count)/(tot_count*2) #Parent 1 allele count
P2_freq<-(P2_count*2 + HE_count)/(tot_count*2) #Parent 2 allele count

if(P1_freq + P2_freq == 1){
  HI<-cbind(specimen,inteHE, P1_freq)
  colnames(HI)<-c("Individual","heterozygosity","hybrid_index")
write.table(HI, file = paste0(specimen,"_HI_HE.csv"), row.names=FALSE, sep="\t", quote=F)
}

