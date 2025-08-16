# bioinformatics
Bioinformatics analysis code
# 01bar_gokegg_term
#!/usr/bin/env Rscript
####################################################################################
#####################################################################################
#####参数获取
#####################################################################################
#suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
    make_option(c("-f", "--infile"),action = "store",type = "character", help = "The Input file"),
    make_option(c("-p", "--group"),action = "store",type = "character",default = "Unbound", help = "The group name"),
    make_option(c("-t", "--outtitle"),action = "store",type = "character",default = "Terms", help = "The title of outimage;default = Term Plot"),
    make_option(c("-n", "--outname"),action = "store",type = "character",default = "term", help = "The name of outimage;default = term.pdf")
)
opt <- parse_args(OptionParser(option_list = option_list))

#opt$infile <- "test.txt"
#opt$infile <- "/users/chengc/dev2016/graphicwork/Macaque/3_pca/All_sample_LncRNA_exp_RPKM_annot_gencode.xls.lncRNA"
#opt$group <- "/users/chengc/dev2016/graphicwork/Macaque/3_pca/Sample_group.txt"
library("ggplot2")
#library("factoextra")
#library("FactoMineR")
#library("stats")
library("ggthemes")

theme_paper <- theme(
    # panel.border = element_rect(fill = NA,colour = "black"),
    # panel.grid.major = element_line(colour = "grey88",
    #                                 size = 0.2),
    # panel.grid.minor = element_line(colour = "grey95",
    #                                 size = 0.5),
    # axis.text.x= element_text(vjust = 1,hjust = 1, angle = 45),
    # legend.position = "top",
    # legend.direction = "horizontal",
    panel.grid.major.x = element_line(colour = "grey70",
                                  size = 0.2),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14,face = "bold"),
    axis.text = element_text(size = 14, colour = "black"))

data = read.table(opt$infile, header=T, com='', quote='',sep="\t", check.names=F)

data <- data[order(data[,6]),]

rlen <- length(rownames(data))

if (rlen == 0){
    print("null data,exit")
    q()
}


#data <- subset(primary_data,Module==opt$group)
data$lp <-  -log(data[,6]+0.00000000001,10)
if (rlen>10){
    data <- data[1:10,]
}

#print(data)

rlen <- length(rownames(data))
clen <- length(colnames(data))

maxlp <-  max(data$lp)
#M1
#data$PT = factor(data$Term,levels = c("Terpenoid backbone biosynthesis","Sulfur relay system","Nicotinate and nicotinamide metabolism","Purine metabolism","DNA replication","Protein export","Glycerophospholipid metabolism","Bacterial secretion system","Homologous recombination","Pyrimidine metabolism","Mismatch repair","Peptidoglycan biosynthesis","Aminoacyl-tRNA biosynthesis","Ribosome"))
##M2
#data$PT = factor(data$Term,levels = c("Pentose phosphate pathway","Microbial metabolism in diverse environments","Alanine, aspartate and glutamate metabolism","Glycerolipid metabolism","Lysine biosynthesis","2-Oxocarboxylic acid metabolism","Nitrogen metabolism","Histidine metabolism","ABC transporters"))
##M3
#data$PT = factor(data$Term,levels = c("Microbial metabolism in diverse environments","Fatty acid biosynthesis","Citrate cycle (TCA cycle)","Biosynthesis of secondary metabolites","Pyruvate metabolism","Metabolic pathways","Amino sugar and nucleotide sugar metabolism","Starch and sucrose metabolism","Propanoate metabolism","Pentose and glucuronate interconversions","Ascorbate and aldarate metabolism","Phosphotransferase system (PTS)","Galactose metabolism"))
##M4
#data$PT = factor(data$Term,levels = c("Alanine, aspartate and glutamate metabolism","Glycine, serine and threonine metabolism","Microbial metabolism in diverse environments","RNA degradation","Purine metabolism","Glycolysis / Gluconeogenesis","Butanoate metabolism"))

##M1_HMP
#data$PT = factor(data$Term,levels = c("Arginine and proline metabolism","Degradation of aromatic compounds","Purine metabolism","Riboflavin metabolism","ABC transporters","Bacterial chemotaxis","Phosphonate and phosphinate metabolism","Bacterial secretion system","Pentose and glucuronate interconversions","Flagellar assembly"))

#Unbound
#data$PT = factor(data$Term,levels = c("RNA degradation","Nicotinate and nicotinamide metabolism","Biosynthesis of amino acids","Homologous recombination","Aminoacyl-tRNA biosynthesis","Pyrimidine metabolism","Amino sugar and nucleotide sugar metabolism","Oxidative phosphorylation","Biosynthesis of secondary metabolites","Metabolic pathways"))


#unUnbound



for(i in 1:rlen){
  data[i,clen+1] <- paste(strwrap(data[i,1],width = 50), collapse = "\n")
  maxlen <- 90


  if (nchar(data[i,clen+1]) > maxlen ){
    while (substr(data[i,clen+1], maxlen, maxlen)!=" " && substr(data[i,clen+1], maxlen, maxlen)!="\n" && nchar(data[i,clen+1]) > maxlen){
        maxlen<- maxlen+1
    }
    data[i,clen+1] <- strtrim(data[i,clen+1],maxlen-1)
    data[i,clen+1] <- paste(data[i,clen+1],"..")
  }

}


data[,clen+1] = factor(data[,clen+1],levels = rev(data[,clen+1]))

ggplot(data, aes(x=data[,clen+1], y=lp)) +
    #
    # geom_bar(stat="identity",fill=rgb(212,93,24, maxColorValue = 255),width=0.5) + theme_minimal() +
    # geom_bar(aes(fill=data$PValue),stat="identity",width=0.5) + theme_minimal() + scale_fill_gradient(low="#198752",high = "#c0f3da")+
    #Unbound
    geom_bar(aes(fill=data[,7]),stat="identity",width=0.5) + scale_fill_gradient(low="#dd4a41",high = "#e98882")+
    #unUnbound
    # geom_bar(aes(fill=data$CP),stat="identity",width=0.5) + scale_fill_gradient(low="#2f8f90",high = "#5bc8ca")+
    #M2
    # geom_bar(aes(fill=data$CP),stat="identity",width=0.5) + theme_minimal() + scale_fill_gradient(low="#00578b",high = "#5994cd")+
    # #M3
    # geom_bar(aes(fill=data$CP),stat="identity",width=0.5) + theme_minimal() + scale_fill_gradient(low="#a53f2a",high = "#f0ccc4")+
    # #M4
    # geom_bar(aes(fill=data$CP),stat="identity",width=0.5) + theme_minimal() + scale_fill_gradient(low="#d8bd00",high = "#ffee73")+
    # #M1 HMP
    # geom_bar(aes(fill=data$CP),stat="identity",width=0.5) + theme_minimal() + scale_fill_gradient(low="#00e0e1",high = "#a3f4f4")+

    coord_flip() + theme_minimal()+theme_paper+
    geom_text(aes(x=data[,clen+1],y=lp + 1.6 * sign(lp),label=paste(data[,4],"/", data[,5],sep = "")),hjust=0.5, size=6.5,
              color="black")+
    ylim(0,maxlp+3)+
    labs(x="Pathway/Term",
         y="-log10 Pvalue",
         title=opt$outname,
         fill="Corrected Pvalue")

ggsave(paste(opt$outname,".pdf",sep = ""), width = 270, height = 130, units = "mm")
ggsave(paste(opt$outname,".png",sep = ""), width = 270, height = 130, units = "mm", dpi=600)
#ggsave(paste(opt$outname,".80.png",sep = ""), width = 270, height = 130, units = "mm", dpi=80)
#ggsave(paste(opt$outname,".300.png",sep = ""), width = 270, height = 130, units = "mm", dpi=300)
#ggsave(paste(opt$outname,".600.png",sep = ""), width = 270, height = 130, units = "mm", dpi=600)


# 02bar_multi_per
#!/usr/bin/env Rscript
#####################################################################################
#####参数获取
#####################################################################################
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))


usage = "Rscript %prog -f -t -n -o

example: Rscript %prog
    -f /users/ablife/ablife-R/Bar/Bar_Rpkm/latest/Mapping_distribution.txt
    -t Mapping_distribution
    -n Mapping_distribution
    -o ./"

option_list <- list(
  make_option(c("-f", "--file"),action = "store",type = "character",
    help = "The Input file"),
  make_option("--log",action = "store_true",default = FALSE,
    help = "handel the data using log"),
  make_option(c("-t", "--title"),action = "store",type = "character",
    help = "The title of outimage"),
  make_option(c("-n", "--filename"),action = "store",type = "character",
    help = "The name of outimage"),
  make_option(c("-X","--Xaxis"),action = "store",type = "character",default = "Type",
    help = "The label of X axis"),
  make_option(c("-Y","--Yaxis"),action = 'store',type = "character",default = "Percentage(%)",
    help = "The label of Y axis"),
  make_option("--ymax",action = "store",type = "integer",default = "100",
    help = "The max of y axis"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
    help = "The outdir")
  )
opt <- parse_args(OptionParser(option_list = option_list,usage=usage)) #####解释器

setwd(opt$outdir)           ##Set the out path

###################################################################################
####加载包 （Load Package）
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(RColorBrewer)
library("ggsci")
library("ggthemes")
###################################################################################
####加载颜色（Load Color）
###################################################################################
colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')


###################################################################################
##数据处理（Deal with Data）
title <- gsub('_',' ',opt$title)
Xaxis <- gsub('_',' ',opt$Xaxis)
Yaxis <- gsub('_',' ',opt$Yaxis)

data <- read.table(opt$file,header=T) #read.table  从外部文本文件读取数据  must with header
colname <- colnames(data) #获取表头名字
dim_data <- dim(data) #获取数据的维度，行列数
ymax <- as.numeric(opt$ymax)
sample <-colname[2:dim_data[2]] #获取部分表头的名字
newdata1 <- data
for (i in 2:dim_data[2])
  {sum <- sum(data[,i])
  for(j in 1:dim_data[1])
    {mean <- data[j,i]/sum
      mean <- round(mean*100,2)
      percentage <- paste(mean,'%',sep='')
      newdata1[j,i]<-percentage
      # mean <- log10(mean)+4  #通过log10处理使得图像好看一些
      data[j,i]<-mean}
  }
newdata <- melt(data,id.vars=colname[1],measure.vars = c(colname[2:dim_data[2]]))
ymax = max(newdata[,3]) + 10
if(ymax > 100){
  ymax =100
}
Percent <- paste(newdata[,3],'%',sep='')

head(newdata)

##################################################################################
###Plot Theme for ABLife 
###theme(),Tha last term without comma
##################################################################################
ablife_theme_bar <- function(base_size = 12){
  library(grid)     ####for using unit function
  theme(
    plot.title = element_text(size=12,lineheight = 10,colour="#000000"),

    axis.title.x = element_text(size=12,colour = "#000000"),
    axis.title.y = element_text(size=12,colour = "#000000"),
    axis.text.x = element_text(size = 12,colour = "#000000"),
    axis.text.y = element_text(size = 12 ,colour = "#000000"),

    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5,"cm"),

    panel.background = element_rect(fill = "white",colour = NA),
    panel.border = element_rect(size = 1,colour = "#000000",fill =NA)
    # panel.border = element_rect(size = 1,colour = "#8B8B8B",fill =NA)
    # panel.grid.major = element_line(size=0.5,colour = "#BFBFBF"),
    # panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
    )
}

theme_paper <- theme(
    # panel.border = element_rect(fill = NA,colour = "black"),
    panel.grid.major = element_line(size = 0.2,linetype = "dashed"),
    # panel.grid.major = element_blank(),
    
    panel.grid.minor = element_line(colour = "grey90",size = 0.2,linetype = "dashed"),
    # axis.text.x= element_text(vjust = 1,hjust = 1, angle = 45),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key = element_rect(fill = 'white', color = 'white'),
    legend.key.size = unit(0.5, "cm"),
    # legend.position = c(0.82, 0.68),
    # legend.position = "top",
    # legend.direction = "vertical",
    legend.title = element_blank(),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(1,1,1,1,unit="cm"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, colour = "black")
)

##################################################################################
###Plot by ggplot2
ggplot(newdata,aes(x = newdata[,1],y=newdata[,3],fill = newdata[,2],stat = "identity"))+
    geom_bar(width = 0.9,stat = "identity",position = "dodge")+labs(x=Xaxis,y= Yaxis)+
    # geom_text(aes(label=Percent),hjust=0,angle = -3,position = position_dodge(.9),size=5)+
    scale_y_continuous(limits=c(0,ymax))+
    coord_flip()+
    # scale_fill_manual(name="Name",values = colour[1:dim_data[2]-1])+
    scale_fill_nejm()+
    # geom_text(aes(label=newdata[,3]),vjust =1.5,position = position_dodge(.9),size=2)+
    theme_bw() + theme_paper
    
#scale_fill_discrete(name="Sample_name")  ##可以独立调整图例的名字，label，breaks

###################################################################################
###Save Plot File 
ggsave(file = paste(opt$filename,"_Bar.pdf",sep=''), width = 180,height = 120,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,"_Bar.png",sep=''), width = 180,height = 120,dpi = 450,units = "mm")

# 03Box_plot
#!/usr/bin/env Rscript

#####################################################################################
#####参数获取
#####################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

usage = "The prog is used to display Box of data,the box display from 1/4 to 3/4 of data,and a line to show the median data.
	if there is whiskers,then it start from the edge of the box and extend to the furthest data point that is within 1.5 times the IQR
	if there are any data points that are past the ends of the whiskers,they are considered outliers and displayed with dots.
Rscript %prog -f -t -n 
example: Rscript %prog
	-f CLIP_Peak_len 
	-t CLIP_Peak_len
	-n CLIP_Peak_len
"

option_list <- list(
	make_option(c("-f", "--file"),action = "store",type = "character",
		help = "The Input file"),
	make_option(c("-t","--title"),action = "store",type = "character",
		help = "The title of the plot"),
	make_option(c("-x","--xaxis"),action = "store",type = "character",default = "Location",
		help = "The name of x-axis "),
	make_option(c("-y","--yaxis"),action = "store",type = "character",default = "Frequency of each base(%)",
		help = "The name of y-axis"),
	make_option(c("-n", "--filename"),action = "store",type = "character",
		help = "The name of outimage"),
	make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
		help = "The outdir")
	)
opt <- parse_args(OptionParser(option_list = option_list)) #####解释器

setwd(opt$outdir)		# Set the Outpath

###################################################################################
####加载包 （Load Package）
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(gtable)
library(grid)
library(RColorBrewer)
###################################################################################
####加载颜色（Load Color）
###################################################################################
colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')
#colour1 <- brewer.pal(8,"Dark2")

#######data format sample###########################################################
#Type    Peaks
#CLIP_overlap    351
#CLIP_overlap    111
#CLIP_overlap    226
#CLIP_overlap    316
#CLIP_overlap    316
#CLIP_overlap    486
#CLIP_sp 346
#CLIP_sp 321
#CLIP_sp 91
#CLIP_sp 291
#CLIP_sp 236
#CLIP_sp 216
#CLIP_sp 191
#####################################################################################

######Instruction for data#######
#first column is the X axis---discrete data
#second column is the Y axis to plot -- continuous data
#the data is to display the summary about X axis data
#################################


###################################################################################
##数据处理（Deal with Data）
data <- read.table(file= opt$file,header = T)#####	must with header
title <- gsub('_',' ',opt$title)
colname <- colnames(data)
dim_data <- dim(data)

Item <- factor(data[,1],levels = unique(data[,1]))
x_axis = gsub('_',' ',opt$xaxis)
y_axis = gsub('_',' ',opt$yaxis)
#######变形数据后直接给variablfe命名即可

##################################################################################
###Plot Theme for ABLife 
###theme(),Tha last term without comma
##################################################################################
ablife_theme_line <- function(base_size = 12){
	library(grid)		####for using unit function
	theme(
			plot.title = element_text(size = 12,lineheight = 100,colour = "black",hjust = 0.5),
			# axis.text.x=element_text(angle=75,hjust=1,size = 5,colour = "black"),
			axis.title.x = element_text(size = 12,colour = "black"),
			axis.text.y = element_text(size = 12 ,colour = "black"),
			axis.title.y = element_text(size = 12,colour = "black"),
			panel.background = element_rect(colour = "black"),
			legend.title = element_text(size = 12),
			legend.text = element_text(size = 12),
			strip.text.x = element_text(colour = "black",size = "8"),
			strip.background = element_rect(colour = "black")
			)
}

##################################################################################
###Plot by ggplot2
ggplot(data)+ ####加载数据
		geom_boxplot(aes(x = Item,y=data[,2],stat = "identity",fill = Type))+									#####基本作图函数
		labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
		# facet_grid(. ~ Sex,scales = "free_x",space="free") +
		theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
			axis.ticks.x = element_blank(),
			legend.position="none"
			# legend.position = c(0.3, 0.92),
			# legend.direction = "horizontal"
			) + 
		ablife_theme_line()+
##scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
		scale_colour_hue(name=opt$legend_name)
#scale_colour_manual("Sample Name",values = colour[1:length(sample)])
#scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks
###################################################################################
###Save Plot File 
ggsave(file = paste(opt$filename,"box.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,"box.png",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")

#ggplot(data)+ ####加载数据
#geom_jitter(alpha=I(1/4),aes(x = Item,y=data[,2],stat = "identity",colour = Type))+									#####基本作图函数
#labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
 		# facet_grid(. ~ Sex,scales = "free_x",space="free") +
 		theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
 			axis.ticks.x = element_blank(),
 			legend.position="none"
 			# legend.position = c(0.3, 0.92),
 			# legend.direction = "horizontal"
 			) + 
 		ablife_theme_line()+
    scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
 		scale_colour_hue(name=opt$legend_name)
    scale_colour_manual("Sample Name",values = colour[1:length(sample)])
    scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks
 ggsave(file = paste(opt$filename,"jitter.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
 #ggsave(file = paste(opt$filename,".pdf",sep=''), width = 12,height = 8,dpi = 450)

#ggplot(data)+ ####加载数据
 		geom_density(alpha = 0.8,aes(x = Item,stat = "identity",fill=Type))+									#####基本作图函数
 		labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
 		theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
 			axis.ticks.x = element_blank(),
 			legend.position="none"
 			# legend.position = c(0.3, 0.92),
 			# legend.direction = "horizontal"
 			) + 
 		ablife_theme_line()+
    #scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
 		scale_colour_hue(name=opt$legend_name)
		scale_colour_manual("Sample Name",values = colour[1:length(sample)])
    scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks
 ggsave(file = paste(opt$filename,"density.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
 ggsave(file = paste(opt$filename,".pdf",sep=''), width = 12,height = 8,dpi = 450)

#ggplot(data)+ ####加载数据
#geom_violin(aes(x = Item,y=data[,2],stat = "identity",fill=Type))+									#####基本作图函数
#labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
#theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
#axis.ticks.x = element_blank(),
#legend.position="none"
##legend.position = c(0.3, 0.92),
##legend.direction = "horizontal"
#) + 
#ablife_theme_line()+
###scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
#scale_colour_hue(name=opt$legend_name)
##scale_colour_manual("Sample Name",values = colour[1:length(sample)])
##scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks
#ggsave(file = paste(opt$filename,"violin.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
##ggsave(file = paste(opt$filename,".pdf",sep=''), width = 12,height = 8,dpi = 450)

#ggplot(data,aes(x=Item,y=data[,2]))+ ####加载数据
#geom_violin(aes(fill=Type))+									#####基本作图函数
#geom_boxplot(width=.1,fill = "black",outlier.colour=NA)+
#stat_summary(fun.y=median,geom="point",fill="white",shape=21,size=2.5)+
#labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
#theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
#axis.ticks.x = element_blank(),
#legend.position="none"
#legend.position = c(0.3, 0.92),
#legend.direction = "horizontal"
#) + 
#ablife_theme_line()+
#scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
#scale_colour_hue(name=opt$legend_name)
#scale_colour_manual("Sample Name",values = colour[1:length(sample)])
##scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks
#ggsave(file = paste(opt$filename,"violin_median.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
#ggsave(file = paste(opt$filename,".pdf",sep=''), width = 12,height = 8,dpi = 450)

# 04heatmapmaker_anno
#!/usr/bin/env Rscript
####################################################################################
#####################################################################################
#####参数获取
#####################################################################################
#suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
  make_option(c("-f", "--file"),action = "store",type = "character",
              help = "The Input file"),
  make_option(c("-t","--type"),action = "store",type = "character",
    help = "heatmap type:
            exp
            exp_norowname
            exp_nox
            exp_nox_log2
            exp_noxy
            exp_noscale
            exp_noscale_nox
            exp_noscale_noy
            exp_noscale_noxy
            exp_noscale_noxy_norowname
            exp_noscale_noxy_log2
            exp_colscale
            exp_colscale_nox
            exp_colscale_noy
            exp_colscale_noxy
            exp_noscale_norowname
            exp_noscale_nox_norowname
            exp_nox_norowname
            exp_nox_norowname_log2
            exp_noy_norowname
            exp_noy_norowname_log2
            exp_noxy_norowname
            exp_noxy_norowname_log2
            samplecor
            samplecor_half_upper
            samplecor_half_lower
            kmeans
            kmeans_nox
    "),
  make_option(c("-c","--colorstyle"),action = "store",type = "character",default="A",
    help = "colorstyle:
            A   :default
            C   :blackandred
            B   :big red
    "),
  make_option(c("-b","--start"),action = "store",type = "integer",default = 1,
              help = "appoint the column of start ; default = 1"),
  make_option(c("-e","--end"),action = "store",type = "integer",default = -1,
              help = "appoint the column of end; default = -1"),
  make_option(c("-x", "--xanno"),action = "store",type = "character",default = "",
              help = "The Anno file for X"),
  make_option(c("-y", "--yanno"),action = "store",type = "character",default = "",
              help = "The Anno file for Y"),
  make_option(c("-w","--width"),action = "store",type = "integer",default = 5,
              help = "width; default = 5"),
  make_option(c("-g","--height"),action = "store",type = "integer",default = 8,
              help = "heigth; default = 8"),
  make_option(c("-n","--filename"),action = "store",type = "character",default="DEG_heatmap",
              help = "The filename of the picture ; default = DEG_heatmap"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir;default = ./")
)
opt <- parse_args(OptionParser(option_list = option_list))
start <- as.numeric(opt$start)
end <- as.numeric(opt$end)
w <- as.numeric(opt$width)
h <- as.numeric(opt$height)

library(pheatmap)
library("RColorBrewer")
library("corrplot")
NO_REUSE = F

 #get the filename to use later
filename <- strsplit(opt$file,"/")[[1]]
filename <- filename[length(filename)]
filename <- sub('.txt','',filename)

 #try to reuse earlier-loaded data if possible

#print('Reading matrix file.')
primary_data = read.table(opt$file, header=T, com='',quote = "", sep="\t", row.names=1, check.names=F)
#primary_data = read.table("Sample_correlation.dat", header=T, com='', sep="\t", row.names=1, check.names=F)


annotation_col = data.frame(
  s = factor(rep(c("s"), length(colnames(primary_data))))
)
rownames(annotation_col) = colnames(primary_data)

if(opt$xanno!=""){
  annotation_col = read.table(opt$xanno, header=T, com='', sep="\t", row.names=1, check.names=F)
}


annotation_row = data.frame(
  s = factor(rep(c("s"), length(rownames(primary_data))))
)
rownames(annotation_row) = rownames(primary_data)

if(opt$yanno!=""){
  annotation_row = read.table(opt$yanno, header=T, com='', sep="\t", row.names=1, check.names=F)
}

if(end> 0){
  primary_data <- primary_data[,start:end]
}
primary_data = as.matrix(primary_data)
data = primary_data

#原花青素
scaleyellowred <- colorRampPalette(c("#08519c","#3182bd","#ffffff","#e6550d","#a63603"),space = "rgb")(500)
if (opt$colorstyle == "B"){
  scaleyellowred <- colorRampPalette(c("#08519c","#08519c","#3182bd","#3182bd","#ffffff","#e6550d","#e6550d","#a63603","#a63603"),space = "rgb")(500)
}
if (opt$colorstyle == "C"){
  scaleyellowred <- colorRampPalette(c("green","green","black","red","red"),space = "rgb")(500)
}
if (opt$colorstyle == "D"){
  scaleyellowred <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
}

annotation_colors <- c("#0072B2","#eb5d46","#F0E442","#009E73","#CC79A7","#5d8ac4","#f39659","#936355","#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#000000","#E69F00","#56B4E9","#D55E00","#CC79A7","#FCFBFD","#EFEDF5","#DADAEB","#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#54278F","#3F007D")
#scaleyellowred <- colorRampPalette(c("green","black","red"),bias=1)(256)

#scaleyellowred <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

##row anno and col anno


if (opt$type == "kmeans_nox"){
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, kmeans_k=4, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,fontsize=15, cellwidth = 15,cellheight = 15,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, kmeans_k=4, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,fontsize=15, cellwidth = 15,cellheight = 15,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_colscale"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp_colscale_nox"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_colscale_noy"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_colscale_noxy"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="column",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_norowname"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_nox_norowname"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_nox"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_noy"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_noxy"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_noxy_log2"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_noxy_norowname"){
  data <- data[rowSums(data) != 0,]
  data <- data[,colSums(data) != 0]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", cluster_cols=F,cluster_rows=F,scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noscale_norowname"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="none",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_log2"){
  data <- data[rowSums(data) != 0,]
  data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_nox"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_nox_log2"){
  data <- data[rowSums(data) != 0,]
  data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noxy"){
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp_noy_norowname"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_noy_norowname_log2"){
  data <- data[rowSums(data) != 0,]
    data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_cols="euclidean",cluster_rows=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp_nox_norowname"){
  data <- data[rowSums(data) != 0,]
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_nox_norowname_log2"){
  data <- data[rowSums(data) != 0,]
    data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}

if (opt$type == "exp_nox_log2"){
  data <- data[rowSums(data) != 0,]
    data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",show_rownames = T,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp_noxy_norowname"){
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "exp_noxy_norowname_log2"){
  data = log2(data+1)
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep='')) #设定格子的尺寸
  pheatmap(data,annotation_col = annotation_col, annotation_row = annotation_row,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_rows=F,cluster_cols=F,scale="row",show_rownames = F,border_color=NA, angle_col = "45",width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep='')) #设定格子的尺寸
}


if (opt$type == "samplecor"){
  sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
  cat(c('Gene\t'),file=paste(opt$filename,".xls",sep=''))
  write.table(sample_cor, file=paste(opt$filename,".xls",sep=''), quote=F, append=T, sep='\t')
  data <- sample_cor
  o = rownames(data)
  sample_dist = dist(t(primary_data), method='euclidean')
  hc = hclust(sample_dist, method='complete')
  data = data[hc$order, hc$order]
  data = data[o, o]
  data

  rowcount <- nrow(data)

  if (rowcount>7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,height=opt$height,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE)

    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,height=opt$height,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE)
  }

  if (rowcount<=7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE,width = 8, height = 7)

    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE,width = 8, height = 7)
  }

}


if (opt$type == "samplecor_half_upper"){
  sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
  cat(c('Gene\t'),file=paste(opt$filename,".xls",sep=''))
  write.table(sample_cor, file=paste(opt$filename,".xls",sep=''), quote=F, append=T, sep='\t')
  data <- sample_cor
  o = rownames(data)
  sample_dist = dist(t(data), method='euclidean')
  hc = hclust(sample_dist, method='complete')
  #hc = hclust(as.dist(1 - data))
  data = data[hc$order, hc$order]
  data[upper.tri(data)] = NA
  data = data[o, o]
  data

  rowcount <- nrow(data)

  if (rowcount>7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE)

    pheatmap(data, cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE)
  }

  if (rowcount<=7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE,width = 8, height = 7)

    pheatmap(data, cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE,width = 8, height = 7)
  }

}


if (opt$type == "samplecor_half_lower"){
  sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
  cat(c('Gene\t'),file=paste(opt$filename,".xls",sep=''))
  write.table(sample_cor, file=paste(opt$filename,".xls",sep=''), quote=F, append=T, sep='\t')
  data <- sample_cor
  o = rownames(data)
  sample_dist = dist(t(data), method='euclidean')
  hc = hclust(sample_dist, method='complete')
  # hc = hclust(as.dist(1 - data))
  data = data[hc$order, hc$order]
  data[lower.tri(data)] = NA
  data = data[o, o]
  data

  rowcount <- nrow(data)

  if (rowcount>7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE)

    pheatmap(data, cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,cellwidth = 15, cellheight = 15, fontsize=15,width=opt$width,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE)
  }

  if (rowcount<=7){
    pheatmap(data, annotation_col = annotation_col, annotation_row = annotation_row,cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".pdf",sep=''),display_numbers = FALSE,width = 8, height = 7)

    pheatmap(data, cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=T,fontsize=20,filename = paste(opt$filename,".png",sep=''),display_numbers = FALSE,width = 8, height = 7)
  }

}

# 04Line_multiSample_multiText

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

usage = "The prog is used to display the distribution of TSS,TTS,startcondon,stopcodon,and there is a break between TTS and TSS.	But this prog can be used to display any multi lines.

option_list <- list(
	make_option(c("-f", "--file"),action = "store",type = "character",
		help = "The Input file"),
	make_option(c("-s", "--samplename"),action = "store",type = "character",
		help = "The Input sample name"),
	make_option(c("-t", "--title"),action = "store",type = "character",
		help = "The title of outimage"),
	make_option(c("-n", "--filename"),action = "store",type = "character",
		help = "The name of outimage"),
	make_option(c("-l","--length"),action = "store",type = "integer",default = "2",
		help = "The length bewteen start to end"),
	make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
		help = "The outdir")
	)
#parser <- OptionParser(usage = "%prog [option] file",option_list =option)

#arguments <- parse_args(parser,positional_arguments = 2)
opt <- parse_args(OptionParser(option_list = option_list,usage=usage))
#######必须有positional_arguments = TRUE才能不设定参数的个数
if(!is.character(opt$filename)){
	opt$filename = gsub(" ","_",opt$title)
}

setwd(opt$outdir)

###################################################################################
####加载包
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(RColorBrewer)
library("ggsci")
library("ggthemes")
###################################################################################
#colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')
colour1 <- brewer.pal(8,"Dark2")
colour2 <- brewer.pal(9,"Set1")
colour3 <- brewer.pal(12,"Paired")
colour<- c(colour1,colour2,colour3)
#colour <-c(colour3[2:12],colour2[8],colour1[4],colour1[7],colour1[8])
#colour <-c("#3C5488CC","#E64B35FF","#4DBBD5FF","#00A087FF",colour1[2],colour1[6:8])

###################################################################################
###数据处理
#多文本，则将其字符串分割处理，然后读取数据
file <- strsplit(opt$file,",")[[1]]
samplename <- strsplit(opt$samplename,",")[[1]]
data <- c()
#######将多个文本的数据（同一种格式），整合在一个文本中。
for(i in 1:length(file)){
	samplename[i] <- gsub('_distance2startcodon','',samplename[i])
	samplename[i] <- gsub('_distance2stopcodon','',samplename[i])
	samplename[i] <- gsub('_distance2tts','',samplename[i])
	samplename[i] <- gsub('_distance2tss','',samplename[i])
	text <- read.table(file = file[i],header = TRUE)
	# if(grepl("stopcodon",file[i]) || grepl("tts",file[i])){
	# # if(grepl("stopcodon",samplename[i]) || grepl("tts",samplename[i])){
	# # if(samplename[i] == "stopcodon" || samplename[i] == "tts"){
	# 	text[,1] <- text[,1]+opt$length*max(text[,1])+500
	# }
	if(grepl("_re_",file[i])){
	# if(grepl("stopcodon",samplename[i]) || grepl("tts",samplename[i])){
	# if(samplename[i] == "stopcodon" || samplename[i] == "tts"){
		text[,2] <- 0-text[,2]
	}
	text[,3] <- samplename[i]
	text[,4] <- colour[i]
	data <- rbind(data,text)
}
Max <- max(data[,2])
samplename
colname <- colnames(text)		
########标示X轴，Y轴名称。
colname[1] <- sub('X.','',colname[1])
title <- gsub('_',' ',opt$title)
#####################################################
#将众多的数据按行加在一起，形成长数据，并且以
#【1】列做横坐标
#【2】列做纵坐标
#【3】列做图例
######################################################

ablife_theme_line <- function(base_size = 12){
	library(grid)		####for using unit function
	theme(
		plot.title = element_text(size=12,lineheight = 10,colour="#000000",vjust = 1),

		axis.title.x = element_text(size=12,colour = "#000000",vjust = 0.5),
		axis.title.y = element_text(size=12,colour = "#000000",vjust = 1),
		axis.text.x = element_text(size = 12,colour = "#000000"),
		axis.text.y = element_text(size = 12 ,colour = "#000000"),
		axis.ticks.length = unit(0.1,"cm"),
		axis.ticks = element_line(colour = "#000000"), 

		# legend.title = element_text(size = 9),
		legend.title = element_blank(),
		legend.text = element_text(size = 9),
		legend.key.size = unit(0.5,"cm"),

		panel.background = element_rect(colour = "black")
		# panel.background = element_rect(fill = "white",colour = NA),
		# # panel.border = element_rect(size = 1,colour = "#8B8B8B",fill =NA),
		# # panel.grid.major = element_line(size=0.1,colour = "#BFBFBF"),
		# # panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
		# panel.border = element_rect(size = 1,colour = "#000000",fill =NA),
		# panel.grid.major.x = element_line(size=0.3,colour = "#000000"),
		# panel.grid.major.y = element_line(size=0.1,colour = "#909090",linetype = "dotted"),
		# panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
		)
}

theme_paper <- theme(
    # panel.border = element_rect(fill = NA,colour = "black"),
    panel.grid.major = element_line(size = 0.2,linetype = "dashed"),
    # panel.grid.major = element_blank(),
    
    panel.grid.minor = element_line(colour = "grey90",size = 0.2,linetype = "dashed"),
    # axis.text.x= element_text(vjust = 1,hjust = 1, angle = 45),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.key = element_rect(fill = 'white', color = 'white'),
    legend.key.size = unit(0.5, "cm"),
    # legend.position = c(0.82, 0.68),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.title = element_blank(),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(1,1,1,1,unit="cm"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, colour = "black")
)

#outname = gsub(" ","_",Image_name)
#png(file=paste(opt$filename,".png",sep=''),pointsize=40,width=1000,height=600)
ggplot(data)+ ####加载数据
		geom_line(aes(x = data[,1],y=data[,2],stat = "identity",group = data[,4],colour = data[,3]),size =1,position = "identity")+									#####基本作图函数
		labs(title = "",x="",y="")+				#####坐标轴调整（图片title，横纵坐标title）
		# ablife_theme_line()+
		# theme(
		# 	# legend.position = c(0.8, 0.78), #adjust for needing
		# 	# legend.position = "right",
		# 	legend.background = element_blank(),
		# 	legend.key = element_blank()
		# 	# legend.direction = "horizontal"
		# 	)+
		# scale_x_continuous(breaks = c(-1000,0,1000,1500,2500,3500),labels = c(-1000,"0\n(5')",1000,-1000,"0\n(3')",1000))+									#####指定x轴的标签
		# scale_y_continuous(limits=c(min(data[,2]),Max+10))+
		# scale_colour_hue(name=filename) ###系统默认
		# scale_x_continuous(trans="log2")+
		# scale_colour_nejm()+
		# scale_color_manual(values=c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00'))+
		# scale_colour_manual("Type",values = colour[1:length(samplename)])+			#########自由调整图形颜色和图例是在画图时自动加载的
		theme_bw() + theme_paper
ggsave(file = paste(opt$filename,".pdf",sep=''), width = 180,height = 90,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,".png",sep=''), width = 180,height = 90,dpi = 450,units = "mm")

# 05Line_multiSample_multiText

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

usage = "The prog is used to display the distribution of TSS,TTS,startcondon,stopcodon,and there is a break between TTS and TSS.	But this prog can be used to display any multi lines.

option_list <- list(
	make_option(c("-f", "--file"),action = "store",type = "character",
		help = "The Input file"),
	make_option(c("-s", "--samplename"),action = "store",type = "character",
		help = "The Input sample name"),
	make_option(c("-t", "--title"),action = "store",type = "character",
		help = "The title of outimage"),
	make_option(c("-n", "--filename"),action = "store",type = "character",
		help = "The name of outimage"),
	make_option(c("-l","--length"),action = "store",type = "integer",default = "2",
		help = "The length bewteen start to end"),
	make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
		help = "The outdir")
	)
# parser <- OptionParser(usage = "%prog [option] file",option_list =option)

# arguments <- parse_args(parser,positional_arguments = 2)
opt <- parse_args(OptionParser(option_list = option_list,usage=usage))
#######必须有positional_arguments = TRUE才能不设定参数的个数
if(!is.character(opt$filename)){
	opt$filename = gsub(" ","_",opt$title)
}

setwd(opt$outdir)

###################################################################################
####加载包
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(RColorBrewer)
library("ggsci")
library("ggthemes")
###################################################################################
#colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')
colour1 <- brewer.pal(8,"Dark2")
colour2 <- brewer.pal(9,"Set1")
colour3 <- brewer.pal(12,"Paired")
colour<- c(colour1,colour2,colour3)
#colour <-c(colour3[2:12],colour2[8],colour1[4],colour1[7],colour1[8])
#colour <-c("#3C5488CC","#E64B35FF","#4DBBD5FF","#00A087FF",colour1[2],colour1[6:8])

###################################################################################
###数据处理
#多文本，则将其字符串分割处理，然后读取数据
file <- strsplit(opt$file,",")[[1]]
samplename <- strsplit(opt$samplename,",")[[1]]
data <- c()
#######将多个文本的数据（同一种格式），整合在一个文本中。
for(i in 1:length(file)){
	samplename[i] <- gsub('_distance2startcodon','',samplename[i])
	samplename[i] <- gsub('_distance2stopcodon','',samplename[i])
	samplename[i] <- gsub('_distance2tts','',samplename[i])
	samplename[i] <- gsub('_distance2tss','',samplename[i])
	text <- read.table(file = file[i],header = TRUE)
	# if(grepl("stopcodon",file[i]) || grepl("tts",file[i])){
	# # if(grepl("stopcodon",samplename[i]) || grepl("tts",samplename[i])){
	# # if(samplename[i] == "stopcodon" || samplename[i] == "tts"){
	# 	text[,1] <- text[,1]+opt$length*max(text[,1])+500
	# }
	if(grepl("_re_",file[i])){
	# if(grepl("stopcodon",samplename[i]) || grepl("tts",samplename[i])){
	# if(samplename[i] == "stopcodon" || samplename[i] == "tts"){
		text[,2] <- 0-text[,2]
	}
	text[,3] <- samplename[i]
	text[,4] <- colour[i]
	data <- rbind(data,text)
}
Max <- max(data[,2])
samplename
colname <- colnames(text)		
########标示X轴，Y轴名称。
colname[1] <- sub('X.','',colname[1])
title <- gsub('_',' ',opt$title)
#####################################################
#将众多的数据按行加在一起，形成长数据，并且以
#【1】列做横坐标
#【2】列做纵坐标
#【3】列做图例
######################################################

ablife_theme_line <- function(base_size = 12){
	library(grid)		####for using unit function
	theme(
		plot.title = element_text(size=12,lineheight = 10,colour="#000000",vjust = 1),

		axis.title.x = element_text(size=12,colour = "#000000",vjust = 0.5),
		axis.title.y = element_text(size=12,colour = "#000000",vjust = 1),
		axis.text.x = element_text(size = 12,colour = "#000000"),
		axis.text.y = element_text(size = 12 ,colour = "#000000"),
		axis.ticks.length = unit(0.1,"cm"),
		axis.ticks = element_line(colour = "#000000"), 

		# legend.title = element_text(size = 9),
		legend.title = element_blank(),
		legend.text = element_text(size = 9),
		legend.key.size = unit(0.5,"cm"),

		panel.background = element_rect(colour = "black")
		# panel.background = element_rect(fill = "white",colour = NA),
		# # panel.border = element_rect(size = 1,colour = "#8B8B8B",fill =NA),
		# # panel.grid.major = element_line(size=0.1,colour = "#BFBFBF"),
		# # panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
		# panel.border = element_rect(size = 1,colour = "#000000",fill =NA),
		# panel.grid.major.x = element_line(size=0.3,colour = "#000000"),
		# panel.grid.major.y = element_line(size=0.1,colour = "#909090",linetype = "dotted"),
		# panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
		)
}

theme_paper <- theme(
    # panel.border = element_rect(fill = NA,colour = "black"),
    panel.grid.major = element_line(size = 0.2,linetype = "dashed"),
    # panel.grid.major = element_blank(),
    
    panel.grid.minor = element_line(colour = "grey90",size = 0.2,linetype = "dashed"),
    # axis.text.x= element_text(vjust = 1,hjust = 1, angle = 45),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.key = element_rect(fill = 'white', color = 'white'),
    legend.key.size = unit(0.5, "cm"),
    # legend.position = c(0.82, 0.68),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.title = element_blank(),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(1,1,1,1,unit="cm"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, colour = "black")
)

#outname = gsub(" ","_",Image_name)
#png(file=paste(opt$filename,".png",sep=''),pointsize=40,width=1000,height=600)
ggplot(data)+ ####加载数据
		geom_line(aes(x = data[,1],y=data[,2],stat = "identity",group = data[,4],colour = data[,3]),size =1,position = "identity")+									#####基本作图函数
		labs(title = "",x="",y="")+				#####坐标轴调整（图片title，横纵坐标title）
		# ablife_theme_line()+
		# theme(
		# 	# legend.position = c(0.8, 0.78), #adjust for needing
		# 	# legend.position = "right",
		# 	legend.background = element_blank(),
		# 	legend.key = element_blank()
		# 	# legend.direction = "horizontal"
		# 	)+
		# scale_x_continuous(breaks = c(-1000,0,1000,1500,2500,3500),labels = c(-1000,"0\n(5')",1000,-1000,"0\n(3')",1000))+									#####指定x轴的标签
		# scale_y_continuous(limits=c(min(data[,2]),Max+10))+
		# scale_colour_hue(name=filename) ###系统默认
		# scale_x_continuous(trans="log2")+
		# scale_colour_nejm()+
		# scale_color_manual(values=c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00'))+
		# scale_colour_manual("Type",values = colour[1:length(samplename)])+			#########自由调整图形颜色和图例是在画图时自动加载的
		theme_bw() + theme_paper
ggsave(file = paste(opt$filename,".pdf",sep=''), width = 180,height = 90,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,".png",sep=''), width = 180,height = 90,dpi = 450,units = "mm")

# 06PCA
#####参数获取
#####################################################################################
#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
  make_option(c("-f", "--infile"),action = "store",type = "character", help = "The Input file"),
  make_option(c("-p", "--group"),action = "store",type = "character", help = "The group file"),
  make_option(c("-e", "--ellipsetype"),action = "store",type = "character", default = "convex",help = "ellipsetype: t,norm,euclid,convex"),
  make_option(c("-c","--column"),action = "store",type = "integer",default = 2,
              help = "the column number for pca group"),
  make_option(c("-t", "--text"),action = "store_true",type = "logical",default = FALSE, help = "annotate for each point,  [default]"),
  make_option(c("-n", "--outname"),action = "store",type = "character",default = "PCA", help = "The name of outimage;default = PCA")
)
opt <- parse_args(OptionParser(option_list = option_list))

#opt$infile <- "/users/chengc/dev2016/graphicwork/Macaque/3_pca/All_sample_merged_RPKM_edSymbol.xls.mRNA"
#opt$infile <- "/users/chengc/dev2016/graphicwork/Macaque/3_pca/All_sample_LncRNA_exp_RPKM_annot_gencode.xls.lncRNA"
#opt$group <- "/users/chengc/dev2016/graphicwork/Macaque/3_pca/Sample_group.txt"

library("factoextra")
library("FactoMineR")
library("stats")
library("ggthemes")
NO_REUSE = F

#opt$file = "Sample_correlation.dat"
#opt$file = "test.txt"

#get the filename to use later

#try to reuse earlier-loaded data if possible

#print('Reading matrix file.')
primary_data = read.table(opt$infile, header=T, com='', sep="\t", row.names=1, check.names=F)
data.class <- read.table(opt$group,header=T,row.names=1,sep='\t')
primary_data <- primary_data[rowSums(primary_data) != 0,]
primary_data <- primary_data[rowMeans(primary_data) != max(primary_data),]
#head(primary_data)

data_matrix <- as.matrix(primary_data)
data_matrix <- t(data_matrix)
data_matrix1 <- data_matrix[,colSums(data_matrix[,1:ncol(data_matrix)]) > 0]

data_matrix1 <- data_matrix1[order(row.names(data_matrix1)),]
data.class <- data.class[order(row.names(data.class)),]

#data.class

#dd <- merge(data_matrix1[,c(1:10)],data.class,by=0,all.x=TRUE)


data.pca <- prcomp(data_matrix1,  scale = TRUE)
data.pca$x
cat(c('sample\t'),file=paste(opt$outname,".PCA.xls",sep=''))
write.table(data.pca$x, file=paste(opt$outname,".PCA.xls",sep=''), quote=F, append=T, sep='\t')


pc <- predict(data.pca)

et <- opt$ellipsetype

a <- as.numeric(opt$column)
a <- a - 1 
data.class[,a]
#data.pca <- PCA(data_matrix1, graph = FALSE)
#f <- fviz_pca_ind(data.pca, label="none", habillage=data.class[,a],addEllipses=TRUE, ellipse.level=0.95,ellipse.type="t",pointsize=4) 
#"t", "norm", "euclid","convex"
p <- fviz_pca_ind(data.pca, label="none", habillage=data.class[,a],addEllipses=TRUE, ellipse.level=0.95,ellipse.type=et,pointsize=2) 

if (opt$text){
  p <- p + 
    geom_text(aes(label=row.names(data.pca$x)),hjust=0, vjust=0)
    }

  p <- p + 
   #scale_shape_manual(values = c(20,18,20,18,20,18,20,18,20)) +
   scale_color_brewer(palette="Set1", direction=-1) +
   scale_fill_brewer(palette="Set1", direction=-1) + 
   theme_bw() +
    theme(
      title = element_text(size = 12),
      # legend.position = c(0.78, 0.7),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.key.height = unit(0.7, "cm"),
      legend.key.width = unit(0.7, "cm"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12))

ggsave(paste(opt$outname,"_0.5.pdf",sep = ""), width = 150, height = 100, units = "mm")

#pdf(file = paste(opt$outname,"_PCA.pdf",sep=''))
#tiff("PCA_plot_time.tif",2400,2400,compression ='lzw',res=300)
#print(p)
#dev.off()



#f <- fviz_pca_biplot(data.pca, 
#habillage = data.class[,1], addEllipses = TRUE,
#col.var = "red", alpha.var ="cos2",
#label = "var") +
#scale_color_brewer(palette="Dark2")+
#theme_minimal()
#f

#07Sample_correlation

rm(list=ls())
library(gplots)
#library('Cairo')

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
	stop('Your input arguments is wrong!\n
		args1:\t the RPKM file\n
		args2:\t column of first sample\n
		args3:\t column of last sample\n
		args4:\t the out directory\n\n')
}
infile <- args[1]
i <- args[2]
j <- args[3]
outdir <- args[4]

data <- read.table(file=infile,header=T,sep='\t',comment.char='')
columnName <- colnames(data)
#Data <- data[,i:j]
Data <- log(data[,i:j]+1,2)
CorData <- cor(data[,i:j])
CorData <- cor(log(data[,i:j]+1,2))
CorData <- CorData^2

setwd(outdir)
cat(c("Sample","\t"),file="Sample_correlation.txt")
write.table(CorData,file="Sample_correlation.txt",
		append=T,quote=F,sep='\t')
png(file='Sample_correlation.png',pointsize =20,width=1600,height=1600,)
par(mar=c(7.1,1.1,0,3.1),omi=c(2.3,0,0.1,2.3))
heatmap.2(CorData,col=cm.colors(256),		#heat.colors bluered(256) greenred(256) cm.colors(256)
		#Rowv=FALSE,
		#Colv=FALSE,#		don't cluster the heatmap
		trace='none',revC=F,#dendrogram='none',
		#hclustfun = hclust,
		na.rm = F,
		#ColSideColors = patientcolors,
		#RowSideColors=RowCol,
		main='Sample Correlation',cex.main=6,
		cexRow=1.5,
		font.lab = 2,font.axis = 2,
		notecex=3,cexCol=1.5,
		#scale = "row",
		scale = 'none',
		#labRow = rownames(X),
		#labCol = colnames(Data),
		symkey=FALSE,key=T,keysize=1,density.info = "none"
)
dev.off()

# 08Scatter_multiSample
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

usage = "The prog is used to display the cumulatived data of FPKM.

example: Rscript /public/pipeline/ablifegraphic/latest/ablife-R/Scatter/latest/Scatter_pas.r 
	-f B246d_pac_50.txt
	-t B246d_pac_50
	-n B246d_pac_50
	"

option_list <- list(
	make_option(c("-f", "--file"),action = "store",type = "character",
		help = "The Input file"),
	make_option(c("-t","--title"),action = "store",type = "character",
		help = "The title of the plot"),
	make_option(c("-n", "--filename"),action = "store",type = "character",
		help = "The name of outimage"),
	make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
		help = "The outdir")
	)
opt <- parse_args(OptionParser(option_list = option_list)) #####解释器

setwd(opt$outdir)						#### Set the outpath

###################################################################################
####加载包 （Load Package）
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(gtable)
library(grid)
library(RColorBrewer)
###################################################################################
####加载颜色（Load Color）
###################################################################################
colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')				###一共14中color
#colour1 <- brewer.pal(8,"Dark2")

######Instruction for data#######
#first column is the X axis
#from second to end is the data for ploting
#first row is the head to declare the data
#################################

###################################################################################
##数据处理（Deal with Data）
title <- gsub('_',' ',opt$title)

data <- read.table(file= opt$file,header = T)#####	must with header
colname <- colnames(data)
dim_data <- dim(data)
newdata <- melt(data,id.vars=colname[1],measure.vars=c(colname[2:dim_data[2]]),variable.name="Base")
#######变形数据后直接给variablfe命名即可

##################################################################################
###Plot Theme for ABLife 
###theme(),Tha last term without comma
##################################################################################
ablife_theme_line <- function(base_size = 12){
	library(grid)		####for using unit function
	theme(
		plot.title = element_text(size=12,colour="#000000",hjust = 0.5),

		axis.title.x = element_text(size=12,colour = "#000000"),
		axis.title.y = element_text(size=12,colour = "#000000"),
		axis.text.x = element_text(size = 12,colour = "#000000"),
		axis.text.y = element_text(size = 12 ,colour = "#000000"),

		legend.title = element_text(size = 12),
		legend.text = element_text(size = 12),
		legend.key.size = unit(1,"cm"),

		panel.background = element_rect(fill = "white",colour = NA),
		panel.border = element_rect(size = 1,colour = "black",fill =NA)
		# panel.border = element_rect(size = 1,colour = "#8B8B8B",fill =NA),
		# panel.grid.major = element_line(size=0.5,colour = "#BFBFBF"),
		# panel.grid.minor = element_line(size=0.1,colour = "#7F7F7F")
		)
}

##################################################################################
###Plot by ggplot2

ggplot(newdata)+ ####加载数据
		geom_point(aes(x = newdata[,1],y=newdata[,3],stat = "identity",group = 1,fill= Base,colour=Base),shape = 2,size =3,position = "identity")+									#####基本作图函数
		labs(title = title,y="Frequency of each base(%)",x="Location")+				#####坐标轴调整（图片title，横纵坐标title）
		ablife_theme_line()

###################################################################################
###Save Plot File 
ggsave(file = paste(opt$filename,"_Scatter.pdf",sep=''), width = 180,height = 120,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,"_Scatter.png",sep=''), width = 180,height = 120,dpi = 450,units = "mm")

# 09Violin_plot
#!/usr/bin/env Rscript

#####################################################################################
#####参数获取
#####################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

usage = "The prog is used to display Violin of data,the Violin display density distribution of data,and a line to show the median data.
	It's used to compare several distributions while they're placed side by side.A violin plot is a kernel density estimate,mirrored so that it forms a symmetrical shape.
Rscript %prog -f -t -n 
example: Rscript %prog
	-f CLIP_Peak_len 
	-t CLIP_Peak_len
	-n CLIP_Peak_len
"

option_list <- list(
	make_option(c("-f", "--file"),action = "store",type = "character",
		help = "The Input file"),
	make_option(c("-t","--title"),action = "store",type = "character",
		help = "The title of the plot"),
	make_option(c("-x","--xaxis"),action = "store",type = "character",default = "Location",
		help = "The name of x-axis "),
	make_option(c("-y","--yaxis"),action = "store",type = "character",default = "Frequency of each base(%)",
		help = "The name of y-axis"),
	make_option(c("-n", "--filename"),action = "store",type = "character",
		help = "The name of outimage"),
	make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
		help = "The outdir")
	)
opt <- parse_args(OptionParser(option_list = option_list)) #####解释器

setwd(opt$outdir)		# Set the Outpath

###################################################################################
####加载包 （Load Package）
###################################################################################
library(ggplot2)
library(reshape2)
library(plotrix)
library(methods)
library(gtable)
library(grid)
library(RColorBrewer)
###################################################################################
####加载颜色（Load Color）
###################################################################################
colour <- c('#85A2EF','#D285EF','#A2EF85','#4682B4','#A0522D','#87CEEB','#6B8E23','#6A5ACD','#E59B95','#EFD285','#B4B643','#2E9AFE','#A1DDBB','#FF8C00')
#colour1 <- brewer.pal(8,"Dark2")

######Instruction for data#######
#first column is the X axis---discrete data
#second column is the Y axis to plot -- continuous data
#the data is to display the summary about X axis data in density.
#################################


###################################################################################
##数据处理（Deal with Data）
data <- read.table(file= opt$file,header = T)#####	must with header
title <- gsub('_',' ',opt$title)
colname <- colnames(data)
dim_data <- dim(data)

Item <- factor(data[,1],levels = data[,1])
x_axis = gsub('_',' ',opt$xaxis)
y_axis = gsub('_',' ',opt$yaxis)
#######变形数据后直接给variablfe命名即可

##################################################################################
###Plot Theme for ABLife 
###theme(),Tha last term without comma
##################################################################################
ablife_theme_line <- function(base_size = 12){
	library(grid)		####for using unit function
	theme(
			plot.title = element_text(size = 12,lineheight = 100,colour = "black",hjust = 0.5),
			# axis.text.x=element_text(angle=75,hjust=1,size = 5,colour = "black"),
			axis.title.x = element_text(size = 12,colour = "black"),
			axis.text.y = element_text(size = 12 ,colour = "black"),
			axis.title.y = element_text(size = 12,colour = "black"),
			panel.background = element_rect(colour = "black"),
			legend.title = element_text(size = 12),
			legend.text = element_text(size = 12),
			strip.text.x = element_text(colour = "black",size = "8"),
			strip.background = element_rect(colour = "black")
			)
}

##################################################################################
###Plot by ggplot2
ggplot(data)+ ####加载数据
		geom_violin(aes(x = Item,y=data[,2],stat = "identity",fill=Type))+									#####基本作图函数
		labs(title = title,y = y_axis , x = x_axis)+				#####坐标轴调整（图片title，横纵坐标title）
		theme(axis.text.x=element_text(angle=60,hjust=1,size = 12,colour = "black"),
			axis.ticks.x = element_blank(),
			legend.position="none"
			# legend.position = c(0.3, 0.92),
			# legend.direction = "horizontal"
			) + 
		ablife_theme_line()+
##scale_y_continuous(limits=c(0,100))+									#####指定y轴的范围
		scale_colour_hue(name=opt$legend_name)
#scale_colour_manual("Sample Name",values = colour[1:length(sample)])
#scale_fill_discrete(name="Sample_name")	##可以独立调整图例的名字，label，breaks

###################################################################################
###Save Plot File 
ggsave(file = paste(opt$filename,"violin.pdf",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")
ggsave(file = paste(opt$filename,"violin.png",sep='_'), width = 105,height =150 ,dpi = 450,units = "mm")

#10soger
#!/usr/bin/env p3
#-*- coding: utf-8 -*-
#gencos
#Utility for visualizing NGS data along gene models
##import re, os, sys, logging, time, datetime

from optparse import OptionParser, OptionGroup

#reload(sys)
#sys.setdefaultencoding('utf-8')

import matplotlib
from numpy import *

#Use PDF backend
matplotlib.use("pdf")

sys.path.insert(1, os.path.split(os.path.realpath(__file__))[0] + "/./")
from lib import Configuration,Tracks,Plot
#from .lib.Plot import plottracks

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


def run(conf_filename, region, output_dir, plot_title=None, plot_label=None):
    """
    plot tracks between region
    """
    plotconfig = Configuration.loadconfiguration(conf_filename, region)
    (tracks,chrome,tx_start, tx_end, graphcoords, graphToGene) = Tracks.maketracks(plotconfig)
    # print(tracks)

    fig = Plot.initFig(plotconfig,tracks)
    Plot.plottracks(fig,plotconfig,tracks,chrome,tx_start, tx_end, graphcoords, graphToGene, plot_title ,plot_label)


def greeting():
    print("Ablife Genome View: Visualize spliced RNA-Seq reads along gene models. ")
    print("See --help for usage.\n")
    print("Manual available at: \n")


def main():
    parser = OptionParser()
    parser.add_option('-c', '--conf', dest='conf', action='store', type='string', help='config file')
    parser.add_option('-r', "--plot-region", dest="plot_region", nargs=1, default=None, help="Plot read densities for a given region. "
                                                                                             "format: chr:start-end")
    parser.add_option('-t', "--plot-title", dest="plot_title", default=None, nargs=1, help="Title of plot: a string that will be displayed at top of plot. Example: "
                                                                                           "--plot-title \"My favorite gene\".")
    parser.add_option('-l', "--plot-label", dest="plot_label", default=None, nargs=1, help="Plot label. If given, plot will be saved in the output directory as "
                                                                                           "the plot label ending in the relevant extension, e.g. <plot_label>.pdf. "
                                                                                           "Example: --plot-label my_gene")
    parser.add_option('-o', "--output-dir", dest="output_dir", nargs=1, default="./", help="Output directory.")
    (options, args) = parser.parse_args()

    if options.conf is None:
        greeting()
        sys.exit(1)

    output_dir = os.path.abspath(os.path.expanduser(options.output_dir))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    plot_title = options.plot_title
    plot_label = options.plot_label

    if options.conf != None:
        region = options.plot_region
        conf_filename = os.path.abspath(os.path.expanduser(options.conf))
        run(conf_filename, region, output_dir, plot_title=plot_title, plot_label=plot_label)


if __name__ == '__main__':
    main()

