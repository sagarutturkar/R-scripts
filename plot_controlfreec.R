library(ggplot2)
library(dplyr)
library(optparse)
library(ggthemes)

# Script to plot control-freec plot - combined chromosomes
# syntax: Rscript --vanilla  plot_freec.R  -i <control-freec-ratio-data>  -p <outpath> -n <Name>")
# Example: Rscript --vanilla  plot_freec.R  -i GRC204500540_sorted.bam_ratio.txt  -p Z:/PCCR/Savran_Cagri/02142020_singlecell_CNV/controlfreec/test_plot/ -n GRC204500540

option_list = list(
  
  make_option(c("-i", "--infile"), type="character", default="NULL", 
              help="ratio file from control-freec", metavar="character"),
  
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path for the output directory", metavar="character"),
  
  make_option(c("-n", "--name"), type="character", default="NULL", 
              help="Name prefix for output", metavar="character"),
  
  make_option(c("-l", "--height"), type="numeric", default="700", 
              help="height for output plot", metavar="numeric"),
  
  make_option(c("-w", "--width"), type="numeric", default="5000", 
              help="width for output plot", metavar="numeric")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#set output directory and read data files
setwd(opt$path)
data = read.table(file = opt$infile, header = TRUE, stringsAsFactors = F, sep = "\t", quote = "")

ploidy = 2
max_level = 3

# function to format data
format_data <- function(data) {
  
  #remove mitochondrial chromosome (mitocondrial)
  data = dplyr::filter(data, !Chromosome == "MT")
  
  # Add new column - "key" denoting copy number call
  data = dplyr::mutate(data, key = case_when(CopyNumber > 2 ~ 'gain',
                                             CopyNumber < 2 ~ 'loss',
                                             TRUE ~ 'neutral'))
  
  # Add new column - "new_ratio" - adjust ratio as in control-freec script
  data = dplyr::mutate(data, new_ratio = case_when(Ratio > max_level ~ max_level,
                                                   TRUE ~ Ratio))
  
  chr_vec = c("1", "2", "3", "4", 
              "5", "6", "7", "8", 
              "9", "10", "11", "12", 
              "13", "14", "15", "16", 
              "17", "18", "19", "20", 
              "21", "22", "X", "Y" )
  
  #factor columns to set the correct order in the plot
  data$Chromosome = factor(data$Chromosome, levels = chr_vec)
  data$key = factor(data$key, levels = c("loss", "neutral", "gain"))
  
  return(data)

}

# function to generate ratio plot
ratio_plot <- function(data) {
  
  #set colors
  col_sample =  c("blue", "green", "red")
  
  p = ggplot(data, aes(x=Start, y=new_ratio*ploidy, color = key)) + 
    geom_point(size= 0.5) +
    coord_cartesian(ylim = c(0, max_level*ploidy)) +
    facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
    #facet_grid(. ~ Chromosome) +
    scale_colour_manual(values=col_sample) +
    labs(title="Normalized copy number Profile", y="Normalized Copy Number", x="Position") +
    theme_clean() +
    theme(panel.border=element_rect(colour="black",size=1, fill=NA))
  
  plot = p + theme(legend.position="bottom",
            plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            panel.spacing = unit(1, "lines"),
            strip.text.x = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=5, stroke=1, alpha = 1)))
  
  return(plot)
  
}

#function to generate log2 plot
log2_plot <- function(data) {
  
  #set colors
  col_sample =  c("blue", "green", "red")
  
  q = ggplot(data, aes(x=Start, y=log2(Ratio), color = key)) + 
    geom_point(size= 0.5) +
    coord_cartesian(ylim = c(-6, 6)) +
    scale_y_continuous(breaks=c(-6, -4, -2, 0, 2, 4, 6)) +
    facet_grid(. ~ Chromosome, scales = "free_x", space = "free_x") +
    #facet_grid(. ~ Chromosome) +
    scale_colour_manual(values=col_sample) +
    labs(title="Normalized copy number Profile (log2)", y="Normalized Copy Number (log2)", x="Position") +
    theme_clean() +
    theme(panel.border=element_rect(colour="black",size=1, fill=NA))
  
  plot = q + theme(legend.position="bottom",
            plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            panel.spacing = unit(1, "lines"),
            strip.text.x = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=5, stroke=1, alpha = 1)))
  
  return(plot)
  
}


data = format_data(data)
filename = paste(opt$name,"plot.png", sep = "_")

#Write regular plot to output
png(filename = filename, res = 200, 
    width = opt$width, height = opt$height, pointsize = 1,  bg = "white")

ratio_plot(data)

dev.off()

filename = paste(opt$name,"plot_log2.png", sep = "_")

#write log2 plot to output
png(filename = filename, res = 200, 
    width = opt$width, height = opt$height, pointsize = 1,  bg = "white")

log2_plot(data)

dev.off()


