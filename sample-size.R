#Calculates scaling factors for input samples to be downsampled to read numbers of ChIP samples
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
work.dir <- args[1]
setwd(work.dir)
read.count.file <- args[2]

df <- read.csv(file=read.count.file, 
               sep=" ", 
               skip=1, 
               header=FALSE, 
               col.names=c("file", "reads"))
df$file <-gsub('.{1}$', '', df$file)

df.non.input <- df[!grepl("input", df$file), ]
df.input <- df[grepl("input", df$file), ]

df.pairs <- read.csv(file="downsample.conf")

df.master <- merge(x=df.non.input, 
                   y=df.pairs,
                   by="file") #matches chip and input file

names(df.input)[names(df)=="file"] <- "input"
df.master <- merge(x=df.master, 
                   y=df.input,
                   by="input")
names(df.master) <- c("input","chip","reads.chip","reads.input")

df.master$`chip/input` <- df.master$reads.chip / df.master$reads.input #scaling factor for input
df.master$factor <- round((df.master$`chip/input`*1000),0) #converts to input for samtools view (seed.factor)
df.master <- df.master[df.master$`chip/input` < 1, ] #selects rows with scaling factor < 1
seed <- 3
df.master$factor <- paste0(seed,".",df.master$factor) #converts to input for samtools view (seed.factor)

df.output <- subset(df.master, select=c("input","chip","factor"))

#prepare output file name after downsampling of bam input file (a unique file name for each chip-input pair is made)
df.output$output.file <- sub("bam/*","downsample/",df.output$input)
df.output$output.file <- sub("-.*","-",df.output$output.file)
df.output$chip.temp <- sub("bam/*","",df.output$chip)
df.output$chip.temp <- sub("-.*","-",df.output$chip.temp)
df.output$extension <- "ds.bam"
df.output$output.file.final <- paste0(df.output$output.file,df.output$chip.temp,df.output$extension)

#select columns for file output
df.output <- subset(df.output, select=c("factor","input","output.file.final","chip"))

#write data frame to file
write.table(df.output, file="scaling_factors.txt",
            sep=",",
            row.names=FALSE)