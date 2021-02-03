library(tidyverse)

setwd("/home/niek/Documents/analyses/ChIP-Seq/SET1B-HIF1-2-H3K4Me3")

df <- read.csv(file="mapped_read_count_dedup.txt", sep=" ", skip=1, header=FALSE, col.names=c("file", "reads"))
df$file <-gsub('.{1}$', '', df$file)

df.non.input <- df[!grepl("input", df$file), ]
df.input <- df[grepl("input", df$file), ]

df.pairs <- read.csv(file="downsample.conf")

df.master <- merge(x=df.non.input, y=df.pairs,by="file") #matches chip and input file

names(df.input)[names(df)=="file"] <- "input"
df.master <- merge(x=df.master, y=df.input,by="input")
names(df.master) <- c("input","chip","reads.chip","reads.input")

df.master$`chip/input` <- df.master$reads.chip / df.master$reads.input
df.master$`input/chip` <- df.master$reads.input / df.master$reads.chip
