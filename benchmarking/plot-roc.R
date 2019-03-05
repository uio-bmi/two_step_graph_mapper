#!/usr/bin/env Rscript
#
# This script a a copy of https://github.com/vgteam/vg/blob/fec51877c7fd21ffbf60cc4c6a157702c112714a/scripts/plot-roc.R
#
# plot-roc.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include>]
#
# plots a pseudo-ROC that allows the comparison of different alignment methods and their mapping quality calculations
# the format is clarified in the map-sim script, and should be a table (tab separated) of:
#      correct   mq    score   aligner
# where "correct" is 0 or 1 depending on whether the alignnment is correct or not and "aligner" labels the mapping method
#
# This is not a true ROC because we are not purely plotting the binary classification performance of
# each of the methods' mapping quality calculation over the same set of candidate alignments.
# Rather, we are mixing both the alignment sensitivity of the method with the MQ classification performance.
# As such we do not ever achieve 100% sensitivity, as we have effectively scaled the y axis (TPR) by the total
# sensitivity of each mapper.


mappercolor <- function(name){
    if (name == "vg"){
        return("#3355FF")
    }
    else if(name == "vg_mitty"){
        return("#79B4FF")
    }
    else if(name == "seven_bridges"){
        return("#BC35C2")
    }
    else if(name == "seven_bridges_mitty"){
        return("#D797DA")
    }
    else if(name == "bwa"){
        return("#B81B1B")
    }
    else if(name == "bwa_untuned"){
        return("#D68B8B")
    }
    else if(name == "two_step_graph_mapper_traversemapped"){
        return("#00768c")
    }
    else if(name == "two_step_graph_mapper_graph_minimap"){
        return("#00598C")
    }
    else if(name == "two_step_graph_mapper_linearmapped"){
        return("#89D2D9")
    }
    else if(name == "two_step_graph_mapper_vg"){
        return("#de9000")
    }
    else if(name == "sb_pe"){
        return("#000000")
    }
    else if(name == "bwa_pe"){
        return("#555555")
    }
    else if(name == "hisat"){
        return("#148046")
    }
    else if(name == "hisat_mitty"){
        return("#5BC38B")
    }
    else{
        stop("Mapper has no color associated with it: ", name)
    }


}

list.of.packages <- c("tidyverse", "ggrepel", "magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=T)
require("tidyverse")
require("ggrepel")
require("scales") # For squish
require("magrittr") # For squish

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name
dat <- read.table(commandArgs(TRUE)[1], header=T)

if (length(commandArgs(TRUE)) > 2) {
    # A set of aligners to plot is specified. Parse it.
    aligner.set <- unlist(strsplit(commandArgs(TRUE)[3], ","))
    # Subset the data to those aligners
    dat <- dat[dat$aligner %in% aligner.set,]
}

# Determine the order of aligners, based on sorting in a dash-separated tag aware manner
aligner.names <- levels(dat$aligner)
print(aligner.names)

mappercolors <- aligner.names
i = 0
for (name in aligner.names){
    i = i + 1
    mappercolors[i] = mappercolor(name)
}

print(mappercolors)

name.lists <- aligner.names %>% (function(name) map(name,  (function(x) as.list(unlist(strsplit(x, "-"))))))
# Transpose name fragments into a list of vectors for each position, with NAs when tag lists end early
max.parts <- max(sapply(name.lists, length))
name.cols <- list()
for (i in 1:max.parts) {
    name.cols[[i]] <- sapply(name.lists, function(x) if (length(x) >= i) { x[[i]] } else { NA })
}
name.order <- do.call(order,name.cols)
dat$aligner <- factor(dat$aligner, levels=aligner.names[name.order])

dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
dat.roc <- dat %>%
    mutate(Positive = correct == 1, Negative = correct == 0) %>%
    group_by(aligner, mq) %>%
    summarise(Positive = sum(Positive), Negative = sum(Negative)) %>%
    arrange(-mq) %>%
    mutate(Total=sum(Positive+Negative)) %>%
    mutate(TPR = cumsum(Positive) / Total, FPR = cumsum(Negative) / Total)

# We want smart scales that know how tiny a rate of things we can care about
total.reads <- max(dat.roc$Total)
min.log10 <- floor(log10(1/total.reads))
max.log10 <- 0
# Work out a set of bounds to draw the plot on
range.log10 <- min.log10 : max.log10
range.unlogged = 10^range.log10

dat.plot <- ggplot(dat.roc, aes( x= FPR, y = TPR, color = aligner, label=mq)) +
    geom_line() + geom_text_repel(data = subset(dat.roc, mq %% 10 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
    geom_point(aes(size=Positive+Negative)) +
    scale_color_manual(values=mappercolors, guide=guide_legend(title=NULL, ncol=2)) +
    scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
    scale_x_log10(limits=c(range.unlogged[1],range.unlogged[length(range.unlogged)]), breaks=range.unlogged, oob=squish) +
    geom_vline(xintercept=1/total.reads) + # vertical line at one wrong read
    theme(legend.position = "none")

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=4)
