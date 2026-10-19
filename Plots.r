library(tidyverse)
library(lubridate)
library(patchwork)
library(ggthemes)
library(ggeffects)
library(brms)
library(cowplot)
library(magick)
library(data.table)
source("hbdseg/theme_simple.R")
library(windowscanr)
library(gt)
library(grid)
library(ggplotify)
library(viridis)
library(scales)
library(ggridges)
library(gt)
library(lme4)
library(broom)
options(scipen=999)
library(RZooRoH)

ppcol <- '#7570b3'
pscol <- '#d95f02'
agecol <- '#e41a1c'
intcol <- c('#dadaeb','#9e9ac8','#3f007d')
intcol2 <- c('#3f007d','#9e9ac8','#bcbddc')
intcol3 <- c('#7f2704','#d94801','#fdae6b')
intcol2.2 <- c('#3f007d','#bcbddc')
intcol3.2 <- c('#7f2704','#fdae6b')
f2_size <- 22
f2_tsize <- 22
tsize <- 28

##################################################################################
##################################################################################
##################################################################################
######## Figure 1:
### SAM + female and pup counts
sams <- read.table('../../../../SAM.txt', h = T)

fems <- read_csv('../Final_models/FurSeal_BI_SpecialStudyBeach_Counts.csv')
fems$Date <- dmy(fems$Date)
yrs <- unique(year(fems$Date))
fems$YEAR<-rep(0,dim(fems)[1])

for (y in yrs){
	for (d in 1:length(fems$Date)){
		if ((year(fems$Date[d]) == y+1) & month(fems$Date[d]) <= 4) {
			fems$YEAR[d] <- y
		} else if ((year(fems$Date[d]) == y) & month(fems$Date[d]) >= 4) {
			fems$YEAR[d] <- y
		}
	}
}

fems <- as.data.frame(fems)
colnames(fems)[4] <- 'NewPups'

adf <- fems |>
	group_by(YEAR) |>
	summarize(Females = max(Females, na.rm = T), Pups = sum(NewPups, na.rm = T)) |>
	ungroup()

trnd <- sams |>
	left_join(adf)

fem_dens <- read.csv('../../../matureFemales19952025.csv')
trnd <- trnd |>
	left_join(fem_dens, by = c('YEAR' = 'year'))
colnames(trnd)[9] <- 'Fem_count'

trnd_bu <- trnd
trnd <- trnd |>
	select(YEAR, Pups, Fem_count) |>
	pivot_longer(cols = !YEAR, names_to = 'age', values_to = 'count')

# Truncated
trnd_bu <- trnd_bu[-c(1:27),]
trnd <- trnd[-c(1:54),]

### lifetime repr. succ.
data <- read.table('PSdataset.txt', h = T)

pups <- data |>
	mutate(Pup_survived = replace(Pup_survived, is.na(Pup_survived), 0)) |>
	mutate(Pup_recruited = replace(Pup_recruited, is.na(Pup_recruited), 0)) |>
	group_by(SampleID) |>
	summarize(PP = sum(Pup_produced), PS = sum(Pup_survived), PR = sum(Pup_recruited))

pups_bu <- pups
pups$PP <- factor(pups$PP, levels = c(1:max(pups$PP)))
pups$PS <- factor(pups$PS, levels = c(0:max(pups$PS)))
pups$PR <- factor(pups$PR, levels = c(0:max(pups$PR)))
ppc <- table(pups$PP)
psc <- table(pups$PS)
prc <- table(pups$PR)
dfpp<-data.frame(Counts=names(ppc),PP=as.numeric(ppc))
dfps<-data.frame(Counts=names(psc),PS=as.numeric(psc))
dfpr<-data.frame(Counts=c(names(prc),as.character(4:12)),PR=c(as.numeric(prc),rep(0,9)))
dfpp$Counts<-factor(dfpp$Counts,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13'))
dfps$Counts<-factor(dfps$Counts,levels=c('0','1','2','3','4','5','6','7','8','9','10','11','12'))
dfpr$Counts<-factor(dfpr$Counts,levels=c('0','1','2','3','4','5','6','7','8','9','10','11','12'))

p1 <- ggplot(trnd_bu, aes(YEAR, ANN)) +
	geom_line(col = 'grey50') +
	geom_point(col = 'grey50', size = 1) +
	geom_smooth(method = 'lm', col = 'black') +
	labs(x='Year',y='Mean SAM', tag = 'B') +
	ylim(-3,NA) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),
	 axis.text.x = element_blank(), 
	 axis.title.x = element_blank(),
	 axis.title.y = element_text(margin = margin(r = 4)),
	 axis.ticks.x=element_blank(), axis.text.y = element_text(color = "black", size = 12),
	 plot.tag = element_text(size = 18))

p2 <- ggplot(trnd, aes(x = YEAR, y = count, color = age)) +
	geom_line() +
	geom_point(size = 1) +
	geom_smooth(method = 'lm',show.legend = FALSE) +
	labs(x='Year',y='Counts', tag = 'C') +
	theme_bw() +
	scale_color_manual(name = '',values=c('grey50', ppcol), labels = c('Breeding females','Pups produced'))+
	scale_fill_manual(values=c('grey50', ppcol), name="",labels=NULL) +
	ylim(NA,1100) +
	theme(axis.line = element_line(colour = "black"),
	legend.position=c(0.2,0.2),
	legend.text = element_text(size = 12),
	legend.key = element_rect(fill = "transparent", linewidth=0, color='transparent'),
	legend.background = element_rect(fill = "transparent"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),
	 axis.text.x = element_text(color = "black", size = 12),
	 axis.text.y = element_text(color = "black", size = 12), 
	 axis.title.x = element_text(vjust = 2),
	 axis.title.y = element_text(margin = margin(r = 4)),
	 plot.tag = element_text(size = 18))

# lft N_pups barplots
pp_hist <- ggplot(data=dfpp, aes(x=Counts)) +
  geom_bar(aes(y=PP), stat="identity", position ="identity", fill=ppcol, color='transparent') +
  labs(x='Produced pups',y='Counts', tag = 'D') +
  ylim(c(0,300)) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     #panel.grid.major = element_blank(),
     #panel.grid.minor = element_blank(),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),
	 axis.text.x = element_text(color = "black", size = 10),
	 axis.text.y = element_text(color = "black", size = 12),
	 axis.title.y = element_text(margin = margin(r = 10)),
	 plot.tag = element_text(size = 18))+
	 expand_limits(x = -0.001, y = 0) +
	 coord_cartesian(expand = FALSE, clip = "off")
	 
ps_hist <- ggplot(data=dfps, aes(x=Counts)) +
  geom_bar(aes(y=PS), stat="identity", position="identity", fill=pscol, color='transparent')+
  labs(x='Survived pups',y='Counts', tag = 'E') +
  ylim(c(0,300)) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     #panel.grid.major = element_blank(),
     #panel.grid.minor = element_blank(),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),
	 axis.text.x = element_text(color = "black", size = 10),
	 axis.text.y = element_text(color = "black", size = 12),
	 axis.title.y = element_text(margin = margin(r = 4)),
	 plot.tag = element_text(size = 18))+
	 expand_limits(x = -0.001, y = 0) +
	 coord_cartesian(expand = FALSE, clip = "off")	 

### Map
img <- image_read('C:/Users/dvendrami/Pictures/2832-world view v2.png')
# or: img <- image_read('C:/Users/dvendrami/Pictures/2832-world view v5 - anomaly trend 2014 to 2024.png')
img_gg <- ggdraw() + draw_image(img) + labs(tag = 'A') + theme_bw() +
		theme(plot.tag = element_text(size = 18), panel.border = element_blank())
		
### Composite figure:
# p3 <- (p1 / plot_spacer() / p2) / (pp_hist | ps_hist) / (pr_hist | Froh) + plot_layout(heights = c(4,-1.05,4,4,4))
p3 <- (p1 / plot_spacer() / p2) + plot_layout(heights = c(4,-0.5,4))

layout <- "
AB
AB
CD
"

Fig1 <- img_gg + p3 + pp_hist + ps_hist + plot_layout(design = layout)  

# Export Figure 1
pdf("No_PR/Figure_1_vF.pdf",width=12, height=8)
Fig1
dev.off()

##################################################################################
##################################################################################
##################################################################################
######## Figure 2 - inbreeding (thanks to Martin Stoffel for code):
data <- read.table('PSdataset.txt', h = T)

# Froh
froh <- data |>
	select(SampleID, Froh) |>
	filter(!duplicated(SampleID)) |>
	mutate(Subset = rep(0, length(unique(data$SampleID))))
 
Froh <- ggplot(data=froh,aes(x=Froh)) +
	geom_histogram(fill='grey50', bins=50) +
	labs(x=expression(italic(F)[ROH]),y='Counts', tag = 'A') +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),axis.text = element_text(color = "black", size = 12),
	 plot.tag = element_text(size = 18)) +
	 scale_x_continuous(limits = c(0,0.12), breaks = seq(0,0.12,0.04)) +
	 scale_y_continuous(limits = c(0,75), breaks = seq(0,75,25))

#### Note: for plotting it is better to use plink output.
# Plink and RZooRoH are highly correlated:
# cor.test(bla$FROH.x,bla$FROH.y)

#         Pearson's product-moment correlation
# 
# data:  bla$FROH.x and bla$FROH.y
# t = 61.571, df = 773, p-value < 0.00000000000000022
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8986432 0.9225973
# sample estimates:
#      cor 
# 0.911389 

roh_lengths <- fread('hbdseg/ZooRoH_max.hom')
chr <- roh_lengths$CHR
chr <- gsub('mscaf_a1_0','',chr)
chr <- gsub('mscaf_a1_','',chr)
roh_lengths$CHR <- as.integer(chr)

# Chr lengths
chr_data <- read_delim("hbdseg/Chromosomes.txt", delim = "\t") %>% 
        rename(size_BP = Length,
               CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:18, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)
sd(num_roh_per_ind$n)

froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))

# ROH length and abundance in the 1% least and most inbred individuals 
num_roh_per_ind %>% 
        left_join(froh) %>% 
       # top_frac(-0.01, FROH) %>%    # top 1% least inbred individuals
        top_n(-7, FROH) %>%      # top 1% most inbred individuals
        summarise(mean(n), mean(KBAVG))

n_roh <- ggplot(data=num_roh_per_ind,aes(x=n)) +
	geom_histogram(fill='grey50', bins=43) +
	labs(x='Number of ROHs per individual',y='Counts', tag = 'B') +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = 12),axis.text = element_text(color = "black", size = 12),
	 plot.tag = element_text(size = 18))
	 #scale_x_continuous(limits = c(0,0.12), breaks = seq(0,0.12,0.04)) +
	 #scale_y_continuous(limits = c(0,75), breaks = seq(0,75,25))

#######

all_roh <- roh_lengths %>% 
        group_by(IID) %>% 
        summarise(sum_roh = sum(KB)) %>% 
        ungroup() %>% 
        arrange(desc(sum_roh))

longest_roh <- all_roh %>% 
        top_n(7)
shortest_roh <- all_roh %>% 
        top_n(-7)
num_ind <- 14

extreme_roh <- rbind(longest_roh, shortest_roh)

df <- roh_lengths %>%
        mutate(POS1 = POS1 / 1e+6,
               POS2 = POS2 / 1e+6,
               MB = KB / 1000)

df <- df %>% filter(IID %in% extreme_roh$IID) %>% 
        mutate(IID = factor(IID, levels = extreme_roh$IID))

yax <- data.frame(IID = fct_inorder(levels(df$IID))) %>%
        mutate(yax = seq(from = 2,
                         to = 2*length(unique(df$IID)),
                         by = 2)) 

df <- left_join(df, yax, by = "IID")

shade <- df %>%
        group_by(CHR) %>%
        summarise(min = min(POS1), max = max(POS2)) %>%
        mutate(min = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                       CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                       CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                               TRUE ~ min)) %>%
        mutate(max = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                       CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                       CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                               TRUE ~ max))

col <- c('grey50', 'grey50')
chr_names <- as.character(1:17)
names(chr_names) <- as.character(1:17)

df %>% 
        filter(MB > 0.01) %>% 
        filter(CHR %in% 1:26) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "#E3E6EB") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                      fill = as.factor(CHR)),  col = "grey", size = 0, alpha = 1) + 
        scale_fill_manual(values = rep(col, 18)) + 
        scale_color_manual(values = rep(col, 18)) +
        scale_y_reverse(expand = c(0, 0)) +
		labs(tag = 'D') +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        theme(#strip.placement = 'outside',
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
                axis.line.x = element_blank(),
                legend.position="none",
                axis.title.x = element_text(margin=margin(t=0)),
                axis.title.y = element_text(margin=margin(r=0)),
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank(),
				plot.tag = element_text(size = 18)) +
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind
ROH_per_ind

# needs to be a grob to display axis labels correctly
pg <- ggplotGrob(ROH_per_ind)

for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
}

ROH_per_ind_grob <- as.ggplot(pg)

#############################################
hom_sum <- fread("hbdseg/ZooRoH_max.hom.summary")
chr <- hom_sum$CHR
chr <- gsub('mscaf_a1_0','',chr)
chr <- gsub('mscaf_a1_','',chr)
hom_sum$CHR <- as.integer(chr)

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))

running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"))

# ROH sharing 1
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/775) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

#cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete.obs")
#plot(running_roh$UNAFF_mean, running_roh$UNAFF_n)

# remove windows without snps
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/775) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

# check distribution of snps
#ggplot(running_roh, aes(win_start, UNAFF_n)) + 
#  geom_point() + geom_smooth(span = 0.1) +
#  facet_wrap(~CHR, scales = "free_x")

# ROH across the genome plot

# Your color palette:
generate_palette <- function(color1, color2, n = 20) {
  palette_func <- colorRampPalette(c(color1, color2))
  palette_func(n)
}


# generate_palette_3 <- function(color1, color2, color3, n = 20) {
#   palette_func <- colorRampPalette(c(color1, color2, color3))
#   palette_func(n)
# }
scales::show_col(my_palette)
 
# fill_cols <- generate_palette(pscol, ppcol)
fill_cols <- generate_palette('#D5D9C3', '#74B7DA') # map colors
#fill_cols <- generate_palette3('#fc8d62', '#66c2a5')

# fill_cols <- generate_palette_3('#fc8d62', '#8da0cb', '#66c2a5')

#fill_cols <- cividis(20)
#fill_cols <- plasma(20)
qn <- scales::rescale(quantile(running_roh_p$UNAFF_mean,
                               probs=seq(0, 1, length.out=length(fill_cols))))

p1 <- ggplot(running_roh_p, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
        #geom_tile(color = "grey", size = 0) +
        geom_tile() +
        theme_simple(base_size = 13, grid_lines = FALSE) + 
        scale_y_continuous(expand = c(0,0))+
        scale_x_continuous(expand = c(0,0), 
                           breaks = seq(0, 300000, by = 50000),
                           labels = as.character(seq(0, 300, 50)))+
        labs(tag = 'C') +
		ylab("Chromosome") +
        scale_fill_gradientn("% of seals with ROH",
                             colors = rev(fill_cols), values = qn,
                             breaks = c(0, 0.05,0.15),
                             labels = c(0, 5, 15),
							 guide = ) +
        facet_grid(CHR~., switch="both") +
        xlab("Position in Mb") +
        theme(
              #axis.title.x = element_text(margin=margin(t=5)),
              #axis.title.y = element_text(margin=margin(r=5)),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.ticks.x = element_line(size = 0.3),
              #axis.title.y = element_blank(),
              #plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = c(0.75,0.65),
              legend.direction = "horizontal",legend.title = element_text(size = 10),
              strip.text.y.left = element_text(size = 10, angle = 0),
            #  axis.line.x = element_line(size = 0.3)
              plot.tag = element_text(size = 18)) +
        guides(fill = guide_colourbar(title.position = "bottom" ,
                                      barwidth = 10.95, barheight = 0.5,
									  title.theme = element_text(margin = margin(b = 8,l=30))))
p1

layout1 <- "
AAAABBBB
AAAABBBB
AAAABBBB
CCDDBBBB
CCDDBBBB
"

layout2 <- "
CCDDBBBB
CCDDBBBB
AAAABBBB
AAAABBBB
AAAABBBB
"

Fig2 <- ROH_per_ind + p1 + Froh + n_roh + plot_layout(design = layout2)

layout3 <- "
AABBCCCC
AABBCCCC
DDDDCCCC
DDDDCCCC
DDDDCCCC
"
Fig3 <- Froh + n_roh + p1 + ROH_per_ind + plot_layout(design = layout3)



# pdf("No_PR/Figure_2_v1.pdf",width=13, height=6)
# Fig2
# dev.off()

# To remove strange vertical lines:
ggsave("No_PR/Figure_2_v5.pdf", Fig3, device = cairo_pdf, width=13, height=6)

##################################################################################
##################################################################################
##################################################################################
######## Figure 3, models: 
data <- read.table('PSdataset.txt', h = T)
data_pp <- read.table('PPdataset.txt', h = T)
load('PP_mod.RData')
load('PS_mod.RData')

sdpr <- function (x){
p.est <- mean(x,na.rm=T)
variance <- (p.est*(1-p.est))/length(x[!is.na(x)])
std.dev <- sqrt(variance)
return(std.dev)
}

mns_pp <- data_pp |>
	group_by(Age_new) |>
	summarize(PP_mean = mean(Pup_produced, na.rm = T), PP_sd = sdpr(Pup_produced)) |>
	rename(Age = Age_new) |>
	ungroup()
mns_pp <- mns_pp[-18,] # NAs

mns_ps <- data |>
	group_by(Age) |>
	summarize(PS_mean = mean(Pup_survived, na.rm = T), PS_sd = sdpr(Pup_survived)) |>
	ungroup()

mns_pr <- data2000 |>
	group_by(Age) |>
	summarize(PR_mean = mean(Pup_Rec, na.rm = T), PR_sd = sdpr(Pup_Rec)) |>
	ungroup()

mns <- left_join(mns_ps, mns_pp) |> left_join(mns_pr)

mns <- as.data.frame(mns)
mns$PP_sd_up <- mns$PP_sd
mns$PP_sd[is.na(mns$PP_sd)] <- 0
mns$PS_sd_up <- mns$PS_sd
mns$PS_sd[is.na(mns$PS_sd)] <- 0
mns$PR_sd_up <- mns$PR_sd
mns$PR_sd[is.na(mns$PR_sd)] <- 0


for (i in 1:length(mns$PP_sd_up)){
 if (!is.na(mns$PP_mean[i]) & mns$PP_mean[i] + mns$PP_sd[i] > 1){
  mns$PP_sd_up[i] <- 1 - mns$PP_mean[i]
 }
 if (!is.na(mns$PS_mean[i]) & mns$PS_mean[i] + mns$PS_sd[i] > 1){
  mns$PS_sd_up[i] <- 1 - mns$PS_mean[i]
 }
 if (!is.na(mns$PR_mean[i]) & mns$PR_mean[i] + mns$PR_sd[i] > 1){
  mns$PR_sd_up[i] <- 1 - mns$PR_mean[i]
 }
}

mns$PP_sd_dn <- mns$PP_sd
mns$PP_sd[is.na(mns$PP_sd)] <- 0
mns$PS_sd_dn <- mns$PS_sd
mns$PS_sd[is.na(mns$PS_sd)] <- 0
mns$PR_sd_dn <- mns$PR_sd
mns$PR_sd[is.na(mns$PR_sd)] <- 0
mns$PR_sd[is.na(mns$PR_sd)] <- 0

for (i in 1:length(mns$PP_sd_dn)){
 if (!is.na(mns$PP_mean[i]) & mns$PP_mean[i] - mns$PP_sd[i] < 0){
  mns$PP_sd_dn[i] <- mns$PP_mean[i]
 }
 if (!is.na(mns$PS_mean[i]) & mns$PS_mean[i] - mns$PS_sd[i] < 0){
  mns$PS_sd_dn[i] <- mns$PS_mean[i]
 }
 if (!is.na(mns$PR_mean[i]) & mns$PR_mean[i] - mns$PR_sd[i] < 0){
  mns$PR_sd_dn[i] <- mns$PR_mean[i]
 }
}

mns$PP_sd_dn95 <- 1.96*mns$PP_sd_dn
mns$PP_sd_up95 <- 1.96*mns$PP_sd_up
mns$PP_sd_up95[18] <- 0.32

#### Big with all
### Complete stick plot
# ES_names <- c('SS','Sex','EQ','PAR','Age',expression(paste("Age"^"2")),'GQ','PYP','Den','EQ*GQ','EQ*Age',expression(paste("EQ*Age"^"2")),'EQ*Sex')
ESnam <- c('Silver spoon','Sex','SAM','PAR','Age',"Age2",'Inbreeding','PPF','Density','SAM*Inbreeding','SAM*Age',"SAM*Age2",'Inbreeding*PPF','Inbreeding*Age','Inbreeding*Age2')
ES_names <- c('Inbreeding*SAM','Inbreeding*PPF',expression(paste("Age"^"2"*"*Inbreeding")),'Age*Inbreeding',expression(paste("Age"^"2"*"*SAM")),'Age*SAM','Sex','Inbreeding',expression(paste("Age"^"2")),'Age','Density','SAM','PPF','PAR','SS')

#load('PPNA_allI_fems_jcat_fin_both_th.RData')
sums <- summary(PP_mod)
ES_pp <- c(sums$fixed$Estimate[2],NA,sums$fixed$Estimate[3:15])
pp_low <- c(sums$fixed[2,3],NA,sums$fixed[3:15,3])
pp_high <- c(sums$fixed[2,4],NA,sums$fixed[3:15,4])
temp <- data.frame(ESnam = ESnam, ES = ES_pp, low = pp_low, high = pp_high)
tbu <- temp
temp[1,] <- tbu[10,]
temp[2,] <- tbu[13,]
temp[3,] <- tbu[15,]
temp[4,] <- tbu[14,]
temp[5,] <- tbu[12,]
temp[6,] <- tbu[11,]
temp[7,] <- tbu[2,]
temp[8,] <- tbu[7,]
temp[9,] <- tbu[6,]
temp[10,] <- tbu[5,]
temp[11,] <- tbu[9,]
temp[12,] <- tbu[3,]
temp[13,] <- tbu[8,]
temp[14,] <- tbu[4,]
temp[15,] <- tbu[1,]
ES_pp <- temp$ES
pp_low <-temp$low
pp_high <- temp$high

#load('PSallI_fems_EQS_jcat_sex_fin_both_th.RData')
sums <- summary(PS_mod)
ES_ps <- sums$fixed$Estimate[2:16]
ps_low <- sums$fixed[2:16,3]
ps_high <- sums$fixed[2:16,4]
temp <- data.frame(ESnam = ESnam, ES = ES_ps, low = ps_low, high = ps_high)
tbu <- temp
temp[1,] <- tbu[10,]
temp[2,] <- tbu[13,]
temp[3,] <- tbu[15,]
temp[4,] <- tbu[14,]
temp[5,] <- tbu[12,]
temp[6,] <- tbu[11,]
temp[7,] <- tbu[2,]
temp[8,] <- tbu[7,]
temp[9,] <- tbu[6,]
temp[10,] <- tbu[5,]
temp[11,] <- tbu[9,]
temp[12,] <- tbu[3,]
temp[13,] <- tbu[8,]
temp[14,] <- tbu[4,]
temp[15,] <- tbu[1,]
ES_ps <- temp$ES
ps_low <-temp$low
ps_high <- temp$high

ESs <- c(ES_ps, ES_pp)
highs <- c(ps_high-ES_ps, pp_high-ES_pp)
lows <- c(ES_ps-ps_low, ES_pp-pp_low)
x <- c(seq(0.8,21.8,1.5), seq(1.2,22.2,1.5))
cols<-c(rep('survived',15),rep('produced',15))
df1 <- data.frame(x=x,y=ESs,Low=lows,High=highs,mod=factor(cols,levels=c('produced','survived')))

ES_2 <- ggplot(df1, aes(x=x, y=y)) +  
   geom_point(aes(col = mod), size = 5) +
   geom_errorbar(aes(ymin=y-lows, ymax=y+highs, col=mod), width=0) +
   labs(x='',y='Scaled effect size', tag = 'A') +
   geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
   ylim(c(-1,1)) + 
   scale_x_continuous(limits=c(0.6,22.4), breaks=c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,17.5,19,20.5,22),labels=ES_names) +
   scale_color_manual(name='', values=c('weight'=pwcol,'produced'=ppcol,'survived'=pscol)) +
   guides(color = guide_legend(override.aes = list(linetype = c(0, 0),
                                                  shape = c(16, 16)))) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
	 legend.position='none',
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize)) +
	 coord_flip()

pp_age <- ggplot(mns, aes(x=Age, y=PP_mean)) +  
   geom_errorbar(aes(ymin=PP_mean-PP_sd_dn95, ymax=PP_mean+PP_sd_up95), width=0, col = ppcol, alpha = .5) +
   stat_smooth(color = 'grey45', fill = 'grey80', method = 'lm', formula = y ~ I(x^2) + x) +
   geom_point(col = ppcol, size = 3) +
   #geom_point(size = 5) +
   #scale_y_continuous(breaks=seq(0,1,0.25)) +
   ylim(0,1) +
   labs(x='Age',y='Pup produced', tag = 'F') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

ps_age <- ggplot(mns, aes(x=Age, y=PS_mean)) +  
   geom_errorbar(aes(ymin=PS_mean-1.96*PS_sd_dn, ymax=PS_mean+1.96*PS_sd_up), width=0, , col = pscol, alpha = .5) +
   stat_smooth(color = 'grey45', fill = 'grey80', method = 'lm', formula = y ~ I(x^2) + x) +
   geom_point(col = pscol, size = 3) +
   ylim(0,1) +
   #geom_point(size = 5) +
   labs(x='Age',y='Pup survived', tag = 'G') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

ps_sam <- ggplot(data, aes(x=Repr_year_SAM_ANN, y=Pup_survived)) +  
   geom_smooth(color = 'grey45', fill = 'grey80', method = 'glm', method.args = list(family="binomial")) +
   geom_point(col = pscol, size = 3, alpha = .02) +
   #geom_point(size = 5) +
   labs(x='SAM',y='Pup survived', tag = 'C') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

# Jcat 
jcatd <- data_pp |>
	group_by(JcatAll) |>
	summarize(PP = mean(Pup_produced), PPSD = sdpr(Pup_produced))
jcatd <- jcatd[-3,]
jcatd$Jcat <- as.character(jcatd$JcatAll)
jcatd$Jcat <- c('0','1')

pp_jcat <- ggplot(jcatd, aes(x=Jcat, y=PP)) +  
   geom_errorbar(aes(ymin=PP-1.96*PPSD, ymax=PP+1.96*PPSD), width=0, , col = ppcol, alpha = .5) +
   #stat_smooth(color = 'grey45', fill = 'grey80', method = 'lm', formula = y ~ I(x^2) + x) +
   geom_point(col = ppcol, size = 3) +
   ylim(0.65,0.85) + 
   #geom_point(size = 5) +
   labs(x='PPF',y='Pup produced', tag = 'B') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

# Density 
pp_dens <- ggplot(data_pp, aes(x=Females_count, y=Pup_produced)) +  
   geom_smooth(aes(y=Pup_produced),color = 'grey45', fill = 'grey80', method = 'glm', method.args = list(family="binomial")) +
   geom_point(aes(y=Pup_produced),col = ppcol, size = 3, alpha = .02) +
   #geom_point(size = 5) +
   labs(x='Density',y='Pup produced', tag = 'D') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

psdata <- data[data$Females_count > 450,]
ps_dens <- ggplot(psdata, aes(x=Females_count)) +  
   geom_smooth(aes(y=Pup_produced),color = 'grey45', fill = 'grey80', method = 'glm', method.args = list(family="binomial")) +
   geom_point(aes(y=Pup_survived),col = pscol, size = 3, alpha = .02) +
   #geom_point(size = 5) +
   labs(x='Density',y='Pup survived', tag = 'E') +
   #coord_cartesian(xlim=c(3,NA), ylim=c(0,NA), expand = TRUE, clip = "off") +
   theme_bw() +
	theme(axis.line = element_line(colour = "black"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

blappp <- blapp
scpp <- scale(data_pp$Repr_year_SAM_ANN)
mupp <- attr(scpp, "scaled:center")
sdpp <- attr(scpp, "scaled:scale")
blappp$x <- blappp$x * sdpp + mupp
inter_pp2 <- ggplot(blappp, aes(x = x, y = predicted, colour = group)) +
  geom_line(size=1) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), alpha = .15, show.legend = FALSE) +
  #facet_wrap(~facet) +
  labs(x='SAM',y='Pup production', tag = 'H') +
  theme_bw() +
  ylim(0.5,1) +
  scale_color_manual(name = '',values=intcol2.2, labels = c('Outbred','Inbred'))+
  scale_fill_manual(values=intcol2, name="",labels=NULL) +
	theme(axis.line = element_line(colour = "black"),
	legend.position=c(0.4,0.25),
	legend.text = element_text(size = 15),
	#legend.title = element_text(size = 12),
	legend.key = element_rect(fill = "transparent", linewidth=0, color='transparent'),
	legend.background = element_rect(fill = "transparent"),
	#legend.key = element_rect(fill = alpha('white',0.99), linewidth=0, color=alpha('white',0.99)),
	#legend.spacing.y = unit(0.01, 'cm'),
	#legend.spacing.x = unit(0.01, 'cm'),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize)) #+
    #guides(fill = guide_legend(byrow = TRUE))

blapss <- blaps
scps <- scale(data$Repr_year_SAM_ANN)
mups <- attr(scps, "scaled:center")
sdps <- attr(scps, "scaled:scale")
blapss$x <- blapss$x * sdps + mups
inter_ps2 <- ggplot(blapss, aes(x = x, y = predicted, colour = group)) +
  geom_line(size=1) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), alpha = .15, show.legend = FALSE) +
  #facet_wrap(~facet) +
  labs(x='SAM',y='Pup survival', tag = 'I') +
  theme_bw() +
  ylim(0.5,1) +
  scale_color_manual(name = '',values=intcol3.2, labels = c('Outbred','Inbred'))+
  scale_fill_manual(values=intcol3, name="",labels=NULL) +
	theme(axis.line = element_line(colour = "black"),
	legend.position=c(0.4,0.25),
	legend.text = element_text(size = 15),
	#legend.title = element_text(size = 12),
	legend.key = element_rect(fill = "transparent", linewidth=0, color='transparent'),
	legend.background = element_rect(fill = "transparent"),
	#legend.key = element_rect(fill = alpha('white',0.99), linewidth=0, color=alpha('white',0.99)),
	#legend.spacing.y = unit(0.01, 'cm'),
	#legend.spacing.x = unit(0.01, 'cm'),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize)) #+
    #guides(fill = guide_legend(byrow = TRUE))

# No legend
ES_2 <- ggplot(df1, aes(x=x, y=y)) +  
   geom_point(aes(col = mod), size = 5) +
   geom_errorbar(aes(ymin=y-lows, ymax=y+highs, col=mod), width=0) +
   labs(x='',y='Scaled effect size', tag = 'A') +
   geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
   ylim(c(-1,1)) + 
   scale_x_continuous(limits=c(0.6,22.4), breaks=c(1,2.5,4,5.5,7,8.5,10,11.5,13,14.5,16,17.5,19,20.5,22),labels=ES_names) +
   scale_color_manual(name='', values=c('weight'=pwcol,'produced'=ppcol,'survived'=pscol)) +
   guides(color = guide_legend(override.aes = list(linetype = c(0, 0),
                                                  shape = c(16, 16)))) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
	 legend.position='none',
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize)) +
	 coord_flip()

layout <- "
AAABBCC
AAADDEE
AAAFFGG
AAAHHII
"

pdf("Figure_3.pdf",width=16, height=15)
ES_2 + pp_jcat + ps_sam + pp_dens + ps_dens + pp_age + ps_age + inter_pp2 + inter_ps2 + plot_layout(design = layout)  
dev.off()

#####################################################################
#####################################################################
################ Figure 4
data <- read.table('PSdataset.txt', h = T)
results <- read.table('results_fig4.txt', h = T)
df_fig4 <- read.table("df_fig4.txt", h = T)
# error bars plotting function									  
mean_ci <- function(x) {
  m  <- mean(x)
  data.frame(
    y    = m,
    ymin = quantile(x,prob=0.025),
    ymax = quantile(x,prob=0.975)
  )
}

scps <- scale(data$Repr_year_SAM_ANN)
mups <- attr(scps, "scaled:center")
sdps <- attr(scps, "scaled:scale")
results$x <- results$Repr_year_SAM_ANN_scaled * sdps + mups


w_fit <- ggplot(results, aes(x = x, y = mean,
                    color = factor(Froh_scaled))) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(Froh_scaled)),
              alpha = 0.2, color = NA,show.legend = FALSE) +
  labs(x = "SAM", y = "Pup weaning posterior probabilities", color = "Inbreeding", fill = "Inbreeding", tag = "A") +
  scale_color_manual(name = '',values=intcol, labels = c('Outbred','Inbred'))+
  scale_fill_manual(values=intcol.2, name="",labels=NULL) +
  theme_bw() +
	theme(axis.line = element_line(colour = "black"),
	legend.position=c(0.4,0.25),
	legend.text = element_text(size = 15),
	legend.key = element_rect(fill = "transparent", linewidth=0, color='transparent'),
	legend.background = element_rect(fill = "transparent"),
     panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize))

betas <- ggplot(df_fig4, aes(x = type, y = value)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_violin(aes(fill = type), trim = FALSE, alpha = 0.4, color = NA) +
  stat_summary(aes(color = type), fun = mean, geom = "point", size = 5) +
  stat_summary(aes(color = type), fun.data = mean_ci,
               geom = "linerange", linewidth = 1) +
  facet_wrap(~env) +
  scale_fill_manual(values = c(
    "PP" = "#8da0cb",
    "PS" = "#fc8d62",
    "PW"  = "#66c2a5")) +
  scale_color_manual(values = c(
    "PP" = "#7570b3",
    "PS" = "#d95f02",
    "PW"  = "#1b9e77")) +
  labs(x = "", y = expression(beta), tag = 'B') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
  strip.background = element_blank(),
  strip.text = element_text(color = "black", size = f2_size),
	 panel.grid.major.y = element_line(linetype = 2),
	 panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
	 axis.title = element_text(vjust = 0, size = f2_size),axis.text.x = element_text(color = "black", size = f2_size),
		,axis.text.y = element_text(color = "black", size = f2_size),
	 plot.tag = element_text(size = tsize),
    legend.position = "none") 


layout <- "
ABB
"

Fig4 <- w_fit + betas + plot_layout(design = layout)  

# Export Figure 1
pdf("Figure_4.pdf",width=15, height=6)
Fig4
dev.off()



