########################################################################################################################################################################################################################################################################################################################################################################
# GENOTYPE OUTPLANTING DATA ANALYSES
# Purpose: Analyzing monitoring data of Acropora cervicornis performance in the nursery and after outplanting, in collaboration with Reef Renewal Foundation Bonaire
# Last updated: 2026 October 26, by Allie Blanchette
########################################################################################################################################################################################################################################################################################################################################################################

##### Publication citation:
# Title: Genotypic variation in outplanted staghorn coral performance and implications for supporting adaptation in restoration
# Authors: Allison Blanchette, Francesca Virdis, Sanne H. Tuijten, Sarah E. Lester, Andrew Rassweiler




#####################################################################################################################################################################################################################################################################
## Load in packages and data
#####################################################################################################################################################################################################################################################################

##### Clear environment
rm(list=ls())  ##specifically clears your global environment


##### Load in relevant packages
library(tidyverse)   # for editing dataframes
library(dplyr)       # for editing dataframes
library(lme4)        # for running GLMEM
library(ggfortify)   # for plotting PCA
library(ggrepel)     # for plotting PCA
library(DHARMa)      # for diagnostics and Moran's I
library(ape)         # for another way to analyze Moran's I
# library(car)       # for Anova (commented out because it masks recode function)
library(viridis)     # for viridis color-blind palette
library(ggplot2)
library(rptR)        # for repeatability results with confidence intervals, citation("rptR")
library(cowplot)     # for ggdraw function
library(lubridate)   # for working with date data types




##### Set working directory
setwd("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\PV data\\Downstream_created")


##### Read in cleaned monitoring data and nursery data from PV_dataframes_250728.R
Mon_wide <- read.csv("Mon_wide.csv") %>% mutate(Geno = as.character(Geno))                        # All monitoring data in wide format = each stressor is a column of 1s and 0s, no redundancy in partial tissue mortality. Row = cluster
Mon_long <- read.csv("Mon_long.csv") %>% mutate(Geno = as.character(Geno))                        # All monitoring data in long format = stressor type is 1 column (variable) and stressor data is 1 column (Stress_pres = 1/0 = presence/absence), there is redundnacy in partial tissue mortality column. Row = stressor variable/data PER cluster
Areas_clusters_all <- read.csv("Areas_clusters_all.csv") %>% mutate(Geno = as.character(Geno))    # All cluster area data per orthomosaic (T1-T3) and calculated growth rate (cm2/day)
Areas_clusters_M03 <- read.csv("Areas_clusters_M03.csv") %>% mutate(Geno = as.character(Geno))    # Month 3 stressors and net change in area, per cluster. NOTE - only use when primary concern is the Net Change in Area. 4 clusters from M03 excluded because they did not have a T2 area measurement: either lost to  Hurricane Ian or T2 measurement manually removed b/c no good (13C03 was covered by a gorg in T2)
Nursery_batches <- read.csv("Nursery_batches.csv") %>% mutate(Geno = as.character(Geno))          # Nursery data, with columns split by Batch



##### Additional wrangling for descriptive stats and plotting

# Make function to calculate confidence intervals around the mean; defaults to 95% CI
ci_mean <- function(y, name, conf.level = 0.95) {
  y <- y[!is.na(y)]
  n <- length(y)
  m <- qt(1 - (1 - conf.level) / 2, df = n - 1)
  se <- sd(y) / sqrt(n)
  mean.y <- mean(y)
  
  ci.lwr <- mean.y - m * se
  ci.upr <- mean.y + m * se
  
  out <- tibble(mean.y, n, se, ci.lwr, ci.upr, conf.level)
  names(out) <- paste(name, c("mean", "n", "se", "ci.lwr", "ci.upr", "conf.level"), sep = "_")
  out
}

# test ci function
y1 <- Mon_wide$DeadTiss_cluster[Mon_wide$Time.Point == "M03" & Mon_wide$Batch == "A" & Mon_wide$Geno == "10"]
ci_mean(y1, name = "y1")         


# Stressor prevalences per genotype
stress_prevs <- Mon_wide %>% 
  select(-notes) %>% 
  group_by(Batch, Time.Point, Geno, Geno_chr) %>% 
  summarise(dis = mean(Disease, na.rm = T),
            cil = mean(Ciliate, na.rm = T),
            cya = mean(Cyanobacteria, na.rm = T),
            bot = mean(Bottom_contact, na.rm = T)) %>% 
  pivot_wider(names_from = Time.Point,
              values_from = c(dis:bot))


# Partial tissue mortality per genotype
pm_avg <- Mon_wide %>% 
  filter(!(Time.Point == "M09" | Time.Point == "M20")) %>%        # remove time points when partial mortality was not collected, to avoid getting warning messages about NaNs
  group_by(Batch, Time.Point, Geno, Geno_chr, Harvest_site) %>% 
  summarise(ci_mean(DeadTiss_cluster, "DeadTiss_cluster"))


# Net change in area per genotype 
Areas_long <- Areas_clusters_all %>% 
  filter(Batch == "A") %>% 
  select(Batch, Geno, Geno_chr, Cluster, GenoCluster, nc_T1T2, nc_T2T3, nc_T1T3) %>% 
  pivot_longer(cols = starts_with("nc"),
               names_to = "duration",
               values_to = "net_change") 

growth_avg_long <- Areas_long %>% 
  group_by(Batch, Geno, Geno_chr, duration) %>% 
  summarise(ci_mean(net_change, "net_change")) %>% 
  mutate(Year = duration,
         Year = recode(Year, "nc_T1T2"="1"),
         Year = recode(Year, "nc_T2T3"="2"),
         Year = recode(Year, "nc_T1T3"="Both Years"))


# Batch correlations 
# - ignore NaNs warnings - they're about M09 and M20 not having any partial mortality data (this is correct, but we want to keep these time points for their stressor data)
batch_corr <- Mon_wide %>% 
  group_by(Batch, Time.Point, Geno, Geno_chr) %>% 
  mutate(pm = DeadTiss_cluster) %>%   # make an abbreviated name for partial mortality before running ci_mean in summarise
  summarise(surv = mean(survived),
            dp = mean(Disease),
            cp = mean(Ciliate),
            fp = mean(Fireworm),
            op = mean(Overgrowth),
            ci_mean(pm, "pm")) %>% 
  pivot_wider(names_from = Batch,
              values_from = c(surv:pm_conf.level))



# Month 3 stressors and net change in area, per genotype
Areas_genos_M03 <- Areas_clusters_M03 %>% 
  filter(Batch != "B") %>%         # removing batch B because it did not have any net change data
  group_by(Batch, Time.Point, Geno, Geno_chr) %>% 
  summarise(
    
    # Summarise stressor prevalences
    dis_prop = mean(Disease),
    cili_prop = mean(Ciliate),
    bottom_prop = mean(Bottom_contact),
    fireworm_prop = mean(Fireworm),
    cyano_prop = mean(Cyanobacteria),
    overgrowth_prop = mean(Overgrowth),
    
    # Summarise net change in coral cover
    ci_mean(nc_T1T2, "nc_T1T2"),    # net change (cm2/day) T1-T2 -- year 1
    ci_mean(nc_T1T3, "nc_T1T3"),    # net change (cm2/day) T1-T3 -- both years
    ci_mean(nc_T2T3, "nc_T2T3")     # net change (cm2/day) T2-T3 -- year 2
    )




###### Join nursery data with different forms of outplant data

# Outplant data: stressors and partial tissue mortality, BOTH BATCHES
Nurs_OP_Monwide_AB <- Mon_wide %>% 
  group_by(Time.Point, Geno, Geno_chr) %>% 
  summarise(dis_OP = mean(Disease, na.rm = T),
            partmort_OP = mean(DeadTiss_cluster, na.rm = T)) %>% 
  left_join(Nursery_batches)

# Outplant data: stressors and partial tissue mortality. BATCH A ONLY
Nurs_OP_Monwide_A <- Mon_wide %>% 
  group_by(Batch, Time.Point, Geno, Geno_chr) %>% 
  summarise(dis_OP = mean(Disease, na.rm = T),
            partmort_OP = mean(DeadTiss_cluster, na.rm = T)) %>% 
  filter(Batch == "A") %>% 
  left_join(Nursery_batches)

# Outplant data: net change in area
Nurs_OP_Areas <- Areas_clusters_all %>% 
  group_by(Geno, Geno_chr) %>% 
  summarise(Avg_nc12 = mean(nc_T1T2, na.rm = T),
            Avg_nc23 = mean(nc_T2T3, na.rm = T),
            Avg_nc13 = mean(nc_T1T3, na.rm = T)) %>% 
  left_join(Nursery_batches)






##### Some useful notes on data subsetting:
# To do multiple subsets within a data frame, use the below syntax
# Don't forget the comma!

# Ex.1: subsets out only rows from Batch A and Month 3
# Mon_wide[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M03",]

# Ex.2: subsets out only rows from either Month 3 OR Month 6
# Mon_wide[Mon_wide$Time.Point == "M03" | Mon_wide$Time.Point == "M06",]




#####################################################################################################################################################################################################################################################################
## Set up plotting aesthetics
#####################################################################################################################################################################################################################################################################

# Typical publication dimensions
# Single-column: 8.5cm wide (x2 = 17)
# Double-column: 18cm wide (x2 = 36)


# Genotype colors: 10-19
scales::show_col(viridis_pal(option = "plasma")(10))
Geno.colors <- c("#F0F921FF", "#FDC926FF", "#FA9E3BFF", "#ED9753FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF", "#0D0887FF")

# Batch colors: A, B
Batch.colors.2 <- c("#005693","#03A4AF")

# Batch colors: A, B, Unknown
Batch.colors.3 <- c("#005693","#03A4AF", "gray40")

# Year colors: 1, 2
Year.colors <- c("#483D8B", "#836fff")



# Plot text sizes
axis.title.sz = 15
axis.txt.sz = 13
axis.txt.geno.sz = 12
leg.txt.sz = 13
annotate.sz = 5

# Plot feature sizes
point.sz = 3.5
error.ln.sz = 0.3
error.wd.sz = 0
regln.sz = 0.75
horizln.sz = 0.6






#####################################################################################################################################################################################################################################################################
## Disease and Ciliate Correlations - Supplement Table S2
#####################################################################################################################################################################################################################################################################


# M03 Batch A - MARGINALLY significant
# Pearson r = 0.634, p = 0.04908
hist(stress_prevs$dis_M03[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$dis_M03[stress_prevs$Batch == "A"])
hist(stress_prevs$cil_M03[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$cil_M03[stress_prevs$Batch == "A"])
cor.test(stress_prevs$dis_M03[stress_prevs$Batch == "A"], stress_prevs$cil_M03[stress_prevs$Batch == "A"], method="pearson")


# M03 Batch B - NOT significant
# Pearson r = 0.556, p = 0.095
# Kendall t = 0.357, p = 0.168
hist(stress_prevs$dis_M03[stress_prevs$Batch == "B"]) # marginally non-normal
shapiro.test(stress_prevs$dis_M03[stress_prevs$Batch == "B"])
hist(stress_prevs$cil_M03[stress_prevs$Batch == "B"]) # normal
shapiro.test(stress_prevs$cil_M03[stress_prevs$Batch == "B"])
cor.test(stress_prevs$dis_M03[stress_prevs$Batch == "B"], stress_prevs$cil_M03[stress_prevs$Batch == "B"], method="pearson")
cor.test(stress_prevs$dis_M03[stress_prevs$Batch == "B"], stress_prevs$cil_M03[stress_prevs$Batch == "B"], method="kendall")


# M06 Batch A - SIGNIFICANT
# Pearson r = 0.753, p = 0.012
hist(stress_prevs$dis_M06[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$dis_M06[stress_prevs$Batch == "A"])
hist(stress_prevs$cil_M06[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$cil_M06[stress_prevs$Batch == "A"])
cor.test(stress_prevs$dis_M06[stress_prevs$Batch == "A"], stress_prevs$cil_M06[stress_prevs$Batch == "A"], method="pearson")


# M06 Batch B - NOT significant
# Pearson r = -0.012, p = 0.974
hist(stress_prevs$dis_M06[stress_prevs$Batch == "B"]) # normal
shapiro.test(stress_prevs$dis_M06[stress_prevs$Batch == "B"])
hist(stress_prevs$cil_M06[stress_prevs$Batch == "B"]) # normal
shapiro.test(stress_prevs$cil_M06[stress_prevs$Batch == "B"])
cor.test(stress_prevs$dis_M06[stress_prevs$Batch == "B"], stress_prevs$cil_M06[stress_prevs$Batch == "B"], method="pearson")

# M09 Batch A - NOT significant
# Kendall t = -0.211, p = 0.478
hist(stress_prevs$dis_M09[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$dis_M09[stress_prevs$Batch == "A"])
hist(stress_prevs$cil_M09[stress_prevs$Batch == "A"]) # extremely non-normal
shapiro.test(stress_prevs$cil_M09[stress_prevs$Batch == "A"])
cor.test(stress_prevs$dis_M09[stress_prevs$Batch == "A"], stress_prevs$cil_M09[stress_prevs$Batch == "A"], method="kendall")


# M20 Batch A - NOT significant
# Pearson r = -0.141, p = 0.698 
# Kendall t = -0.225, p = 0.3995
hist(stress_prevs$dis_M20[stress_prevs$Batch == "A"]) # normal
shapiro.test(stress_prevs$dis_M20[stress_prevs$Batch == "A"])
hist(stress_prevs$cil_M20[stress_prevs$Batch == "A"]) # non-normal
shapiro.test(stress_prevs$cil_M20[stress_prevs$Batch == "A"])
cor.test(stress_prevs$dis_M20[stress_prevs$Batch == "A"], stress_prevs$cil_M20[stress_prevs$Batch == "A"], method="pearson")
cor.test(stress_prevs$dis_M20[stress_prevs$Batch == "A"], stress_prevs$cil_M20[stress_prevs$Batch == "A"], method="kendall")




#####################################################################################################################################################################################################################################################################
## Summary Figures and Tables - Plotting Figure 2, Supplement Table S3
#####################################################################################################################################################################################################################################################################

##### FIGURE 2 - stressor prevalences across time points, averaged to batch

# load in calendar codes
Batch_Month_Calendar_codes <- read.csv("Batch_Month_Calendar_codes.csv") %>% 
  mutate(Date_ymd = as.Date(Date_ymd)) 
  

# CORAL COVER

# load in unknown and missing coral cover data 
unknown_cover <- read.csv("unknown_cover.csv")

# summarize coral cover per orthomosaic time point
summary_cover <- Areas_clusters_all %>% 
  select(-TagLab.Centroid.x, -TagLab.Centroid.y) %>% 
  pivot_longer(cols = starts_with("T"),
               names_to = "Time.Point",    
               values_to = "cover_cm2") %>% 
  group_by(Batch, Time.Point) %>% 
  summarise(total_cover_cm2 = sum(cover_cm2, na.rm = T)) %>% 
  mutate(total_cover_m2 = total_cover_cm2/10000) %>% 
  full_join(unknown_cover) %>% 
  left_join(Batch_Month_Calendar_codes)


# STRESSORS 

# summarize stressor prevalence data per monitoring time point
summary_stress <- Mon_long %>% 
  mutate(variable=recode(variable, Competition="Overgrowth"),
         variable=recode(variable, Cyanobacteria="Overgrowth")) %>% 
  group_by(Batch, Time.Point, variable) %>% 
  summarise(prevalence = mean(Stress_pres, na.rm = T)) %>%  # Need to remove NAs that arise from clusters that are still listed in the dataset but had actyually died (and so their stressor data is NA)
  filter(variable!= "Bottom_contact",
         variable!= "Unknown",
         variable!= "Snail") %>%
  left_join(Batch_Month_Calendar_codes)

summary_stress_table <- summary_stress %>% 
  ungroup() %>% 
  select(Batch, Calendar, variable, prevalence) %>% 
  spread(variable, prevalence)

write.csv(summary_stress_table, "summary_stress_table.csv", row.names = F)



  
  


# FIGURE 2 

# 2A: coral cover  
fig2A <- 
  ggplot(data= summary_cover, aes(x = Date_ymd, y = total_cover_m2, fill = Batch))+
    geom_bar(color="black", linewidth = 0.1, width = 100, stat="identity", aes(color = Batch), position = position_dodge())+
    scale_fill_manual(values = Batch.colors.3)+
    geom_vline(xintercept = as.Date("2021-09-16"), linetype = "dashed", color = "#005693",  linewidth = 1)+  # Batch A was actually outplanted October 29, but I shifted the dotted line so that it would not overlap bar graph
    geom_vline(xintercept = as.Date("2023-06-12"), linetype = "dashed", color = "#03A4AF",  linewidth = 1)+  # Batch B outplanting
    geom_vline(xintercept = as.Date("2022-09-23"), linetype = "dashed", color = "orange",  linewidth = 1)+   # Hurricane Ian
    geom_vline(xintercept = as.Date("2024-02-08"), linetype = "dashed", color = "orange",  linewidth = 1)+   # wind reversal
    scale_y_continuous(limits = c(0,16),
                       breaks = seq(0,16, by = 3))+
    scale_x_date(limits = c(as.Date("2021-09-15"),
                            as.Date("2024-09-20")),
                 breaks = seq.Date(from = as.Date("2021-11-01"),
                                   to = as.Date("2024-09-20"),
                                   by = "6 months"),
                 date_labels = "%m/%y")+
    theme_classic()+
    xlab("")+
    ylab(expression(paste("Coral cover (",m^{2},")")))+
    theme_classic()+
    theme(axis.title.y = element_text(color="black", size = axis.title.sz),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(color="black",size = axis.txt.sz),
          legend.title=element_blank(),
          legend.direction = "horizontal",
          legend.text=element_text(color="black",size = leg.txt.sz),
          legend.key.size = unit(0.3, 'cm'),
          legend.position = c(0.14, 0.93),
          legend.key.spacing.x = unit(1.5, "cm"))

# 2B: stressors
fig2B <- 
  ggplot(data= summary_stress, aes(x = Date_ymd, y = prevalence, fill = reorder(variable, -prevalence)))+
  geom_bar(color="black", linewidth = 0.1, width = 66, stat="identity", aes(color = reorder(variable, -prevalence)), position = position_dodge())+
  scale_fill_manual(values=c("gray10", "darkcyan", "coral2", "violetred4"), breaks = c('Disease','Ciliate', 'Fireworm', 'Overgrowth'))+
  xlab("")+
  ylab("Stressor prevalence")+
  ylim(0,0.7)+
  scale_x_date(limits = c(as.Date("2021-09-15"),
                          as.Date("2024-09-20")),
               breaks = seq.Date(from = as.Date("2021-11-01"),
                                 to = as.Date("2024-09-20"),
                                 by = "6 months"),
               date_labels = "%b %Y")+
  theme_classic()+
  theme(axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black",size = axis.title.sz, angle=30, vjust=0.5),
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.title=element_blank(),
        legend.direction = "horizontal",
        legend.text=element_text(color="black",size = leg.txt.sz),
        legend.key.size = unit(0.3, 'cm'),
        legend.position = c(0.21, 0.95),
        legend.key.spacing.x = unit(0.8, "cm"))


# Combine  
fig2 <-
  plot_grid(
  fig2A,
  fig2B,
  ncol = 1,
  align = "v",
  labels = c("A", "B"),  # optional labels
  label_size = 15,
  label_y = c(1.0, 1.1),
  rel_heights = c(0.9, 1.1)
)

# Save!
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch2_Punt Vierkant outplanting experiment\\Figures\\fig2\\fig2_graph.tiff",
       fig2,
       bg = "white",
       dpi=600,
       units="cm",
       width = 34,
       height = 14)







##### SUPPLEMENT TABLE S3 - large table of averages per genotype and across site of survivorship, partial tissue mortality, areas, and stressor prevalences

# summary of stressors
summ_stress_geno <- Mon_long %>% 
  group_by(variable, Batch, Time.Point, Geno) %>% 
  summarise(Geno_avg = mean(Stress_pres, na.rm = T)) %>% 
  pivot_wider(names_from = Geno,
              values_from = Geno_avg) %>% 
  filter(variable != "Unknown")

summ_stress_df <- read.csv("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\PV data\\summ_stress_df.csv") %>% 
  pivot_wider(names_from = Geno,
              values_from = Geno_avg)

summ_stress_site <- Mon_long %>% 
  group_by(variable, Batch, Time.Point) %>% 
  summarise(site_avg = mean(Stress_pres, na.rm = T)) %>% 
  filter(variable != "Unknown") %>% 
  left_join(summ_stress_geno) %>% 
  full_join(summ_stress_df)

# summary of partial mortality
summ_pm_geno <- Mon_wide %>% 
  group_by(Batch, Time.Point, Geno) %>% 
  summarise(Geno_avg = mean(DeadTiss_cluster)) %>% 
  pivot_wider(names_from = Geno,
              values_from = Geno_avg) %>% 
  mutate(variable = "Partial tissue mortality") %>% 
  filter(Time.Point %in% c("M03", "M06"))

summ_pm_site <- Mon_wide %>% 
  group_by(Batch, Time.Point) %>% 
  summarise(site_avg = mean(DeadTiss_cluster)) %>% 
  mutate(variable = "Partial tissue mortality") %>% 
  filter(Time.Point %in% c("M03", "M06")) %>% 
  left_join(summ_pm_geno)

# summary of cluster areas
summ_areas <- Areas_clusters_all %>% 
  select(Batch, Geno, T1_cm2, T2_cm2, T3_cm2) %>% 
  pivot_longer(cols = starts_with("T"),
               names_to = "Time.Point",
               values_to = "area_cm2") %>% 
  filter( !(Batch == "B" & Time.Point %in% c("T1_cm2","T2_cm2")))

summ_areas_geno <- summ_areas %>% 
  group_by(Batch, Time.Point, Geno) %>% 
  summarise(Geno_avg = mean(area_cm2, na.rm = T)) %>% 
  pivot_wider(names_from = Geno,
              values_from = Geno_avg) %>% 
  mutate(variable = "Cluster area (cm2)")

summ_areas_site <- summ_areas %>% 
  group_by(Batch, Time.Point) %>% 
  summarise(site_avg = mean(area_cm2, na.rm = T)) %>% 
  left_join(summ_areas_geno)

# summary of survivorship - STILL NEED TO DO
summ_surv_geno <- Mon_wide %>% 
  group_by(Batch, Time.Point, Geno) %>% 
  summarise(Geno_avg = mean(survived)) %>% 
  pivot_wider(names_from = Geno,
              values_from = Geno_avg)

summ_surv_site <- Mon_wide %>% 
  group_by(Batch, Time.Point) %>% 
  summarise(site_avg = mean(survived)) %>% 
  left_join(summ_surv_geno) %>% 
  mutate(variable = "Cluster survivorship")
  


# combine all summaries for mega summary table
summary_large <- summ_surv_site %>% 
  full_join(summ_pm_site) %>% 
  full_join(summ_areas_site) %>% 
  full_join(summ_stress_site) %>% 
  mutate(variable = recode(variable, "Bottom_contact" = "Bottom contact prevalence"),
         variable = recode(variable, "Ciliate" = "Ciliate infection prevalence"),
         variable = recode(variable, "Cyanobacteria" = "Cyanobacteria prevalence"),
         variable = recode(variable, "Disease" = "Disease prevalence"),
         variable = recode(variable, "Fireworm" = "Fireworm predation prevalence"),
         variable = recode(variable, "Overgrowth" = "Overgrowth prevalence"),
         variable = recode(variable, "Snail" = "Snail prevalence"),
         Time.Point = recode(Time.Point, "T1_cm2" = "T1"),
         Time.Point = recode(Time.Point, "T2_cm2" = "T2"),
         Time.Point = recode(Time.Point, "T3_cm2" = "T3"))

write.csv(summary_large, "summary_large_table.csv", row.names = F)








#####################################################################################################################################################################################################################################################################
## Variance Partitioning - Supplement Table S4
#####################################################################################################################################################################################################################################################################

##### Notes:
# Block types: Block_Cat (Inshore or Offshore), Actual.space (Full or Open [space])
# using rptR package (Stoffel et al. 2017) as a means of variance partitioning. See https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html 
# summary(model$mod) to get variance values
# summary(model) to get bootstrap and permutation results for repeatability (ie % variance explained)



##### Disease prevalence

# Months 3 and 6 have 3 factors: (1|Geno) + (1|Block_Cat) + (1|Batch)
rpt_d_03 <- rpt(Disease ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                grname = c("Geno", "Block_Cat", "Batch", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M03",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_d_03$mod)
summary(rpt_d_03)

rpt_d_06 <- rpt(Disease ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                grname = c("Geno", "Block_Cat", "Batch", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M06",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_d_06$mod)
summary(rpt_d_06)

# Months 9 and 20 have 2 factors: (1|Geno) + (1|Block_Cat)
rpt_d_09 <- rpt(Disease ~ (1|Geno) + (1|Block_Cat), 
                grname = c("Geno", "Block_Cat", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M09",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_d_09$mod)
summary(rpt_d_09)

rpt_d_20 <- rpt(Disease ~ (1|Geno) + (1|Block_Cat), 
                grname = c("Geno", "Block_Cat", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M20",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_d_20$mod)
summary(rpt_d_20)



##### Ciliate prevalence

# Months 3 and 6 have 3 factors: (1|Geno) + (1|Block_Cat) + (1|Batch)
rpt_c_03 <- rpt(Ciliate ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                grname = c("Geno", "Block_Cat", "Batch", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M03",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_c_03$mod)
summary(rpt_c_03)

rpt_c_06 <- rpt(Ciliate ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                grname = c("Geno", "Block_Cat", "Batch", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M06",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_c_06$mod)
summary(rpt_c_06)

# Months 9 and 20 have 2 factors: (1|Geno) + (1|Block_Cat)
rpt_c_09 <- rpt(Ciliate ~ (1|Geno) + (1|Block_Cat), 
                grname = c("Geno", "Block_Cat", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M09",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_c_09$mod)
summary(rpt_c_09)

rpt_c_20 <- rpt(Ciliate ~ (1|Geno) + (1|Block_Cat), 
                grname = c("Geno", "Block_Cat", "Residual"), 
                data = Mon_wide[Mon_wide$Time.Point == "M20",], 
                datatype = "Binary", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_c_20$mod)
summary(rpt_c_20)



##### Partial mortality

# Months 3 and 6 have 3 factors: (1|Geno) + (1|Block_Cat) + (1|Batch)
rpt_pm_03 <- rpt(DeadTiss_arcsin ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                grname = c("Geno", "Block_Cat", "Batch", "Residual"),
                data = Mon_wide[Mon_wide$Time.Point == "M03",], 
                datatype = "Gaussian", 
                nboot = 1000, 
                npermut = 1000)

summary(rpt_pm_03$mod)
summary(rpt_pm_03)

rpt_pm_06 <- rpt(DeadTiss_arcsin ~ (1|Geno) + (1|Block_Cat) + (1|Batch), 
                 grname = c("Geno", "Block_Cat", "Batch", "Residual"),
                 data = Mon_wide[Mon_wide$Time.Point == "M06",], 
                 datatype = "Gaussian", 
                 nboot = 1000, 
                 npermut = 1000)

summary(rpt_pm_06$mod)
summary(rpt_pm_06)



##### Net change in area

# Year 1
rpt_nca_12 <- rpt(nc_T1T2 ~ (1|Geno) + (1|Block_Cat), 
                 grname = c("Geno", "Block_Cat", "Residual"),
                 data = Areas_clusters_all, 
                 datatype = "Gaussian", 
                 nboot = 1000, 
                 npermut = 1000)

summary(rpt_nca_12$mod)
summary(rpt_nca_12)

# Year 2
rpt_nca_23 <- rpt(nc_T2T3 ~ (1|Geno) + (1|Block_Cat), 
                 grname = c("Geno", "Block_Cat", "Residual"),
                 data = Areas_clusters_all, 
                 datatype = "Gaussian", 
                 nboot = 1000, 
                 npermut = 1000)

summary(rpt_nca_23$mod)
summary(rpt_nca_23)

# Both Years
rpt_nca_13 <- rpt(nc_T1T3 ~ (1|Geno) + (1|Block_Cat), 
                 grname = c("Geno", "Block_Cat", "Residual"),
                 data = Areas_clusters_all, 
                 datatype = "Gaussian", 
                 nboot = 1000, 
                 npermut = 1000)

summary(rpt_nca_13$mod)
summary(rpt_nca_13)




##### Diagnostics
# Paste rpt model from above in glmmTMB format
# Partial tissue mortality example:
rpt_diag <-glmmTMB(DeadTiss_arcsin ~ 1 + (1|Geno) + (1|Block_Cat) +  (1|Batch),
                      family = gaussian,
                      data = Mon_wide[Mon_wide$Time.Point == "M03",])

par(mfrow=c(1,2))
simulationOutput <- simulateResiduals(fittedModel = rpt_pm_diag, n = 250)
testZeroInflation(simulationOutput) # look for red line in center and not sig pvalue
testDispersion(simulationOutput)    # look for red line in center and not sig pvalue
plot(simulationOutput)              # look for points to match on line on the left graph and for the lines to be horizontal on the right graphs



#####################################################################################################################################################################################################################################################################
## Genotype Variation: Stressor Prevalence - Table 1 stats
#####################################################################################################################################################################################################################################################################


##### Modify Batch and time point subsetting as needed
# Time.Point == "M03" or "M06" or "M09" or "M20"
# Batch == "A" or "B"

##### Stressor variables measured:
# Disease
# Ciliate
# Cyanobacteria
# Bottom_contact
# Fireworm
# Overgrowth
# Snail
# Unknown


##### Mean (Range)
# Modify Batch and time point subsetting as needed

# Disease
mean(stress_prevs$dis_M03[stress_prevs$Batch == "A"])
min(stress_prevs$dis_M03[stress_prevs$Batch == "A"])
max(stress_prevs$dis_M03[stress_prevs$Batch == "A"])

# Ciliate
mean(stress_prevs$cil_M03[stress_prevs$Batch == "A"])
min(stress_prevs$cil_M03[stress_prevs$Batch == "A"])
max(stress_prevs$cil_M03[stress_prevs$Batch == "A"])




##### Generalized Linear Model -- Disease & Ciliate

# Months 3 and 6 (Batch as fixed effect)
sp <- glm(Ciliate ~ Geno + Batch,
              family = binomial,
              data = Mon_wide[Mon_wide$Time.Point == "M03",])
sp_null <- glm(Ciliate ~ Batch,
                   family = binomial,
                   data = Mon_wide[Mon_wide$Time.Point == "M03",])
anova(sp, sp_null, test="Chisq") 


# Months 9 and 20 (Only Batch A)
sp <- glm(Disease ~ Geno,
          family = binomial,
          data = Mon_wide[Mon_wide$Time.Point == "M09" & Mon_wide$Batch == "A",])
sp_null <- glm(Disease ~ 1,
               family = binomial,
               data = Mon_wide[Mon_wide$Time.Point == "M09" & Mon_wide$Batch == "A",])
anova(sp, sp_null, test="Chisq") 




##### Diagnostics
# Binomial (logistic regression) assumptions:
# Response data is a vector of 0s and 1s 
# In logistic regression, overdispersion is not possible


##### Moran's I to check for signal of spatial autocorrelation, using DHARMa
# run on response variables that exhibited genotypic variation: Disease M03, Ciliate M03, Ciliate M06

# test batches separately because xy coordinates were derived in separate unliked orthomosaics (T1, T3), so they're not in the exact same xy coordinate space (eg T1 is portrait orientation, T3 is landscape)
Mon_wide_sub <- Mon_wide %>% 
  filter(Batch == "B")

# run test model to calculate Moran's I
test_model <- glm(Disease ~ Geno,
          family = binomial,
          data = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",])

# simulate residuals
sim <- simulateResiduals(test_model)

# test for spatial autocorrelation
testSpatialAutocorrelation(sim,
                           x = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",]$TagLab.Centroid.x,
                           y = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",]$TagLab.Centroid.y)









#####################################################################################################################################################################################################################################################################
## Genotype Variation: Partial Tissue Mortality - Table 1 stats
#####################################################################################################################################################################################################################################################################


##### Modify Batch and time point as needed
# Time.Point == "M03" or "M06"
# Batch == "A" or "B"


##### Mean (Range)
# Use non-transformed data
mean(pm_avg$DeadTiss_avg[pm_avg$Batch == "B" & pm_avg$Time.Point == "M03"])
min(pm_avg$DeadTiss_avg[pm_avg$Batch == "B" & pm_avg$Time.Point == "M03"])
max(pm_avg$DeadTiss_avg[pm_avg$Batch == "B" & pm_avg$Time.Point == "M03"])


##### Generalized Linear Model
# use arcsin-transformed data
pm <- glm(DeadTiss_arcsin ~ Geno + Batch,
              family = gaussian,
              data = Mon_wide[Mon_wide$Time.Point == "M03",])
pm_null <- glm(DeadTiss_arcsin ~ Batch,
                   family = gaussian,
                   data = Mon_wide[Mon_wide$Time.Point == "M03",])
anova(pm, pm_null, test="Chisq") 




##### Diagnostics
# use arcsin-transformed data

# Data distributions
hist(Mon_wide$DeadTiss_arcsin[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M06"])
shapiro.test(Mon_wide$DeadTiss_arcsin[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M06"]) # p-value should be > 0.05 to be considered NOT significantly different from normal distribution

# Base R
pm_diag <- lm(DeadTiss_arcsin ~ Geno, data = Mon_wide[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M06",]) # need to subset within batches, can't run on mixed effects model
par(mfrow=c(2,2))
plot(pm_diag)

# DHARMa
simulationOutput <- simulateResiduals(fittedModel = pm, n = 250)
par(mfrow=c(1,2))
testZeroInflation(simulationOutput) # look for red line in center and not sig pvalue
testDispersion(simulationOutput)  # look for red line in center and not sig pvalue
plot(simulationOutput) # look for points to match on line on the left graph and for the lines to be horizontal on the right graphs




##### Moran's I to check for signal of spatial autocorrelation, using DHARMa
# run on response variables that exhibited genotypic variation: partial tissue mortality M03 nd M06

# test batches separately because xy coordinates were derived in separate unliked orthomosaics (T1, T3), so they're not in the exact same xy coordinate space (eg T1 is portrait orientation, T3 is landscape)
Mon_wide_sub <- Mon_wide %>% 
  filter(Batch == "B")

# run test model to calculate Moran's I
test_model <- glm(DeadTiss_arcsin ~ Geno,
                  family = gaussian,
                  data = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",])

# simulate residuals
sim <- simulateResiduals(test_model)

# test for spatial autocorrelation
testSpatialAutocorrelation(sim,
                           x = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",]$TagLab.Centroid.x,
                           y = Mon_wide_sub[Mon_wide_sub$Time.Point == "M03",]$TagLab.Centroid.y)





#####################################################################################################################################################################################################################################################################
## Genotype Variation: Net Change in Area - Table 1 stats
#####################################################################################################################################################################################################################################################################

##### Modify Year as needed
# Year == "1" or "2" or "Both Years"

##### Response variables
# nc_T1T2 = Year 1
# nc_T2T3 = Year 2
# nc_T1T3 = Both years


##### Mean (Range)
mean(growth_avg_long$avg_nc[growth_avg_long$Year == "Both Years"])
min(growth_avg_long$avg_nc[growth_avg_long$Year == "Both Years"])
max(growth_avg_long$avg_nc[growth_avg_long$Year == "Both Years"])



##### Generalized Linear Model
nca <- glm(nc_T2T3 ~ Geno,
               family = gaussian,
               data = Areas_clusters_all[Areas_clusters_all$Batch == "A",])

nca_null <- glm(nc_T2T3 ~ 1,
                   family = gaussian,
                   data = Areas_clusters_all[Areas_clusters_all$Batch == "A",])

anova(nca, nca_null, test="Chisq") 




##### Diagnostics

# Data distributions
hist(Areas_clusters_all$nc_T1T2)
shapiro.test(Areas_clusters_all$nc_T1T2) # p-value should be > 0.05 to be considered NOT significantly different from normal distribution

# Base R
nca_diag <- lm(nc_T1T2 ~ Geno, data = Areas_clusters_all)
par(mfrow=c(2,2))
plot(nca_diag)

# DHARMa
simulationOutput <- simulateResiduals(fittedModel = nca, n = 250)
par(mfrow=c(1,2))
testZeroInflation(simulationOutput) ## look for red line in center and not sig pvalue
testDispersion(simulationOutput)  ## look for red line in center and not sig pvalue
plot(simulationOutput) # look for points to match on line on the left graph and for the lines to be horizontal on the right graphs




##### Moran's I to check for signal of spatial autocorrelation, using DHARMa
# run on response variables that exhibited genotypic variation: net change in coral cover Y1, Y2, Both years

# run test model to calculate Moran's I
test_model <- glm(nc_T2T3 ~ Geno,
                  family = gaussian,
                  data = Areas_clusters_all[Areas_clusters_all$Batch == "A",])

# deal with differences in dimensions between residuals and xy coords: some nc values are NA for clusters that didn't survive until later timepoints - glm drops them on its own, but we need to manually ID and subset them out for testSpatialAutocorrelation()
keep <- as.numeric(rownames(model.frame(test_model))) 

# simulate residuals
sim <- simulateResiduals(test_model)

# test for spatial autocorrelation
testSpatialAutocorrelation(sim,
                           x = Areas_clusters_all$TagLab.Centroid.x[keep],
                           y = Areas_clusters_all$TagLab.Centroid.y[keep])








#####################################################################################################################################################################################################################################################################
## Batch-Batch and Year-Year Correlations - Table 2 stats
#####################################################################################################################################################################################################################################################################

##### Disease prevalence
# Diagnostics - not normal distribution, use "kendall"
hist(batch_corr$dp_A[batch_corr$Time.Point == "M06"])
shapiro.test(batch_corr$dp_A[batch_corr$Time.Point == "M06"])
hist(batch_corr$dp_B[batch_corr$Time.Point == "M06"])
shapiro.test(batch_corr$dp_B[batch_corr$Time.Point == "M06"])

# kendall correlation test
cor.test(batch_corr$dp_A[batch_corr$Time.Point == "M06"], batch_corr$dp_B[batch_corr$Time.Point == "M06"], method = "kendall")



##### Ciliate prevalence
# Diagnostics - normal distribution, use "pearson"
hist(batch_corr$cp_A[batch_corr$Time.Point == "M03"])
shapiro.test(batch_corr$cp_A[batch_corr$Time.Point == "M03"])
hist(batch_corr$cp_B[batch_corr$Time.Point == "M03"])
shapiro.test(batch_corr$cp_B[batch_corr$Time.Point == "M03"])

# pearson correlation test
cor.test(batch_corr$cp_A[batch_corr$Time.Point == "M03"], batch_corr$cp_B[batch_corr$Time.Point == "M03"], method = "pearson")



##### Partial tissue mortality
# Diagnostics - normal distribution, use "pearson"
hist(batch_corr$pm_A[batch_corr$Time.Point == "M06"])
shapiro.test(batch_corr$pm_A[batch_corr$Time.Point == "M06"])
hist(batch_corr$pm_B[batch_corr$Time.Point == "M06"])
shapiro.test(batch_corr$pm_B[batch_corr$Time.Point == "M06"])

# pearson correlation test
cor.test(batch_corr$pm_A[batch_corr$Time.Point == "M06"], batch_corr$pm_B[batch_corr$Time.Point == "M06"], method = "pearson")



##### Net change in area
# linear regression
years <- lm(avg_nc_T2T3 ~ avg_nc_T1T2, data = Areas_genos_M03)
summary(years)

##### Diagnostics

# Data distributions
hist(Areas_genos_M03$avg_nc_T1T2)
shapiro.test(Areas_genos_M03$avg_nc_T1T2)
hist(Areas_genos_M03$avg_nc_T2T3)
shapiro.test(Areas_genos_M03$avg_nc_T2T3)

# Base R
par(mfrow=c(2,2))
plot(years)

##### Diagnostics -- DHARMa
simulationOutput <- simulateResiduals(fittedModel = years, n = 250)
par(mfrow=c(1,2))
testZeroInflation(simulationOutput) ## look for red line in center and not sig pvalue
testDispersion(simulationOutput)  ## look for red line in center and not sig pvalue
plot(simulationOutput) # look for points to match on line on the left graph and for the lines to be horizontal on the right graphs





#####################################################################################################################################################################################################################################################################
## Disease and Net Change in Area Association - Stats associated with Figure 4
#####################################################################################################################################################################################################################################################################


##### Geno 14 Batch A - disease prevalence and net change in area metrics
(stress_prevs$dis_M03[stress_prevs$Batch == "A" & stress_prevs$Geno == 14])              # G14 disease prevalence across all clusters at M03 
max(stress_prevs$dis_M03[stress_prevs$Batch == "A"])                                     # Maximum disease prevalence across all clusters at M03 
(Areas_genos_M03$dis_prop[Areas_genos_M03$Batch == "A" & Areas_genos_M03$Geno == 14])    # G14 disease prevalence across clusters at M03 that survived to T2
max(Areas_genos_M03$dis_prop[Areas_genos_M03$Batch == "A"])                              # Maximum disease prevalence across clusters at M03 that survived to T2
(Areas_genos_M03$avg_nc_T1T2[Areas_genos_M03$Batch == "A" & Areas_genos_M03$Geno == 14]) # G14 net change in area after 1 year
min(Areas_genos_M03$avg_nc_T1T2[Areas_genos_M03$Batch == "A"])                           # Maximum net change in area after 1 year




##### Regression - Disease @M03 & Net Change in Area @Y1 

# One year, full dataset
reg1 <- lm(avg_nc_T1T2 ~ dis_prop, data = Areas_genos_M03)
summary(reg1)

# One year, exclude G14
reg2 <- lm(avg_nc_T1T2 ~ dis_prop, data = Areas_genos_M03[Areas_genos_M03$Geno != "14",])
summary(reg2)

# Two years, full dataset
reg3 <- lm(avg_nc_T1T3 ~ dis_prop, data = Areas_genos_M03)
summary(reg3)


# Diagnostics -- data distributions
hist(Areas_genos_M03$avg_nc_T1T2)
shapiro.test(Areas_genos_M03$avg_nc_T1T2) # p-value should be > 0.05 to be considered NOT significantly different from normal distribution

# Diagnostics -- Base R
par(mfrow=c(2,2))
plot(reg1)

# Diagnostics -- DHARMa
simulationOutput <- simulateResiduals(fittedModel = reg1, n = 250)
par(mfrow=c(1,2))
testZeroInflation(simulationOutput) ## look for red line in center and not sig pvalue
testDispersion(simulationOutput)  ## look for red line in center and not sig pvalue
plot(simulationOutput) # look for points to match on line on the left graph and for the lines to be horizontal on the right graphs












#####################################################################################################################################################################################################################################################################
## Nursery Analyses - Table 3 stats
#####################################################################################################################################################################################################################################################################

##### Definitions and y-axis metrics

# ProductionCycle_ID ---------- Unique identifier for each Production cycle

# CycleStartDate -------------- The date the production cycle starts; initiated by hanging initial stock on the Unit
# CycleEndDate ---------------- Date the production cycle ends; initiated manually normally when no more corals are outplanted
# CycleDuration --------------- Days = CycleStartDate - CycleEndDate 

# Initial_Broodstock ---------- number of corals that are left after outplanting and are used to start a new production cycle. These fragments can differ in size
# SumOfTotal_Coral_Produced --- initial stock + corals produced
# OutputRatio ----------------- SumOfOutplanted / SumOfTotal_Coral_Produced *100%
# Net_Produced ---------------- total produced minus initial broodstock

# Production.Capacity.Days ---- The number of days it takes for a unit to reach 95% capacity from the CycleStartDate
# Ready.To.OP.Days.Calculated - The number of days it takes for a unit to reach 95% capacity from the CycleStartDate AND reach 80% Ready to Outplant
# Ready.To.Outplant ----------- number of fragments that are size ready to be outplanted

# newfrgs_yr: new fragments per, calculated based on the number of days (Production capacity days) it took to fill 95% of tree (95 slots - initial broodstock): (95-IB)/PCD
# newfrgs_r_lg: same as above, but logged each of the numerator values first: log(95)-log(IB)/PCD

# Stressor pcts: percents calculated as number of fragments with the stressor divided by number of total number fragments (initial + produced) * 100
# Mortality.Pct
# Survivor.Pct
# Disease.Pct



##### Correlations between nursery and outplant performance:

# Nursery disease-based mortality (Batch A cycle): Disease.Pct_A
# Outplant disease prevalence: dis_OP (must subset Time.Point)

# Nursery total natural mortality (Batch A cycle): Mort.Natural.Pct_A
# Outplant partial tissue mortality: partmort_OP (must subset Time.Point)

# Nursery growth rate (Batch A cycle): newfrgs_yr_lg_A
# Outplant net change in area: Avg_nc12, Avg_nc23, Avg_nc13



########## BATCH A

# Choose time point to focus on
Month_subset_A <- Nurs_OP_Monwide_A$Time.Point == "M03"
Month_subset_A <- Nurs_OP_Monwide_A$Time.Point == "M06"
Month_subset_A <- Nurs_OP_Monwide_A$Time.Point == "M09"
Month_subset_A <- Nurs_OP_Monwide_A$Time.Point == "M20"



##### Disease correlations
# Check data distributions for nursery and outplant metrics -- Nursery data is not normal, Outplant data is normal
hist(Nurs_OP_Monwide_A$Disease.Pct_A[Month_subset_A])
hist(Nurs_OP_Monwide_A$dis_OP[Month_subset_A])
shapiro.test(Nurs_OP_Monwide_A$Disease.Pct_A[Month_subset_A])
shapiro.test(Nurs_OP_Monwide_A$dis_OP[Month_subset_A])
# Correlation test
cor.test(Nurs_OP_Monwide_A$Disease.Pct_A[Month_subset_A], Nurs_OP_Monwide_A$dis_OP[Month_subset_A], method ="kendall")



##### Mortality correlations
# Check data distributions for nursery and outplant metrics -- Nursery data is not normal, Outplant data is normal
hist(Nurs_OP_Monwide_A$Mort.Natural.Pct_A[Month_subset_A])
hist(Nurs_OP_Monwide_A$partmort_OP[Month_subset_A])
shapiro.test(Nurs_OP_Monwide_A$Mort.Natural.Pct_A[Month_subset_A])
shapiro.test(Nurs_OP_Monwide_A$partmort_OP[Month_subset_A])
# Correlation test
cor.test(Nurs_OP_Monwide_A$Mort.Natural.Pct_A[Month_subset_A], Nurs_OP_Monwide_A$partmort_OP[Month_subset_A], method ="kendall")



##### Growth correlations 
# Check data distributions for nursery and outplant metrics -- Nursery data is not normal, Outplant data is normal
hist(Nurs_OP_Areas$newfrgs_yr_lg_A)
hist(Nurs_OP_Areas$Avg_nc12)
hist(Nurs_OP_Areas$Avg_nc23)
hist(Nurs_OP_Areas$Avg_nc13)
shapiro.test(Nurs_OP_Areas$newfrgs_yr_lg_A)
shapiro.test(Nurs_OP_Areas$Avg_nc12)
shapiro.test(Nurs_OP_Areas$Avg_nc23)
shapiro.test(Nurs_OP_Areas$Avg_nc13)
# Correlation test
cor.test(Nurs_OP_Areas$newfrgs_yr_lg_A, Nurs_OP_Areas$Avg_nc12, method ="kendall")
cor.test(Nurs_OP_Areas$newfrgs_yr_lg_A, Nurs_OP_Areas$Avg_nc23, method ="kendall")
cor.test(Nurs_OP_Areas$newfrgs_yr_lg_A, Nurs_OP_Areas$Avg_nc13, method ="kendall")




########## BATCHES A AND B

# Choose time point to focus on
Month_subset_AB <- Nurs_OP_Monwide_AB$Time.Point == "M03"
Month_subset_AB <- Nurs_OP_Monwide_AB$Time.Point == "M06"
Month_subset_AB <- Nurs_OP_Monwide_AB$Time.Point == "M09"
Month_subset_AB <- Nurs_OP_Monwide_AB$Time.Point == "M20"



##### Disease correlations
# Check data distributions for nursery and outplant metrics -- Nursery data is not normal, Outplant data is normal
hist(Nurs_OP_Monwide_AB$Disease.Pct_A[Month_subset_AB])
hist(Nurs_OP_Monwide_AB$dis_OP[Month_subset_AB])
shapiro.test(Nurs_OP_Monwide_AB$Disease.Pct_A[Month_subset_AB])
shapiro.test(Nurs_OP_Monwide_AB$dis_OP[Month_subset_AB])
# Correlation test
cor.test(Nurs_OP_Monwide_AB$Disease.Pct_A[Month_subset_AB], Nurs_OP_Monwide_AB$dis_OP[Month_subset_AB], method ="kendall")



##### Mortality correlations
# Check data distributions for nursery and outplant metrics -- Nursery data is not normal, Outplant data is normal
hist(Nurs_OP_Monwide_AB$Mort.Natural.Pct_A[Month_subset_AB])
hist(Nurs_OP_Monwide_AB$partmort_OP[Month_subset_AB])
shapiro.test(Nurs_OP_Monwide_AB$Mort.Natural.Pct_A[Month_subset_AB])
shapiro.test(Nurs_OP_Monwide_AB$partmort_OP[Month_subset_AB])
# Correlation test
cor.test(Nurs_OP_Monwide_AB$Mort.Natural.Pct_A[Month_subset_AB], Nurs_OP_Monwide_AB$partmort_OP[Month_subset_AB], method ="kendall")



##### Growth correlations 
# Not relevant, net change in area is only Batch A


# Calculations in PV_dataframes_241231.R




#####################################################################################################################################################################################################################################################################
## Other In-Text Reported Stats
#####################################################################################################################################################################################################################################################################


##### Percent of Batch A outplanted clusters that experienced disease at least once by month 20
# Batch A only, all time points
dhist_BatchA <- Mon_wide %>% 
  filter(Batch == "A") %>% 
  group_by(Batch, Geno, Geno_chr, Cluster, GenoCluster) %>%               # per GenoCluster (the rest of the grouping variables are just useful metadata)
  summarise(dis_occurrences = sum(Disease, na.rm = T),                               # across all documented time points, how many was that cluster sick. Each row is per cluster and should be 185 total rows to account for every outplanted cluster
            n_Time.Points = sum(survived)) %>%                            # how many time points are you summing across
  mutate(dis_history = ifelse(dis_occurrences > 0, 1, 0),                 # turning disease occurrences into a binary sick/healthy disease history. 1 = got sick at least 1x. 0 = never ever got sick <3
         diseased_life_proportion = dis_occurrences/n_Time.Points) #%>%    # calculating proportion of life spent sick 
#filter(n_Time.Points == "4")                                              # can filter out clusters that didn't make it to the final Time.Point 

mean(dhist_BatchA$dis_history, na.rm = T)
  


#####################################################################################################################################################################################################################################################################
## Plotting Figure 3
#####################################################################################################################################################################################################################################################################


##### 3A - disease across genos
Dis_AB <- 
  ggplot(data = stress_prevs, aes(x = Geno, y = dis_M03, group = Batch))+
  geom_hline(yintercept= mean(Mon_wide$Disease[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M03"]), linetype="dashed", color = Batch.colors.2[1], linewidth = horizln.sz)+  
  geom_hline(yintercept= mean(Mon_wide$Disease[Mon_wide$Batch == "B" & Mon_wide$Time.Point == "M03"]), linetype="dashed", color = Batch.colors.2[2], linewidth = horizln.sz)+
  geom_point(size = point.sz, aes(color = Batch))+
  scale_color_manual(values= Batch.colors.2)+
  ylim(0,1)+
  ylab("Disease prevalence - M3")+
  xlab("Genotype")+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black",size = axis.txt.geno.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.title = element_text(size = leg.txt.sz),
        legend.text = element_text(size = leg.txt.sz))



##### 3B - partial mortality across genos
PartMort_AB <- 
  ggplot(data= pm_avg[pm_avg$Time.Point == "M06",], aes(x = Geno, y = DeadTiss_cluster_mean, group = Batch))+
  geom_hline(yintercept= mean(Mon_wide$DeadTiss_cluster[Mon_wide$Batch == "A" & Mon_wide$Time.Point == "M06"]), linetype="dashed", color = Batch.colors.2[1], linewidth = horizln.sz)+
  geom_hline(yintercept= mean(Mon_wide$DeadTiss_cluster[Mon_wide$Batch == "B" & Mon_wide$Time.Point == "M06"]), linetype="dashed", color = Batch.colors.2[2], linewidth = horizln.sz)+
  geom_point(size = point.sz, aes(color = Batch))+
  geom_errorbar(aes(ymin=DeadTiss_cluster_ci.lwr, ymax=DeadTiss_cluster_ci.upr, color = Batch), width = error.wd.sz, linewidth = error.ln.sz)+
  scale_color_manual(values=Batch.colors.2)+
  ylim(0,100)+
  ylab("Partial mortality - M6")+
  xlab("Genotype")+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black",size = axis.txt.geno.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.title = element_text(size = leg.txt.sz),
        legend.text = element_text(size = leg.txt.sz))



##### 3C - net change in area across genos
NCA_12 <- 
  ggplot(data = growth_avg_long[growth_avg_long$Year != "Both Years",], aes(x = Geno, y = net_change_mean, group = Year))+
  geom_hline(yintercept = mean(Areas_clusters_all$nc_T1T2, na.rm = T), linetype="dashed", color = Year.colors[1], linewidth = horizln.sz)+
  geom_hline(yintercept = mean(Areas_clusters_all$nc_T2T3, na.rm = T), linetype="dashed", color = Year.colors[2], linewidth = horizln.sz)+
  geom_errorbar(aes(ymin = net_change_ci.lwr, ymax = net_change_ci.upr, color = Year), width = error.wd.sz, linewidth = error.ln.sz)+
  geom_point(size = point.sz, aes(color = Year))+
  scale_color_manual(values = Year.colors)+
  ylab("Net change in coral cover")+
  #ylab(expression(paste("Net change coral cover (",cm^{2},"/day)")))+
  xlab("Genotype")+
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.geno.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.title = element_text(size = leg.txt.sz),
        legend.text = element_text(size = leg.txt.sz))



##### 3D - disease batch correlations
Batchcorr_Dis <-
  ggplot(data = batch_corr[batch_corr$Time.Point == "M03",], aes(x = dp_A, y = dp_B))+
  geom_point(size = point.sz, aes(color = Geno))+
  geom_smooth(method= "lm", color = "black", se=F, linewidth = regln.sz)+
  scale_colour_manual(values = Geno.colors)+
  xlab("Batch A disease prevalence")+
  ylab("Batch B disease prevalence")+
  xlim(0,1)+
  ylim(0,1)+
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.position = "none")



##### 3E - partial mortality batch correlations
Batchcorr_PM <-
  ggplot(data = batch_corr[batch_corr$Time.Point == "M06",], aes(x = pm_mean_A, y = pm_mean_B))+
  geom_point(size = point.sz, aes(color = Geno))+
  geom_errorbar(aes(ymin = pm_ci.lwr_B, ymax = pm_ci.upr_B, color = Geno), width = error.wd.sz, linewidth = error.ln.sz)+
  geom_errorbar(aes(xmin = pm_ci.lwr_A, xmax = pm_ci.upr_A, color = Geno), height = error.wd.sz, linewidth = error.ln.sz, orientation = "y")+  
  geom_smooth(method= "lm", color = "black", se=F, linewidth = regln.sz)+
  scale_colour_manual(values = Geno.colors)+
  xlab("Batch A partial mortality")+
  ylab("Batch B partial mortality")+
  xlim(0,100)+
  ylim(0,100)+
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.position = "none")



##### 3F - net change in area year correlations
Yearcorr_NCA <-
  ggplot(data = Areas_genos_M03[Areas_genos_M03$Batch == "A",], aes(x = nc_T1T2_mean, y = nc_T2T3_mean))+
  geom_point(size = point.sz, aes(color = Geno))+
  geom_errorbar(aes(ymin = nc_T2T3_ci.lwr, ymax = nc_T2T3_ci.upr, color = Geno), width = error.wd.sz, linewidth = error.ln.sz)+
  geom_errorbar(aes(xmin = nc_T1T2_ci.lwr, xmax = nc_T1T2_ci.upr, color = Geno), width = error.wd.sz, linewidth = error.ln.sz, orientation = "y")+
  geom_smooth(method= "lm", color = "black", se=F, linewidth = regln.sz)+
  scale_colour_manual(values = Geno.colors)+
  xlab(expression(paste("Year 1 net change coral cover")))+
  ylab(expression(paste("Year 2 net change coral cover")))+
  guides(color=guide_legend(title="Genotype"))+
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.title = element_text(size = leg.txt.sz),
        legend.text = element_text(size = leg.txt.sz))


##### Put all the pieces together

fig3<-
  ggdraw()+
  draw_plot(Dis_AB, x = 0, y = 0.51, width = 0.32, height = 0.43)+
  draw_plot(PartMort_AB, x = 0.33, y = 0.51, width = 0.32, height = 0.43)+
  draw_plot(NCA_12, x = 0.66, y = 0.51, width = 0.32, height = 0.43)+
  draw_plot(Batchcorr_Dis, x = 0, y = 0, width = 0.3, height = 0.43)+
  draw_plot(Batchcorr_PM, x = 0.33, y = 0, width = 0.3, height = 0.43)+
  draw_plot(Yearcorr_NCA, x = 0.66, y = 0, width = 0.32, height = 0.43)+
  draw_plot_label(label= c("A", "B", "C",
                           "D", "E", "F"), size=16,
                  x= c(0.06, 0.385, 0.7,
                       0.06, 0.385, 0.705), 
                  y = c(0.94, 0.94, 0.94, 
                        0.43, 0.43, 0.43))+
  draw_label(expression(bold("Variation in performance among genotypes")), x = 0.5, y = 0.96, vjust = 0, size = axis.title.sz) +
  draw_label(expression(bold("Consistency in performance between batches and years")), x = 0.5, y = 0.45, vjust = 0, size = axis.title.sz)
  


ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch2_Punt Vierkant outplanting experiment\\Figures\\fig3.tiff", fig3, dpi=600, units="cm", width = 36, height = 22, bg = "white")


#####################################################################################################################################################################################################################################################################
## Plotting Figure 4
#####################################################################################################################################################################################################################################################################

##### Figure 4 -- disease/net change regression 
# Areas_genos_M03 made in STRESSOR TOLERANCE section (sub-section Regression) above, lines ~1180

DisEffReg <- 
  ggplot(data= Areas_genos_M03[Areas_genos_M03$Batch == "A",], aes(x = dis_prop, y = nc_T1T2_mean))+
  geom_point(size = point.sz, aes(color = Geno))+
  geom_errorbar(aes(ymin = nc_T1T2_ci.lwr, ymax = nc_T1T2_ci.upr, color = Geno), width = error.wd.sz, linewidth = error.ln.sz)+
  geom_smooth(method="lm", color = "black", se=F, linewidth = regln.sz)+
  scale_colour_manual(values = Geno.colors)+
  xlim(0,1)+
  xlab("Disease prevalence - M3")+
  ylab(expression(paste("Net change in coral cover (",cm^{2},"/day) - Year 1")))+
  theme_classic()+
  guides(color=guide_legend(title="Genotype"))+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black",size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.title = element_text(size = leg.txt.sz),
        legend.text = element_text(size = leg.txt.sz))+
  annotate("text", x=.8, y=1.5, 
           label = 'paste(R^2, "= 0.53, p = 0.017")',
           size = annotate.sz,
           parse =TRUE)

ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch2_Punt Vierkant outplanting experiment\\Figures\\fig4.tiff", DisEffReg, dpi=600, units="cm", width = 14, height = 12, bg = "white")










#####################################################################################################################################################################################################################################################################
## Plotting Figure 5
#####################################################################################################################################################################################################################################################################

nurs_dis <-
  ggplot(data = Nurs_OP_Monwide_A[Nurs_OP_Monwide_A$Time.Point == "M03",],
         aes(x = Disease.Tot.Pct_A, y = dis_OP))+
  geom_point(size = point.sz, aes(color = Geno), position=position_jitter(0.2))+
  scale_colour_manual(values = Geno.colors)+
  xlab("Disease-based mortality (%)")+
  ylab("Disease prevalence - M3")+
  xlim(-2,100)+
  ylim(0,1)+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black",size = axis.txt.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.position = "none")+
  annotate("text", x=72, y=0.95, label= "tau = -0.05, p = 0.859", size = annotate.sz)


nurs_mort <-
  ggplot(data = Nurs_OP_Monwide_A[Nurs_OP_Monwide_A$Time.Point == "M03",],
         aes(x = Mort.Natural.Pct_A, y = partmort_OP))+
  geom_point(size = point.sz+1, aes(color = Geno), position=position_jitter(0.2))+
  scale_colour_manual(values = Geno.colors)+
  xlab("Total natural mortality (%)")+
  ylab("Partial mortality (%) - M3")+
  xlim(-2,100)+
  ylim(0,100)+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black",size = axis.txt.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.position = "none")+
  annotate("text", x=75, y=95, label= "tau = 0.46, p = 0.082", size = annotate.sz)


nurs_ncowth <-
  ggplot(data = Nurs_OP_Areas, aes(x = newfrgs_yr_lg_A, y = Avg_nc12))+
  geom_point(size = point.sz+1, aes(color = Geno), position=position_jitter(0.2))+
  scale_colour_manual(values = Geno.colors)+
  guides(color=guide_legend(title="Genotype"))+
  xlab("Growth rate (new frag. per frag. per yr)")+
  ylab(expression(paste("Net change in coral cover (",cm^{2},"/day) - Y1")))+
  ylim(-2.2, 1.3)+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size = axis.title.sz), 
        axis.title.y = element_text(color="black", size = axis.title.sz-1.5),
        axis.text.x = element_text(color="black",size = axis.txt.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        legend.title = element_text(color = "black", size = leg.txt.sz),
        legend.text=element_text(color="black",size = leg.txt.sz))+
  annotate("text", x=7.2, y=1.2, label= "tau = 0.37, p = 0.163", size = annotate.sz)


fig5<-
  ggdraw()+
  draw_plot(nurs_dis, x = 0.05, y = 0.1, width = 0.29, height = 0.9)+
  draw_plot(nurs_mort, x = 0.35, y = 0.1, width = 0.29, height = 0.9)+
  draw_plot(nurs_ncowth, x = 0.65, y = 0.1, width = 0.35, height = 0.9)+
  draw_plot_label(label= c("A", "B", "C"), size=18,
                  x= c(0.11, 0.4, 0.705), 
                  y = c(0.99, 0.99, 0.99))+
  draw_label(expression(bold("Nursery performance")), x = 0.5, y = 0.035, vjust = 0.5, angle = 0, size = axis.title.sz+3)+
  draw_label(expression(bold("Outplant performance")), x = 0.01, y = 0.5, vjust = 0.5, angle = 90, size = axis.title.sz+3)


ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch2_Punt Vierkant outplanting experiment\\Figures\\fig5.tiff", fig5, dpi=600, units="cm", width = 36, height = 11, bg = "white")








