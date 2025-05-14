
#clear data & gc
rm(list=ls(all=TRUE));gc() 

#load packages
library(tidyverse)
library(stringi)
library(data.table)
library(respirometry)
library(respR)
library(FishResp)
library(dplyr)
library(purrr)
library(dplyr)
library(rlist)
library(lubridate)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggsignif)
library(ggtext)
library(ggplot2)
library(sjPlot) # table functions
library(sjmisc) # sample data
library(lme4) # fitting models
library(stargazer)
library(broom.mixed)
library(knitr)
library(nlme)
library(lme4)
library(nlme)
library(broom.mixed)
library(gt)
library(emmeans)



setwd("")
############
#  SET UP  #
############

# load mmr 
mmr <- read.delim("Laukaa_2023_MMR.txt", stringsAsFactors=TRUE)

# load smr
smr <- read.csv("Laukaa_Summer_2023_SMR_data_SUMMARY.csv", stringsAsFactors=TRUE)

# sampling data
chamberdat <- read.csv("combinedchamberdata.csv")
chamberdat$Chamber.No <- paste0("CH", chamberdat$chamber)
chamberdat$Chamber_ID <- as.factor(chamberdat$Chamber_ID)
# 

# family data
fam <- read.csv("Respirometry chamber plan.csv", stringsAsFactors=TRUE)
fam$Chamber <- as.factor(fam$Chamber)
fam$Family <- ifelse(fam$Family == "yes", NA, fam$Family)

chamberdat <- left_join(fam, chamberdat, by = c("Chamber" = "Chamber_ID", "Batch"))

chamberdat <- drop_na(chamberdat, "Family")
chamberdat <- chamberdat %>% filter(Batch != "B10" | Batch != "B9" | Dead != "dead" | Dead != "Dead") %>% dplyr:: select(-WHO, X.1, X.2)


#vgll3 <- read.csv("Complete_Data_Renamed/Laukaa_2023_vgll3.csv", stringsAsFactors=TRUE)
#vgll3 <- vgll3 %>% filter(Batch != "B10" | Batch != "B9" | Batch != "B9.2" ) %>% dplyr::select(-c(WHO, X, X.1))
#chamberdat <- right_join(chamberdat, vgll3, by = c("Chamber" = "Chamber.ID", "Batch", "Notes", "DNA_well", "Sampling.date"))


data <- left_join(smr, chamberdat, by = c("Ind"="Fish_ID", "Chamber.No"))
data <- left_join(data, mmr, by =c("Chamber.No", "pumpresp", "Batch"))

data$Temperature <- as.factor(data$Temperature)
data$Batch <- as.factor(data$Batch)
data$pumpresp <- as.factor(data$pumpresp)
data$Chamber.No <- as.factor(data$Chamber.No)
data$Family <- as.factor(data$Family)
data$Week <- recode(data$Batch, 
                    B1 = 1,
                    B2 = 1,
                    B3 = 2,
                    B4 = 2,
                    B5 = 3,
                    B6 = 3,
                    B7 = 4,
                    B8 = 4,
                    B9 = 4,
                    B9.2 = 4,
                    B10 = 4)
data$Week <- as.factor(data$Week)

#DATA IS PREPARED#

# create massadj AS, SMR and MMR
data <- data %>%
  mutate(SMRmassadj = resid(lm(log10(data$MLND) ~ log10(data$Weight_g))),
         MMRmassadj = resid(lm(log10(data$MR.abs) ~ log10(data$Weight_g))))
data <- data %>% 
  mutate(ASabs =  data$MR.abs-data$MLND)
#NAs produced
data$ASmass = resid(lm(log10(data$ASabs) ~ log10(data$Weight_g), na.action=na.exclude))

#Data visually inspected using hist, PR2.B1.2.3 removed from MR analysis.

#Points considered for removal:
#row 12, MR.abs outlier, PR2.B1.2.3
#row 12, MR.mass outlier, PR2.B1.2.3
#row 22 ASabs outlier PR3.B1.3.3
#row 26 ASabs outlier, PR3.B2.3.3
#row 27 ASmass outlier PR3.B4.3.1
data <- filter(data, Ind != "PR2.B1.2.3")

# refit data without the outliers
data <- data %>%
  mutate(SMRmassadj = resid(lm(log10(data$MLND) ~ log10(data$Weight_g))),
         MMRmassadj = resid(lm(log10(data$MR.abs) ~ log10(data$Weight_g))))
data <- data %>% 
  mutate(ASabs =  data$MR.abs-data$MLND)
#NAs produced
data$ASmass = resid(lm(log10(data$ASabs) ~ log10(data$Weight_g), na.action=na.exclude))



##################
#  PANEL GRAPHS  #
##################


lukePalette <- c("#FF8200", "#00B5E2", "#78BE20", "#54585A", "#004152", "#00442B", "#000000",
                     "#7F3F98", "#0033A0", "#CC0D82")


# ggplots: 
#mass vs SMR, MMR and AS
#Jenni has approved
plot1 <- ggplot(data, aes(x = log10(Weight_g), y = log10(MLND), colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = "Log Mass (g)",
       y = bquote(bold(atop("Log SMR (mg O"[2]*"h"^-1*")")))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)") +
  ylim(-1.5,1)

plot2 <- ggplot(data, aes(x = log10(Weight_g), y = log10(MR.abs), colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = "Log Mass (g)",
       y = bquote(bold(atop("Log MMR (mg O"[2]*"h"^-1*")")))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
  geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
  labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)")+ 
  ylim(c(-1.5, 1))

plot3 <- ggplot(data, aes(x = log10(Weight_g), y = log10(ASabs), colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = "Log Mass (g)",
       y = bquote(bold(atop("Log Aerobic Scope (mg O"[2]*"h"^-1*")")))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
  geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
  labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)") +
 ylim(c(-1.5, 1))

ggarrange(plot1, plot2, plot3, ncol = 3, common.legend = TRUE, legend="right")
require(grid)   # for the textGrob() function

figure1 <- ggarrange(plot1 + rremove("xlab"), plot2 + rremove("xlab"),plot3 + rremove("xlab"),# remove axis labels from plots
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "right",
                    align = "hv", font.label = list(size = 22)
                    ) 

figure1 <- annotate_figure(figure1, bottom = textGrob("Log Mass (g)", gp = gpar(fontface = "bold", cex = 1.3)))


ggsave("scatters.png", plot = figure1, width = 16, height = 8, units = "in", dpi = 300)




options(scipen = 999, digits=4)


#temp vs SMR, MMR and AS#temp vs SMR, MMR and ASTemperature

# Combine all variables into a long format

plot4 <- ggplot(data, aes(x = Temperature, y = SMRmassadj, colour = Temperature, fill = Temperature)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(aes(shape = Week), width = 0.2, height = 0, size = 3) +
  scale_fill_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
  scale_colour_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
    labs(x = "Temperature (\u00B0C)",
         y = bquote(bold("Mass-Adjusted SMR"))) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  )+ 
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.6, 0.6, 0.2), 
                     labels = scales::label_number(accuracy = 0.1)) +
geom_signif(test = t.test, comparisons = list(c("16", "22")), annotations = "**",
            map_signif_level=TRUE, y_position = 0.35,
            textsize = 8, colour = "black") +
  geom_signif(test = t.test, comparisons = list(c("16", "19")), annotations = "*",
              map_signif_level=TRUE, y_position = 0.25,
              textsize = 8, colour = "black")

plot4

plot5 <- ggplot(data, aes(x = Temperature, y = MMRmassadj, colour = Temperature, fill = Temperature)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(aes(shape = Week), width = 0.2, height = 0, size = 3) +
  scale_fill_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
  scale_colour_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
  labs(x = "Temperature (\u00B0C)",
       y = bquote(bold("Mass-Adjusted MMR"))) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  )+ 
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.8, 0.8, 0.2)) 

plot6 <- ggplot(data, aes(x = Temperature, y = ASmass, colour = Temperature, fill = Temperature)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(aes(shape = Week), width = 0.2, height = 0, size = 3) +
  scale_fill_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
  scale_colour_manual(guide = "none", values = c("#0033A0", "#7F3F98", "#CC0D82")) +
  labs(x = "Temperature (\u00B0C)",
       y = bquote(bold("Mass-Adjusted Aerobic Scope"))) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  )  +
  scale_y_continuous(limits = c(-1.4, 0.8), breaks = seq(-2, 2, 0.2)) 
 # geom_signif(test = t.test, comparisons = list(c("16", "22")), annotations = "*",
   #           map_signif_level=TRUE, y_position = 0.55,
    #          textsize = 5, colour = "black")


ggarrange(plot4, plot5, plot6, ncol = 3, common.legend = TRUE, legend="right")

require(grid)   # for the textGrob() function

figure2 <- ggarrange(plot4 + rremove("xlab"), plot5 + rremove("xlab"),plot6 + rremove("xlab"),# remove axis labels from plots
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "right",
                    align = "hv", font.label = list(size = 22))

figure2 <- annotate_figure(figure2, bottom = textGrob("Temperature (\u00B0C)", gp = gpar(fontface = "bold", cex = 1.3)))


ggsave("boxplots.png", plot = figure2, width = 16, height = 8, units = "in", dpi = 300)




plot9 <- ggplot(data, 
               aes(x = MMRmassadj, y = ASmass, colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = bquote(bold("Mass-Adjusted MMR")),
       y = bquote(bold("Mass-Adjusted Aerobic Scope"))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
  geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
  labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)") 
  #ylim(-0.7, 0.7)

plot8 <- ggplot(data, 
       aes(x = SMRmassadj, y = ASmass, colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = bquote(bold("Mass-Adjusted SMR")),
       y = bquote(bold("Mass-Adjusted Aerobic Scope"))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
  geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
  labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)") 

plot7 <- ggplot(data, 
                aes(x = SMRmassadj, y = MMRmassadj, colour = Temperature, fill = Temperature)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = T, alpha = 0.2) +  # Add linear regression line
  labs(x = bquote(bold("Mass-Adjusted SMR")),
       y = bquote(bold("Mass-Adjusted MMR"))) +
  scale_fill_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200"), labels=c("16", "19", "22", "All Data")) +
  theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_text(face = "bold", size = 20)
  ) +
  geom_smooth(method = lm, aes(group=1, colour = "#FF8200", fill = "#FF8200"), alpha = 0.3)+
  labs(color = "Temperature (\u00B0C)", fill = "Temperature (\u00B0C)") 


ggarrange(plot7, plot8, plot9, ncol = 3, common.legend = TRUE, legend="right")
require(grid)   # for the textGrob() function

figure3 <- ggarrange(plot7, plot8,plot9,# remove axis labels from plots
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "right",
                    align = "hv", font.label = list(size = 22))

ggsave("scatters2.png", plot = figure3, width = 16, height = 8, units = "in", dpi = 300)


#################
#   STATISTICS  #
#################

extract_random_effects <- function(model, data) {
  re_var <- as.data.frame(VarCorr(model))
  num_groups <- length(unique(data$Batch))
  
  intercept_var <- round(re_var$vcov[1], 4)
  intercept_sd <- round(sqrt(re_var$vcov[1]), 4)
  residual_var <- round(re_var$vcov[2], 4)
  residual_sd <- round(sqrt(re_var$vcov[2]), 4)
  
  return(c(intercept_var, intercept_sd, residual_var, residual_sd, num_groups))
}

data1 <- subset(data, !is.na(Family))
data1 <- data1 %>% select(-contains("X"))
data1$logSMR <- log10(data1$MLND)
data1$logMMR <- log10(data1$MR.abs)
data1$logAS <- log10(data1$ASabs)
data1$logweight <- log10(data1$Weight_g)


SMRmodel <- lmer(log10(MLND) ~ log10(Weight_g) + Temperature + (1|Batch) , data = data1)
summary(SMRmodel) #significance
SMRmodel.emm <- emmeans::emmeans(SMRmodel,  "Temperature")
pairs <- pairs(SMRmodel.emm)

# Format pairwise comparisons as a table
pairs_SMRmodel_emm_tidy <- as.data.frame(pairs)
# Use knitr to create an HTML table
html_table <- kable(pairs_SMRmodel_emm_tidy, 
                    caption = "Pairwise Comparisons for Temperature in SMR Model",
                    digits = 3, 
                    format = "html", 
                    col.names = c("Comparison", "Estimate", "SE", "df", "t-value", "p-value"))

# Save the table as HTML
writeLines(html_table, "SMRmodel_emmeans_pairs_table.html")
re_model1 <- extract_random_effects(SMRmodel, data1)
class(SMRmodel) <- "lmerMod"


MMRmodel <- lmer(log10(MR.abs) ~ log10(Weight_g) + Temperature + (1|Batch) , data = data1)
re_model2 <- extract_random_effects(MMRmodel, data1)
summary(MMRmodel) #no significance
class(MMRmodel) <- "lmerMod"

#lmer singular result, needs lme
ASmodel <- lme(log10(ASabs) ~ log10(Weight_g) + Temperature, random = ~1|Batch, data = data1, method = "REML")
ASmodel.emm <- emmeans::emmeans(ASmodel,  "Temperature")
pairs(ASmodel.emm) #no signifcance

#this is the table in the paper
tab_model(SMRmodel, MMRmodel, ASmodel,
          show.re.var = TRUE,
          dv.labels = c("SMR", "MMR", "AS"),
          file = "allMod_tabmodel.html",
          p.style = "stars")

######################
# descriptive tables #
######################


desc_stats <- data1 %>% 
  dplyr::group_by(Week, Temperature) %>% 
  dplyr::arrange(Week, Temperature) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Mean_Mass = mean(Weight_g, na.rm = TRUE),
    SD_Mass = sd(Weight_g, na.rm = TRUE),
    Min_Mass = min(Weight_g, na.rm = TRUE),
    Max_Mass = max(Weight_g, na.rm = TRUE),
    Mean_Length = mean(Total_length_mm, na.rm = TRUE),
    SD_Length = sd(Total_length_mm, na.rm = TRUE),
    Min_Length = min(Total_length_mm, na.rm = TRUE),
    Max_Length = max(Total_length_mm, na.rm = TRUE))

desc_stats2 <- data1 %>% dplyr::summarise(
  n = dplyr::n(),
  Mean_Mass = mean(Weight_g, na.rm = TRUE),
  SD_Mass = sd(Weight_g, na.rm = TRUE),
  Min_Mass = min(Weight_g, na.rm = TRUE),
  Max_Mass = max(Weight_g, na.rm = TRUE),
  Mean_Length = mean(Total_length_mm, na.rm = TRUE),
  SD_Length = sd(Total_length_mm, na.rm = TRUE),
  Min_Length = min(Total_length_mm, na.rm = TRUE),
  Max_Length = max(Total_length_mm, na.rm = TRUE))   

desc <- dplyr::bind_rows(desc_stats,desc_stats2)

desc$Week <- as.factor(desc$Week)
desc$Week <- ifelse(!is.na(desc$Week),desc$Week, "Total")
desc$Temperature <- as.factor(desc$Temperature)
desc$Temperature <- ifelse(!is.na(desc$Temperature),desc$Temperature, "Total")

## MR descriptive stats
desc_stats3 <- data1 %>% 
  dplyr::group_by(Week, Temperature) %>% 
  dplyr::arrange(Week, Temperature) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Mean_SMR = mean(MLND, na.rm = TRUE),
    SD_SMR = sd(MLND, na.rm = TRUE),
    Min_SMR = min(MLND, na.rm = TRUE),
    Max_SMR = max(MLND, na.rm = TRUE),
    Mean_MMR = mean(MR.abs, na.rm = TRUE),
    SD_MMR = sd(MR.abs, na.rm = TRUE),
    Min_MMR = min(MR.abs, na.rm = TRUE),
    Max_MMR = max(MR.abs, na.rm = TRUE),
    Mean_AS = mean(ASabs, na.rm = TRUE),
    SD_AS = sd(ASabs, na.rm = TRUE),
    Min_AS = min(ASabs, na.rm = TRUE),
    Max_AS = max(ASabs, na.rm = TRUE))


desc_stats4 <- data1 %>% dplyr::summarise(
  n = dplyr::n(),   
  Mean_SMR = mean(MLND, na.rm = TRUE),
  SD_SMR = sd(MLND, na.rm = TRUE),
  Min_SMR = min(MLND, na.rm = TRUE),
  Max_SMR = max(MLND, na.rm = TRUE),
  Mean_MMR = mean(MR.abs, na.rm = TRUE),
  SD_MMR = sd(MR.abs, na.rm = TRUE),
  Min_MMR = min(MR.abs, na.rm = TRUE),
  Max_MMR = max(MR.abs, na.rm = TRUE),
  Mean_AS = mean(ASabs, na.rm = TRUE),
  SD_AS = sd(ASabs, na.rm = TRUE),
  Min_AS = min(ASabs, na.rm = TRUE),
  Max_AS = max(ASabs, na.rm = TRUE))

desc2 <- rbind(desc_stats3,desc_stats4)

desc %>%
  gt() %>%
  cols_label(  
    Week = "Week",
    Temperature = "Temperature",
    Mean_Mass = "Mean Mass",
    SD_Mass = "SD Mass",
    Min_Mass = "Min Mass",
    Max_Mass = "Max Mass",
    Mean_Length = "Mean Length",
    SD_Length = "SD Length",
    Min_Length = "Min Length",
    Max_Length = "Max Length") %>%
  fmt_number(columns = c(3:9), decimals = 2) %>%
  tab_header(title = "Table 1: Descriptive Statistics for Mass (g) and Length (mm)")

desc2 %>%
  gt() %>%
  cols_label(
    Week = "Week",
    Temperature = "Temperature",
    Mean_SMR = "Mean SMR",
    SD_SMR = "SD SMR",
    Min_SMR = "Min SMR",
    Max_SMR = "Max SMR",
    Mean_MMR = "Mean MMR",
    SD_MMR = "SD MMR",
    Min_MMR = "Min MMR",
    Max_MMR = "Max MMR",
    Mean_AS = "Mean AS",
    SD_AS = "SD AS",
    Min_AS = "Min AS",
    Max_AS = "Max AS") %>%
  fmt_number(columns = c(4:14), decimals = 2) %>%
  tab_header(title = "Table 1: Descriptive Statistics for MR")


# regression equations for mass

SMRlm <- lm(log10(data$MLND) ~ log10(data$Weight_g))
MMRlm <- lm(log10(data$MR.abs) ~ log10(data$Weight_g))
ASlm <- lm(log10(data$ASabs) ~ log10(data$Weight_g))

model <- SMRlm

extract_lm_info <- function(model) {
  # Extract model summary
  model_summary <- summary(model)
  
  # Extract coefficients
  log_a <- model_summary$coefficients[1, 1]  # Intercept (loga)
  b <- model_summary$coefficients[2, 1]      # Slope (b)
  
  # Extract confidence intervals
  conf_int <- confint(model)  # 95% CI by default
  log_a_lower <- conf_int[1, 1]
  log_a_upper <- conf_int[1, 2]
  b_lower <- conf_int[2, 1]
  b_upper <- conf_int[2, 2]
  
  # Extract R-squared value
  r2 <- model_summary$r.squared
  
  # Extract p-value for the slope (b)
  p_value <- model_summary$coefficients[2, 4]
  
  # Create a table with the extracted values
  result_table <- data.frame(
    loga = log_a, 
    loga_lower = log_a_lower,
    loga_upper = log_a_upper,
    b = b, 
    b_lower = b_lower,
    b_upper = b_upper,
    R2 = r2, 
    P = p_value
  )
  
  return(result_table)
}

# Regression equation of the linear model of SMR and weight
regeq <- rbind(extract_lm_info(SMRlm),  # Apply the function to extract information
               extract_lm_info(MMRlm),  # Apply the function to extract information
               extract_lm_info(ASlm))  # Apply the function to extract information


####################
# predicted weight #
####################

#In some table, or in the text, can you add predicted means and CI for an
#average size fish for SMR at different temperatures? If you have the model as above, 
#where mass in the in the model rather than adjusted before, you can get predictions with predict(model) or ggpredict(), 
#but it doesn’t always  adjust for log-transformation automatically so you have to back-transform them.

library(ggeffects)
SMRmodel <- lmer(log10(MLND) ~ log10(Weight_g) + Temperature + (1|Batch), data = data1) #model
avg_weight <- mean(data1$Weight_g, na.rm = TRUE) # mean weight

newdata <- data.frame(
  Weight_g = avg_weight,
  Temperature = unique(data1$Temperature)) 

newdata$log10_pred <- predict(SMRmodel, newdata = newdata, re.form = NA)  # fixed effects only
newdata$SMR_pred <- 10^newdata$log10_pred #back transform log10
newdata$Batch <- data1$Batch[1]  # or any valid level from data1$Batch

library(merTools)

# Get prediction intervals on the log10 scale
preds_ci <- predictInterval(SMRmodel, newdata = newdata, level = 0.95, n.sims = 1000, stat = "mean", type = "linear.prediction", which = "fixed")

# Add and back-transform
newdata$log10_lower <- preds_ci$lwr
newdata$log10_upper <- preds_ci$upr
newdata$SMR_lower <- 10^newdata$log10_lower
newdata$SMR_upper <- 10^newdata$log10_upper
newdata[, c("Temperature", "SMR_pred", "SMR_lower", "SMR_upper")]


##################
# Big temp Graph #
##################

#Load Temperature data
T4 <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/T4.csv")
T11 <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/T11 F9.csv")
T12 <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/T12 F1.csv")
T13 <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/T13 F4.csv")
Teggs <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/Teggs.csv")
TAcc <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/TAcclimation.csv")
TBase <- read.csv("C:/Users/03243918/Valtion/ThermoEvo - General/Laukaa Summer 2023/Experiment_temperature_loggers/GraphFiles/TBase.csv")


#filter and format data so it all can be joined
Teggs$Date <- as.Date(Teggs$Date, format = "%d/%m/%Y")
tex <- filter(Teggs, Date > "2023-05-25")
Temps <- rbind(T11, T12, T13)
Temps$Date <- as.Date(Temps$Date, format = "%d/%m/%Y")
Temps$DateTime <- as.POSIXct(as.character(paste(Temps$Date, Temps$Time)), format="%Y-%m-%d %H:%M:%S")
Temps$Tank <- as.factor(Temps$Tank)
Temps <- subset(Temps,  "2023-08-17" > Date )
TAcc$Date <- as.Date(TAcc$Date, format = "%d/%m/%Y")
TAcc$DateTime <- as.POSIXct(as.character(paste(TAcc$Date, TAcc$Time)), format="%Y-%m-%d %H:%M:%S")
T4$Date <- as.Date(T4$Date, format = "%d/%m/%Y")
T4$DateTime <- as.POSIXct(as.character(paste(T4$Date, T4$Time)), format="%Y-%m-%d %H:%M:%S")
TAcc$Tank <- as.factor(TAcc$Tank)
TBase$Date <- as.Date(TBase$Date, format = "%d/%m/%Y")
TBase$DateTime <- as.POSIXct(as.character(paste(TBase$Date, TBase$Time)), format="%Y-%m-%d %H:%M:%S")
TBase$Tank <- as.factor(TBase$Tank)


TAcc <- subset(TAcc, Date > "2023-07-27")
TAcc <- subset(TAcc,  "2023-08-01" > Date )
TAcc <- rbind(TAcc, T4, TBase)
TAcc$Experiment <-  "Control"
df <- as_tibble(seq(1, 19, 1))
df$Date <- seq.Date(as.Date("2023-07-27", format="%Y-%m-%d"),as.Date("2023-08-14", format="%Y-%m-%d"),1)
df$Time <- (rep("00:01:00", 19))
df$DateTime <- as.POSIXct(as.character(paste(df$Date, df$Time)), format="%Y-%m-%d %H:%M:%S")
df$Experiment <- "Control"
df <- df[1:19,c(2:5)]
df$Tank <- "Control"
df$Experiment <- "Control"
df$Temp <- "NA"
df <- df[, c(1,2,6, 5, 3, 4)]
df <- rbind(df, TAcc)
df$movTemp <- 0

Temps$Experiment <- Temps$Tank
levels(Temps$Experiment) <- c('Heatwave Tank 1', 'Heatwave Tank 2', 'Heatwave Tank 3', "Control", "Control")
Temps <- plyr::rbind.fill(Temps, df)
Temps <- Temps %>% arrange(DateTime)
Temps <- subset(Temps, Date > "2023-07-27")
Temps <- subset(Temps,  "2023-08-22" > Date )
Temps1 <- Temps
Temps <- Temps1
Temps$Temp <- as.numeric(Temps$Temp)
Temps$Tank <- as.factor(Temps$Tank)

library(zoo)

#smooth control data and filter right parts
Temps <- Temps %>% 
  group_by(Tank) %>% 
  arrange(DateTime) %>%
  mutate(movTemp = ifelse(Experiment != "Control", 
                          zoo::rollmean(Temp, k=40, fill = NA, align = "right"),
                          Temp))

Temps$Experiment <- as.factor(Temps$Experiment)
Temps <- Temps %>% filter(!(Experiment == "Heatwave Tank 1"& "2023-08-18" > Date & Date > "2023-08-10"  & 21.6 > movTemp))


#Plot here
bigtemp <- ggplot(subset(Temps, 22.1 > movTemp | movTemp > 15), 
       aes(x = DateTime, 
           y = Temp, 
           colour = Experiment, 
           fill = Experiment)) +
  geom_point(alpha = 0.05,na.rm=F, size = 3) +
  geom_line(data = subset(Temps, Experiment != "Control" ),aes(x = DateTime, 
                 y = movTemp),linewidth = 1) +
  geom_line(data = subset(Temps, Tank == "Acc" ),aes(x = DateTime, 
                                                                y = zoo::rollmean(Temp, k=40, fill = NA, align = "right")),linewidth = 1) +
  geom_line(data = subset(Temps, Tank == "T4" ),aes(x = DateTime, 
                                                     y = zoo::rollmean(Temp, k=40, fill = NA, align = "right")),linewidth = 1) +
# Add linear regression line
labs(x = "Date",
     y = "Temperature (\u00B0C)") +
  scale_fill_manual( values = c("#6c9bff", "#cca7db", "#f88bce", "#ffcb94")) +
  scale_colour_manual( values = c("#0033A0", "#7F3F98", "#CC0D82", "#FF8200")) +
  ggpubr::theme_classic2()+
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 18)
  ) +
  labs(color = "Tank", fill = "Tank") +
  geom_segment(aes(x=(as.POSIXct("2023-08-02 09:00:00")), y = 16, yend = 22.8), colour = "#0033A0", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-03 09:00:00")), y = 16, yend = 22.8), colour = "#0033A0", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-09 09:00:00")), y = 19, yend = 22.8), colour = "#54585A", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-10 09:00:00")), y = 19, yend = 22.8), colour = "#54585A", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-16 09:00:00")), y = 22, yend = 22.8), colour = "#54585A", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-17 09:00:00")), y = 22, yend = 22.8), colour = "#54585A", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-21 09:00:00")), y = 16, yend = 22.8), colour = "#0033A0", linetype = 2, linewidth=1) + 
  geom_segment(aes(x=(as.POSIXct("2023-08-22 09:00:00")), y = 16, yend = 22.8), colour = "#0033A0", linetype = 2, linewidth=1) + 
  annotate("text", x=as.POSIXct("2023-08-02 09:00:00"), y=23.5, label="B1", colour = "#0033A0")+
  annotate("text", x=as.POSIXct("2023-08-03 09:00:00"), y=23.5, label="B2", colour = "#0033A0")+
  annotate("text", x=as.POSIXct("2023-08-09 09:00:00"), y=23.5, label="B3")+
  annotate("text", x=as.POSIXct("2023-08-10 09:00:00"), y=23.5, label="B4")+
  annotate("text", x=as.POSIXct("2023-08-16 09:00:00"), y=23.5, label="B5")+
  annotate("text", x=as.POSIXct("2023-08-17 09:00:00"), y=23.5, label="B6")+
  annotate("text", x=as.POSIXct("2023-08-21 09:00:00"), y=23.5, label="B7", colour = "#0033A0")+
  annotate("text", x=as.POSIXct("2023-08-22 09:00:00"), y=23.5, label="B8", colour = "#0033A0")+
  annotate("text", x=as.POSIXct("2023-08-02 21:00:00"), y=24.5, label="16(\u00B0C)", colour = "#0033A0")+
  annotate("text", x=as.POSIXct("2023-08-09 21:00:00"), y=24.5, label="19(\u00B0C)")+
  annotate("text", x=as.POSIXct("2023-08-16 21:00:00"), y=24.5, label="22(\u00B0C)")+
  annotate("text", x=as.POSIXct("2023-08-21 21:00:00"), y=24.5, label="16(\u00B0C)", colour = "#0033A0")

ggsave("allTemps.png", plot = bigtemp, width = 16, height = 8, units = "in", dpi = 300)






#################
#   CORR TABLE  #
#################
library(psych)

# Assuming your data has a column for MR values (MR.abs) and a Temperature factor
data16 <- subset(data1, Temperature == 16)
data16 <- data16[, c("MLND", "MR.abs", "ASabs", "Weight_g", "SMRmassadj", "MMRmassadj","ASmass")]
# Custom names for the variables
custom_names <- c("SMR","MMR","AS",  "Mass")
vars <- c("MLND", "MR.abs", "ASabs", "Weight_g", "SMRmassadj", "MMRmassadj","ASmass")
cor_results <- data.frame(var1 = character(),
                          var2 = character(),
                          r = numeric(),
                          p = numeric(),
                          temp = numeric(),
                          stringsAsFactors = FALSE)

for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    test <- cor.test(data16[[vars[i]]], data16[[vars[j]]], use = "complete.obs")
    cor_results <- rbind(cor_results, data.frame(
      var1 = vars[i],
      var2 = vars[j],
      r = round(test$estimate, 2),
      p = test$p.value,
      temp = 16
    ))
  }
}

cor16 <- cor_results
cor_results1 <- data.frame(var1 = character(),
                          var2 = character(),
                          r = numeric(),
                          p = numeric(),
                          temp = numeric(),
                          stringsAsFactors = FALSE)
data19 <- subset(data1, Temperature == 19)
data19 <- data19[, c("MLND", "MR.abs", "ASabs", "Weight_g", "SMRmassadj", "MMRmassadj","ASmass")]
rm(i)
rm(j)
for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    test <- cor.test(data19[[vars[i]]], data19[[vars[j]]], use = "complete.obs")
    cor_results1 <- rbind(cor_results1, data.frame(
      var1 = vars[i],
      var2 = vars[j],
      r = round(test$estimate, 2),
      p = test$p.value,
      temp = 19
    ))
  }
}

cor19 <- cor_results1


data22 <- subset(data1, Temperature == 22)
data22 <- data22[, c("MLND", "MR.abs", "ASabs", "Weight_g", "SMRmassadj", "MMRmassadj","ASmass")]
rm(i)
rm(j)
cor_results2 <- data.frame(var1 = character(),
                          var2 = character(),
                          r = numeric(),
                          p = numeric(),
                          temp = numeric(),
                          stringsAsFactors = FALSE)
for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    test <- cor.test(data22[[vars[i]]], data22[[vars[j]]], use = "complete.obs")
    cor_results2 <- rbind(cor_results2, data.frame(
      var1 = vars[i],
      var2 = vars[j],
      r = round(test$estimate, 2),
      p = test$p.value,
      temp = 22
    ))
  }
}

cor22 <- cor_results2

cordata <- rbind(cor16, cor19, cor22)
#### this table in the paper
