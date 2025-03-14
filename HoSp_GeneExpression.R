#Data analysis House Sparrow Gene Expression

#Load Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotrix)

#for circadian rhythm 
library(circacompare)


#Load working directory 
#Click on more (bottom right panel)
#Click "Set As Working Directory"
setwd("")


#Load data
hs <- read.csv("")

#look at data
str(hs)

#separate by gene
bmal <- subset(hs, Gene == "BMAL")
cry = subset(hs, Gene == "CRY1")
per2 <- subset(hs, Gene == "PER2")
per3 <- subset(hs, Gene == "PER3")


#view Data
# Plotting
ggplot(bmal, aes(x = ZT, y = Expression, group = site, color = site)) + 
  geom_point() + # Add points
  theme_classic() + # Use a minimal theme
  labs(title = "Expression by ZT and Site", # Add a title
       x = "Zeitgeber Time (ZT)", # Label the x-axis
       y = "Expression") + # Label the y-axis
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) # Customize colors

ggplot(cry, aes(x = ZT, y = Expression, group = site, color = site)) + 
  geom_point() + # Add points
  theme_classic() + # Use a minimal theme
  labs(title = "Expression by ZT and Site", # Add a title
       x = "Zeitgeber Time (ZT)", # Label the x-axis
       y = "Expression") + # Label the y-axis
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) # Customize colors

ggplot(per2, aes(x = ZT, y = Expression, group = site, color = site)) + 
  geom_point() + # Add points
  theme_classic() + # Use a minimal theme
  labs(title = "Expression by ZT and Site", # Add a title
       x = "Zeitgeber Time (ZT)", # Label the x-axis
       y = "Expression") + # Label the y-axis
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) # Customize colors

ggplot(per3, aes(x = ZT, y = Expression, group = site, color = site)) + 
  geom_point() + # Add points
  theme_classic() + # Use a minimal theme
  labs(title = "Expression by ZT and Site", # Add a title
       x = "Zeitgeber Time (ZT)", # Label the x-axis
       y = "Expression") + # Label the y-axis
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) # Customize colors

#Get stats
bmal_sum <- bmal %>%
  group_by(site, ZT) %>%
  summarize(Mean = mean(Expression, na.rm=TRUE), 
            SEM=std.error(Expression, na.rm=TRUE), n = n())


# Summarize for cry
cry_sum <- cry %>%
  group_by(site, ZT) %>%
  summarize(Mean = mean(Expression, na.rm = TRUE), 
            SEM = std.error(Expression, na.rm = TRUE), n = n())

# Summarize for per2
per2_sum <- per2 %>%
  group_by(site, ZT) %>%
  summarize(Mean = mean(Expression, na.rm = TRUE), 
            SEM = std.error(Expression, na.rm = TRUE), n = n())

# Summarize for per3
per3_sum <- per3 %>%
  group_by(site, ZT) %>%
  summarize(Mean = mean(Expression, na.rm = TRUE), 
            SEM = std.error(Expression, na.rm = TRUE), n = n())

###
#view again
###

# Function to plot each gene
plot_gene_expression <- function(data, gene_name) {
  ggplot(data, aes(x = ZT, y = Mean, color = site, group = site)) +
    geom_line() + # Add lines to connect points
    geom_point() + # Add points for Mean
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) + # Add error bars for SE
    scale_color_brewer(palette = "Set1") + # Use a color palette
    theme_minimal() + # Use a minimal theme for aesthetics
    labs(title = paste("Expression of", gene_name, "by ZT and Site"),
         x = "Zeitgeber Time (ZT)",
         y = "Mean Expression ± SE",
         color = "Site") + 
    theme(legend.position = "bottom") # Adjust legend position
}

# Plot for bmal
bmal_plot <- plot_gene_expression(bmal_sum, "BMAL")
print(bmal_plot)


# Plot for cry
cry_plot <- plot_gene_expression(cry_sum, "CRY")
print(cry_plot)

# Plot for per2
per2_plot <- plot_gene_expression(per2_sum, "PER2")
print(per2_plot)

# Plot for per3
per3_plot <- plot_gene_expression(per3_sum, "PER3")
print(per3_plot)



#Rhythmic analysis for BMAL1

b_greens <- subset(bmal, site != "campus")

res <- circacompare(x=b_greens, col_time="ZT", col_group = "site", col_outcome = "Expression")
res



b_campus <- subset(bmal, site != "greens")

res <- circacompare(x=b_campus, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


b_farm <- subset(bmal, site != "farm")

res <- circacompare(x=b_farm, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


############

#BMAL1
#Anova for each time point
tp = 13
dat = subset(bmal, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 1
dat = subset(bmal, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 7
dat = subset(bmal, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 19
dat = subset(bmal, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


######################


#Rhythmic Analysis for cry1


c_greens <- subset(cry, site != "campus")

res <- circacompare(x=c_greens, col_time="ZT", col_group = "site", col_outcome = "Expression")
res

#alpha_threshold needs to be higher than p=0.18 to work

c_greens <- subset(cry, site != "campus")

res <- circacompare(x=c_greens, col_time="ZT", col_group = "site", col_outcome = "Expression", alpha_threshold = 0.2)
res


c_campus <- subset(cry, site != "greens")

res <- circacompare(x=c_campus, col_time="ZT", col_group = "site", col_outcome = "Expression")
res

c_farm <- subset(cry, site != "farm")

res <- circacompare(x=c_farm, col_time="ZT", col_group = "site", col_outcome = "Expression")
res



#cry1
#ANOVA for each time point

tp = 1
dat = subset(cry, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 7
dat = subset(cry, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

tp = 13
dat = subset(cry, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

tp = 19
dat = subset(cry, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


############################

#Rhythmic Analysis for per2

p_greens <- subset(per2, site != "campus")

res <- circacompare(x=p_greens, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


p_campus <- subset(per2, site != "greens")

res <- circacompare(x=p_campus, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


p_farm <- subset(per2, site != "farm")

res <- circacompare(x=c_farm, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


#per2
#ANOVA for each time point

tp = 1
dat = subset(per2, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

tp = 7
dat = subset(per2, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 13
dat = subset(per2, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 19
dat = subset(per2, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

#########################

#Rhythmic Analysis of per3

p3_greens <- subset(per3, site != "campus")

res <- circacompare(x=p_greens, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


p3_campus <- subset(per3, site != "greens")

res <- circacompare(x=p_campus, col_time="ZT", col_group = "site", col_outcome = "Expression")
res


p3_farm <- subset(per3, site != "farm")

res <- circacompare(x=p3_farm, col_time="ZT", col_group = "site", col_outcome = "Expression")
res

#per3
#ANOVA for each time point

tp = 1
dat = subset(per3, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

tp = 7
dat = subset(per3, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

tp = 13
dat = subset(per3, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc


tp = 19
dat = subset(per3, ZT == tp)
rest = aov(data = dat, Expression ~ site)
summary(rest)
posthoc <- TukeyHSD(rest)
posthoc

#####################################
#Plotting figures 

Fig1 <- ggplot(data = bmal, aes(x = ZT, y = Expression)) +
  geom_rect(aes(xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf), fill ="grey", alpha = .7)+
  geom_point(position = position_jitterdodge(), aes(shape = site, color = site), alpha = 0.5, size = 2) + # Individual points
  geom_errorbar(data = bmal_sum, aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, size = 1) + # Error bars
  geom_line(data = bmal_sum, aes(y = Mean, group = site, color = site), size = 1) + # Connect means with lines
  geom_point(data = bmal_sum, aes(y = Mean, fill = site, shape = site, color = site), size = 3) + # Mean points with filled color and black outline
  scale_fill_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5f5953")) +
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) + # Outline color for treatment
  scale_shape_manual(values = c("campus" = 17, "greens" = 16, "farm" = 15)) + # Shapes for different sites
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  theme_classic() +
  labs(x = "ZT",
       y = expression(italic(bmal)~"Expression")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=15, family= "Times", color = "black"),
    axis.title.x = element_text( size=17, family= "Times", vjust=-0.2),
    axis.title.y = element_text( size=17, family= "Times", vjust=2),
    axis.text.y  = element_text( size=15, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), # Adjusted to show legend text
    legend.position = "right" # Adjusted to show the legend
  )
Fig1


Fig2 <- ggplot(data = cry, aes(x = ZT, y = Expression)) +
  geom_rect(aes(xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf), fill ="grey", alpha = .7)+
  geom_point(position = position_jitterdodge(), aes(shape = site, color = site), alpha = 0.5, size = 2) + # Individual points
  geom_errorbar(data = cry_sum, aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, size = 1) + # Error bars
  geom_line(data = cry_sum, aes(y = Mean, group = site, color = site), size = 1) + # Connect means with lines
  geom_point(data = cry_sum, aes(y = Mean, fill = site, shape = site, color = site), size = 3) + # Mean points with filled color and black outline
  scale_fill_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5f5953")) +
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) + # Outline color for treatment
  scale_shape_manual(values = c("campus" = 17, "greens" = 16, "farm" = 15)) + # Shapes for different sites
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  theme_classic() +
  labs(x = "ZT",
       y = expression(italic(cry1)~"Expression")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=15, family= "Times", color = "black"),
    axis.title.x = element_text( size=17, family= "Times", vjust=-0.2),
    axis.title.y = element_text( size=17, family= "Times", vjust=2),
    axis.text.y  = element_text( size=15, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), # Adjusted to show legend text
    legend.position = "right" # Adjusted to show the legend
  )
Fig2



Fig3 <- ggplot(data = per2, aes(x = ZT, y = Expression)) +
  geom_rect(aes(xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf), fill ="grey", alpha = .7)+
  geom_point(position = position_jitterdodge(), aes(shape = site, color = site), alpha = 0.5, size = 2) + # Individual points
  geom_errorbar(data = per2_sum, aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, size = 1) + # Error bars
  geom_line(data = per2_sum, aes(y = Mean, group = site, color = site), size = 1) + # Connect means with lines
  geom_point(data = per2_sum, aes(y = Mean, fill = site, shape = site, color = site), size = 3) + # Mean points with filled color and black outline
  scale_fill_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5f5953")) +
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) + # Outline color for treatment
  scale_shape_manual(values = c("campus" = 17, "greens" = 16, "farm" = 15)) + # Shapes for different sites
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  theme_classic() +
  labs(x = "ZT",
       y = expression(italic(per2)~"Expression")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=15, family= "Times", color = "black"),
    axis.title.x = element_text( size=17, family= "Times", vjust=-0.2),
    axis.title.y = element_text( size=17, family= "Times", vjust=2),
    axis.text.y  = element_text( size=15, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), # Adjusted to show legend text
    legend.position = "right" # Adjusted to show the legend
  )
Fig3


Fig4 <- ggplot(data = per3, aes(x = ZT, y = Expression)) +
  geom_rect(aes(xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf), fill ="grey", alpha = .7)+
  geom_point(position = position_jitterdodge(), aes(shape = site, color = site), alpha = 0.5, size = 2) + # Individual points
  geom_errorbar(data = per3_sum, aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, size = 1) + # Error bars
  geom_line(data = per3_sum, aes(y = Mean, group = site, color = site), size = 1) + # Connect means with lines
  geom_point(data = per3_sum, aes(y = Mean, fill = site, shape = site, color = site), size = 3) + # Mean points with filled color and black outline
  scale_fill_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5f5953")) +
  scale_color_manual(values = c("campus" = "#F2EE05", "greens" = "#1A25FF", "farm" = "#5F5953")) + # Outline color for treatment
  scale_shape_manual(values = c("campus" = 17, "greens" = 16, "farm" = 15)) + # Shapes for different sites
  scale_x_continuous(breaks=c(1,7,13,19), limits=c(0,24))+
  theme_classic() +
  labs(x = "ZT",
       y = expression(italic(per3)~"Expression")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=15, family= "Times", color = "black"),
    axis.title.x = element_text( size=17, family= "Times", vjust=-0.2),
    axis.title.y = element_text( size=17, family= "Times", vjust=2),
    axis.text.y  = element_text( size=15, family= "Times"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), # Adjusted to show legend text
    legend.position = "right" # Adjusted to show the legend
  )
Fig4


