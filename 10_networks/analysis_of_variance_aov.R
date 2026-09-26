# turn off scientific notation
options(scipen=999)

# install packages
#install.packages("ggpubr")
#install.packages("multcomp")

# load libraries
library(ggpubr)
library(multcomp)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[11], plotColors[6], plotColors[4], plotColors[5])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NA_22Sep2026/unsigned_power2/analysis_of_variance"
dir.create(workingDir)
setwd(workingDir)

# import expression data
inputTable <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NA_22Sep2026/unsigned_power2/eigengeneExpression.csv", row.names="gene")

# transpose expression data
inputTable <- as.data.frame(t(inputTable))

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# setup data frame
expData <- merge(targets, inputTable, by = 'row.names') 

# convert mouseline to a factor
expData$mouseline <- factor(expData$mouseline,
                            levels = c("CAST", "SARA", "MANF"),
                            labels = c("CAST", "SARA", "MANF"))

# convert treatment to a factor
expData$treatment <- factor(expData$treatment,
                            levels = c("Ad_lib", "Food_Restriction"),
                            labels = c("Ad_lib", "Food_Restriction"))

# convert tissue to a factor
expData$tissue <- factor(expData$tissue,
                           levels = c("BAT", "WAT", "HYP", "LIV"),
                           labels = c("BAT", "WAT", "HYP", "LIV"))

# convert sex to a factor
expData$sex <- factor(expData$sex,
                              levels = c("Male", "Female"),
                              labels = c("Male", "Female"))

# set row names
rownames(expData) <- expData$Row.names

# remove the Row.names column
#expData <- expData[,2:6]
expData <- expData[,-1]

# remove the individual column
expData <- expData[,-1]

# check the structure of the expression data
#str(expData)

# make frequency tables
#table(expData$treatment, expData$genotype, expData$tolerance, expData$replicate)

# loop over each module
expData_updated <- expData
for (curMod in 5:ncol(expData)) {
  # test
  #curMod <- 5
  
  # retrieve module column name
  modName <- colnames(expData_updated)[curMod]
  
  # update module column name
  colnames(expData_updated)[curMod] <- "expression"
  
  # compute two-way anova
  affy.aov <- aov(expression ~ treatment + sex + tissue + tissue * mouseline + treatment * mouseline + treatment * tissue, data = expData_updated)
  
  # output summary statistics
  sumAffy <- summary(affy.aov)
  
  # write summary statistics to a file
  exportFile <- paste(modName, "ANOVA_summary.csv", sep="_")
  capture.output(sumAffy, file=exportFile)
  
  # create a colored box plot
  exportFile <- paste(modName, "mouseline_tissue_coloredBoxPlot.png", sep="_")
  png(file=exportFile)
  ggboxplot(data=expData_updated, x="mouseline", y="expression", color="tissue",
            palette = plotColorSubset)
  dev.off()
  
  # create a colored box plot
  exportFile <- paste(modName, "treatment_tissue_coloredBoxPlot.png", sep="_")
  png(file=exportFile)
  ggboxplot(data=expData_updated, x="treatment", y="expression", color="tissue",
            palette = plotColorSubset)
  dev.off()
  
  # create a colored box plot
  exportFile <- paste(modName, "treatment_mouseline_coloredBoxPlot.png", sep="_")
  png(file=exportFile)
  ggboxplot(data=expData_updated, x="treatment", y="expression", color="mouseline",
            palette = plotColorSubset)
  dev.off()
  
  # box plot with two variable factors
  #exportFile <- paste(modName, "twoVariableFactorBoxPlot.png", sep="_")
  #png(file=exportFile)
  #boxplot(expression ~ treatment + sex + tissue + tissue * mouseline + treatment * mouseline + treatment * tissue, data=expData_updated, frame=FALSE,
  #        col = plotColorSubset)
  #dev.off()
  
  # two way interaction plot
  exportFile <- paste(modName, "twoWayInteractionPlot.png", sep="_")
  png(file=exportFile)
  interaction.plot(x.factor = expData_updated$treatment, trace.factor = expData_updated$tissue,
                   response = expData_updated$expression, fun = mean,
                   type = "b", legend = TRUE,
                   xlab = "Treatment", ylab="Expression",
                   pch=c(1,19), col = plotColorSubset)
  dev.off()
  
  
  ## check the validity of ANOVA assumptions
  # The data must be regularly distributed, and the variation between groups must be homogeneous
  
  ## examine the assumption of homogeneity of variance
  # the residuals versus fits graphic are used to assess for variance homogeneity
  
  # plot homogeneity of variances
  exportFile <- paste(modName, "homogeneityPlot.png", sep="_")
  png(file=exportFile)
  plot(affy.aov, 1)
  dev.off()
  
  # examine the assumption of normality
  # in a residuals’ normality plot the residuals quantiles are displayed against the normal distribution quantiles
  # the residuals’ normal probability plot is used to confirm that the residuals are normally distributed
  # the residuals’ normal probability plot should roughly follow a straight line
  
  # normality plot
  exportFile <- paste(modName, "normalityPlot.png", sep="_")
  png(file=exportFile)
  plot(affy.aov, 2)
  dev.off()
  
  # extract the residuals
  aovRes <- residuals(object = affy.aov)
  
  # run Shapiro-Wilk test
  swTest <- shapiro.test(x = aovRes)
  
  # write test statistics to a file
  exportFile <- paste(modName, "shapiroTest_summary.csv", sep="_")
  capture.output(swTest, file=exportFile)
  
  # update module column name
  colnames(expData_updated)[curMod] <- modName
}

# for i in /Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NA_22Sep2026/unsigned_power2/analysis_of_variance/*_ANOVA_summary.csv; do echo $i; cat $i | tail -n+2 | head -7 | tr -s " " | sed "s/< /</g" | tr " " "," | cut -d"," -f1,6- | sed "s/,$//g"; done > /Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NA_22Sep2026/unsigned_power2/analysis_of_variance/ANOVA_summary.csv
