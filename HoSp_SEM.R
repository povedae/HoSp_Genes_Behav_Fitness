## PATH ANALYSIS ##
library(tidyr)
library(lavaan)
library(plotrix)

SEM_df<-read.csv("~/HoSp Data/prepared datasets/MEGA_df.csv",header=TRUE,na.strings=c("","NA"))

  SEM_df$farm <-as.numeric(SEM_df$farm)
  SEM_df$campus <-as.numeric(SEM_df$campus)
  SEM_df$greens <-as.numeric(SEM_df$greens)
  SEM_df$onset_rel <- as.numeric(SEM_df$onset_rel)
  SEM_df$offset_rel <- as.numeric(SEM_df$offset_rel)
  SEM_df$avgmassyoung <- as.numeric(SEM_df$avgmassyoung)
  SEM_df$indv_feedrate <- as.numeric(SEM_df$indv_feedrate)
  
  SEM_df$onset_rel<-SEM_df$onset_rel/60
  SEM_df$offset_rel<-SEM_df$offset_rel/60
  
  
SEM_df2<-subset(SEM_df,select=c(-tagID,-sex,-male,-smi,-site,-lux,-nest,-nest_feedrate))


CVridge_ons <- cv.glmnet(as.matrix(SEM_df2[, c("farm","campus","greens")]),SEM_df2$onset_rel,alpha=0)
CVlambda_ons <- CVridge$lambda.min

ridgemodel_ons <- glmnet(as.matrix(SEM_df2[, c("farm","campus","greens")]),SEM_df2$onset_rel,alpha=0)

ridgecoef_ons <- coef(ridgemodel_ons, s=CVlambda_ons)


SEM_df2$onset_rel<-scale(SEM_df2$onset_rel,center=TRUE,scale=TRUE)
SEM_df2$farm<-scale(SEM_df2$farm,center=TRUE,scale=TRUE)
SEM_df2$campus<-scale(SEM_df2$campus,center=TRUE,scale=TRUE)
SEM_df2$greens<-scale(SEM_df2$greens,center=TRUE,scale=TRUE)

SEM_df2$avgmassyoung<-scale(SEM_df2$avgmassyoung,center=TRUE,scale=TRUE)

ALANmodel<- "
    activity =~ onset_rel + indvfeedrate 
    activity ~ farm + campus + greens
    
"

ALAN_SEM <- sem(ALANmodel,data=SEM_df2)
cor(SEM_df2[,c(1:3,7)], use="pairwise", method="spearman")^2


semPower.getDF(ALANmodel)

semPower::semPower.getDf(ALANmodel)



# one call could do all of this but my attempts were wrong so I'm using an inefficient method
SEM_df1$onset_rel<-scale(SEM_df1$onset_rel,center=TRUE,scale=TRUE)
SEM_df1$offset_rel<-scale(SEM_df1$offset_rel,center=TRUE,scale=TRUE)
SEM_df1$noctact<-scale(SEM_df1$noctact,center=TRUE,scale=TRUE)
SEM_df1$avgmassyoung<-scale(SEM_df1$avgmassyoung,center=TRUE,scale=TRUE)
SEM_df1$indv_feedrate<-scale(SEM_df1$indv_feedrate,center=TRUE,scale=TRUE)



SEM_df3<-subset(SEM_df2,select=c(-offset_rel,-noctact,-avgmassyoung,-indv_feedrate))
cor(SEM_df3, use="pairwise", method="spearman")^2


tot_model1<-lm(onset_rel~site+sex+smi+offset_rel+noctact+avgmassyoung+indv_feedrate,data=SEM_df)
cor()

tot_model2<-lm(offset_rel~site+sex+smi+noctact+avgmassyoung+indv_feedrate,data=SEM_df)


tot_model3<-lm(noctact~site+sex+smi+avgmassyoung+indv_feedrate,data=SEM_df)


tot_model4<-lm(avgmassyoung~site+sex+smi+indv_feedrate,data=SEM_df)


tot_model5<-lm(indv_feedrate~site+sex+smi,data=SEM_df)


SEM_df1<-subset(SEM_df,select=c(-tagID,-sex,-male,-smi,-lux,-nest,-nest_feedrate))
SEM_df1$site<-as.factor(SEM_df1$site)

cor(SEM_df1[,5:9], use="pairwise", method="spearman")^2


ALANmodel<- " 

    avgmassyoung ~ onset_rel + offset_rel + noctact
    avgmassyoung ~ indv_feedrate
    indv_feedrate ~ onset_rel + offset_rel + noctact

  "

ALAN_SEM<-sem(ALANmodel,data=SEM_df1)

ALAN_SEM<-cfa(ALANmodel,data=SEM_df1, group="site")
summary(ALAN_SEM, standardized=T, fit.measures=T, rsq=T)


library(semPlot)
semPaths(ALAN_SEM,'std',layout='tree',edge.label.cex = 1.5,label.cex=1.5)





ALANmodel2<- " 

    avgmassyoung ~ noctact
    avgmassyoung ~ indv_feedrate
    indv_feedrate ~ noctact
    
"

ALAN_SEM2<-cfa(ALANmodel2,data=SEM_df1, group="site")
summary(ALAN_SEM2, standardized=T, fit.measures=T, rsq=T)


library(semPlot)
semPaths(ALAN_SEM2,'std',layout='tree',edge.label.cex = 1.5,label.cex=1.5)

ggplot ()



ALANmodel3<- " 

    avgmassyoung ~ indv_feedrate + farm + campus + greens

"
ALAN_SEM3<-cfa(ALANmodel3,data=SEM_df1, group="site")
summary(ALAN_SEM3, standardized=T, fit.measures=T, rsq=T)
