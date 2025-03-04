## MODELS ##
library(dplyr)
library(plotrix)
library(lme4)
library(rstatix)
library(emmeans)
library(nlme)
library(tidyr)

# ONSET -------------------------------------------------------------------

# ONSET: SEX DIFFS AT LOW ALAN
onset_low<-onsetoffset %>% filter(site=='farm')  #isolate data from one site
ONSETlow_summ<-group_by(onset_low,sex) %>% summarize(meanmin=mean(avg_onset_rel/60,na.rm=TRUE), 
  sd=sd(avg_onset_rel/60,na.rm=TRUE),se=std.error(avg_onset_rel/60,na.rm=TRUE))  #summarize
ONSETlow_summ<-as.data.frame(ONSETlow_summ);      ONSETlow_summ  #view

wilcox.test(avg_onset_rel~sex,data=onset_low)

# ONSET: SEX DIFFS AT MED ALAN
onset_med<-onsetoffset %>% filter(site=='campus')
ONSETmed_summ<-group_by(onset_med,sex) %>% summarize(meanmin=mean(avg_onset_rel/60,na.rm=TRUE),
  sd=sd(avg_onset_rel/60,na.rm=TRUE),se=std.error(avg_onset_rel/60,na.rm=TRUE))
ONSETmed_summ<-as.data.frame(ONSETmed_summ);      ONSETmed_summ

wilcox.test(avg_onset_rel~sex,data=onset_med)

# ONSET: SEX DIFFS AT HIGH ALAN
onset_high<-onsetoffset %>% filter(site=='greens')
ONSEThigh_summ<-group_by(onset_high,sex) %>% summarize(meanmin=mean(avg_onset_rel/60,na.rm=TRUE),
  sd=sd(avg_onset_rel/60,na.rm=TRUE),se=std.error(avg_onset_rel/60,na.rm=TRUE))
ONSEThigh_summ<-as.data.frame(ONSEThigh_summ);      ONSEThigh_summ

wilcox.test(avg_onset_rel~sex,data=onset_high)

# ONSET: smi DIFFS AT LOW ALAN
cor.test(x=onset_low$avg_onset_rel,y=onset_low$smi,method='kendall')

# ONSET: smi DIFFS AT MED ALAN
cor.test(x=onset_med$avg_onset_rel,y=onset_med$smi,method='kendall')

# ONSET: smi DIFFS AT HIGH ALAN
cor.test(x=onset_high$avg_onset_rel,y=onset_high$smi,method='kendall')

# *** ONSET ACROSS SITES
ONSET_SITEsumm<-group_by(onsetoffset,site) %>% summarize(meanmin=mean(avg_onset_rel,na.rm=TRUE),
  sd=sd(avg_onset_rel,na.rm=TRUE),se=std.error(avg_onset_rel,na.rm=TRUE)) #summary with mean, sd, se
ONSET_SITEsumm<-as.data.frame(ONSET_SITEsumm);  ONSET_SITEsumm  #convert to df and view

KTonset<-kruskal.test(avg_onset_rel~site,data=onsetoffset); KTonset  #Kruskal-Wallis rank sum test
as.data.frame(dunn_test(onsetoffset,avg_onset_rel~site,p.adjust.method="holm"))  #Dunn test

# OFFSET ------------------------------------------------------------------

# OFFSET: SEX DIFFS AT LOW ALAN
offset_low<-onsetoffset %>% filter(site=='farm')
OFFSETlow_summ<-group_by(offset_low,sex) %>% summarize(meanmin=mean(avg_offset_rel/60,na.rm=TRUE),
  sd=sd(avg_offset_rel/60,na.rm=TRUE),se=std.error(avg_offset_rel/60,na.rm=TRUE)) 
OFFSETlow_summ<-as.data.frame(OFFSETlow_summ);      OFFSETlow_summ  

wilcox.test(avg_offset_rel~sex,data=offset_low)

# OFFSET: SEX DIFFS AT MED ALAN
offset_med<-onsetoffset %>% filter(site=='campus')
OFFSETmed_summ<-group_by(offset_med,sex) %>% summarize(meanmin=mean(avg_offset_rel/60,na.rm=TRUE),
  sd=sd(avg_offset_rel/60,na.rm=TRUE),se=std.error(avg_offset_rel/60,na.rm=TRUE))
OFFSETmed_summ<-as.data.frame(OFFSETmed_summ);      OFFSETmed_summ

wilcox.test(avg_offset_rel~sex,data=offset_med)

# OFFSET: SEX DIFFS AT HIGH ALAN
offset_high<-onsetoffset %>% filter(site=='greens')
OFFSEThigh_summ<-group_by(offset_high,sex) %>% summarize(meanmin=mean(avg_offset_rel/60,na.rm=TRUE),
  sd=sd(avg_offset_rel/60,na.rm=TRUE),se=std.error(avg_offset_rel/60,na.rm=TRUE))
OFFSEThigh_summ<-as.data.frame(OFFSEThigh_summ);      OFFSEThigh_summ

wilcox.test(avg_offset_rel~sex,data=offset_high)

# OFFSET: smi DIFFS AT LOW ALAN
cor.test(x=offset_low$avg_offset_rel,y=offset_low$smi,method='kendall')

# OFFSET: smi DIFFS AT MED ALAN
cor.test(x=offset_med$avg_offset_rel,y=offset_med$smi,method='kendall')

# OFFSET: smi DIFFS AT HIGH ALAN
cor.test(x=offset_high$avg_offset_rel,y=offset_high$smi,method='kendall')


# *** offset across sites
OFFSET_SITEsumm<-group_by(onsetoffset,site) %>% summarize(meanmin=mean(avg_offset_rel,na.rm=TRUE),
  sd=sd(avg_offset_rel,na.rm=TRUE),se=std.error(avg_offset_rel,na.rm=TRUE)) 
OFFSET_SITEsumm<-as.data.frame(OFFSET_SITEsumm);  OFFSET_SITEsumm 

kruskal.test(avg_offset_rel~site,data=onsetoffset)
as.data.frame(dunn_test(onsetoffset,avg_offset_rel~site,p.adjust.method="holm"))  


# NOCTURNAL ACTIVITY ------------------------------------------------------

# NOCTURNAL ACTIVITY: SEX DIFFS AT LOW ALAN
noctact_low<-nightprop %>% filter(site=='farm')
NOCTACTlow_summ<-as.data.frame(group_by(noctact_low,sex) %>% summarize(meanprop=mean(avg_prop,na.rm=TRUE),
  sd=sd(avg_prop,na.rm=TRUE),se=std.error(avg_prop,na.rm=TRUE))); NOCTACTlow_summ

wilcox.test(avg_prop~sex,data=noctact_low)

# NOCTURNAL ACTIVITY: SEX DIFFS AT MED ALAN
noctact_med<-nightprop %>% filter(site=='campus')
NOCTACTmed_summ<-as.data.frame(group_by(noctact_med,sex) %>% summarize(meanprop=mean(avg_prop,na.rm=TRUE),
  sd=sd(avg_prop,na.rm=TRUE),se=std.error(avg_prop,na.rm=TRUE))); NOCTACTmed_summ

wilcox.test(avg_prop~sex,data=noctact_med)

# NOCTURNAL ACTIVITY: SEX DIFFS AT HIGH ALAN
noctact_high<-nightprop %>% filter(site=='greens')
NOCTACThigh_summ<-as.data.frame(group_by(noctact_high,sex) %>% summarize(meanprop=mean(avg_prop,na.rm=TRUE),
  sd=sd(avg_prop,na.rm=TRUE),se=std.error(avg_prop,na.rm=TRUE))); NOCTACThigh_summ

wilcox.test(avg_prop~sex,data=noctact_high)

# NOCTURNAL ACTIVITY: smi DIFFS AT LOW ALAN
  plot(proportion~smi,data=noctact_low)
  
  NAsmilow_model<-glmer(proportion~smi+(1|tagID),family=binomial(link="logit"),
      weights=tot_nightdur,data=noctact_low); summary(NAsmilow_model)

# NOCTURNAL ACTIVITY: smi DIFFS AT MED ALAN
  plot(proportion~smi,data=noctact_med)
  
  NAsmimed_model<-glmer(proportion~smi+(1|tagID),family=binomial(link="logit"),
      weights=tot_nightdur,data=noctact_med); summary(NAsmimed_model)

# NOCTURNAL ACTIVITY: smi DIFFS AT HIGH ALAN
  plot(proportion~smi,data=noctact_high)

  NAsmihigh_model<-glmer(proportion~smi+(1|tagID),family=binomial(link="logit"),
      weights=tot_nightdur,data=noctact_high); summary(NAsmihigh_model)
  
  # nocturnal activity across sites
  NOCTACT_SITEsumm<-as.data.frame(group_by(nightprop,site) %>% summarize(meanprop=mean(proportion,na.rm=TRUE),
    sd=sd(proportion,na.rm=TRUE),se=std.error(proportion,na.rm=TRUE)));  NOCTACT_SITEsumm

  modelNA<-glmer(proportion~site+sex+smi+(1|tagID),family=binomial(link="logit"),
        weights=tot_nightdur,data=nightprop);    summary(modelNA)

  emmeans(modelNA, list(pairwise~site), adjust='tukey')


##__ FEEDING RATES __
#************#    
# __TOTAL FEED RATE model -------------------------------------------------
as.data.frame(group_by(reprodata,location) %>% summarize(meanfeed=mean(feedtot,na.rm=TRUE),
  sd=sd(feedtot,na.rm=TRUE),se=std.error(feedtot,na.rm=TRUE)))

modelFR<-glm(feedtot~location+smi1+smi2,family=poisson(link="log"), data=reprodata);summary(modelFR)
modelFR.aov<-aov(modelFR); summary(modelFR.aov)
TukeyHSD(modelFR.aov,which='location')

cor.test(x=reprodata$feedtot,y=reprodata$smi1,method='pearson')
cor.test(x=reprodata$feedtot,y=reprodata$smi2,method='pearson')


fitness_high<-reprodata %>% filter(location=='greens')
cor.test(x=fitness_high$feedtot,y=fitness_high$smi2,method='kendall')

fitness_med<-reprodata %>% filter(location=='campus')
cor.test(x=fitness_med$feedtot,y=fitness_med$smi1,method='kendall')

fitness_low<-reprodata %>% filter(location=='farm')
cor.test(x=fitness_low$feedtot,y=fitness_low$smi1,method='kendall')

as.data.frame(group_by(reprodata,location) %>% summarize(meanfeed=mean(feedtot,na.rm=TRUE),
  sd=sd(feedtot,na.rm=TRUE),se=std.error(feedtot,na.rm=TRUE)))


feeddata<-reprodata %>% pivot_longer(cols=c('mfeed','ffeed'),names_to='sex',values_to='feedrate')
feeddata2<-data.frame(site=feeddata$location,sex=feeddata$sex,feedrate=feeddata$feedrate)
feeddata2$sex<-as.factor(feeddata2$sex)

FR_low<-feeddata2 %>% filter(site=='farm')
wilcox.test(feedrate~sex,data=FR_low)
t.test(feedrate~sex,data=FR_low)

FR_med<-feeddata2 %>% filter(site=='campus')
wilcox.test(feedrate~sex,data=FR_med)
t.test(feedrate~sex,data=FR_med)

FR_high<-feeddata2 %>% filter(site=='greens')
wilcox.test(feedrate~sex,data=FR_high)
t.test(feedrate~sex,data=FR_high)

FRlow_summ<-group_by(FRlow,location) %>% summarize(meanMale=mean(mfeed,na.rm=TRUE), 
    meanFem=mean(ffeed,na.rm=TRUE),   sd=sd(mfeed,na.rm=TRUE),se=std.error(mfeed,na.rm=TRUE))  #summarize
  FRlow_summ<-as.data.frame(FRlow_summ);      FRlow_summ  #view
  

feedtable<-data.frame(site=reprodata$location,mfeed=reprodata$mfeed,ffeed=reprodata$ffeed)
feedtable<-feedtable %>% group_by(site) %>% summarise(meanMale=mean(mfeed),meanFem=mean(ffeed))

FRlow<-feeddata2 %>% filter(site=='farm')




##*******##

#************# 
# __OFFSPRING MASS model --------------------------------------------------
as.data.frame(group_by(reprodata,location) %>% summarize(meanmass=mean(avgmassyoung,na.rm=TRUE),
  sd=sd(avgmassyoung,na.rm=TRUE),se=std.error(avgmassyoung,na.rm=TRUE)))

modelB<-glm(avgmassyoung~location+smi1+smi2, data=reprodata)
summary(modelB)
modelB.aov<-aov(modelB); summary(modelB.aov)

TukeyHSD(modelB.aov,which='location')
cor.test(x=reprodata$avgmassyoung,y=reprodata$smi1,method='pearson')
cor.test(x=reprodata$avgmassyoung,y=reprodata$smi2,method='pearson')


cor.test(x=fitness_high$avgmassyoung,y=fitness_high$smi2,method='kendall')

cor.test(x=fitness_med$avgmassyoung,y=fitness_med$smi2,method='kendall')

cor.test(x=fitness_low$avgmassyoung,y=fitness_low$smi2,method='kendall')

#************# 
# __BROOD SIZE model ------------------------------------------------------
modelC1<-glm(broodsize~location+smi1+smi2,family=poisson(link="log"),data=reprodata);summary(modelC1)
modelC1.aov<-aov(modelC1); summary(modelC1.aov)
TukeyHSD(modelC1.aov,which="location")


