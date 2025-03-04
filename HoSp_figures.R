## FIGURES ##
library(ggplot2)
library(plotrix)
library(dplyr)
library(ggbeeswarm)

  ## ONSET -------------------------------------------------------------------
    # LOAD "onsetoffset" from "HoSp_processing.R"
    # summarize relative onset (sd, se, mean) across sites
    ONS_summary <- onsetoffset %>% dplyr::group_by(site) %>% dplyr::summarise(sd = sd(avg_onset_rel,
        na.rm =TRUE), se=std.error(avg_onset_rel,na.rm =TRUE),avg_onset_rel = mean(avg_onset_rel,
        na.rm =TRUE)); print.data.frame(ONS_summary)
  
    ggplot(data=onsetoffset, aes(x=site,y=avg_onset_rel,na.rm=TRUE)) +  
        # adds geom_jitter but decreases point overlap
      geom_quasirandom(size=3.0,alpha=0.3,width=0.2) + 
        # line to indicate sunrise
      geom_hline(yintercept=0, linewidth=1,linetype="dashed")  +
        # plot mean
      geom_point(aes(y=avg_onset_rel),data = ONS_summary, size=5.0) +
        # add error bars using se
      geom_errorbar(aes(ymin = avg_onset_rel-se, ymax = avg_onset_rel+se),data = ONS_summary,
          linewidth=1.2,width=0.07) +
      ylab("Mean relative activity onset (min)") + 
      theme_classic() + theme(axis.line=element_line(linewidth=0.5),axis.title.x=element_blank(),
          axis.text=element_text(size=24),axis.title.y=element_text(size=26), legend.position="none",
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))  +
      scale_x_discrete(labels=c("Low ALAN","Medium ALAN",
          "High ALAN")) +  theme (legend.position = "none") + ylim(-60,20)


  ## OFFSET ------------------------------------------------------------------
    # LOAD "onsetoffset" from "HoSp_processing.R"
    # summarize relative offset (sd, se, mean) across sites
    OFFS_summary <- onsetoffset %>% dplyr::group_by(site) %>% dplyr::summarise(sd = sd(avg_offset_rel,
        na.rm=TRUE),se=std.error(avg_offset_rel,na.rm =TRUE), avg_offset_rel = mean(avg_offset_rel,
        na.rm =TRUE));  print.data.frame(OFFS_summary)
    
    ggplot(data=onsetoffset, aes(x=site,y=avg_offset_rel,na.rm=TRUE)) +  
      geom_quasirandom(size=3.0,alpha=0.3,width=0.2) + 
      geom_hline(yintercept=0, linewidth=1,linetype="dashed") +
      geom_point(aes(y=avg_offset_rel),data = OFFS_summary,  size=5.0) +
      geom_errorbar(aes(ymin = avg_offset_rel-se, ymax = avg_offset_rel+se),data = OFFS_summary,
        linewidth=1.2,width=0.07) +
      ylab("Mean relative activity offset (min)") + 
      theme_classic() + theme(axis.line=element_line(linewidth=0.5),axis.title.x=element_blank(),
        axis.text=element_text(size=24),axis.title.y=element_text(size=26), legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) + 
      scale_x_discrete(labels=c("Low ALAN","Medium ALAN","High ALAN")) +  ylim (-23,60)
    
  ## NOCTURNAL ACTIVITY ------------------------------------------------------
    # 3 bar plots
    NA_summary <- nightprop %>% dplyr::group_by(site) %>% dplyr::summarise(sd = sd(avg_prop,
        na.rm=TRUE),se=std.error(avg_prop,na.rm =TRUE), proportion = mean(avg_prop, na.rm=TRUE))
        print.data.frame(NA_summary)
        
    ggplot(data=NA_summary, aes(x=site,y=proportion)) +
      geom_bar(stat="identity",position="stack")+
      geom_errorbar(aes(ymin = proportion-se, ymax = proportion+se),data=NA_summary,
          linewidth=1.2,width=0.07)+
      ylab("Proportion of night spent active") + 
      theme_classic() + theme(axis.line=element_line(linewidth=0.5),axis.title.x=element_blank(),
          axis.text=element_text(size=24),axis.title.y=element_text(size=26),
          legend.title=element_blank(),legend.text=element_blank(), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))  +
      scale_x_discrete(labels=c("Low ALAN","Medium ALAN","High ALAN")) +  
      theme(legend.position = "none") +
      ylim (0,0.20)
    
    

  ## OFFSPRING MASS ----------------------------------------------------------
    df.summaryA <- reprodata %>% dplyr::group_by(location) %>% dplyr::summarise(
      sd = sd(avgmassyoung,na.rm=TRUE),se=std.error(avgmassyoung,na.rm =TRUE), 
      avgmassyoung = mean(avgmassyoung,na.rm =TRUE)); print.data.frame(df.summaryA)
        
    ggplot(data=reprodata, aes(x=location,y=avgmassyoung)) +
      geom_quasirandom(size=3.0,alpha=0.3,width=0.2)+
      geom_point(aes(y = avgmassyoung),data=df.summaryA,size=5.0)+ 
      geom_errorbar(aes(ymin = avgmassyoung-se, ymax = avgmassyoung+se),data=df.summaryA,
          linewidth=1.2,width=0.06)+ 
      ylab("Mean offspring mass (grams)")+ 
      theme_classic() + theme(axis.line=element_line(linewidth=0.5),axis.title.x=element_blank(),
          axis.text=element_text(size=24),legend.title=element_blank(),legend.text=element_blank(), 
          axis.title.y=element_text(size=26),plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))  +
      scale_x_discrete(labels=c("Low ALAN","Medium ALAN",
          "High ALAN")) +  theme (legend.position = "none")

  ## FEEDING RATE ------------------------------------------------------------
    feedlong<-reprodata%>%pivot_longer(cols=c("mfeed","ffeed"),names_to="parent",values_to="feedrate")
    df.summaryB1 <- reprodata %>% dplyr::group_by(location) %>% dplyr::summarise(
      sd = sd(mfeed,na.rm=TRUE),se=std.error(mfeed,na.rm =TRUE), mfeed = mean(mfeed,na.rm =TRUE))
      print.data.frame(df.summaryB1)
    
    df.summaryB2 <- reprodata %>% dplyr::group_by(location) %>% dplyr::summarise(
      sd = sd(ffeed,na.rm=TRUE),se=std.error(ffeed,na.rm =TRUE), ffeed = mean(ffeed,na.rm =TRUE))
      print.data.frame(df.summaryB2)
    
    df.summary.feed<-merge(df.summaryB1,df.summaryB2,by="location")
    df.summary.feed<-df.summary.feed%>%pivot_longer(cols=c("mfeed","ffeed"),names_to="parent",
      values_to="feedrate")
    df.summary.feed$parent<-factor(df.summary.feed$parent)
    df.summary.feed<-df.summary.feed%>%reframe(location=location,
      se=ifelse(parent=="mfeed",se.x,se.y), parent=parent,feedrate=feedrate)
        
    ggplot(data=feedlong, aes(x=location,y=feedrate)) +
      geom_quasirandom(aes(color=parent,shape=parent),size=3.0,alpha=0.3,width=0.12,
          dodge.width=0.3) +
      geom_point(aes(color=parent,shape=parent),data=df.summary.feed,size=5.0,
          position=position_dodge(width=0.3)) + geom_errorbar(aes(ymin=feedrate-se,
          ymax=feedrate+se,color=parent),data=df.summary.feed, linewidth=1.2,width=0.1,
          position=position_dodge(width=0.3)) +
      scale_color_manual(values=c("red","blue"),labels=c("Female","Male")) +
      scale_shape_manual(values=c(16,17),labels=c("Female","Male")) + 
      labs (color=expression(underline("Parent")),shape=expression(underline("Parent"))) +
      ylab("Parental feeding rate (visits/hr)") + 
      theme_classic() + theme(axis.line=element_line(linewidth=0.5),axis.title.x=element_blank(),
          axis.text=element_text(size=24),legend.title=element_text(size=19),
          legend.text=element_text(size=18), axis.title.y=element_text(size=26),
          plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),legend.position="inside",
          legend.position.inside=c(0.82,0.85))  +
      scale_x_discrete(labels=c("Low ALAN","Medium ALAN",
          "High ALAN")) + ylim(c(4,15))

    ggplot( data=SEM_df, aes(x=noctact,y=avgmassyoung,color=site)) + geom_point()
    
    ggplot( data=SEM_df, aes(x=indv_feedrate,y=avgmassyoung,color=site)) + geom_point()
    
    ggplot( data=SEM_df, aes(x=noctact,y=indv_feedrate,color=site)) + geom_point()
    
    model<-glm(avgmassyoung~onset_rel*site,data=SEM_df);summary(model)
    ggplot( data=SEM_df, aes(x=onset_rel,y=avgmassyoung,color=site)) + geom_point()
    
    
    
    
    
    
    