
## PREPARE TEST-READY DATA FROM RAW DATA ##-----------------------------------
library(dplyr)
library(data.table)
library(lubridate)
library(hms)

  # ADD ACTIVITY ("0"/"1") AND "DAY"/"NIGHT" CATEGORIZATION TO RAW DATA --------
    
    ## Create lists of files and names for tag datsets
      csvlist<-list.files("~/HoSp Data/raw datasets/tagdata")
      tagnum<-gsub(".csv","",csvlist)
      tagnum<-substr(tagnum,5,7)
      csvnames<-paste("tag",tagnum,sep="_")
      setwd("~/HoSp Data/raw datasets/tagdata")
      for (k in 1:length(csvlist)) {
        assign(csvnames[k],read.csv(csvlist[k],header=TRUE))  }
      taglist<-mget(ls(pattern="tag_"))

    ## Sunrise and sunset times from timeanddate.com
      sunrisetimes<-as.POSIXct(c("2023-05-14 5:46","2023-05-15 5:45",
        "2023-05-16 5:44","2023-05-17 5:43","2023-05-18 5:42",
        "2023-05-19 5:41","2023-05-20 5:41","2023-05-21 5:40",
        "2023-05-22 5:39","2023-05-23 5:38","2023-05-24 5:38",
        "2023-05-25 5:37","2023-05-26 5:37","2023-05-27 5:36",
        "2023-05-28 5:35","2023-05-29 5:35","2023-05-30 5:34",
        "2023-05-31 5:34","2023-06-01 5:34"),format="%Y-%m-%d %H:%M")
      sunsettimes<-as.POSIXct(c("2023-05-14 20:05","2023-05-15 20:06",
        "2023-05-16 20:07","2023-05-17 20:08","2023-05-18 20:09",
        "2023-05-19 20:10","2023-05-20 20:11","2023-05-21 20:11",
        "2023-05-22 20:12","2023-05-23 20:13","2023-05-24 20:14",
        "2023-05-25 20:15","2023-05-26 20:16","2023-05-27 20:16",
        "2023-05-28 20:17","2023-05-29 20:18","2023-05-30 20:19",
        "2023-05-31 20:19","2023-06-01 20:20"),format="%Y-%m-%d %H:%M")

        ###____ custom function to add Julian day, activity status (0/1), & day/night to raw datasets.
        full_data<-function(dataset){
          #categorizes values in "Date" column as a date. Used when calculating day number.
          dataset$Date<-as.Date(dataset$Date,format="%m/%d/%Y",tz="Etc/GMT+7")
          dataset$Time<-as.character(dataset$Time)
          #omit values of "NA"
          dataset<-na.omit(dataset)
          #create empty columns equal in length to the entire dataset
          daynum<-rep(0,length(dataset$Date))
          dbdiff<-rep(0,length(dataset$Date))
          activity<-rep(0,length(dataset$Date))
          night_day<-rep(0,length(dataset$Date))
          
          datetime<-strptime(paste(dataset$Date,dataset$Time),format="%Y-%m-%d %H:%M:%S")
          
          #determine day number by subtracting the date of day 1 from "current" date.
          for (a in 1:length(dataset$Date)) {
            daynum[a]<-as.numeric((dataset$Date[a])-(as.Date("12/31/2022", 
                format="%m/%d/%Y",tz="Etc/GMT+7")))                }
          #calculate difference between subsequent RSSI readings.
          for (b in 2:length(dataset$RSSI)) {
            dbdiff[b]<-abs(dataset$RSSI[b]-dataset$RSSI[b-1])  }
          #if difference in RSSI >/= 10, activity="1"; if < 10, activity="0".
          for (c in 2:length(dataset$dfdiff)) {
            activity[c]<-ifelse(dbdiff[c]>=10,1,0)    }
          
          #if time in datetime[d] lies between sunrise and sunset, then "day"- otherwise "night"
          #subsets narrow comparison down to sunset/rise dates that match the date of datetime[d]
          #date matters because sunrise and sunset times are different every day
          for (d in 1:length(dataset$Date)) {
            night_day[d]<-ifelse(
              datetime[d] > sunrisetimes[as.Date(sunrisetimes,format="%Y-%m-%d",
                  tz="Etc/GMT+7") %in% dataset$Date[d]]  &&
              datetime[d] < sunsettimes[as.Date(sunsettimes,format="%Y-%m-%d",
                  tz="Etc/GMT+7") %in% dataset$Date[d]] ,"day","night")
          }
          #creates new data frame from raw data that includes day, activity value, and day/night.
          fulldata<-data.frame(date=dataset$Date,time=dataset$Time,day=daynum,
              ID=dataset$TagID.BPM,db=dataset$RSSI,dbdiff=dbdiff,activity=activity,day.night=night_day)
      }

        ## To create one dataset:
          tag_1_full<-full_data(tag_1)
            ## Skip to "load data" section to quickly load all tag datasets.
          
 
        



   
  # SMI ---------------------------------------------------------------------
    ## function that calculates scaled mass index from weight and tarsus length, per Peig & Green 2009.
    smi<-function(M,L){   plot(log(M)~log(L))
      {
        if(require(smatr)){
          SMA<-sma(log(M)~log(L))
          bSMA<-coef(SMA)[2]  }
        else { 
          OLS<-lm(log(M)~log(L))
          bOLS<-coef(OLS)[2]
          r<-cor.test(x=log(M),y=log(L),method='pearson')$estimate
          outliers<-which(abs(rstandard(OLS))>3)
          bSMA<-bOLS/r
        }   }   
      L0<-median(L,na.rm=T)
      SMi<-M*((L0/L)^bSMA)
      return(SMi) 
      }
      
      captdata<-read.csv("~/HoSp Data/raw datasets/capture data.csv",header=TRUE,na.strings=c("","NA"))
      smiframe<-data.frame(site=captdata$location,sex=captdata$sex,ID=captdata$radiotag,
        mass=captdata$weight,tarsus=captdata$tarsus,smi=smi(captdata$weight,captdata$tarsus))
      smiframe$site<-as.factor(smiframe$site)
      smiframe$site<-factor(smiframe$site, levels=c("farm","campus", "greens"))
        smiID<-smiframe
        write.csv(smiID, "~/HoSp Data/prepared datasets/smiID.csv", row.names=FALSE)

  # NOCTURNAL ACTIVITY -----------------------------------------------------
    # load all tag datasets first.
    ## create a dataframe with ALL tags
    alltaglist<-list(tag_8_full,tag_5_full,tag_11_full,tag_9_full,tag_14_full,tag_15_full,tag_29_full,
      tag_30_full,tag_51_full,tag_52_full,tag_59_full,tag_60_full,tag_33_full,tag_34_full,tag_39_full,
      tag_40_full,tag_46_full,tag_47_full,tag_55_full,tag_56_full,tag_50_full,tag_10_full,tag_16_full,
      tag_27_full,tag_28_full,tag_62_full,tag_1_full,tag_23_full,tag_44_full,tag_45_full,tag_35_full,
      tag_36_full,tag_17_full,tag_18_full,tag_21_full,tag_22_full,tag_26_full,tag_19_full,tag_20_full,
      tag_24_full,tag_25_full,tag_31_full,tag_32_full,tag_37_full,tag_38_full,tag_48_full,tag_49_full,
      tag_6_full,tag_7_full,tag_57_full,tag_58_full,tag_12_full,tag_13_full,tag_42_full,tag_43_full,
      tag_53_full,tag_54_full)
    alltags<-rbindlist(alltaglist)
    
    write.csv(alltags, "~/HoSp Data/prepared datasets/alltags.csv", row.names=FALSE)
    
      #convert 'day' to numeric
      alltags$day<-as.character(alltags$day)
      alltags$day<-as.numeric(alltags$day)
      
      #filter all data by entries where it is nighttime.
      nightfreq<-alltags %>% filter(day.night=='night')
      
      nightfreq$time<-lubridate::hms(nightfreq$time)

      #then sum activity values to get frequency of nighttime bouts for each half of the night
      nightfreq<-nightfreq %>% group_by(day,ID) %>% summarize(freq_dawn=sum(activity[hour(time)<12]), 
        freq_dusk=sum(activity[hour(time)>12])); colnames(nightfreq)<-c("day","tagID","freq_dawn",
          "freq_dusk")
      #create new columns that convert frequency to activity duration in mins
      nightfreq<-nightfreq %>% mutate(dur_dawn=(freq_dawn*10)/60, dur_dusk=(freq_dusk*10)/60)
    
    #difftime sunset/rise
    sunrisetimes<-ymd_hm(c("2023-05-14 5:46","2023-05-15 5:45",
        "2023-05-16 5:44","2023-05-17 5:43","2023-05-18 5:42",
        "2023-05-19 5:41","2023-05-20 5:41","2023-05-21 5:40",
        "2023-05-22 5:39","2023-05-23 5:38","2023-05-24 5:38",
        "2023-05-25 5:37","2023-05-26 5:37","2023-05-27 5:36",
        "2023-05-28 5:35","2023-05-29 5:35","2023-05-30 5:34",
        "2023-05-31 5:34","2023-06-01 5:34"), tz="US/Pacific")
      sunrise1<-data.frame(date=as_date(sunrisetimes),datetime=sunrisetimes)
      
    sunsettimes<-ymd_hm(c("2023-05-14 20:05","2023-05-15 20:06",
        "2023-05-16 20:07","2023-05-17 20:08","2023-05-18 20:09",
        "2023-05-19 20:10","2023-05-20 20:11","2023-05-21 20:11",
        "2023-05-22 20:12","2023-05-23 20:13","2023-05-24 20:14",
        "2023-05-25 20:15","2023-05-26 20:16","2023-05-27 20:16",
        "2023-05-28 20:17","2023-05-29 20:18","2023-05-30 20:19",
        "2023-05-31 20:19","2023-06-01 20:20"),tz="US/Pacific")
      sunset1<-data.frame(date=as_date(sunsettimes),datetime=sunsettimes)
      
      #create data frame with date, sunrise time, sunset time
      sunlist1<-inner_join(sunrise1,sunset1,by="date"); colnames(sunlist1)<-c("date","sunrise","sunset")
      
      #add column for julian day and difftimes between 0:00-sunrise and sunset-24:00
      sunlist1<-sunlist1%>%mutate(day=as.numeric((sunlist1$date)-(mdy("12/31/2022"))),
        nightdur_1=as.numeric(difftime(as_hms(sunrise),as_hms("00:00:00"),units="mins")),
        nightdur_2=as.numeric(difftime(as_hms("24:00:00"),as_hms(sunset),units="mins")))

      #create dataframe 
      nightprop<-nightfreq%>%mutate(nightdur_1=sunlist1$nightdur_1[sunlist1$day %in% day], 
        nightdur_2=sunlist1$nightdur_2[sunlist1$day %in% day])
        
        
        
        
        prop_1=(dur_dawn/sunlist1$nightdur_1[sunlist1$day %in% day])
        prop_2=(dur_dusk/sunlist1$nightdur_2[sunlist1$day %in% day])
      
      write.csv(nightprop, "~/HoSp Data/prepared datasets/nightpropX.csv", row.names=FALSE)
      
      #open 'nightprop' in Excel and add columns for site, sex, and smi using XLOOKUP and 'smiID'
      # manually removed values for first half of first day and second half of last day for which
          # no activity data exist
      # manually calculated proportion for each day and ID and average proportion per ID
      
    
  # REPRO DATA ---------------------------------------
    reprodata<-read.csv("~/HoSp Data/raw datasets/reprod_data.csv",header=TRUE,na.strings=c("","NA"))
    reprodata<-reprodata %>% mutate(reprodata,feedtot=reprodata$mfeed+reprodata$ffeed)
    ## write as csv and edit in excel:
      # add smi of parent 1 and parent 2 using XLOOKUP from 'smiID'
      write.csv(reprodata, "~/HoSp Data/prepared datasets/reprod_data.csv", row.names=FALSE)


  # SEM DF ------------------------------------------------------------------
    ## MEGA_df.csv was created in Excel by amalgamating capture data.csv, reprod_data.csv, onsetoffset.csv,
      # smiID, and nightprop2.csv using Lookup functions. Load from "HoSp_SEM" script.

      
      
      
      
## OR LOAD TEST-READY DATASETS -------------------------------------------------
  # LOAD: SMI ---------------------------------------------------------------
    smiID<-read.csv("~/HoSp Data/prepared datasets/smiID.csv",header=TRUE,na.strings=c("","NA"))

  # LOAD: all tag datasets -----------------------------------------------------
      fulltaglist<-list.files("~/HoSp Data/prepared datasets/full tags")
    
      for (o in 1:length(fulltaglist)) {
        setwd("~/HoSp Data/prepared datasets/full tags")
        load(fulltaglist[o])
      }

  # LOAD: ONSET & OFFSET ----------------------------------------------------
    onsetoffset<-read.csv("~/HoSp Data/prepared datasets/onsetoffset.csv",header=TRUE,na.strings=c("","NA"))
      onsetoffset$site<-as.factor(onsetoffset$site)
      onsetoffset$site<-factor(onsetoffset$site, levels=c("farm", "campus", "greens"))
      onsetoffset$sex<-as.factor(onsetoffset$sex)
      onsetoffset$tagID<-as.factor(onsetoffset$tagID)
      onsetoffset$smi<-as.numeric(onsetoffset$smi)
      
      onsetoffset$avg_onset_rel<-as.numeric(onsetoffset$avg_onset_rel)
      onsetoffset$avg_offset_rel<-as.numeric(onsetoffset$avg_offset_rel)
      
      # average on/offsets are in seconds. convert to minutes.
      onsetoffset$avg_onset_rel<-onsetoffset$avg_onset_rel/60
      onsetoffset$avg_offset_rel<-onsetoffset$avg_offset_rel/60
      

  # LOAD: NOCTURNAL ACTIVITY ------------------------------------------------
      nightprop<-read.csv("~/HoSp Data/prepared datasets/nightprop2.csv",header=TRUE,na.strings=c("","NA"))
      nightprop$site<-as.factor(nightprop$site)
      nightprop$site<-factor(nightprop$site, levels=c("farm","campus", "greens"))
      nightprop$sex<-as.factor(nightprop$sex)
      nightprop$tagID<-as.factor(nightprop$tagID)
      nightprop$smi<-as.numeric(nightprop$smi)
      
      nightprop$proportion<-as.numeric(nightprop$proportion)
      nightprop$avg_prop<-as.numeric(nightprop$avg_prop)
      nightprop$tot_duration<-as.numeric(nightprop$tot_duration)
      nightprop$tot_nightdur<-as.numeric(nightprop$tot_nightdur)


  # LOAD: FITNESS DATA ------------------------------------------------------
    reprodata<-read.csv("~/HoSp Data/prepared datasets/reprod_data2.csv",header=TRUE,na.strings=c("","NA"))
      reprodata$location<-as.factor(reprodata$location)
      reprodata$location<-factor(reprodata$location, levels=c("farm","campus","greens"))
      reprodata$broodsize<-as.numeric(reprodata$broodsize)
      reprodata$avgmassyoung<-as.numeric(reprodata$avgmassyoung)
      reprodata$mfeed<-as.numeric(reprodata$mfeed)
      reprodata$ffeed<-as.numeric(reprodata$ffeed)
      reprodata$feedtot<-as.numeric(reprodata$feedtot)
      
      reprodata$smi1<-as.numeric(reprodata$smi1)
      reprodata$smi2<-as.numeric(reprodata$smi2)

  
          
          