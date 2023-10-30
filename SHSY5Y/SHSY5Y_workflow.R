list.of.packages <- c("dplyr","EnvStats","stringr","gtools","ggplot2","ggprism","yaml","ggpubr","patchwork","tcpl","drc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(EnvStats)
library(stringr)
library(gtools) 
library(ggplot2) #needed for plotting
library(ggprism) #needed for Prism themes
library(ggpubr) #needed for plot alignment
library(patchwork) #needed to combine plots
library(drc)
library(tcpl)
library(yaml)

#Functions
Rosner_outlier_removal = function(x){
  k = as.integer(0.1*length(x))
  if(k<3){
    k = 3
  }
  Rosner = tryCatch(rosnerTest(x,k),warning=function(w){Rosner=NA})
  if(!is.na(Rosner)){
    Resulta = Rosner$all.stats
    Outliers = Resulta$Outlier == T
    Outliers = dplyr::filter(Resulta,Outliers==T)
    for(k in 1:length(Outliers$Outlier)){
      r = Outliers$Obs.Num[k]
      x[r] = NA
    }
  }
  return(x)
}

Split_collapsed = function(x){
  y = stringr::str_split(x,";")
  return(as.numeric(unlist(y)))
}

create_names = function(x){
  names = rep(x$Measurement,each=2)
  names1 = rep(x$Sample_ID,each=2)
  what = rep(c("Cytotoxicity","Effect"),times=nrow(x))
  final = rep("c",length(names))
  for(u in 1:length(names)){
    final[u] = paste0(names[u],"_",names1[u],"_",what[u])
  }
  return(as.vector(final))
}

asteriks = function(x,y){
  for(z in 1:length(y)){
    if(is.na(y[z])){
      x[z] = paste0(x[z],"*")
    }
  }
  return(x)
}

calculate_SE_ac10 <- function(tp,bt,y,ac50, SE_ac50, slope, SE_slope) {
  # Calculate partial derivatives
  d_ac50 <- 1
  d_slope <- -log10((tp - bt) / (y - bt) - 1) / (slope^2)
  
  # Calculate SE_ac10 using error propagation formula
  SE_ac10 <- sqrt(
    (d_ac50^2) * (SE_ac50^2) +
      (d_slope^2) * (SE_slope^2)
  )
  
  return(SE_ac10)
}

makeExcelreadable = function(x){
  return(stringr::str_replace_all(x,"e-","E-"))
}
################################################################################

#Sets the current project as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Gets project name
project = dirname(rstudioapi::getActiveDocumentContext()$path)
assays = list.files(project)
assays = paste0(assays[!grepl("Summary",assays)],collapse = "|")
project = unlist(stringr::str_split(project,"/"))[length(unlist(stringr::str_split(project,"/")))-1]

#Summary
Summary_path = list.files("../","Summary",full.names = T)
Summary_path = Summary_path[!grepl("filled|~",Summary_path)]
Summary_list = openxlsx::read.xlsx(Summary_path)
Summary_results = openxlsx::read.xlsx(Summary_path,2)

if(!file.exists("./Plots")){
  dir.create("./Plots")
}

#settings
settings = yaml::yaml.load_file(paste0(list.files(pattern="yaml",full.names = T)))
list2env(settings,envir = .GlobalEnv)

#Search for each sample i - merged data will be included 
for(i in 1:nrow(Summary_list)){
  if(is.na(Summary_list$SHSY5Y[i])){
    next
  }
  Sample_ID = Summary_list$Sample_ID[i]
  #Unified_IDs = as.character(sapply(Summary_list$Sample_ID,Unify_IDs))
  Unified_IDs = Summary_list$Sample_ID
  frame = Summary_list[(Unified_IDs==Sample_ID)&!is.na(Summary_list$SHSY5Y),]
  if(nrow(frame)==0){
    next
  }
  Nr_of_samples = nrow(Summary_list[(Unified_IDs==Sample_ID)&!is.na(Summary_list$SHSY5Y),])
  workframe = as.data.frame(matrix(nrow = Nr_of_samples,ncol = 10))
  colnames(workframe)=c("Sample_ID","Measurement","plate","Column","Unit","Concentration","Effect","Cytotox","Effect_Control","Cytotox_Control")
  workframe$Sample_ID = frame$Sample_ID
  workframe$Measurement = frame$Measurement
  workframe$Unit = frame$Unit[1]
  
  if(length(workframe$Measurement[!is.na(workframe$Measurement)])==nrow(workframe)){
    next
  }
  
  Measurements = list.files(full.names = T)
  Measurements = Measurements[!grepl("tcpl|Plots",Measurements)]
  
  #-----------------------------------------------------------------------------
  
  #Set up workframe for analysis - get all coordinates of samples over all measurements
  Starti = 1
  o = 1
  max = length(Measurements)
  for(j in 1:nrow(workframe)){
    Dosingsheet = list.files(Measurements[o],pattern = "dosing",full.names = T)
    Dosingsheet = Dosingsheet[!grepl("~",Dosingsheet)]
    Samplelist = openxlsx::read.xlsx(Dosingsheet)
    matchi = match(workframe$Sample_ID[j],Samplelist$Sample_ID[Starti:length(Samplelist$Sample_ID)])+Starti-1
    while(is.na(matchi)|matchi>nrow(Samplelist)){
      Starti = 1
      o = o+1
      if(o>max){
        break
      }
      Dosingsheet = list.files(Measurements[o],pattern = "dosing",full.names = T)
      Dosingsheet = Dosingsheet[!grepl("~",Dosingsheet)]
      Samplelist = openxlsx::read.xlsx(Dosingsheet)
      matchi = match(workframe$Sample_ID[j],Samplelist$Sample_ID[Starti:length(Samplelist$Sample_ID)])+Starti-1
    }
    Starti=matchi+1
    workframe$Measurement[j] = Measurements[o]
    workframe$plate[j] = Samplelist$Plate[matchi]
    workframe$Column[j] = Samplelist$Col.[matchi]
    
    Dosing_type = Samplelist$Dosing[match(workframe$plate[j],Samplelist$Plate)]
    if(Dosing_type=="serial"){
      n = 2
    } else {
      n = 3
    }
    n = n+(workframe$plate[j]-1)*2
    Concentrations = openxlsx::read.xlsx(Dosingsheet,n,rowNames = T)
    Row = match(workframe$Column[j],rownames(Concentrations))
    workframe$Concentration[j] = paste0(Concentrations[Row,1:22],collapse = ";")
    #Raw data
    Rawdata = list.files(Measurements[o],paste0("p",workframe$plate[j]),full.names = T)
    Neurite_file = Rawdata[grepl("t48",Rawdata)]
    Green_file = Rawdata[grepl("Green",Rawdata)]
    Red_file = Rawdata[grepl("Red",Rawdata)]
    
    Neurite = read.table(Neurite_file,dec = ",",skip = 9)
    Green = read.table(Green_file,dec = ",",skip = 9)
    Green[Green==0] = 1e-30
    Red = read.table(Red_file,dec = ",",skip = 9)
    
    if(F%in%apply(Neurite,c(1,2),is.numeric)){
      Neurite=read.table(Neurite_file,dec = ".",skip = 9)
      Green=read.table(Green_file,dec = ".",skip = 9)
      Green[Green==0] = 1e-30
      Red=read.table(Red_file,dec = ".",skip = 9)
    }
    
    #Neurite outgrowth
    #Unexposed_mean = mean(Rosner_outlier_removal(unlist(c(Neurite[,23:24]))),na.rm = T)
    Unexposed_mean = mean(outliers::rm.outlier(as.numeric(unlist(Neurite[,23:24]))),na.rm=T)
    Neurite_outgrowth = 100-(Neurite/Unexposed_mean*100)
    workframe$Effect[j] = paste0(Neurite_outgrowth[Row,1:22],collapse = ";")
    workframe$Effect_Control = paste0(unlist(c(Neurite_outgrowth[,23:24])),collapse = ";")
    
    #Viability
    Values = 100*(1-Red/Green)
    Values[Values<0] = 0
    #Unexposed_mean = mean(Rosner_outlier_removal(unlist(c(Values[,23:24]))),na.rm = T)
    Unexposed_mean = mean(outliers::rm.outlier(as.numeric(unlist(Values[,23:24]))),na.rm=T)
    Cytotoxicity = 100-(Values/Unexposed_mean*100)
    workframe$Cytotox[j] = paste0(Cytotoxicity[Row,1:22],collapse = ";")
    workframe$Cytotox_Control = paste0(unlist(c(Cytotoxicity[,23:24])),collapse = ";")
  }
  
  #-----------------------------------------------------------------------------
  
  if(combine_repeats==F&!grepl(Reference,Sample_ID)){
    indices = which(Unified_IDs==Sample_ID)
    indices = indices[!is.na(Summary_list$SHSY5Y[indices])]
    #If sample is used multiple times in summary but not called by the respective assay
    if(!(i%in%indices)){
      next
    }
    
    #Update summary
    Summary_list$Measurement[indices] = workframe$Measurement
    Summary_list$Plate[indices] = workframe$plate
    Summary_list$Position[indices] = workframe$Column
    
    workframe = workframe[match(i,indices),]
    level = match(i,indices)
    if(level>1){
      workframe$Sample_ID = paste0(workframe$Sample_ID,"_repeat",(level-1))
      Sample_ID = workframe$Sample_ID
    }
    
  } else {
    indices = which(Unified_IDs==Sample_ID)
    indices = indices[!is.na(Summary_list$SHSY5Y[indices])]
    
    #If sample is used multiple times in summary but not called by the respective assay
    if(!(i%in%indices)){
      next
    }
    
    #Update summary
    Summary_list$Measurement[indices] = workframe$Measurement
    Summary_list$Plate[indices] = workframe$plate
    Summary_list$Position[indices] = workframe$Column
    
    level = match(i,indices)
    if(level>1){
      next
    }
  }
  
  #-----------------------------------------------------------------------------
  
  #Fit data from workframe, create plots, set up input for Prism, and calculate EC10 and IC10 
  
  #Setup
  Concentrations_linear = as.vector(sapply(workframe$Concentration,Split_collapsed))
  Concentrations_log = log10(Concentrations_linear)
  
  Effect = as.vector(sapply(workframe$Effect,Split_collapsed))
  Effect[Effect==Inf|Effect==-Inf] = NA
  Effect_all = Effect
  Effect_linear = Effect
  Effect_prism = Effect
  Effect_linear_prism = Effect
  Effect_control = as.vector(sapply(workframe$Effect_Control,Split_collapsed))
  Effect_control_prism = Effect_control
  Removed = outliers::rm.outlier(Effect_control)
  Effect_control[!(Effect_control%in%Removed)] = NA
  
  Cytotoxicity = as.vector(sapply(workframe$Cytotox,Split_collapsed))
  Cytotoxicity[Cytotoxicity==Inf|Cytotoxicity==-Inf] = NA
  Cytotoxicity_all = Cytotoxicity
  Cytotoxicity_prism = Cytotoxicity
  Cytotoxicity_linear = Cytotoxicity
  Cytotoxicity_linear_prism = Cytotoxicity
  Cytotoxicity_control = as.vector(sapply(workframe$Cytotox_Control,Split_collapsed))
  Cytotoxicity_control_prism = Cytotoxicity_control
  Cytotoxcicity_control_all = Cytotoxicity_control
  Removed = outliers::rm.outlier(Cytotoxicity_control)
  Cytotoxicity_control[!(Cytotoxicity_control%in%Removed)] = NA
  #Cytotoxicity_control = Rosner_outlier_removal(Cytotoxicity_control)
  
  #Set up control frames for plots
  Linear_all_control_frame = data.frame("Concentration"=10^log_concentration_unexposed_controls,
                                        "Name"=rep(workframe$Measurement,each=32),
                                        "Cytotox"=Cytotoxicity_control,
                                        "Effect"=Effect_control)
  Linear_all_control_frame = Linear_all_control_frame[!is.na(Linear_all_control_frame$Cytotox)&!is.na(Linear_all_control_frame$Effect),]
  Log_all_control_frame = Linear_all_control_frame
  Log_all_control_frame[,1] = log10(Linear_all_control_frame[,1])
  
  
  Linear_cyt_control_frame = data.frame("x"=10^log_concentration_unexposed_controls,
                                        "y"=Cytotoxicity_control,
                                        "Name"=rep(workframe$Measurement,each=32))
  Linear_cyt_control_frame = Linear_cyt_control_frame[Linear_cyt_control_frame$y<linear_cutoff_cytotoxicity,]
  
  Linear_effect_control_frame = data.frame("x"=10^log_concentration_unexposed_controls,
                                           "y"=Effect_control,
                                           "Name"=rep(workframe$Measurement,each=32))
  Linear_effect_control_frame = Linear_effect_control_frame[Linear_effect_control_frame$y<linear_cutoff_effect,]
  
  Log_cyt_control_frame = data.frame("x"=log_concentration_unexposed_controls,
                                     "y"=Cytotoxicity_control,
                                     "Name"=rep(workframe$Measurement,each=32))
  Log_effect_control_frame = data.frame("x"=log_concentration_unexposed_controls,
                                        "y"=Effect_control,
                                        "Name"=rep(workframe$Measurement,each=32))
  
  #-----------------------------------------------------------------------------
  #Cytotoxicity
  #Test if Cytotoxicity is different from control
  #t_test = t.test(Effect,Effect_control,alternative = "greater")
  upper_control = quantile(Cytotoxicity_control,probs=c(0.95),na.rm=T)
  #lower_control = quantile(Effect_control,probs=c(0.05),na.rm=T)
  replicates = as.data.frame(table(Concentrations_linear))[1,2]
  control_threshold = replicates*2
  above = length(Cytotoxicity[Cytotoxicity>upper_control])
  #below = length(Effect[Effect<lower_control])
  f_test = var.test(Cytotoxicity,Cytotoxicity_control)
  f_test_p = f_test$p.value
  #below = length(Effect[Effect<lower_control])
  
  if(above<control_threshold|f_test_p>0.05){
    Cytotoxicity_dev_Control = F
  } else {
    Cytotoxicity_dev_Control = T
  }
  
  #Linear Regression
  cytotoxicity_linear_split = as.data.frame(matrix(Cytotoxicity_linear,nrow=22))
  
  #use the 1 % quantile of the control as threshold, but round up to 1% (outliers could be problematic)
  precipitation_threshold = plyr::round_any(quantile(Cytotoxicity_control,probs = c(0.01),na.rm = T),1,f=ceiling)
  if(precipitation_threshold>-10){
    precipitation_threshold = -10
  }
  cytotoxicity_linear_split = cytotoxicity_linear_split<(precipitation_threshold)
  for(c in 1:ncol(cytotoxicity_linear_split)){
    if(!is.na(match(T,cytotoxicity_linear_split[,c]))&match(T,cytotoxicity_linear_split[,c])<=4){
      k = match(F,cytotoxicity_linear_split[match(T,cytotoxicity_linear_split[,c]):length(cytotoxicity_linear_split[,c]),c])-1
      cytotoxicity_linear_split[1:k,c] = T
    } else {
      cytotoxicity_linear_split[,c] = F
    }
  }
  Precipitation_filter_linear = c(cytotoxicity_linear_split)
  
  Cytotoxicity_linear[Cytotoxicity_linear>linear_cutoff_cytotoxicity] = NA
  Cytotoxicity_linear[Precipitation_filter_linear] = NA
  
  Linear_cyt_frame = data.frame("x"=Concentrations_linear,
                                "y"=Cytotoxicity_linear,
                                "Name"=rep(workframe$Measurement,each=22))
  frame = as.data.frame(cbind(Concentrations_linear,Cytotoxicity_linear))
  Linear_cyt_frame = Linear_cyt_frame[!is.na(frame$Cytotoxicity_linear),]
  Cytotoxicity_linear[is.na(frame$Cytotoxicity_linear)] = NA
  frame = frame[!is.na(frame$Cytotoxicity_linear),]
  colnames(frame) = c("x","y")
  
  linear_regression_cytotox = tryCatch(lm(data=frame,y~0+x),error=function(e){linear_regression_cytotox=NA})
  influentials = c()
  if(!is.na(linear_regression_cytotox)){
    r2 = round(summary(linear_regression_cytotox)$adj.r.squared,digits = 2)
    Concis = frame[!is.na(frame[,2]),1]
    Concis = unique(Concis)
    if(r2 < 0.8 & length(Concis)>5){
      redo_calib = T
      r2_initial = r2
      while(redo_calib == T){
        if(length(Concis)==4){
          break
        }
        cooksD = cooks.distance(linear_regression_cytotox)
        influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
        Distances = cooksD[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
        if(length(influential)==0){
          influential = NaN
        }
        if(!is.nan(influential)){
          influential = influential[match(max(Distances),Distances)]
          frame_new = frame
          frame_new$y[match(max(influential),frame_new$y)] = NaN
          frame_new$y[frame_new$y%in%influential] = NaN
          Filter_frame = !is.nan(frame_new$y)
          frame_new = frame_new[!is.nan(frame_new$y),]
          model_new = lm(formula=y~0+x,data = frame_new)
          r2 = round(summary(model_new)$adj.r.squared,digits = 2)
          if(r2<r2_initial){
            redo_calib = F
            r2 = r2_initial
          } else {
            linear_regression_cytotox = model_new
            if(r2 >= 0.8){
              redo_calib = F
            }
            frame = frame_new
            Linear_cyt_frame = dplyr::filter(Linear_cyt_frame,Filter_frame)
            Cytotoxicity_linear[!Filter_frame] = NA
            Concis = frame[!is.na(frame[,2]),1]
            Concis = unique(Concis)
            linear_regression_cytotox = model_new
            r2_initial = r2
            influentials = append(influentials,influential)
          }
        } else {
          redo_calib = F
        }
      }
    }
    
    slope = as.numeric(linear_regression_cytotox$coefficients[1])
    SE_slope = sqrt(diag(vcov(linear_regression_cytotox)))/nobs(linear_regression_cytotox)
    
    plot(frame$y~frame$x,type="p")
    abline(linear_regression_cytotox)
    
  } else {
    Cytotoxicity_linear_prism = asteriks(Cytotoxicity_linear_prism,Cytotoxicity_linear)
  }
  if(!is.na(linear_regression_cytotox)){
    #Detect plateau and remove since this is start of log-logistic part of the fit
    unique_concentrations = nrow(as.data.frame(table(frame$x)))
    frame_plateau = as.data.frame(matrix(nrow=unique_concentrations,ncol=2))
    colnames(frame_plateau) = colnames(frame)
    frame_plateau$x = as.data.frame(table(frame$x))$Var1
    for(f in 1:nrow(frame_plateau)){
      frame_plateau$y[f] = mean(frame$y[frame$x==frame_plateau$x[f]],na.rm=T)
    }
    if(length(frame_plateau$y[frame_plateau$y>10])>2){
      frame_plateau = frame_plateau[order(frame_plateau$x,decreasing = F),]
      frame_plateau$x = as.numeric(as.character(frame_plateau$x))
      model = tryCatch(drc::drm(y~x,data=frame_plateau,fct=drc::MM.2()),
                       error=function(e){model=NA})
      if(!is.na(model)){
        model_linear = lm(y~0+x,data=frame_plateau)
        predi = predict(model,newdata = frame_plateau)
        predi_linear = predict(model_linear)
        plot(frame_plateau$y~frame_plateau$x,type="p")
        lines(frame_plateau$x,predi,col="red")
        lines(frame_plateau$x,predi_linear,col="blue")
        MM_AIC = round(as.numeric(2*3-2*(logLik(model))),digits = 0)
        Linear_AIC = round(as.numeric(2*1-2*(logLik(model_linear))),digits = 0)
        while(MM_AIC<(Linear_AIC-4)&length(frame_plateau$y)>4){
          frame_plateau = frame_plateau[-nrow(frame_plateau),]
          model = tryCatch(drc::drm(y~x,data=frame_plateau,fct=drc::MM.2()),
                           error=function(e){model=NA})
          if(is.na(model)){
            break
          }
          model_linear = lm(y~0+x,data=frame_plateau)
          predi = predict(model,newdata = frame_plateau)
          predi_linear = predict(model_linear)
          plot(frame_plateau$y~frame_plateau$x,type="p")
          lines(frame_plateau$x,predi,col="red")
          lines(frame_plateau$x,predi_linear,col="blue")
          MM_AIC = round(as.numeric(2*3-2*(logLik(model))),digits = 0)
          Linear_AIC = round(as.numeric(2*1-2*(logLik(model_linear))),digits = 0)
        }
        highest_conc = max(frame_plateau$x,na.rm=T)
      } else {
        highest_conc = max(frame$x,na.rm=T)
      }
    } else {
      highest_conc = max(frame$x,na.rm=T)
    }
    
    Cytotoxicity_linear[frame$x>=highest_conc] = NA
    Linear_cyt_frame = Linear_cyt_frame[frame$x<=highest_conc,]
    frame = frame[frame$x<=highest_conc,]
    
    linear_regression_cytotox = tryCatch(lm(data=frame,y~0+x),error=function(e){linear_regression_cytotox=NA})
    influentials = c()
    if(!is.na(linear_regression_cytotox)){
      r2 = round(summary(linear_regression_cytotox)$adj.r.squared,digits = 2)
      Concis = frame[!is.na(frame[,2]),1]
      Concis = unique(Concis)
      if(r2 < 0.8 & length(Concis)>5){
        redo_calib = T
        r2_initial = r2
        while(redo_calib == T){
          if(length(Concis)==4){
            break
          }
          cooksD = cooks.distance(linear_regression_cytotox)
          influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          Distances = cooksD[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          if(length(influential)==0){
            influential = NaN
          }
          if(!is.nan(influential)){
            influential = influential[match(max(Distances),Distances)]
            frame_new = frame
            frame_new$y[match(max(influential),frame_new$y)] = NaN
            frame_new$y[frame_new$y%in%influential] = NaN
            Filter_frame = !is.nan(frame_new$y)
            frame_new = frame_new[!is.nan(frame_new$y),]
            model_new = lm(formula=y~0+x,data = frame_new)
            r2 = round(summary(model_new)$adj.r.squared,digits = 2)
            if(r2<r2_initial){
              redo_calib = F
              r2 = r2_initial
            } else {
              linear_regression_cytotox = model_new
              if(r2 >= 0.8){
                redo_calib = F
              }
              frame = frame_new
              Linear_cyt_frame = dplyr::filter(Linear_cyt_frame,Filter_frame)
              Cytotoxicity_linear[!Filter_frame] = NA
              Concis = frame[!is.na(frame[,2]),1]
              Concis = unique(Concis)
              linear_regression_cytotox = model_new
              r2_initial = r2
              influentials = append(influentials,influential)
            }
          } else {
            redo_calib = F
          }
        }
      }
      
      Cytotoxicity_linear_prism = asteriks(Cytotoxicity_linear_prism,Cytotoxicity_linear)
      
      slope = as.numeric(linear_regression_cytotox$coefficients[1])
      SE_slope = sqrt(diag(vcov(linear_regression_cytotox)))/nobs(linear_regression_cytotox)
      
      plot(frame$y~frame$x,type="p")
      abline(linear_regression_cytotox)
      
    }
  }
  
  if(is.na(linear_regression_cytotox)|max(frame$x,na.rm=T)<(10/slope)|slope<0){
    IC10_linear = "no cytotoxicity"
    SE_IC10_linear = ""
    R2_IC10_linear = ""
  } else {
    IC10_linear = 10/slope
    SE_IC10_linear = 10/slope^2*SE_slope
    R2_IC10_linear = r2
    
    #Combine repeats
    replicates = as.data.frame(table(frame$x))[1,2]
    dose_new = frame$x[seq(1,length(frame$x),by=replicates)]
    
    if(length(dose_new[dose_new>IC10_linear])<2){
      #IC10_linear = "not enough datapoints > threshold"
      #SE_IC10_linear = ""
      #R2_IC10_linear = ""
      confidence_factor_cyt_linear = 0.5
    } else {
      confidence_factor_cyt_linear = 1
    }
    
  }
  
  #Log-logistic regression
  cytotoxicity_split = as.data.frame(matrix(Cytotoxicity,nrow=22))
  #use the 1 % quantile of the control as threshold, but round up to 1% (outliers could be problematic)
  precipitation_threshold = plyr::round_any(quantile(Cytotoxicity_control,probs = c(0.01),na.rm = T),1,f=ceiling)
  if(precipitation_threshold>-10){
    precipitation_threshold = -10
  }
  cytotoxicity_split = cytotoxicity_split<(precipitation_threshold)
  for(c in 1:ncol(cytotoxicity_split)){
    if(!is.na(match(T,cytotoxicity_split[,c]))&match(T,cytotoxicity_split[,c])<=4){
      k = match(F,cytotoxicity_split[match(T,cytotoxicity_split[,c]):length(cytotoxicity_split[,c]),c])-1
      cytotoxicity_split[1:k,c] = T
    } else {
      cytotoxicity_split[,c] = F
    }
  }
  Precipitation_filter_log = c(cytotoxicity_split)
  Cytotoxicity[Precipitation_filter_log] = NA
  Log_cyt_frame = data.frame("x"=Concentrations_log,
                             "y"=Cytotoxicity,
                             "Name"=rep(workframe$Measurement,each=22))
  frame = as.data.frame(cbind(Concentrations_log,Cytotoxicity))
  colnames(frame) = c("x","y")
  Log_cyt_frame = Log_cyt_frame[order(Log_cyt_frame$x,decreasing = T),]
  Log_cyt_frame = Log_cyt_frame[(!is.na(Log_cyt_frame$y)&Log_cyt_frame$y<=log_y_axis_limits_cytotoxicity),]
  Cytotoxicity[is.na(frame$y)&frame$y>log_y_axis_limits_cytotoxicity] = NA
  order_coordinates = order(frame$x,decreasing = T)
  frame = frame[(!is.na(frame$y)&frame$y<=log_y_axis_limits_cytotoxicity),]
  Cytotoxicity = Cytotoxicity[order_coordinates]
  frame = frame[order(frame$x,decreasing = T),]
  
  #bmad is the absolute deviation in assay - I use the control as reference
  bmad = sum(Cytotoxicity_control-mean(Cytotoxicity_control,na.rm = T),na.rm = T)/length(Cytotoxicity_control[!is.na(Cytotoxicity_control)])
  loglogistic_regression_cytotoxicity = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
  
  #Check AICs
  AICs = c(loglogistic_regression_cytotoxicity$cnst_aic,
           loglogistic_regression_cytotoxicity$gnls_aic,
           loglogistic_regression_cytotoxicity$hill_aic)
  
  model_names = c("cnst",
                  "gnls",
                  "hill")
  
  model = match(min(AICs,na.rm = T),AICs)
  
  if(is.na(model)){
    model = 1
  }
  
  if(model==2){
    model = 3
  }
  
  if(model!=1){
    if(abs(AICs[model]-AICs[1])<3){
      model = 1
    }
  }
  
  final_model_cyt = model_names[model]
  
  if(final_model_cyt=="cnst"){
    IC10_log = "no effect"
    SE_IC10_log = ""
    R2_IC10_log = ""
  } else if(final_model_cyt=="gnls"){
    maximum = max(frame$y,na.rm=T)
    index = match(maximum,frame$y)
    if(index<nrow(frame)){
      #frame is from highest to lowest concentration
      frame[c(frame$y[1:index]<(0.9*maximum),rep(F,index+1,length(frame$y))),] = NA
      Cytotoxicity[is.na(frame$x)] = NA
      Log_cyt_frame = Log_cyt_frame[!is.na(frame$x),]
      frame = frame[!is.na(frame$x),]
    }
    loglogistic_regression_cytotoxicity = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
    Differences = abs(frame$y-loglogistic_regression_cytotoxicity$hill_modl)
    if(length(Differences)==0){
      Differences = 0
    }
    if(length(Differences)>=5){
      Differences = Rosner_outlier_removal(Differences)
      frame = frame[!is.na(Differences),]
      Cytotoxicity[is.na(Differences)] = NA
      Log_cyt_frame = Log_cyt_frame[!is.na(Differences),]
      loglogistic_regression_cytotoxicity = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
      Top = loglogistic_regression_cytotoxicity$hill_tp
      Bottom = 0
      log_Slope = loglogistic_regression_cytotoxicity$hill_gw
      IC50_log = loglogistic_regression_cytotoxicity$hill_ga
      SEM_Slope = loglogistic_regression_cytotoxicity$hill_gw_sd/sqrt(nrow(Log_cyt_frame))
      SEM_IC50 = loglogistic_regression_cytotoxicity$hill_ga_sd/sqrt(nrow(Log_cyt_frame))
      SEM_IC10 = calculate_SE_ac10(Top,Bottom,10,IC50_log,SEM_IC50,log_Slope,SEM_Slope)
      
      if(is.nan(SEM_IC10)){
        fit = tryCatch(drc::drm(y~x, data=Log_cyt_frame, fct=LL.4(fixed = c(NA,0,NA,NA),names = c("Slope", "Lower Limit", "Upper Limit", "AC50")), logDose = 10),error=function(e){fit=NA})
        if(is.na(fit)){
          SEM_IC10 = "no error calculation possible"
        } else {
          SEM_IC10 = as.numeric(drc::ED(fit,respLev=10)[2])
        }
      }
      
      L0 = exp((AICs[1]-(2*1 + (2*1*(1+1)/(nrow(Log_cyt_frame)-1-1))))/(-2))
      Lm = exp((AICs[3]-(2*4 + (2*4*(4+1)/(nrow(Log_cyt_frame)-4-1))))/(-2))
      
      if(L0==0){
        L0 = 1e-1000
      }
      CoxSnellR2 = 1-((L0/Lm)^(2/nrow(frame)))
      upper_bond = 1-L0^(2/nrow(frame))
      corrected_CoxSnellR2 = CoxSnellR2/upper_bond
      
      plot(frame$y~frame$x,type="p")
      lines(frame$x,loglogistic_regression_cytotoxicity$hill_modl,type="l",col="red") 
      IC10_log = tcpl::tcplHillACXX(10/Top*100,tp=Top,ga=IC50_log,gw=log_Slope)
    } else {
      IC10_log = "no cytotoxicity"
      SE_IC10_log = ""
      R2_IC10_log = ""
    }
    
    
  } else {
    Differences = abs(frame$y-loglogistic_regression_cytotoxicity$hill_modl)
    if(length(Differences)==0){
      Differences = 0
    }
    if(length(Differences)>=5){
      Differences = Rosner_outlier_removal(Differences)
      frame = frame[!is.na(Differences),]
      Cytotoxicity[is.na(Differences)] = NA
      Log_cyt_frame = Log_cyt_frame[!is.na(Differences),]
      loglogistic_regression_cytotoxicity = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
      Top = loglogistic_regression_cytotoxicity$hill_tp
      Bottom = 0
      log_Slope = loglogistic_regression_cytotoxicity$hill_gw
      IC50_log = loglogistic_regression_cytotoxicity$hill_ga
      SEM_Slope = loglogistic_regression_cytotoxicity$hill_gw_sd/sqrt(nrow(Log_cyt_frame))
      SEM_IC50 = loglogistic_regression_cytotoxicity$hill_ga_sd/sqrt(nrow(Log_cyt_frame))
      SEM_IC10 = calculate_SE_ac10(Top,Bottom,10,IC50_log,SEM_IC50,log_Slope,SEM_Slope)
      
      if(is.nan(SEM_IC10)){
        fit = tryCatch(drc::drm(y~x, data=Log_cyt_frame, fct=LL.4(fixed = c(NA,0,NA,NA),names = c("Slope", "Lower Limit", "Upper Limit", "AC50")), logDose = 10),error=function(e){fit=NA})
        if(is.na(fit)){
          SEM_IC10 = "no error calculation possible"
        } else {
          SEM_IC10 = as.numeric(drc::ED(fit,respLev=10)[2])
        }
      }
      
      L0 = exp((AICs[1]-(2*1 + (2*1*(1+1)/(nrow(Log_cyt_frame)-1-1))))/(-2))
      Lm = exp((AICs[3]-(2*4 + (2*4*(4+1)/(nrow(Log_cyt_frame)-4-1))))/(-2))
      
      if(L0==0){
        L0 = 1e-1000
      }
      CoxSnellR2 = 1-((L0/Lm)^(2/nrow(frame)))
      upper_bond = 1-L0^(2/nrow(frame))
      corrected_CoxSnellR2 = CoxSnellR2/upper_bond
      
      plot(frame$y~frame$x,type="p")
      lines(frame$x,loglogistic_regression_cytotoxicity$hill_modl,type="l",col="red")
      
      IC10_log = tcpl::tcplHillACXX(10/Top*100,tp=Top,ga=IC50_log,gw=log_Slope)
    } else {
      IC10_log = "no cytotoxicity"
      SE_IC10_log = ""
      R2_IC10_log = ""
    }
    
  }
  
  Cytotoxicity_prism = asteriks(Cytotoxicity_prism,Cytotoxicity)
  Cytotoxicity_copy = Cytotoxicity
  for(v in 1:length(order_coordinates)){
    Cytotoxicity[order_coordinates[v]] = Cytotoxicity_copy[v]
  }
  
  if(is.nan(IC10_log)|!is.numeric(IC10_log)|max(frame$x,na.rm=T)<IC10_log){
    IC10_log = "no cytotoxicity"
    SE_IC10_log = ""
    R2_IC10_log = ""
  } else {
    IC10_log = IC10_log
    SE_IC10_log = SEM_IC10
    R2_IC10_log = round(corrected_CoxSnellR2,digits = 2)
    
    #Combine repeats
    replicates = as.data.frame(table(frame$x))[1,2]
    dose_new = frame$x[seq(1,length(frame$x),by=replicates)]
    
    if(length(dose_new[dose_new>IC10_log])<2){
      #IC10_log = "not enough datapoints > threshold"
      #SE_IC10_log = ""
      #R2_IC10_log = ""
      confidence_factor_cyt_log = 0.5
    } else {
      confidence_factor_cyt_log = 1
    }
  }
  
  #Effect
  #t_test = t.test(Effect,Effect_control,alternative = "greater")
  upper_control = quantile(Effect_control,probs=c(1),na.rm=T)
  #lower_control = quantile(Effect_control,probs=c(0.05),na.rm=T)
  replicates = as.data.frame(table(Concentrations_linear))[1,2]
  control_threshold = replicates*2
  above = length(Effect[Effect>upper_control])
  
  f_test = var.test(Effect,Effect_control)
  f_test_p = f_test$p.value
  #below = length(Effect[Effect<lower_control])
  
  if(above<control_threshold|f_test_p>0.05){
    Effect_dev_Control = F
  } else {
    Effect_dev_Control = T
  }
  #Linear Regression
  Effect_linear[Effect_linear>linear_cutoff_effect] = NA
  Effect_linear[Precipitation_filter_linear] = NA
  frame = as.data.frame(cbind(Concentrations_linear,Effect_linear))
  Linear_effect_frame = data.frame("x"=Concentrations_linear,
                                   "y"=Effect_linear,
                                   "Name"=rep(workframe$Measurement,each=22))
  Linear_effect_frame = Linear_effect_frame[!is.na(frame$Effect_linear),]
  Effect_linear[is.na(frame$Effect_linear)] = NA
  frame = frame[!is.na(frame$Effect_linear),]
  colnames(frame) = c("x","y")
  if(!is.numeric(IC10_linear)|IC10_linear<0){
    IC10_test = Inf
  } else {
    IC10_test = IC10_linear
  }
  if(cytotoxicity_cutoff==T){
    Linear_effect_frame = Linear_effect_frame[frame$x<=IC10_test,]
    Effect_linear[frame$x>=IC10_test] = NA
    frame = frame[frame$x<=IC10_test,]
  }
  if(nrow(frame)>5){
    linear_regression_effect = tryCatch(lm(data=frame,y~0+x),error=function(e){linear_regression_effect=NA})
    influentials = c()
    if(!is.na(linear_regression_effect)){
      r2 = round(summary(linear_regression_effect)$adj.r.squared,digits = 2)
      Concis = frame[!is.na(frame[,2]),1]
      Concis = unique(Concis)
      if(r2 < 0.8 & length(Concis)>5){
        redo_calib = T
        r2_initial = r2
        while(redo_calib == T){
          if(length(Concis)==4){
            break
          }
          cooksD = cooks.distance(linear_regression_effect)
          influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          Distances = cooksD[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
          if(length(influential)==0){
            influential = NaN
          }
          if(!is.nan(influential)){
            influential = influential[match(max(Distances),Distances)]
            frame_new = frame
            frame_new$y[match(max(influential),frame_new$y)] = NaN
            frame_new$y[frame_new$y%in%influential] = NaN
            Filter_frame = !is.nan(frame_new$y)
            frame_new = frame_new[!is.nan(frame_new$y),]
            model_new = lm(formula=y~0+x,data = frame_new)
            r2 = round(summary(model_new)$adj.r.squared,digits = 2)
            if(r2<r2_initial){
              redo_calib = F
              r2 = r2_initial
            } else {
              linear_regression_effect = model_new
              if(r2 >= 0.8){
                redo_calib = F
              }
              frame = frame_new
              Linear_effect_frame = dplyr::filter(Linear_effect_frame,Filter_frame)
              Effect_linear[!Filter_frame] = NA
              Concis = frame[!is.na(frame[,2]),1]
              Concis = unique(Concis)
              linear_regression_effect = model_new
              r2_initial = r2
              influentials = append(influentials,influential)
            }
          } else {
            redo_calib = F
          }
        }
      }
      
      #Effect_linear[Effect_linear%in%influentials] = NA
      #Effect_linear_prism = asteriks(Effect_linear_prism,Effect_linear)
      
      slope = as.numeric(linear_regression_effect$coefficients[1])
      SE_slope = sqrt(diag(vcov(linear_regression_effect)))/nobs(linear_regression_effect)
      
      plot(frame$y~frame$x,type="p")
      abline(linear_regression_effect)
      
    } else {
      Effect_linear_prism = asteriks(Effect_linear_prism,Effect_linear)
    }
    
    if(!is.na(linear_regression_effect)){
      #Detect plateaus and cut since its start of log-logistic fit
      unique_concentrations = nrow(as.data.frame(table(frame$x)))
      frame_plateau = as.data.frame(matrix(nrow=unique_concentrations,ncol=2))
      colnames(frame_plateau) = colnames(frame)
      frame_plateau$x = as.data.frame(table(frame$x))$Var1
      for(f in 1:nrow(frame_plateau)){
        frame_plateau$y[f] = mean(frame$y[frame$x==frame_plateau$x[f]],na.rm=T)
      }
      if(length(frame_plateau$y[frame_plateau$y>10])>2){
        frame_plateau = frame_plateau[order(frame_plateau$x,decreasing = F),]
        frame_plateau$x = as.numeric(as.character(frame_plateau$x))
        model = tryCatch(drc::drm(y~x,data=frame_plateau,fct=drc::MM.2()),
                         error=function(e){model=NA})
        if(!is.na(model)){
          model_linear = lm(y~0+x,data=frame_plateau)
          predi = predict(model,newdata = frame_plateau)
          predi_linear = predict(model_linear)
          plot(frame_plateau$y~frame_plateau$x,type="p")
          lines(frame_plateau$x,predi,col="red")
          lines(frame_plateau$x,predi_linear,col="blue")
          MM_AIC = round(as.numeric(2*3-2*(logLik(model))),digits = 0)
          Linear_AIC = round(as.numeric(2*1-2*(logLik(model_linear))),digits = 0)
          while(MM_AIC<(Linear_AIC-4)&length(frame_plateau$y)>4){
            frame_plateau = frame_plateau[-nrow(frame_plateau),]
            model = tryCatch(drc::drm(y~x,data=frame_plateau,fct=drc::MM.2()),
                             error=function(e){model=NA})
            if(is.na(model)){
              break
            }
            model_linear = lm(y~0+x,data=frame_plateau)
            predi = predict(model,newdata = frame_plateau)
            predi_linear = predict(model_linear)
            plot(frame_plateau$y~frame_plateau$x,type="p")
            lines(frame_plateau$x,predi,col="red")
            lines(frame_plateau$x,predi_linear,col="blue")
            MM_AIC = round(as.numeric(2*3-2*(logLik(model))),digits = 0)
            Linear_AIC = round(as.numeric(2*1-2*(logLik(model_linear))),digits = 0)
          }
          highest_conc = max(frame_plateau$x,na.rm=T)
        } else {
          highest_conc = max(frame$x,na.rm=T)
        }
      } else {
        highest_conc = max(frame$x,na.rm=T)
      }
      
      Effect_linear[frame$x>=highest_conc] = NA
      Linear_effect_frame = Linear_effect_frame[frame$x<=highest_conc,]
      frame = frame[frame$x<=highest_conc,]
      
      linear_regression_effect = tryCatch(lm(data=frame,y~0+x,weights = x^2),error=function(e){linear_regression_effect=NA})
      influentials = c()
      if(!is.na(linear_regression_effect)){
        r2 = round(summary(linear_regression_effect)$adj.r.squared,digits = 2)
        Concis = frame[!is.na(frame[,2]),1]
        Concis = unique(Concis)
        if(r2 < 0.8 & length(Concis)>5){
          redo_calib = T
          r2_initial = r2
          while(redo_calib == T){
            if(length(Concis)==4){
              break
            }
            cooksD = cooks.distance(linear_regression_effect)
            influential <- frame$y[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
            Distances = cooksD[cooksD > (3 * mean(cooksD, na.rm = TRUE))]
            if(length(influential)==0){
              influential = NaN
            }
            if(!is.nan(influential)){
              influential = influential[match(max(Distances),Distances)]
              frame_new = frame
              frame_new$y[match(max(influential),frame_new$y)] = NaN
              frame_new$y[frame_new$y%in%influential] = NaN
              Filter_frame = !is.nan(frame_new$y)
              frame_new = frame_new[!is.nan(frame_new$y),]
              model_new = lm(formula=y~0+x,data = frame_new)
              r2 = round(summary(model_new)$adj.r.squared,digits = 2)
              if(r2<r2_initial){
                redo_calib = F
                r2 = r2_initial
              } else {
                linear_regression_effect = model_new
                if(r2 >= 0.8){
                  redo_calib = F
                }
                frame = frame_new
                Linear_effect_frame = dplyr::filter(Linear_effect_frame,Filter_frame)
                Effect_linear[!Filter_frame] = NA
                Concis = frame[!is.na(frame[,2]),1]
                Concis = unique(Concis)
                linear_regression_effect = model_new
                r2_initial = r2
                influentials = append(influentials,influential)
              }
            } else {
              redo_calib = F
            }
          }
        }
        
        #Effect_linear[Effect_linear%in%influentials] = NA
        Effect_linear_prism = asteriks(Effect_linear_prism,Effect_linear)
        
        slope = as.numeric(linear_regression_effect$coefficients[1])
        SE_slope = sqrt(diag(vcov(linear_regression_effect)))/nobs(linear_regression_effect)
        
        plot(frame$y~frame$x,type="p")
        abline(linear_regression_effect)
        
      }
      
      if(is.na(linear_regression_effect)|max(frame$x,na.rm=T)<(10/slope)|slope<0){
        EC10_linear = "no effect"
        SE_EC10_linear = ""
        R2_EC10_linear = ""
      } else {
        EC10_linear = 10/slope
        SE_EC10_linear = 10/slope^2*SE_slope
        R2_EC10_linear = r2
        
        if(is.numeric(IC10_linear)&IC10_linear<EC10_linear){
          EC10_linear = "masked by cytotoxicity"
          SE_EC10_linear = ""
          R2_EC10_linear = ""
        }
      }
      #Combine repeats
      replicates = as.data.frame(table(frame$x))[1,2]
      dose_new = frame$x[seq(1,length(frame$x),by=replicates)]
      
      if(length(dose_new[dose_new>=EC10_linear])<2){
        #EC10_linear = "not enough datapoints > threshold"
        #SE_EC10_linear = ""
        #R2_EC10_linear = ""
        confidence_factor_effect_linear = 0.5
      } else {
        confidence_factor_effect_linear = 1
      }
    }
  } else {
    EC10_linear = "masked by cytotoxicity"
    SE_EC10_linear = ""
    R2_EC10_linear = ""
  }
  
  #Log-logistic regression
  Effect[Precipitation_filter_log] = NA
  frame = as.data.frame(cbind(Concentrations_log,Effect))
  colnames(frame) = c("x","y")
  Log_effect_frame = data.frame("x"=Concentrations_log,
                                "y"=Effect,
                                "Name"=rep(workframe$Measurement,each=22))
  if(is.numeric(IC10_log)){
    if((10^IC10_log)<0){
      IC10_test = Inf
    } else {
      IC10_test = IC10_log
    }
  } else {
    IC10_test = Inf
  }
  
  if(cytotoxicity_cutoff==T){
    Log_effect_frame = Log_effect_frame[frame$x<=IC10_test,]
    Effect[frame$x>=IC10_test] = NA
    order_coordinates = order(frame$x,decreasing = T)
    frame = frame[frame$x<=IC10_test,]
  }
  
  if(nrow(frame)>5){
    Log_effect_frame = Log_effect_frame[(!is.na(Log_effect_frame$y)&Log_effect_frame$y<=log_y_axis_limits_effects),]
    Effect[(is.na(frame$y)&frame$y>log_y_axis_limits_effects)] = NA
    frame = frame[(!is.na(frame$y)&frame$y<=log_y_axis_limits_effects),]
    Log_effect_frame = Log_effect_frame[order(frame$x,decreasing = T),]
    Effect = Effect[order_coordinates]
    frame = frame[order(frame$x,decreasing = T),]
    
    #bmad is the absolute deviation in assay - I use the control as reference
    bmad = sum(Effect_control-mean(Effect_control,na.rm = T),na.rm = T)/length(Effect_control[!is.na(Effect_control)])
    loglogistic_regression_Effect = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
    
    #Check AICs
    AICs = c(loglogistic_regression_Effect$cnst_aic,
             loglogistic_regression_Effect$gnls_aic,
             loglogistic_regression_Effect$hill_aic)
    
    model_names = c("cnst",
                    "gnls",
                    "hill")
    
    model = match(min(AICs,na.rm = T),AICs)
    
    if(is.na(model)){
      model = 1
    }
    
    if(model==2){
      if(abs(AICs[model]-AICs[3])<3){
        model = 3
      }
    }
    
    if(model!=1){
      if(abs(AICs[model]-AICs[1])<3){
        model = 1
      }
    }
    
    final_model_effect = model_names[model]
    
    if(final_model_effect=="cnst"){
      EC10_log = "no effect"
      SE_EC10_log = ""
      R2_EC10_log = ""
    } else if(final_model_effect=="gnls"){
      maximum = max(frame$y,na.rm=T)
      index = match(maximum,frame$y)
      if(index<nrow(frame)){
        #frame is from highest to lowest concentration
        frame[c(frame$y[1:index]<(0.9*maximum),rep(F,index+1,length(frame$y))),] = NA
        Effect[is.na(frame$x)] = NA
        Log_effect_frame = Log_effect_frame[!is.na(frame$x),]
        frame = frame[!is.na(frame$x),]
      }
      loglogistic_regression_Effect = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
      Differences = abs(frame$y-loglogistic_regression_Effect$hill_modl)
      if(length(Differences)==0){
        Differences = 0
      }
      if(length(Differences)>=5){
        Differences = Rosner_outlier_removal(Differences)
        frame = frame[!is.na(Differences),]
        Effect[is.na(Differences)] = NA
        Log_effect_frame = Log_effect_frame[!is.na(Differences),]
        loglogistic_regression_Effect = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
        Top = loglogistic_regression_Effect$hill_tp
        Bottom = 0
        log_Slope = loglogistic_regression_Effect$hill_gw
        EC50_log = loglogistic_regression_Effect$hill_ga
        SEM_Slope = loglogistic_regression_Effect$hill_gw_sd/sqrt(nrow(Log_effect_frame))
        SEM_EC50 = loglogistic_regression_Effect$hill_ga_sd/sqrt(nrow(Log_effect_frame))
        SEM_EC10 = calculate_SE_ac10(Top,Bottom,10,EC50_log,SEM_EC50,log_Slope,SEM_Slope)
        
        if(is.nan(SEM_EC10)){
          fit = tryCatch(drc::drm(y~x, data=Log_effect_frame, fct=LL.4(fixed = c(NA,0,NA,NA),names = c("Slope", "Lower Limit", "Upper Limit", "AC50")), logDose = 10),error=function(e){fit=NA})
          if(is.na(fit)){
            SEM_EC10 = "no error calculation possible"
          } else {
            SEM_EC10 = as.numeric(drc::ED(fit,respLev=10)[2])
          }
        }
        
        L0 = exp((AICs[1]-(2*1 + (2*1*(1+1)/(nrow(Log_effect_frame)-1-1))))/(-2))
        Lm = exp((AICs[3]-(2*4 + (2*4*(4+1)/(nrow(Log_effect_frame)-4-1))))/(-2))
        
        if(L0==0){
          L0 = 1e-1000
        }
        CoxSnellR2 = 1-((L0/Lm)^(2/nrow(frame)))
        upper_bond = 1-L0^(2/nrow(frame))
        corrected_CoxSnellR2 = CoxSnellR2/upper_bond
        
        plot(frame$y~frame$x,type="p")
        lines(frame$x,loglogistic_regression_Effect$hill_modl,type="l",col="red") 
        EC10_log = tcpl::tcplHillACXX(10/Top*100,tp=Top,ga=EC50_log,gw=log_Slope)
      } else {
        EC10_log = "no effect"
        SE_EC10_log = ""
        R2_EC10_log = ""
      }
      
      
    } else {
      Differences = abs(frame$y-loglogistic_regression_Effect$hill_modl)
      if(length(Differences)==0){
        Differences = 0
      }
      if(length(Differences)>=5){
        Differences = Rosner_outlier_removal(Differences)
        frame = frame[!is.na(Differences),]
        Effect[is.na(Differences)] = NA
        Log_effect_frame = Log_effect_frame[!is.na(Differences),]
        loglogistic_regression_Effect = tcpl::tcplFit(logc=frame$x,resp=frame$y,bmad=bmad,force.fit = T)
        Top = loglogistic_regression_Effect$hill_tp
        Bottom = 0
        log_Slope = loglogistic_regression_Effect$hill_gw
        EC50_log = loglogistic_regression_Effect$hill_ga
        
        #We have the standard deviation but want to calculate the standard error
        SEM_Slope = loglogistic_regression_Effect$hill_gw_sd/sqrt(nrow(frame))
        SEM_EC50 = loglogistic_regression_Effect$hill_ga_sd/sqrt(nrow(frame))
        SEM_EC10 = calculate_SE_ac10(Top,Bottom,10,EC50_log,SEM_EC50,log_Slope,SEM_Slope)
        
        if(is.nan(SEM_EC10)){
          SEM_EC10 = "no error calculation possible"
        }
        
        L0 = exp((AICs[1]-(2*1 + (2*1*(1+1)/(nrow(Log_effect_frame)-1-1))))/(-2))
        Lm = exp((AICs[3]-(2*4 + (2*4*(4+1)/(nrow(Log_effect_frame)-4-1))))/(-2))
        
        if(L0==0){
          L0 = 1e-1000
        }
        CoxSnellR2 = 1-((L0/Lm)^(2/nrow(frame)))
        upper_bond = 1-L0^(2/nrow(frame))
        corrected_CoxSnellR2 = CoxSnellR2/upper_bond
        
        plot(frame$y~frame$x,type="p")
        lines(frame$x,loglogistic_regression_Effect$hill_modl,type="l",col="red")
        
        EC10_log = tcpl::tcplHillACXX(10/Top*100,tp=Top,ga=EC50_log,gw=log_Slope)
        abline(v=EC10_log)
        abline(h=10)
        
      } else {
        EC10_log = "no effect"
        EC50_log = "no effect"
        SE_EC10_log = ""
        R2_EC10_log = ""
      }
      
    }
    
    Effect_prism = asteriks(Effect_prism,Effect)
    Effect_copy = Effect
    for(v in 1:length(order_coordinates)){
      Effect[order_coordinates[v]] = Effect_copy[v]
    }
    
    if(is.nan(EC10_log)|!is.numeric(EC10_log)|10^(max(frame$x,na.rm=T))<EC10_log){
      EC10_log = "no effect"
      SE_EC10_log = ""
      R2_EC10_log = ""
    } else {
      EC10_log = EC10_log
      SE_EC10_log = SEM_EC10
      R2_EC10_log = round(corrected_CoxSnellR2,digits = 2)
      if(is.numeric(IC10_log)&IC10_log<EC10_log){
        EC10_log = "masked by cytotoxicity"
        SE_EC10_log = ""
        R2_EC10_log = ""
      }
    }
    #Combine repeats
    replicates = as.data.frame(table(frame$x))[1,2]
    dose_new = frame$x[seq(1,length(frame$x),by=replicates)]
    
    if(length(dose_new[dose_new>EC10_log])<2){
      #EC10_log = "not enough data points > threshold"
      #SE_EC10_log = ""
      #R2_EC10_log = ""
      confidence_factor_effect_log = 0.5
    } else {
      confidence_factor_effect_log = 1
    }
  } else {
    EC10_log = "masked by cytotoxicity"
    SE_EC10_log = ""
    R2_EC10_log = ""
  }

  #Results in Summary
  d = NA
  Summary_results = Summary_results
  Final_ID = Sample_ID
  #Final_ID = unlist(stringr::str_split(Sample_ID,"-"))[1]
  #For PDFs
  Final_ID = unlist(stringr::str_replace(Final_ID,"/","_"))
  Final_ID = unlist(stringr::str_replace(Final_ID,":","_"))
  if(Final_ID%in%Summary_results$Sample_ID){
    d = match(Final_ID,Summary_list$Sample_ID)
    Results_frame = Summary_results[d,]
  }
  
  Results_frame = as.data.frame(matrix(nrow=1,ncol = ncol(Summary_results)))
  colnames(Results_frame) = colnames(Summary_results)
  
  Results_frame$Sample_ID = Final_ID
  Results_frame$Unit = workframe$Unit[1]
  
  #Effects
  Results_frame$SHSY5Y_Effect_dev_Control = Effect_dev_Control
  Results_frame$SHSY5Y_linmodel_EC10 = EC10_linear
  Results_frame$SE_SHSY5Y_linmodel_EC10 = SE_EC10_linear
  Results_frame$R2_SHSY5Y_linmodel_EC10 = R2_EC10_linear
  if(!is.na(as.numeric(EC10_linear))&!is.na(as.numeric(SE_EC10_linear))){
    Results_frame$SHSY5Y_linmodel_confidence_EC10 = plyr::round_any((1-abs(as.numeric(SE_EC10_linear/as.numeric(EC10_linear))))*confidence_factor_effect_linear*R2_EC10_linear,0.01,ceiling)
  } else {
    Results_frame$SHSY5Y_linmodel_confidence_EC10 = NA
  }
  Results_frame$SHSY5Y_logmodel_EC10 = EC10_log
  Results_frame$SE_SHSY5Y_logmodel_EC10 = SE_EC10_log
  Results_frame$R2_SHSY5Y_logmodel_EC10 = R2_EC10_log
  if(!is.na(as.numeric(EC10_log))&!is.na(as.numeric(SE_EC10_log))){
    Results_frame$SHSY5Y_logmodel_confidence_EC10 = plyr::round_any((1-abs(as.numeric(SE_EC10_log/as.numeric(EC10_log))))*confidence_factor_effect_log*R2_EC10_log,0.01,ceiling)
  } else {
    Results_frame$SHSY5Y_logmodel_confidence_EC10 = NA
  }
  Results_frame$SHSY5Y_logmodel_linear_EC10 = 10^(as.numeric(EC10_log))
  
  #Cytotoxicity
  Results_frame$SHSY5Y_Cytotoxicity_dev_Control = Cytotoxicity_dev_Control
  Results_frame$SHSY5Y_linmodel_IC10 = IC10_linear
  Results_frame$SE_SHSY5Y_linmodel_IC10 = SE_IC10_linear
  Results_frame$R2_SHSY5Y_linmodel_IC10 = R2_IC10_linear
  if(!is.na(as.numeric(IC10_linear))&!is.na(as.numeric(SE_IC10_linear))){
    Results_frame$SHSY5Y_linmodel_confidence_IC10 = plyr::round_any((1-abs(as.numeric(SE_IC10_linear/as.numeric(IC10_linear))))*confidence_factor_cyt_linear*R2_IC10_linear,0.01,ceiling)
  } else {
    Results_frame$SHSY5Y_linmodel_confidence_IC10 = NA
  }
  Results_frame$SHSY5Y_logmodel_IC10 = IC10_log
  Results_frame$SE_SHSY5Y_logmodel_IC10 = SE_IC10_log
  Results_frame$R2_SHSY5Y_logmodel_IC10 = R2_IC10_log
  Results_frame$SHSY5Y_logmodel_linear_IC10 = 10^(as.numeric(IC10_log))
  if(!is.na(as.numeric(IC10_log))&!is.na(as.numeric(SE_IC10_log))){
    Results_frame$SHSY5Y_logmodel_confidence_IC10 = plyr::round_any((1-abs(as.numeric(SE_IC10_log/as.numeric(IC10_log))))*confidence_factor_cyt_log*R2_IC10_log,0.01,ceiling)
  } else {
    Results_frame$SHSY5Y_logmodel_confidence_IC10 = NA
  }
  
  if(!is.na(d)){
    Summary_results[d,] = Results_frame
  } else {
    Summary_results = rbind(Summary_results,Results_frame)
  }
  
  #-----------------------------------------------------------------------------
  
  graph_path = "./Plots"
  
  Prism = openxlsx::createWorkbook()
  
  #Prism
  #Sheet 1
  Linear_all = as.data.frame(matrix(nrow = length(Effect_all),ncol = (nrow(workframe)*2+1)))
  colnames(Linear_all) = c(paste0("Concentration (",workframe$Unit[1],")"),create_names(workframe))
  Linear_all[,1] = Concentrations_linear
  Start = 1
  Interval = 22
  Starts = seq(1,10000,22)
  Ends = seq(22,10000,22)
  Intervals = as.data.frame(matrix(nrow = nrow(workframe)*2,ncol = 2))
  colnames(Intervals) = c("Start","End")
  Intervals$Start = rep(Starts[1:nrow(workframe)],each=2)
  Intervals$End = rep(Ends[1:nrow(workframe)],each=2)
  for(h in 1:nrow(Intervals)){
    if(!even(h)){
      Linear_all[Intervals$Start[h]:Intervals$End[h],(h+1)] = Cytotoxicity_all[Intervals$Start[h]:Intervals$End[h]]
    } else {
      Linear_all[Intervals$Start[h]:Intervals$End[h],(h+1)] = Effect_all[Intervals$Start[h]:Intervals$End[h]]
    }
  }
  openxlsx::addWorksheet(Prism,"Linear_all")
  openxlsx::writeData(Prism,"Linear_all",Linear_all)
  
  #Sheet 2
  Log_all = Linear_all
  Log_all[,1] = Concentrations_log
  openxlsx::addWorksheet(Prism,"Log_all")
  openxlsx::writeData(Prism,"Log_all",Log_all)
  
  #Sheet3
  Linear_cutoff = as.data.frame(matrix(nrow = length(Effect_linear),ncol = (nrow(workframe)*2+1)))
  colnames(Linear_cutoff) = c(paste0("Concentration (",workframe$Unit[1],")"),create_names(workframe))
  Linear_cutoff[,1] = Concentrations_linear
  Start = 1
  Interval = 22
  Starts = seq(1,10000,22)
  Ends = seq(22,10000,22)
  Intervals = as.data.frame(matrix(nrow = nrow(workframe)*2,ncol = 2))
  colnames(Intervals) = c("Start","End")
  Intervals$Start = rep(Starts[1:nrow(workframe)],each=2)
  Intervals$End = rep(Ends[1:nrow(workframe)],each=2)
  for(h in 1:nrow(Intervals)){
    if(!even(h)){
      Linear_cutoff[Intervals$Start[h]:Intervals$End[h],(h+1)] = Cytotoxicity_linear_prism[Intervals$Start[h]:Intervals$End[h]]
    } else {
      Linear_cutoff[Intervals$Start[h]:Intervals$End[h],(h+1)] = Effect_linear_prism[Intervals$Start[h]:Intervals$End[h]]
    }
  }
  openxlsx::addWorksheet(Prism,"Linear_cutoff")
  openxlsx::writeData(Prism,"Linear_cutoff",Linear_cutoff)
  
  #Sheet 4
  Log_w_excluded = as.data.frame(matrix(nrow = length(Effect),ncol = (nrow(workframe)*2+1)))
  colnames(Log_w_excluded) = c(paste0("Concentration (",workframe$Unit[1],")"),create_names(workframe))
  Log_w_excluded[,1] = Concentrations_log
  Start = 1
  Interval = 22
  Starts = seq(1,10000,22)
  Ends = seq(22,10000,22)
  Intervals = as.data.frame(matrix(nrow = nrow(workframe)*2,ncol = 2))
  colnames(Intervals) = c("Start","End")
  Intervals$Start = rep(Starts[1:nrow(workframe)],each=2)
  Intervals$End = rep(Ends[1:nrow(workframe)],each=2)
  for(h in 1:nrow(Intervals)){
    if(!even(h)){
      Log_w_excluded[Intervals$Start[h]:Intervals$End[h],(h+1)] = Cytotoxicity_prism[Intervals$Start[h]:Intervals$End[h]]
    } else {
      Log_w_excluded[Intervals$Start[h]:Intervals$End[h],(h+1)] = Effect_prism[Intervals$Start[h]:Intervals$End[h]]
    }
  }
  openxlsx::addWorksheet(Prism,"Log_w_excluded")
  openxlsx::writeData(Prism,"Log_w_excluded",Log_w_excluded)
  
  #Controls
  Effect_control_prism = asteriks(Effect_control_prism,Effect_control)
  Cytotoxicity_control_prism = asteriks(Cytotoxicity_control_prism,Cytotoxicity_control)
  
  #Sheet 5
  Controls = as.data.frame(matrix(nrow = length(Cytotoxicity_control),ncol = (nrow(workframe)*2+1)))
  colnames(Controls) = c(paste0("Concentration (",workframe$Unit[1],")"),create_names(workframe))
  Controls[,1] = 10^log_concentration_unexposed_controls
  Start = 1
  Inverval = 32 #one control is set for 2*16 samples
  Intervals = as.data.frame(matrix(nrow = nrow(workframe)*2,ncol = 2))
  colnames(Intervals) = c("Start","End")
  Starts = seq(1,10000,32)
  Ends = seq(32,10000,32)
  Intervals$Start = rep(Starts[1:nrow(workframe)],each=2)
  Intervals$End = rep(Ends[1:nrow(workframe)],each=2)
  for(h in 1:nrow(Intervals)){
    if(!even(h)){
      Controls[Intervals$Start[h]:Intervals$End[h],(h+1)] = Cytotoxicity_control_prism[Intervals$Start[h]:Intervals$End[h]]
    } else {
      Controls[Intervals$Start[h]:Intervals$End[h],(h+1)] = Effect_control_prism[Intervals$Start[h]:Intervals$End[h]]
    }
  }
  openxlsx::addWorksheet(Prism,"Controls")
  openxlsx::writeData(Prism,"Controls",Controls)
  
  #Sheet 6
  Controls_log = Controls
  colnames(Controls_log) = c(paste0("log Concentration (",workframe$Unit[1],")"),create_names(workframe))
  Controls[,1] = log_concentration_unexposed_controls
  openxlsx::addWorksheet(Prism,"Controls_log")
  openxlsx::writeData(Prism,"Controls_log",Controls_log)
  
  openxlsx::saveWorkbook(Prism,paste0(graph_path,"/",Final_ID,".xlsx"),overwrite = T)
  
  #-----------------------------------------------------------------------------
  
  #Visualization
  Effect_color = rgb(206/255,6/255,101/255,1)
  Cyto_color = rgb(7/255,126/255,151/255,1)
  #Linear regression
  Linear_all_frame = as.data.frame(matrix(nrow=length(Concentrations_linear),ncol=4))
  colnames(Linear_all_frame) = c("Concentration","Name","Cytotox","Effect")
  Linear_all_frame[,1] = Concentrations_linear
  Measis = as.data.frame(table(workframe$Measurement))
  
  Linear_all_frame[,2] = rep(workframe$Measurement,each=22)
  Linear_all_frame[,3] = Cytotoxicity_all
  Linear_all_frame[,4] = Effect_all
  
  p_left = ggplot(Linear_all_control_frame,aes(x=Concentration)) + 
    geom_point(data=Linear_all_control_frame,aes(y = Cytotox,shape = Name),fill=Cyto_color,size=8) + 
    geom_point(data=Linear_all_control_frame,aes(y = Effect,shape= Name),fill=Effect_color,size=8) +
    scale_y_continuous(
      name = "Cytotoxicity [%]", limits = c(min_y_axis_limits_cytotoxicity,log_y_axis_limits_cytotoxicity),position = "left") + theme(axis.ticks=element_blank()) + xlim((10^log_concentration_unexposed_controls)*0.99999,(10^log_concentration_unexposed_controls)*1.00001)
  
  p_right = ggplot(Linear_all_frame,aes(x=Concentration)) + 
    geom_point(data=Linear_all_frame,aes(y = Cytotox,shape = Name),fill=Cyto_color,size=8) + 
    geom_point(data=Linear_all_frame,aes(y = Effect,shape= Name),fill=Effect_color,size=8) +
    scale_y_continuous(
      name = "Neurite outgrowth inhibition [%]", limits = c(min_y_axis_limits_effects,log_y_axis_limits_effects),position="right")
  
  Linear_all_plot_left = p_left + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=32,
                                      color=Cyto_color,
                                      face="bold",
                                      angle=90)) + 
    scale_shape_prism(palette="filled")
  
  Linear_all_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.title.y.right = element_text(size=32,
                                            color=Effect_color,
                                            face="bold",
                                            angle=270)) + 
    xlab(paste0("Concentration (",workframe$Unit[1],")")) + 
    scale_shape_prism(palette="filled")
  
  Linear_all_plot = (Linear_all_plot_left|Linear_all_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Linear_all_plot
  
  #Fit
  if(is.na(linear_regression_cytotox)|nrow(Linear_cyt_frame)<3){
    linear_cyt_pred = 0
    Linear_cyt_frame = data.frame("x"=Concentrations_linear,
                                  "y"=Cytotoxicity_linear,
                                  "Name"=rep(workframe$Measurement,each=22))
  } else {
    linear_cyt_pred = predict(linear_regression_cytotox,Linear_cyt_frame)
  }
  
  p_left = ggplot(Linear_cyt_control_frame,aes(x=x)) + 
    geom_point(data=Linear_cyt_control_frame,aes(y = y,shape = Name),fill=Cyto_color,size=8) + 
    geom_line(aes(y=10),linetype="dashed",color=Cyto_color,linewidth=2) + 
    scale_y_continuous(
      name = "Cytotoxicity [%]",position = "left") + theme(axis.ticks=element_blank()) + xlim((10^log_concentration_unexposed_controls)*0.99999,(10^log_concentration_unexposed_controls)*1.00001) +
    coord_cartesian(ylim=c(min_y_axis_limits_cytotoxicity,linear_y_axis_limits_cytotoxicity),clip = "on")
  
  p_right = ggplot(Linear_cyt_frame,aes(x=x)) + 
    geom_point(data=Linear_cyt_frame,aes(y = y,shape = Name),fill=Cyto_color,size=8) + labs(y=NULL) +
    geom_line(aes(y=10),linetype="dashed",color=Cyto_color,linewidth=2) + 
    geom_line(aes(y = linear_cyt_pred),color=Cyto_color,linewidth=3) + coord_cartesian(ylim=c(min_y_axis_limits_cytotoxicity,linear_y_axis_limits_cytotoxicity),clip = "on")
  
  
  Linear_cyt_plot_left = p_left + theme_prism(base_size = 30) + ylab("Cytotoxicity [%]") + 
    theme(legend.position = "none",
          axis.title.y = element_text(size=32,
                                      color=Cyto_color,
                                      face="bold",
                                      angle=90),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_shape_prism(palette="filled")
  
  Linear_cyt_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none") + 
    xlab(paste0("Concentration (",workframe$Unit[1],")")) +
    scale_shape_prism(palette="filled") + labs(y=NULL) + guides(y = "none")
  
  Linear_cyt_plot = (Linear_cyt_plot_left|Linear_cyt_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Linear_cyt_plot
  #------------------
  if(is.na(linear_regression_effect)|nrow(Linear_effect_frame)<3){
    linear_effect_pred = 0
    Linear_effect_frame = data.frame("x"=Concentrations_linear,
                                     "y"=Effect_linear,
                                     "Name"=rep(workframe$Measurement,each=22))
  } else {
    linear_effect_pred = predict(linear_regression_effect,Linear_effect_frame)
  }
  
  p_left = ggplot(Linear_effect_control_frame,aes(x=x)) + 
    geom_point(data=Linear_effect_control_frame,aes(y = y,shape = Name),fill=Effect_color,size=8) + 
    geom_line(aes(y=10),linetype="dashed",color=Effect_color,linewidth=2) + 
    labs(y=NULL) + ylim((min(c(Linear_effect_frame$y,Linear_effect_control_frame$y),na.rm=T)-5),linear_y_axis_limits_cytotoxicity) + xlim((10^log_concentration_unexposed_controls)*0.99999,(10^log_concentration_unexposed_controls)*1.00001)
  
  p_right = ggplot(Linear_effect_frame,aes(x=x)) + 
    geom_point(data=Linear_effect_frame,aes(y = y,shape = Name),fill=Effect_color,size=8) +
    geom_line(aes(y=10),linetype="dashed",color=Effect_color,linewidth=2) + 
    geom_line(aes(y = linear_effect_pred),color=Effect_color,linewidth=3) +
    scale_y_continuous(
      name = "Neurite outgrowth inhibition [%]",position = "right") + coord_cartesian(ylim=c(min_y_axis_limits_effects,linear_y_axis_limits_cytotoxicity),clip = "on")
  
  
  Linear_effect_plot_left = p_left + theme_prism(base_size = 30) + labs(y=NULL) + guides(y="none") + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_shape_prism(palette="filled")
  
  Linear_effect_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.title.y.right = element_text(size=32,
                                            color=Effect_color,
                                            face="bold",
                                            angle=90)) + 
    xlab(paste0("Concentration (",workframe$Unit[1],")")) +
    scale_shape_prism(palette="filled")
  
  Linear_effect_plot = (Linear_effect_plot_left|Linear_effect_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Linear_effect_plot
  
  #Log-logistic
  #--------------------------
  Log_all_frame = Linear_all_frame
  Log_all_frame[,1] = Concentrations_log
  
  p_left = ggplot(Log_all_control_frame,aes(x=Concentration)) + 
    geom_point(data=Log_all_control_frame,aes(y = Cytotox,shape = Name),fill=Cyto_color,size=8) + 
    geom_point(data=Log_all_control_frame,aes(y = Effect,shape= Name),fill=Effect_color,size=8) +
    scale_y_continuous(
      name = "Cytotoxicity [%]", limits = c(min_y_axis_limits_cytotoxicity,log_y_axis_limits_cytotoxicity),position = "left") + theme(axis.ticks=element_blank()) + xlim((log_concentration_unexposed_controls)*0.99999,(log_concentration_unexposed_controls)*1.00001)
  
  p_right = ggplot(Log_all_frame,aes(x=Concentration)) + 
    geom_point(data=Log_all_frame,aes(y = Cytotox,shape = Name),fill=Cyto_color,size=8) + 
    geom_point(data=Log_all_frame,aes(y = Effect,shape= Name),fill=Effect_color,size=8) +
    scale_y_continuous(
      name = "Neurite outgrowth inhibition [%]", limits = c(min_y_axis_limits_effects,log_y_axis_limits_effects),position="right")
  
  Log_all_plot_left = p_left + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=32,
                                      color=Cyto_color,
                                      face="bold",
                                      angle=90)) + 
    scale_shape_prism(palette="filled")
  
  Log_all_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.title.y.right = element_text(size=32,
                                            color=Effect_color,
                                            face="bold",
                                            angle=270)) + 
    xlab(paste0("log Concentration (",workframe$Unit[1],")")) + 
    scale_shape_prism(palette="filled")
  
  Log_all_plot = (Log_all_plot_left|Log_all_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Log_all_plot
  
  #Fit
  if(is.na(loglogistic_regression_cytotoxicity)|final_model_cyt=="cnst"|nrow(Log_cyt_frame)<3){
    log_cyt_pred = 0
    Log_cyt_frame = data.frame("x"=Concentrations_log,
                               "y"=Cytotoxicity,
                               "Name"=rep(workframe$Measurement,each=22))
  } else {
    log_cyt_pred = loglogistic_regression_cytotoxicity$hill_modl
  }
  
  if(is.null(log_cyt_pred)){
    log_cyt_pred = 0
    Log_cyt_frame = data.frame("x"=Concentrations_log,
                               "y"=Cytotoxicity,
                               "Name"=rep(workframe$Measurement,each=22))
  }
  
  p_left = ggplot(Log_cyt_control_frame,aes(x=x)) + 
    geom_point(data=Log_cyt_control_frame,aes(y = y,shape = Name),fill=Cyto_color,size=8) + 
    geom_line(aes(y=10),linetype="dashed",color=Cyto_color,linewidth=2) + 
    labs(y=NULL) + ylim((min(c(Log_cyt_frame$y,Log_cyt_control_frame$y),na.rm=T)-5),log_y_axis_limits_cytotoxicity) + xlim((log_concentration_unexposed_controls)*0.99999,(log_concentration_unexposed_controls)*1.00001)
  
  
  p_right = ggplot(Log_cyt_frame,aes(x=x)) + 
    geom_point(data=Log_cyt_frame,aes(y = y,shape = Name),fill=Cyto_color,size=8) + labs(y=NULL) +
    geom_line(aes(y=10),linetype="dashed",color=Cyto_color,linewidth=2) + 
    geom_line(aes(y = log_cyt_pred),color=Cyto_color,linewidth=3) + coord_cartesian(ylim=c(min_y_axis_limits_cytotoxicity,log_y_axis_limits_cytotoxicity),clip = "on")
  
  
  Log_cyt_plot_left = p_left + theme_prism(base_size = 30) + ylab("Cytotoxicity [%]") + 
    theme(legend.position = "none",
          axis.title.y = element_text(size=32,
                                      color=Cyto_color,
                                      face="bold",
                                      angle=90),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_shape_prism(palette="filled")
  
  Log_cyt_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none") + 
    xlab(paste0("log Concentration (",workframe$Unit[1],")")) +
    scale_shape_prism(palette="filled") + labs(y=NULL) + guides(y = "none")
  
  Log_cyt_plot = (Log_cyt_plot_left|Log_cyt_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Log_cyt_plot
  
  #---------------------
  
  if(is.na(loglogistic_regression_Effect)|final_model_effect=="cnst"|nrow(Log_effect_frame)<3){
    log_effect_pred = 0
    Log_effect_frame = data.frame("x"=Concentrations_log,
                                  "y"=Effect,
                                  "Name"=rep(workframe$Measurement,each=22))
  } else {
    log_effect_pred = loglogistic_regression_Effect$hill_modl
  }
  
  if(is.null(log_effect_pred)){
    log_effect_pred = 0
    Log_effect_frame = data.frame("x"=Concentrations_log,
                                  "y"=Effect,
                                  "Name"=rep(workframe$Measurement,each=22))
  }
  
  p_left = ggplot(Log_effect_control_frame,aes(x=x)) + 
    geom_point(data=Log_effect_control_frame,aes(y = y,shape = Name),fill=Effect_color,size=8) + 
    geom_line(aes(y=10),linetype="dashed",color=Effect_color,linewidth=2) + 
    labs(y=NULL) + ylim((min(c(Log_effect_frame$y,Log_effect_control_frame$y),na.rm=T)-5),log_y_axis_limits_effects) + xlim((log_concentration_unexposed_controls)*0.99999,(log_concentration_unexposed_controls)*1.00001)
  
  p_right = ggplot(Log_effect_frame,aes(x=x)) + 
    geom_point(data=Log_effect_frame,aes(y = y,shape = Name),fill=Effect_color,size=8) +
    geom_line(aes(y=10),linetype="dashed",color=Effect_color,linewidth=2) + 
    geom_line(aes(y = log_effect_pred),color=Effect_color,linewidth=3) +
    scale_y_continuous(
      name = "Neurite outgrowth inhibition [%]",position = "right") + coord_cartesian(ylim=c(min_y_axis_limits_effects,log_y_axis_limits_effects),clip = "on")
  
  
  Log_effect_plot_left = p_left + theme_prism(base_size = 30) + labs(y=NULL) + guides(y="none") + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    scale_shape_prism(palette="filled")
  
  Log_effect_plot_right = p_right + theme_prism(base_size = 30) + 
    theme(legend.position = "none",
          axis.title.y.right = element_text(size=32,
                                            color=Effect_color,
                                            face="bold",
                                            angle=90)) + 
    xlab(paste0("log Concentration (",workframe$Unit[1],")")) +
    scale_shape_prism(palette="filled")
  
  Log_effect_plot = (Log_effect_plot_left|Log_effect_plot_right) + plot_layout(ncol=2,widths=c(0.1,1))
  Log_effect_plot
  
  #Theme for all
  if(use_letters_for_subplots==T){
    if(which_plots=="all"){
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 15)
      plot = ggarrange(Linear_all_plot,Linear_cyt_plot,Linear_effect_plot,
                       Log_all_plot,Log_cyt_plot,Log_effect_plot,labels = c("A","B","C","D","E","F"),ncol = 3,nrow = 2,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 50)))
      dev.off()
    } else if(which_plots=="log"){
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 9)
      plot = ggarrange(Log_all_plot,Log_cyt_plot,Log_effect_plot,labels = c("A","B","C"),ncol = 3,nrow = 1,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 40)))
      dev.off()
    } else {
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 9)
      plot = ggarrange(Linear_all_plot,Linear_cyt_plot,Linear_effect_plot,labels = c("A","B","C"),ncol = 3,nrow = 1,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 40)))
      dev.off()
    }
  } else {
    if(which_plots=="all"){
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 15)
      plot = ggarrange(Linear_all_plot,Linear_cyt_plot,Linear_effect_plot,
                       Log_all_plot,Log_cyt_plot,Log_effect_plot,ncol = 3,nrow = 2,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 50)))
      dev.off()
    } else if(which_plots=="log"){
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 9)
      plot = ggarrange(Log_all_plot,Log_cyt_plot,Log_effect_plot,ncol = 3,nrow = 1,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 40)))
      dev.off()
    } else {
      pdf(file=paste0(graph_path,paste0("/",Final_ID,".pdf")),width = 30,height = 9)
      plot = ggarrange(Linear_all_plot,Linear_cyt_plot,Linear_effect_plot,ncol = 3,nrow = 1,
                       font.label = list(size = 35, color = "black", face = "bold", family = NULL),
                       align = "h")
      print(annotate_figure(plot, top = text_grob(Final_ID, 
                                                  color = "black", face = "bold", size = 40)))
      dev.off()
    }
  }
  
  #Get rid of fits so that next samples does not have it
  linear_regression_cytotox = NA
  linear_regression_effect = NA
  loglogistic_regression_cytotoxicity = NA
  loglogistic_regression_Effect = NA
  
}

#Update Summary file with results
Output = openxlsx::createWorkbook()
openxlsx::addWorksheet(Output,"Sample_List")
openxlsx::writeData(Output,"Sample_List",Summary_list)
openxlsx::addWorksheet(Output,"Results")
Summary_results = as.data.frame(apply(Summary_results,c(1,2),makeExcelreadable))
openxlsx::writeData(Output,"Results",Summary_results)


openxlsx::addWorksheet(Output,"Summary")

Final_Summary_Effects = as.data.frame(matrix(nrow=nrow(Summary_results),ncol=5))
colnames(Final_Summary_Effects) = c("Sample","Unit","EC10","SE_EC10","Confident_EC10?")
Final_Summary_Cytotoxicity = as.data.frame(matrix(nrow=nrow(Summary_results),ncol=3))
colnames(Final_Summary_Cytotoxicity) = c("IC10","SE_IC10","Confident_IC10?")

Final_Summary_Effects$Sample = Summary_results$Sample_ID
Final_Summary_Effects$Unit = Summary_results$Unit

#Effects
if(Summary_model_effects=="linear"){
  for(i in 1:nrow(Summary_results)){
    if(is.na(as.numeric(Summary_results$SHSY5Y_linmodel_EC10[i]))){
      if(grepl("masked",Summary_results$SHSY5Y_linmodel_EC10[i])){
        Final_Summary_Effects[i,3] = "masked by cytotoxicity"
        Final_Summary_Effects[i,4] = NA
        Final_Summary_Effects[i,5] = FALSE
      } else {
        Final_Summary_Effects[i,3] = "no effect"
        Final_Summary_Effects[i,4] = NA
        Final_Summary_Effects[i,5] = FALSE
      }
    } else if(is.na(as.numeric(Summary_results$SHSY5Y_linmodel_confidence_EC10[i]))|as.numeric(Summary_results$SHSY5Y_linmodel_confidence_EC10[i])<0.70){
      Final_Summary_Effects[i,3] = Summary_results$SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,4] = Summary_results$SE_SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,5] = FALSE
    } else {
      Final_Summary_Effects[i,3] = Summary_results$SHSY5Y_linmodel_EC10[i]
      Final_Summary_Effects[i,4] = Summary_results$SE_SHSY5Y_linmodel_EC10[i]
      Final_Summary_Effects[i,5] = TRUE
    }
  }
} else {
  colnames(Final_Summary_Effects) = c("Sample","Unit","log_EC10","log_SE_EC10","Confident_log_EC10?")
  for(i in 1:nrow(Summary_results)){
    if(is.na(as.numeric(Summary_results$SHSY5Y_logmodel_EC10[i]))){
      if(grepl("masked",Summary_results$SHSY5Y_logmodel_EC10[i])){
        Final_Summary_Effects[i,3] = "masked by cytotoxicity"
        Final_Summary_Effects[i,4] = NA
        Final_Summary_Effects[i,5] = FALSE
      } else {
        Final_Summary_Effects[i,3] = "no effect"
        Final_Summary_Effects[i,4] = NA
        Final_Summary_Effects[i,5] = FALSE
      }
    } else if(is.na(as.numeric(Summary_results$SHSY5Y_logmodel_confidence_EC10[i]))|as.numeric(Summary_results$SHSY5Y_logmodel_confidence_EC10[i])<0.70){
      Final_Summary_Effects[i,3] = Summary_results$SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,4] = Summary_results$SE_SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,5] = FALSE
    } else {
      Final_Summary_Effects[i,3] = Summary_results$SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,4] = Summary_results$SE_SHSY5Y_logmodel_EC10[i]
      Final_Summary_Effects[i,5] = TRUE
    }
  }
}

#Cytotoxicity
if(Summary_model_cytotoxicity=="linear"){
  for(i in 1:nrow(Summary_results)){
    if(is.na(as.numeric(Summary_results$SHSY5Y_linmodel_IC10[i]))){
      Final_Summary_Cytotoxicity[i,1] = "no cytotoxicity"
      Final_Summary_Cytotoxicity[i,2] = NA
      Final_Summary_Cytotoxicity[i,3] = FALSE
    } else if(is.na(as.numeric(Summary_results$SHSY5Y_linmodel_confidence_IC10[i]))|as.numeric(Summary_results$SHSY5Y_linmodel_confidence_IC10[i])<0.70){
      Final_Summary_Cytotoxicity[i,1] = Summary_results$SHSY5Y_linmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,2] = Summary_results$SE_SHSY5Y_linmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,3] = FALSE
    } else {
      Final_Summary_Cytotoxicity[i,1] = Summary_results$SHSY5Y_linmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,2] = Summary_results$SE_SHSY5Y_linmodel_IC10[i]+
      Final_Summary_Cytotoxicity[i,3] = TRUE
    }
  }
} else {
  colnames(Final_Summary_Cytotoxicity) = c("log_IC10","log_SE_IC10","Confident_log_IC10?")
  for(i in 1:nrow(Summary_results)){
    if(is.na(as.numeric(Summary_results$SHSY5Y_logmodel_IC10[i]))){
      Final_Summary_Cytotoxicity[i,1] = "no cytotoxicity"
      Final_Summary_Cytotoxicity[i,2] = NA
      Final_Summary_Cytotoxicity[i,3] = FALSE
    } else if(is.na(as.numeric(Summary_results$SHSY5Y_logmodel_confidence_IC10[i]))|as.numeric(Summary_results$SHSY5Y_logmodel_confidence_IC10[i])<0.70){
      Final_Summary_Cytotoxicity[i,1] = Summary_results$SHSY5Y_logmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,2] = Summary_results$SE_SHSY5Y_logmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,3] = FALSE
    } else {
      Final_Summary_Cytotoxicity[i,1] = Summary_results$SHSY5Y_logmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,2] = Summary_results$SE_SHSY5Y_logmodel_IC10[i]
      Final_Summary_Cytotoxicity[i,3] = TRUE
    }
  }
}

Final_Summary = cbind(Final_Summary_Effects,Final_Summary_Cytotoxicity)

openxlsx::writeData(Output,"Summary",Final_Summary)

Output_path = paste0(stringr::str_split(Summary_path,"/",simplify = T))
Output_path = paste0(Output_path[1:(length(Output_path)-1)],collapse = "/")
openxlsx::saveWorkbook(Output,paste0(Output_path,"/Summary_",project,"_SHSY5Y_filled.xlsx"),overwrite = T)

