
# Author: Lyron Juan Winderbaum, email: lyron.winderbaum@student.adelaide.edu.au

library(stringr)
library(base)
library(data.table)
library(reshape2)
library(ggplot2)
library(plyr)
library(grid)

# This finds unique elements in a vector and counts their occurances
cU <- function(x){ data.table(x)[, .N, keyby= x] }

# These are functions for extracting the X,Y coordinates and Region Numbers from peaklist files.
Xcoord   <- function(x) as.numeric(substring(str_extract(x,"X\\d{3,4}"),2))
Ycoord   <- function(x) as.numeric(substring(str_extract(x,"Y\\d{3,4}"),2))
RegionNo <- function(x) as.numeric(substring(str_extract(x,"R\\d{2,3}"),2))
Peaklist_ID <- function(x) str_extract(x,"R\\d{2,3}X\\d{3,4}Y\\d{3,4}")

# For reading in peaklist files into a comprehensive peaklist data.frame in R.
readPeaklists <- function(dataset_name, peaklist_folder_name = "peaklists", parent_folder_name = ".", data_folder = "./data"){
  # If the folder with the data where not in the current folder you give the folder 
  # it is in in parent_folder_name
  # If the subfolder containing the peaklists is not named "peaklists", for example
  # if it where named "peaklists_quad" or something, you would put 
  # peaklist_folder_name = "peaklists_quad".
  
  ###################################################
  # Find the peaklist files and extract their names #
  ###################################################
  peaklist_folder_path <- paste(parent_folder_name,dataset_name,peaklist_folder_name,sep="/")
  # Matches files to the regular expression for peaklist files.
  peaklist_file_names <- list.files(path = peaklist_folder_path,"R\\d{2,3}X\\d{3,4}Y\\d{3,4}.txt")
  if (length(peaklist_file_names) == 0){
    print("ERROR: No peaklist files found, aborting readPeaklists command")
  } else {
    #############################
    # Generates fExists and LXY #
    #############################
    # Extracts the Region numbers and X,Y coordinates for all the spectra.
    LXY <- data.frame(fname = peaklist_file_names,
                      Peaklist = sapply(peaklist_file_names,Peaklist_ID,USE.NAMES=FALSE),
                      Region = sapply(peaklist_file_names,RegionNo,USE.NAMES=FALSE),
                      X = sapply(peaklist_file_names,Xcoord,USE.NAMES=FALSE),
                      Y = sapply(peaklist_file_names,Ycoord,USE.NAMES=FALSE),
                      Acquisition = rep(0,length(peaklist_file_names))
                      )
    minX <- min(LXY$X)
    minY <- min(LXY$Y)
    X <- minX:max(LXY$X)
    Y <- minY:max(LXY$Y)
    fExists <- matrix(rep(0,(length(X)+2)*(length(Y)+2)),nrow=length(X)+2)
    count <- 0
    for (region in sort(unique(LXY$Region))){
      for (y in sort(unique(LXY[LXY$Region==region,]$Y))){
        for (x in sort(unique(LXY[LXY$Region==region & LXY$Y==y,]$X))){
          count <- count + 1
          fExists[x-minX+2,y-minY+2] <- count
          if (sum(LXY$Region==region & LXY$Y==y & LXY$X==x) == 1){
            LXY[LXY$Region==region & LXY$Y==y & LXY$X==x,]$Acquisition = count
          } else {
            print("WARNING: THIS SHOULD NOT HAPPEN.")
          }
        }
      }
    }
    write.table(fExists,file=paste(data_folder, "/", dataset_name,"_fExists.txt",sep=""),
                sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(LXY,file=paste(data_folder, "/", dataset_name,"_LXY.txt",sep=""),
                sep="\t",row.names=FALSE,col.names=TRUE)
    
    ################################################################
    # Read the peaklist files and compile a comprehensive peaklist #
    ################################################################
    peaklist_not_empty = logical(length=length(peaklist_file_names))
    peaklist_all_does_not_exist <-TRUE
    peaklist_temp_does_not_exist <-TRUE
    count <- 0
    # For each peaklist file
    for (spec_idx in 1:length(peaklist_file_names)){ 
      fname <- peaklist_file_names[spec_idx]
      # Read the peaklist file
      peaklist_cur <- read.table(paste(peaklist_folder_path,fname,sep="/"), header=TRUE)  
      # Check that it is not empty (no peaks)
      if (nrow(peaklist_cur) > 0){
        # It's not empty!
        peaklist_not_empty[spec_idx] <- TRUE
        # Annotate peaks with the name of their parent peaklist and acquisition number.
        peaklist_cur <- transform(peaklist_cur, Peaklist = Peaklist_ID(fname))
        peaklist_cur <- transform(peaklist_cur, Acquisition = LXY[LXY$fname == fname,]$Acquisition)
        # Add the peaks just read in to a temorary caching variable.
        if (peaklist_temp_does_not_exist) {
          peaklist_temp <- peaklist_cur
          peaklist_temp_does_not_exist <- FALSE
        } else {
          peaklist_temp <- rbind(peaklist_temp,peaklist_cur)
        }
        count <- count + 1
      }
      # Add the cached peaks onto the comprehensive peaklist
      # In order to fine-tune performance you can tweak "x" in (count == "x"),
      # depending on how large you expect your peaklist files to be.
      # This number determines how many peaklist files it will cache before
      # adding them to the comprehensive peaklist.
      if (count == 1000){
        if (peaklist_all_does_not_exist) {
          peaklist_all <- peaklist_temp
          peaklist_all_does_not_exist <- FALSE
        } else {
          peaklist_all <- rbind(peaklist_all,peaklist_temp)
        }    
        count <- 0
        peaklist_temp_does_not_exist <- TRUE
        # Diagnostic, each time we combine cached peaks we print the peaklist 
        # we most recently read to allow you to see progress.
        print(fname)
      }
    }
    # add any remaining cached peaks to the comprehensive peaklist.
    if (peaklist_all_does_not_exist) {
      peaklist_all <- peaklist_temp
    } else {
      peaklist_all <- rbind(peaklist_all,peaklist_temp)
    }
    # Diagnostic: Lets you know how many empty peaklists where found total.
    print(paste(toString(sum(!peaklist_not_empty)),"Empty Peaklists."))
    # could save peaklist_not_empty here, if you needed it for anything later on.
    write.table(peaklist_all,file=paste(data_folder, "/",dataset_name,"_comprehensive_peaklist.txt",sep=""),
                sep="\t",row.names=FALSE)
    return(peaklist_all)
  }
}

load_peaklist <- function(dataset_name,data_folder="./data"){
  peaklist_all <- read.table(paste(data_folder, "/", dataset_name, "_comprehensive_peaklist.txt",sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  return(peaklist_all)
}

load_LXY <- function(dataset_name,data_folder="./data"){
  LXY <- read.table(paste(data_folder, "/", dataset_name, "_LXY.txt",sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  return(LXY)
}

load_fExists <- function(dataset_name,data_folder="./data"){
  fExists <- read.table(paste(data_folder, "/", dataset_name, "_fExists.txt",sep=""), sep="\t", header=FALSE)
  return(fExists)
}



# This is for taking a subset of a peaklist centered around certain known m/z values, 
# for example calibrants. fixed Da bins (use_ppm = FALSE) or ppm based tolerances
# (use_ppm = TRUE) are both supported.
mzMatch <- function(peaklist_in,mzList,binMargin=0.3,use_ppm=FALSE) {
  peaklist_subset_does_not_exist = TRUE
  if (use_ppm){
    binMargin_ppm <- binMargin
  }
  for (i in 1:length(mzList)){
    if (use_ppm){
      binMargin <- binMargin_ppm*mzList[i]/1000000
    }
    if (peaklist_subset_does_not_exist){
      peaklist_subset <- transform(peaklist_in[which(abs(peaklist_in$m.z - mzList[i])<binMargin),],
                                   PeakGroup=mzList[i])
      peaklist_subset_does_not_exist = FALSE
    } else {
      peaklist_subset <- rbind(peaklist_subset,transform(peaklist_in[which(abs(peaklist_in$m.z - mzList[i])<binMargin),],
                                                         PeakGroup=mzList[i]))
    }
  }
  return(peaklist_subset)
}


# Does a peak-grouping, and annotates peaks by peakgroup, in case you want that.
# If you set a non-zero minGroupSize here any peaks not allocated to groups will be 
# annotates peakgroup zero. This can always be done later though, with the cU function 
# in the localFunctions.R file for example, so I reccomend leaving minGroupSize = 0 here.
# You could modify tol (the tolerance used) if you wish however.
groupPeaks <- function(peaklist_in,tol = 0.1, minGroupSize = 0) {
  peaklist_in <- peaklist_in[order(peaklist_in$m.z),]
  nPeaks <- nrow(peaklist_in)
  peaklist_in <- transform(peaklist_in,PeakGroup = 1)
  for (i in which(peaklist_in[2:nPeaks,]$m.z - peaklist_in[1:(nPeaks-1),]$m.z > tol)){
    peaklist_in$PeakGroup[(i+1):nPeaks] <- peaklist_in$PeakGroup[(i+1):nPeaks] + 1 
  }
  for (p in which(cU(peaklist_in$PeakGroup)$N < minGroupSize)){
    peaklist_in[which(peaklist_in$PeakGroup == p),]$PeakGroup <- 0
  }
  return(peaklist_in)
}

dbscan_lw <- function(peaklist_in,eps=0.05,mnpts=100,cvar="m.z",pp=TRUE){
  # Implements DBSCAN* as in section 3 of:
  #
  # Campello, Ricardo JGB, Davoud Moulavi, and Joerg Sander. 
  # "Density-based clustering based on hierarchical density estimates." 
  # In Advances in Knowledge Discovery and Data Mining, pp. 160-172. 
  # Springer Berlin Heidelberg, 2013.
  #
  # For one-dimensional objects only. Intended for clustering
  # peaks by m/z location in MALDI Imaging.
  
  if(sum(cvar==names(peaklist_in))==0){
    error("Non-existent variable selected for clustering.")
  }
  
  n <- nrow(peaklist_in)
  # sort
  if(pp){
    print("Sorting...")
  }
  peaklist_in <- peaklist_in[order(peaklist_in[,cvar]),]
  if(pp){
    print("Done")
  }
  p_locs <- peaklist_in[,cvar]
  p_core <- rep(0,n)
  if(pp){
    print("Counting neighbours")
  }
  for(i in 1:(mnpts+1)){
    temp = (p_locs[(1+i):n] - p_locs[1:(n-i)]) <= eps
    p_core[(1+i):n] = p_core[(1+i):n] + temp
    p_core[1:(n-i)] = p_core[1:(n-i)] + temp
    if(pp){
      print(paste(toString(i-1),"/",toString(mnpts)))
    }
  }
  # Identify core points
  p_core <- p_core >= mnpts
  n_core <- sum(p_core)
  if(n_core == 0){
    print("No core points")
    return(FALSE)
  }
  # Check there is more than one core point
  if(n_core == 1){
    peaklist_in$PeakGroup <- as.numeric(p_core)
    return(peaklist_in)
  }
  # calc pairwise (adjacent) distances between core points (only)
  p_locs <- p_locs[p_core]
  d_pair <- p_locs[2:n_core] - p_locs[1:(n_core-1)] 
  clus <- c(1,1+cumsum(d_pair > eps))
  
  peaklist_in$PeakGroup = 0
  peaklist_in[p_core,"PeakGroup"] = clus
  return(peaklist_in)
}









DIPPS <- function(pl_in){
  # peaklist_all should be a data.frame with at least 
  # three variables:
  #  - Acquisition (identifying unique spectra)
  #  - PeakGroup (identifying peakgroups)
  #  - Group (identifying regions to be compared by 
  #     DIPPS), and should take three values:
  #      - 1 coding for the `downregulated' group, and
  #      - 2 coding for the `upregulated group.
  
  # Commented out -- meaning assume input is unique
  # Check for multiple peaks and remove.
  #   pl_in = unique(pl_in)
  
  nSpec_d = length(unique(subset(pl_in,Group == 1)$Acquisition))
  nSpec_u = length(unique(subset(pl_in,Group == 2)$Acquisition))
  
  prop = ddply(pl_in,
               c("PeakGroup","Group"),
               summarise,
               p = length(Acquisition)
  )
  prop[prop$Group == 1,"p"] = prop[prop$Group == 1,"p"]/nSpec_d
  prop[prop$Group == 2,"p"] = prop[prop$Group == 2,"p"]/nSpec_u
  
  prop <- reshape(prop,
                  timevar = "Group",
                  idvar="PeakGroup",
                  direction="wide")
  names(prop) = c("PeakGroup","p.d","p.u")
  
  prop = replace(prop,is.na(prop),0)
  prop$d = prop$p.u - prop$p.d
  
  return(prop)
}




dippsHeur <- function(pl_in,dsum_in){
  # Takes the output of DIPPS as input, and calculates
  # a heuristic cutoff for the `optimal' number of 
  # variables with highest DIPPS.
  
  u.m = dcast(subset(pl_in,Group==2),
              Acquisition~PeakGroup,
              value.var="Group")
  acq = u.m[,1]
  u.m = as.matrix(u.m[,-1])
  u.m[!is.na(u.m)] = 1
  u.m[is.na(u.m)] = 0
  acq = data.frame(Acquisition = acq, 
                   nPeaks = rowSums(u.m))
  u.m = u.m/sqrt(acq$nPeaks)
  c.u = colMeans(u.m)
  c.u = c.u/norm(c.u,"2")
  
  temp = match(colnames(u.m),dsum_in$PeakGroup)
  dsum_in$c.u = 0
  dsum_in[temp,"c.u"] = c.u
  
  # Find data-driven cutoff according to the DIPPS method using cosine distance.
  curMinCosD = 10
  sortedDIPPS = sort(dsum_in$d,index.return=TRUE)
  # vN = 1:floor(nrow(Summary_merged)/2)
  vN = 1:nrow(dsum_in)
  cosD = vN
  for (n in vN){
    dsum_in$t = 0
    dsum_in[tail(sortedDIPPS$ix,n),]$t = 1/sqrt(n)
    cosD[n] = 1 - sum(dsum_in$t * dsum_in$c.u)
  }
  nStar = vN[which.min(cosD)]
  
  return(nStar)
}










# Plots a spatial image
spatialPlot <- function(peaklist_in,fExists_in,
                        plot_var = "intensity",
                        plot_var_transform = "none",
                        plot_var_type = "continuous",
                        mult_peaks = "average",
                        save_plot = FALSE,
                        plot_name_in = "",
                        minX_in = 1,
                        minY_in = 1,
                        display_pixel_borders = FALSE,
                        display_legend = TRUE,
                        return_mI.m = FALSE
                        ){
  
  if (plot_var_type == "continuous"){
    if (is.na(match(plot_var,c("m.z","SN","QualityFactor","Resolution","intensity","area","count")))){
      print(paste("Plotting Variable",plot_var,"is not currently supported."))
      print("Defaulting to intensity")
      plot_var <- "intensity"
      print("In order to plot non-standard variables the spatialPlot() function will need to be modified.")
    } else {
      if (plot_var == "count" & mult_peaks != "sum") {
        mult_peaks <- "sum"
      }
    }
    if (is.na(match(mult_peaks,c("average","sum","max")))){
      print(paste("Method for reducing multiple peaks",mult_peaks,"is not currently supported."))
      print("Defaulting to averaging")
      mult_peaks <- "average"
      print("In order to use non-standard methods the spatialPlot() function will need to be modified.")    
    }
  } else if (plot_var_type != "categorical") {
    print(paste("Variable type",plot_var_type,"is not currently supported."))
    print("Defaulting to categorical")
    plot_var_type <- "categorical"
    print("In order to use non-standard methods the spatialPlot() function will need to be modified.")          
  }
  
  temp = cU(peaklist_in$Acquisition)
  # Deal with multiple peaks.
  
  if (plot_var_type == "continuous"){
    if (plot_var == "count") {
      peaklist_in$count = 1
    }
    if (sum(temp$N>1) > 0) {
      peaklist_in <- switch(mult_peaks,
                            average = switch(plot_var,
                                             m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = mean(m.z)),
                                             SN            = ddply(peaklist_in,"Acquisition",summarise,SN = mean(SN)),
                                             QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = mean(QualityFactor)),
                                             Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = mean(Resolution)),
                                             intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = mean(intensity)),
                                             area          = ddply(peaklist_in,"Acquisition",summarise,area = mean(area))
                            ),
                            sum = switch(plot_var,
                                             m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = sum(m.z)),
                                             SN            = ddply(peaklist_in,"Acquisition",summarise,SN = sum(SN)),
                                             QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = sum(QualityFactor)),
                                             Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = sum(Resolution)),
                                             intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = sum(intensity)),
                                             area          = ddply(peaklist_in,"Acquisition",summarise,area = sum(area)),
                                             count         = ddply(peaklist_in,"Acquisition",summarise,count = sum(count))
                            ),
                            max = switch(plot_var,
                                             m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = max(m.z)),
                                             SN            = ddply(peaklist_in,"Acquisition",summarise,SN = max(SN)),
                                             QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = max(QualityFactor)),
                                             Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = max(Resolution)),
                                             intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = max(intensity)),
                                             area          = ddply(peaklist_in,"Acquisition",summarise,area = max(area))
                            )
      )
    }
  } else if (sum(temp$N>1) > 0) {
    stop("In order to plot categorical variables spectra must be uniquely specified.")
  }
  
  if (plot_var_type == "continuous"){
    if (is.na(match(plot_var_transform,c("none","log")))){
      print(paste("Transformation",plot_var_transform,"is not currently supported."))
      print("Defaulting to none")
      plot_var_transform = "none"
      print("In order to use non-standard transformations the spatialPlot() function will need to be modified.")
    }
    if (plot_var_transform == "log"){
      min_plot_var = min(round(log(1 + peaklist_in[,plot_var]),5))
      max_plot_var = max(round(log(1 + peaklist_in[,plot_var]),5))
    } else {
      min_plot_var = min(peaklist_in[,plot_var])
      max_plot_var = max(peaklist_in[,plot_var])
    }
  } else {
    plot_var_transform <- "none"
    if (is.factor(peaklist_in[,plot_var])) {
      peaklist_in[,plot_var] = levels(peaklist_in[,plot_var])[as.numeric(peaklist_in[,plot_var])]
    }
  }
  
  if (is.data.frame(fExists_in)){
    fExists_in = list(fExists_in)
  }
  mI.m_out <- as.list(1:length(fExists_in))
  for (region_idx in 1:length(fExists_in)){
    fExists = fExists_in[[region_idx]]
    if (length(fExists_in) != 1){
      if (length(fExists_in) == length(plot_name_in)){
        plot_name = plot_name_in[region_idx]
      } else {
        plot_name = toString(region_idx)
      }
      if (length(fExists_in) == length(minX_in)){
        minX = minX_in[region_idx]
      } else {
        minX = 1
      }
      if (length(fExists_in) == length(minY_in)){
        minY = minY_in[region_idx]
      } else {
        minY = 1
      }
    } else {
      plot_name = plot_name_in
      minX = minX_in
      minY = minY_in
    }
    
    # Generate an image matrix mI
    mI <- fExists
    
    mI[mI > 0] <- match(mI[mI > 0],peaklist_in$Acquisition)
    subset_of_peaklist = mI[mI > 0 & !is.na(mI)]
    mI[mI > 0 & !is.na(mI)] <- peaklist_in[subset_of_peaklist,plot_var]
        
    mI$X = (1:nrow(mI)) + minX - 2
    mI.m = melt(mI,id.var="X")
    mI.m$Y = (as.numeric(substring(mI.m$variable,2))) + minY - 2
    
    if (plot_var_transform == "log"){
      mI.m$value = round(log(1 + mI.m$value),5)
    }
    
    # Make a new variable empty for acquisition regions that have no peaks, 
    #    and for non-acquisition regions set value to NA     
    mI.m$empty = FALSE
    if (sum(is.na(mI.m$value)) > 0) {
      mI.m[is.na(mI.m$value),]$empty = TRUE
    }
    if (sum(mI.m[!is.na(mI.m$value),"value"]==0) > 0){
      mI.m[replace(mI.m$value,is.na(mI.m$value),1)==0,]$value = NA
    }
    
    if (plot_var_type == "categorical"){
      mI.m$value <- factor(mI.m$value)
    }
    
    if (!return_mI.m){
      p = ggplot(mI.m,aes(X,Y))
      if (display_pixel_borders){
        p = p + geom_tile(data = mI.m,aes(fill=value,alpha=as.numeric(!is.na(value))),colour="grey")
      } else {
        p = p + geom_tile(data = mI.m,aes(fill=value,alpha=as.numeric(!is.na(value))),colour=NA)  
      }
      if (display_legend){
        p = p + guides(alpha = FALSE)
      } else {
        p = p + guides(alpha = FALSE,fill=FALSE)
      }
      p = p + geom_tile(data = mI.m,alpha=0.5*as.numeric(mI.m$empty))
      if (plot_var_type == "continuous"){
        p = p + scale_fill_gradient(name=plot_var, limits = c(min_plot_var,max_plot_var))
      }
      #     p = p + ggtitle(plot_name)
      p = p + coord_fixed()
      p = p + scale_y_reverse()
      
      if (plot_name != ""){
        plot_name = paste(plot_name,"_",sep="")
      }
      
      if (save_plot){
        ggsave(paste(plot_name,paste(mult_peaks,plot_var,"transformation",plot_var_transform,sep="_"),".png",sep=""),p)
      }
      
      if (region_idx == 1){
        p_list = list(p)
      } else {
        p_list[[region_idx]] = p
      }
    } else {
      mI.m_out[[region_idx]] = mI.m
    }
  }
  if (return_mI.m){
    if (length(fExists_in)==1){
      return(mI.m_out[[1]])
    } else {
      return(mI.m_out)
    }
  } else {
    if (length(fExists_in)==1){
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  }
}

# Plots an acquisition plot (with acquisition order on the x-axis.)
acquisitionPlot <- function(peaklist_in,
                            plot_var = "intensity",
                            plot_var_transform = "none",
                            mult_peaks = "average",
                            save_plot = FALSE,
                            plot_name = ""){
  
  if (is.na(match(plot_var,c("m.z","SN","QualityFactor","Resolution","intensity","area","count")))){
    print(paste("Plotting Variable",plot_var,"is not currently supported."))
    print("Defaulting to intensity")
    plot_var <- "intensity"
    print("In order to plot non-standard variables the spatialPlot() function will need to be modified.")
  } else {
    if (plot_var == "count" & mult_peaks != "sum") {
      mult_peaks <- "sum"
    }
  }
  if (is.na(match(mult_peaks,c("average","sum","max")))){
    print(paste("Method for reducing multiple peaks",mult_peaks,"is not currently supported."))
    print("Defaulting to averaging")
    mult_peaks <- "average"
    print("In order to use non-standard methods the spatialPlot() function will need to be modified.")    
  }
  
  # Generate an image matrix mI
  temp = cU(peaklist_in$Acquisition)
  # Deal with multiple peaks.
  if (plot_var == "count") {
    peaklist_in$count = 1
  }
  if (sum(temp$N>1) > 0) {
    peaklist_in <- switch(mult_peaks,
                          average = switch(plot_var,
                                           m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = mean(m.z)),
                                           SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = mean(SN)),
                                           QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = mean(QualityFactor)),
                                           Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = mean(Resolution)),
                                           intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = mean(intensity)),
                                           area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = mean(area))
                          ),
                          sum = switch(plot_var,
                                       m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = sum(m.z)),
                                       SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = sum(SN)),
                                       QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = sum(QualityFactor)),
                                       Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = sum(Resolution)),
                                       intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = sum(intensity)),
                                       area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = sum(area)),
                                       count         = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,count = sum(count))
                          ),
                          max = switch(plot_var,
                                       m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = max(m.z)),
                                       SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = max(SN)),
                                       QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = max(QualityFactor)),
                                       Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = max(Resolution)),
                                       intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = max(intensity)),
                                       area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = max(area))
                          )
    )
  }
  
  if (is.na(match(plot_var_transform,c("none","log")))){
    print(paste("Transformation",plot_var_transform,"is not currently supported."))
    print("Defaulting to none")
    plot_var_transform = "none"
    print("In order to use non-standard transformations the spatialPlot() function will need to be modified.")
  }
  if (plot_var_transform == "log"){
    peaklist_in[,plot_var] = round(log(1 + peaklist_in[,plot_var]),5)
  }
  
  
  p <- switch(plot_var,
              m.z           = ggplot(peaklist_in,aes(x=Acquisition,y=m.z)),
              SN            = ggplot(peaklist_in,aes(x=Acquisition,y=SN)),
              QualityFactor = ggplot(peaklist_in,aes(x=Acquisition,y=QualityFactor)),
              Resolution    = ggplot(peaklist_in,aes(x=Acquisition,y=Resolution)),
              intensity     = ggplot(peaklist_in,aes(x=Acquisition,y=intensity)),
              area          = ggplot(peaklist_in,aes(x=Acquisition,y=area)),
              count         = ggplot(peaklist_in,aes(x=Acquisition,y=count))
  )
  p <- p + layer(geom = "point",alpha=I(1/12))
  p <- p + layer(geom="smooth",method="gam",formula=y~s(x, bs="cs"),size=I(2),level=0.99)
  p <- p + ggtitle(plot_name)
  
  
  if (plot_name != ""){
    plot_name = paste(plot_name,"_",sep="")
  }
  
  if (save_plot){
    ggsave(paste(plot_name,paste(mult_peaks,plot_var,"transformation",plot_var_transform,sep="_"),".png",sep=""),p)
  }
  
  return(p)
  
}















































