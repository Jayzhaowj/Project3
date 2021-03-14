library(R.matlab)
root_dir <- getwd()
cond_types <- c("stand_pull", "stand_rotate", "walk_pull", "walk_rotate")


###### load the cluster ######
cluster <- readMat("./data/cluster_1222018.mat")

###### labels #######
labels <- c("left occipital", "right occipital", 
            "left sensorimotor", "anterior cingulate", 
            "right sensorimotor", "posterior parietal", 
            "supplementary motor", "anterior parietal")

index <- 1:33
index_not_in_list <- c(14, 18, 21)
set_index <- index[-index_not_in_list]

##### index of subjects ######
subjects_ind <- matrix(NA, nrow = 30, ncol = 8)
for(i in 3:10){
  for(j in 1:30){
    subjects_ind[j, i-2] <- j %in% as.numeric(names(table(cluster$cluster[, , i]$sets)))
  }
}

###### subjects index #######
subjects <- sapply(1:30, function(x) labels[subjects_ind[x, ]])


###### clean data function #######
clean_data <- function(file_name, subjects_index){
  eeg_data <- rep(list(list(data = NA, names = NA, times = NA)), 4)
  for(i in 1:4){
    path_set <- paste0(root_dir, "/data/", cond_types[i], "/", file_name)
    EEG <- readMat(path_set)
    var_names <- dimnames(EEG$EEG)[[1]]
    
    ### number of clusters (the last eight are EMG)
    n_chans <- EEG$EEG[[which(var_names == "nbchan")]] 
    ### number of epochs
    n_trials <- EEG$EEG[[which(var_names == "trials")]]
    ### number of time points
    times <- EEG$EEG[[which(var_names == "times")]]
    n_times <- length(times)
    ### load fdt file
    fdt_name <- unlist(EEG$EEG[[which(var_names=="datfile")]])
    path_fdt <- paste0(root_dir, "/data/", cond_types[i], "/", fdt_name)
    fdt_file <- file(path_fdt, "rb")
    signals <- readBin(fdt_file,
                       "double",
                       n = n_chans * n_trials * n_times,
                       size = 4,
                       endian = "little")
    close(fdt_file)

    if(n_chans != length(subjects[[subjects_index]])+8){
      print(paste0("Current condition: ", cond_types[i]))
      print(paste0("number of clusters in set: ", n_chans))
      print(paste0("number of clusters in clusters: ", length(subjects[[subjects_index]])+8))
      stop("differernt cluster number")
    }
    dim(signals) <- c(n_chans, n_times, max(n_trials, 1))
    eeg_data[[i]]$data <- signals[(1:length(subjects[[subjects_index]])), ,]
    eeg_data[[i]]$names <- subjects[[subjects_index]]
    eeg_data[[i]]$times <- times
    cat("\n The Perturbation condition: ", cond_types[i], "\n")
    cat("\n Number of cluster: ", n_chans, "\n")
    cat("\n Number of time points: ", n_times, "\n")
    cat("\n Number of epochs: ", n_trials, "\n")
  }
  names(eeg_data) <- cond_types
  return(eeg_data)
}

####### retrieve from subject 1 to 30 except for subject 8. #######
sink(file = paste0(getwd(), "/data/log.txt"))
iter_index <- c(1:7,9:30)
for(i in iter_index){
  file_name <- paste0("S", set_index[i], ".set")
  EEG_data <- clean_data(file_name = file_name, subjects_index = i)
  save("EEG_data", file = paste0(root_dir, "/data/S", i, ".RData"))
  cat("\n The subject ", i, " has been retrieved. \n")
  cat("\n ========================================== \n")
}
sink()




