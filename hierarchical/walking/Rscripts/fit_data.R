library(PARCOR)
######## load dataset #########
root_dir <- getwd()

######## subject id #######
subject_id <- 6
load(file = paste0(getwd(), "/data/S", subject_id,".RData"))


####### name of cluster area #######
cluster_area <- EEG_data$stand_pull$names

####### potential model order ########
P <- 10


####### sample size ##########
sample_size <- 1000
####### construct discount factor #########
delta <- seq(0.995, 0.999, by = 0.001)
delta_matrix <- as.matrix(expand.grid(delta, delta))
result <- hparcor(yt=EEG_data$stand_pull$data[1, , ], P=P, delta = delta_matrix, 
                  sample_size=sample_size, chains=1, DIC = TRUE, uncertainty = FALSE)
