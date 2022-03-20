####Soundscape networks#####

library(tuneR)
library(seewave)
library(reshape2)
library(plyr)
library(vegan)
library(Rtsne)
library(fpc)
library(ggplot2)
library(igraph)
library(parallel)

####Version 2. For multiple files from a single folder####
##Including parallelization. Oh shit. Let's go##

##Where are the files
setwd("~/Desktop/Demasiado grandes para dropbox/datos_escorial/ESCORIAL_2/")
directory <- "~/Desktop/Demasiado grandes para dropbox/datos_escorial/ESCORIAL_2/"

#wav_files <- as.list(dir(path = directory, pattern = "wav$", ignore.case = TRUE))
#audio_files <- vector("list", length(wav_files))
#names(audio_files) <- as.character(wav_files)

wav_files <- list.files(pattern=".wav", recursive=F)

##Parallelization parameters##
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

###THE BIG LOOP####

result <- parLapply(cl, wav_files, function(x){
  library(tuneR)
  library(seewave)
  library(reshape2)
  library(plyr)
  library(vegan)
  library(Rtsne)
  library(fpc)
  library(igraph)
  library(parallel)

#for(j in 1:length(wav_files)) {

  ####Cut file into 1s chunks and put into a list####
  #After noise removal, that I have already done#
  
  duration = 1
  start_times = seq(0, 420, by = duration) #esto son los 8 minutos, para evitar los audios cortos.
  tr1 <- vector("list", sum(start_times == start_times))
  
  for (i in seq_along(start_times)) {
    tr1[[i]] = readWave(filename = x,
                        from = start_times[i],
                        to = start_times[i] + duration,
                        units = 'seconds')
    #b <- start_times[i]
    #a <- tr1[[i]]
    #writeWave(a, filename=paste("201229_0155_", as.character(b), ".wav", sep = ""))
  }
  names(tr1)<- paste(x, start_times, sep="_")
  
  
  ####Estimate frequency spectrum for these chunks####
  tr2 <- lapply(tr1, function(x){spec(x, f=48000, channel=1, wl=512, PMF=TRUE, fftw=TRUE, plot=FALSE, from=0)})
  tr2 <- lapply(tr2, function(x){mutate(as.data.frame(x), rounded=as.data.frame(x)[,1])})
  for(i in 1:length(tr2)) {
    tr2[[i]]$rounded <- round(tr2[[i]]$rounded, digits=2)
  }
  tr2 <- lapply(tr2, function(x){aggregate(x$y~x$rounded, FUN=mean)}) #este último paso tarda un pelín#
   
  ####Reshape list of frequencies to a standard dataframe####
  spectrum_results <- data.frame(matrix(NA, nrow = nrow(as.data.frame(tr2[[1]]$`x$rounded`)), ncol = length(tr2)))
  
  for(i in 1:length(tr2)) {
    spectrum_results$frecuencia <- as.data.frame(tr2[[1]])[,1]
    row.names(spectrum_results) <- as.data.frame(tr2[[1]])[,1]
    spectrum_results[[i]] <- as.data.frame(tr2[[i]])[,2]
  }
  
  names(spectrum_results) <- c(names(tr2), "frecuencia")
  names(spectrum_results)
  
  plot(spectrum_results[[3]], log="y", type="l")
  
  ####Calculate chunk intensity, shannon, and filtering metrics: I may not use them later####
  frequency_filter <- subset(spectrum_results, spectrum_results$frecuencia > 0.5)
  frequency_filter <- subset(frequency_filter, spectrum_results$frecuencia < 17) #uso 17 xq he puesto eso en script_escorial_preliminar_analysus
  chunk_intensity <- colSums(frequency_filter[,1:ncol(frequency_filter)-1])
  plot(chunk_intensity)
  #chunk_intensity <- 1-chunk_intensity
  
  ####Standardize: logaritmic + relative####
  frequency_filter$frecuencia <- NULL
  frequency_filter_log <- 100-abs(log(frequency_filter))
  spectrum.decostand <- decostand(frequency_filter_log, method="total", MARGIN=2)
  #spectrum.decostand <- decostand(frequency_filter, method="total", MARGIN=2)
  colSums(spectrum.decostand)
  spectrum.decostand <- t(spectrum.decostand)
  #chunk_shannon <-  diversity(spectrum.decostand, index="shannon")
  #plot(chunk_shannon)
  
  ####Remove chunks that do not have calls or are very silent####
  #Perhaps with an index, equitativity or total intensity (filter those with the lowest intensity) or something like that#
  #or quantiles based on intensity thresohlds#
  
  #tofilter_intensity <- chunk_intensity[chunk_intensity > quantile(chunk_intensity, probs=0.8)]
  #tofilter_intensity <- flatness[flatness < quantile(flatness, probs=0.8)]
  #spectrum_subset <- spectrum.decostand[names(tofilter_intensity),]
  
  ####Find similar calls (tsne) and make clusters (pamk)####
  #spectrum.dist <- vegdist(spectrum_subset, method="bray")
  spectrum.dist <- vegdist(spectrum.decostand, method="bray")
  rtsne_out <- Rtsne(as.matrix(spectrum.dist), is_distance=T, pca = T, dims=2, verbose = TRUE, theta=0.5, perplexity=nrow(spectrum.decostand)^(1/2), pca_scale=T, pca_center=T, max_iter=2000)
  rtsne_plot <- as.data.frame(cbind(rtsne_out$Y[,1], rtsne_out$Y[,2]))
  plot(rtsne_plot)
  row.names(rtsne_plot) <- labels(spectrum.dist)
  
  #Intermediate stept to visually examine
  spectrum.tsne.points <- merge(rtsne_plot, chunk_intensity, by.x="row.names", by.y="row.names")
  
  order.pamk <- pamk((scale(spectrum.tsne.points[,c(2:4)])), krange=10:20, usepam=T, criterion="ch") 
  #order.pamk$pamobject[3] #if clara with multiasw is used
  order.pamk <- as.data.frame(order.pamk$pamobject[3])
  
  spectrum.tsne.points.clust <- cbind.data.frame(spectrum.tsne.points, order.pamk)
  #spectrum.tsne.points.clust.bis <- subset(spectrum.tsne.points, spectrum.tsne.points$y > quantile(chunk_intensity, probs=0.80))
  
  ####Matrix: time-distances####
  ab <- strsplit(spectrum.tsne.points$Row.names, split="_")
  names(ab) <- spectrum.tsne.points$Row.names
  for(i in 1:length(ab)) {
    ab[[i]] <- ab[[i]][4] #OJO, en escorial1 aquí poner un CINCO.
  }
  ab <- melt(ab)
  ab$value <- as.numeric(ab$value)
  row.names(ab) <- ab$value
  
  timedist.matrix <- as.matrix(vegdist(ab$value, method="euclidean"))
  timedist.matrix[lower.tri(timedist.matrix, diag=T)] <- NA
  timedist.matrix <- as.data.frame(as.matrix(timedist.matrix))
  
  names(timedist.matrix) <- ab[,2]
  row.names(timedist.matrix) <- ab[,2]
  
  ####Matrix: distribution of clusters vs. time####
  clusters.time.matrix <- melt(as.matrix(timedist.matrix))
  clusters.time.matrix <- merge(clusters.time.matrix, spectrum.tsne.points.clust[,c(1,5)], by.x="Var1", by.y="Row.names")
  clusters.time.matrix <- merge(clusters.time.matrix, spectrum.tsne.points.clust[,c(1,5)], by.x="Var2", by.y="Row.names")
  
  ####Filter distances above above a threshold. For example, farther than 2 seconds####
  clusters.time.matrix <- subset(clusters.time.matrix, clusters.time.matrix$value < 3)
  clusters.time.matrix <- subset(clusters.time.matrix, clusters.time.matrix$clustering.x != clusters.time.matrix$clustering.y)
  clusters.time.matrix$value <- 3-clusters.time.matrix$value
  
  ####Average by cluster, to use only those for clustering####
  clusters.time.matrix.2 <- aggregate(clusters.time.matrix$value~clusters.time.matrix$clustering.x+clusters.time.matrix$clustering.y, FUN=sum) #como una suma ponderada. Las veces que aparece, multiplicado por el peso (1 o 2 segundos de distancia)
  names(clusters.time.matrix.2)[3] <- "weight"
  
  #Filter those weights that are too low. Ideally a probability. For now, below the sample mean.
  clusters.time.matrix.2 <- subset(clusters.time.matrix.2, clusters.time.matrix.2$weight >= mean(clusters.time.matrix.2$weight))
  
  ####Transform edgelist into a network, and estimate modularity####
  network <- graph_from_data_frame(clusters.time.matrix.2, directed=F)
  is_weighted(network)
  plot(network)
  module.value <- modularity(cluster_walktrap(network)) #sí que está cogiendo los pesos.
  nr.modules <- max(membership(cluster_walktrap(network)))
  Sys.time()
  
  ####Store modularity, number of modules, clustering coefficient, and number of pamkclusters####
  nr.clusters <- max(V(network))
  nr.components <- components(network)$no
  aveplen <- mean_distance(network, directed = "FALSE")
  
    ####Create matrix to store info for each element of the loop
  Table <- data.frame(filename=x, nr.clusters = nr.clusters, nr.components=nr.components, module.value=module.value, nr.modules=nr.modules, aveplen=aveplen)
  #audio_files[[j]] <- cbind(nr.clusters, nr.components, module.value, nr.modules, aveplen)
  return(Table)
})
stopCluster(cl)

library(dplyr)
properties <- bind_rows(result)
#prop96 <- properties
#prop1 <- properties
#prop2 <- properties

#properties <- melt(audio_files, level=1)
#properties <- reshape(properties[,c(2:4)], direction="wide", timevar = "Var2", idvar = "L1")
