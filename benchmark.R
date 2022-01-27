library(data.table)
library(pheatmap)
library(icarmc)
library(stringr)
library(cluster)
library(mclust)
library(WGCNA)
library(flashClust)
library(fastICA)
library(clValid)
library(ggplot2)

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)

#source("./methods/utils.R")
setwd('C:/Users/molot/Desktop/Projects/RProjects')
# ------------------------------------------------------------------------------
# Generate sample for benchmark
# ------------------------------------------------------------------------------

mPaths <- list.files("./maayan_creeds_subset",
                     recursive = TRUE, all.files = TRUE,
                     full.names = TRUE, pattern = "rda")
listGseSignature <- readRDS("gse_signature_list.rds")
listSpecies <- readRDS('metaDt.rds')
set.seed(1)
mPaths <- sample(mPaths, 100)




internal_metric <- function(cl, dist_mat) {
    eval <- silhouette(dist=dist_mat, x=cl)
    return(summary(eval)$avg.width)
}

external_metric <- function(df, pred_cl, real_cl) {
    max_fscore <- -1
    max_ind <- -1
    best_prec <- -1
    best_rec <- -1
    for (cl in 1:length(unique(pred_cl))) {
        genes <- rownames(df)[pred_cl == cl]
        prec <- overlap(genes, real_cl) / length(genes)
        rec <- overlap(genes, real_cl) / length(real_cl)
        fscore <- 2 * (prec * rec) / (prec + rec)
        if (is.na(fscore)) {
            fscore <- -10
        }
        if (max_fscore < fscore) {
            max_fscore <- fscore
            max_ind <- cl
            best_prec <- prec
            best_rec <- rec
        }
    }
    return(list('fscore'=max_fscore, 'cl'=max_ind, 
                'precision'=best_prec, 'recall'=best_rec))
}

overlap <- function(v1, v2) {
    cnt <- 0
    for (v in v1) {
        cnt <- cnt + ifelse(v %in% v2, 1, 0)
    }
    return(cnt)
}


convert_to_entrezid <- function(geneset, species) {
    if (species == 'mouse') {
        return(mapIds(org.Mm.eg.db, keys=geneset, 
                    column="ENTREZID", keytype="SYMBOL"))
    } else {
        return(mapIds(org.Hs.eg.db, keys=geneset, 
                    column="ENTREZID", keytype="SYMBOL"))
    }
}
# ------------------------------------------------------------------------------
# MODELS
# ------------------------------------------------------------------------------

get_kmeans_metrics <- function(df, real_geneset, diss_mat) {
    kmeans_internal <- c()
    kmeans_prec <- c()
    kmeans_rec <- c()
    kmeans_fscore <- c()
    
    clusters <- c(4,5,6,8,10)
    for (ncluster in clusters) {
        kmeans_res <- kmeans(df, centers=ncluster)
        kmeans_cl <- kmeans_res$cluster
        kmeans_internal <- c(kmeans_internal, 
                             internal_metric(kmeans_cl, diss_mat))
        
        ext <- external_metric(df, kmeans_cl, real_geneset)
        kmeans_prec <- c(kmeans_prec, ext$precision)
        kmeans_rec <- c(kmeans_rec, ext$recall)
        kmeans_fscore <- c(kmeans_fscore, ext$fscore)
    }
    int <- max(kmeans_internal)
    ext <- max(kmeans_fscore)
    best_int <- clusters[which.max(kmeans_internal)]
    best_ext <- clusters[which.max(kmeans_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}

get_kmed_metrics <- function(df, real_geneset, diss_mat) {
    kmeds_internal <- c()
    kmeds_prec <- c()
    kmeds_rec <- c()
    kmeds_fscore <- c()

    clusters <- c(4,5,6,8,10)
    for (ncluster in clusters) {
        kmeds_res <- pam(x=df, diss=diss_mat, k=ncluster)
        kmeds_cl <- kmeds_res$clustering
        kmeds_internal <- c(kmeds_internal, 
                            internal_metric(kmeds_cl, diss_mat))
        ext <- external_metric(df, kmeds_cl, real_geneset)
        kmeds_prec <- c(kmeds_prec, ext$precision)
        kmeds_rec <- c(kmeds_rec, ext$recall)
        kmeds_fscore <- c(kmeds_fscore, ext$fscore)
    }

    int <- max(kmeds_internal)
    ext <- max(kmeds_fscore)
    best_int <- clusters[which.max(kmeds_internal)]
    best_ext <- clusters[which.max(kmeds_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}

get_GMM_metrics <- function(df, real_geneset, diss_mat) {
    GMM_internal <- c()
    GMM_prec <- c()
    GMM_rec <- c()
    GMM_fscore <- c()

    clusters <- c(4,5,6,8,10)
    for (ncluster in clusters) {
        GMM_res <- Mclust(df, G=ncluster)
        GMM_cl <- GMM_res$classification
        GMM_internal <- c(GMM_internal, 
                          internal_metric(GMM_cl, diss_mat))
        ext <- external_metric(df, GMM_cl, real_geneset)
        GMM_prec <- c(GMM_prec, ext$precision)
        GMM_rec <- c(GMM_rec, ext$recall)
        GMM_fscore <- c(GMM_fscore, ext$fscore)
    }

    int <- max(GMM_internal)
    ext <- max(GMM_fscore)
    best_int <- clusters[which.max(GMM_internal)]
    best_ext <- clusters[which.max(GMM_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}


get_WGCNA_metrics <- function(df, real_geneset, diss_mat) {
    
    datExpr <- t(df)
    softPower <- 6
    
    adjacency <- adjacency(datExpr, power = softPower, type = "signed")
    
    TOM <- TOMsimilarity(adjacency, TOMType="signed") # specify network type
    dissTOM <- 1-TOM
    geneTree <- flashClust(d = as.dist(dissTOM), method='complete')
    
    WGCNA_internal <- c()
    WGCNA_prec <- c()
    WGCNA_rec <- c()
    WGCNA_fscore <- c()

    clusters <- c(4,5,6,8,10)
    for (ncluster in clusters) {
        WGCNA_cl <- cutree(geneTree, k=ncluster)
        WGCNA_internal <- c(WGCNA_internal, 
                            internal_metric(WGCNA_cl, diss_mat))
        ext <- external_metric(df, WGCNA_cl, real_geneset)
        WGCNA_prec <- c(WGCNA_prec, ext$precision)
        WGCNA_rec <- c(WGCNA_rec, ext$recall)
        WGCNA_fscore <- c(WGCNA_fscore, ext$fscore)
    }

    int <- max(WGCNA_internal)
    ext <- max(WGCNA_fscore)
    best_int <- clusters[which.max(WGCNA_internal)]
    best_ext <- clusters[which.max(WGCNA_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}

get_hclust_metrics <- function(df, real_geneset, diss_mat) {
    geneTree <- flashClust(d = as.dist(diss_mat), method='complete')
    
    hclust_internal <- c()
    hclust_prec <- c()
    hclust_rec <- c()
    hclust_fscore <- c()

    clusters <- c(4,5,6,8,10)
    for (ncluster in clusters) {
        hclust_cl <- cutree(geneTree, k=ncluster)
        hclust_internal <- c(hclust_internal, 
                             internal_metric(hclust_cl, diss_mat))
        ext <- external_metric(df, hclust_cl, real_geneset)
        hclust_prec <- c(hclust_prec, ext$precision)
        hclust_rec <- c(hclust_rec, ext$recall)
        hclust_fscore <- c(hclust_fscore, ext$fscore)
    }

    int <- max(hclust_internal)
    ext <- max(hclust_fscore)
    best_int <- clusters[which.max(hclust_internal)]
    best_ext <- clusters[which.max(hclust_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}


nrun_ICA <- function(df, ncomp, nrun) {
    opt_metric <- 10**10
    opt_res <- NA
    for (i in 1:nrun) {
        res <- fastICA(df, ncomp)
        approx <- res$S %*% res$A
        if (opt_metric > mean((approx - df)^2)) {
            opt_metric <- mean((approx - df)^2)
            opt_res <- res
        }
    }
    return(opt_res)
}
only_outlliers <- function(v) {
    return(ifelse(abs(v - mean(v)) > 2*sd(v), v/abs(v), 0))
}

get_ICA_metrics <- function(df, real_geneset, diss_mat) {
    res <- nrun_ICA(df, getNC(df), 50)
    B <- res$S
    W <- res$A

    # delete the intercept vector
    ind = -1
    max_range = -1
    for (i in 1:5) {
        range <- quantile(B[, i])[4] - quantile(B[, i])[2]
        if (max_range < range) {
            max_range <- range
            ind <- i
        }  
    }
    B <- B[, -c(ind)]
    
    # only leave outliers with +1 / -1 coefficients
    B <- apply(B, 2, only_outlliers)
    
    ICA_internal <- c()
    ICA_prec <- c()
    ICA_rec <- c()
    ICA_fscore <- c()
    
    clusters <- 4:10#c(4,5,6,8,10)
    for (ncluster in clusters) {
        ICA_cl <- kmeans(B, ncluster, nstart=5)$cluster
        ICA_internal <- c(ICA_internal, 
                          internal_metric(ICA_cl, diss_mat))
        ext <- external_metric(df, ICA_cl, real_geneset)
        ICA_prec <- c(ICA_prec, ext$precision)
        ICA_rec <- c(ICA_rec, ext$recall)
        ICA_fscore <- c(ICA_fscore, ext$fscore)
    }

    int <- max(ICA_internal)
    ext <- max(ICA_fscore)
    best_int <- clusters[which.max(ICA_internal)]
    best_ext <- clusters[which.max(ICA_fscore)]
    return(list(intMetric=int, extMetric=ext, 
                intNC=best_int, extNC=best_ext))
}
    
# ------------------------------------------------------------------------------
# Calculate models
# ------------------------------------------------------------------------------



#######KMEANS##########
resdt <- rbindlist(lapply(mPaths[1:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_kmeans_metrics(E, geneset, dist(E))

        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'km')
        
        return(res)
        
        
        
    }, error = function(e){
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
write.csv(resdt, file='metrics/metrics_kmeans.txt')


#######ICA##########
resdt_ica <- rbindlist(lapply(mPaths[1:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_ICA_metrics(E, geneset, dist(E))
        
        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'ica')
        
        return(res)
        
        
        
    }, error = function(e){
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
write.csv(resdt_ica, file='metrics/metrics_ICA3.txt')


#######HCLUST##########
resdt_hc <- rbindlist(lapply(mPaths[1:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_hclust_metrics(E, geneset, dist(E))
        
        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'hc')
        
        return(res)
        
        
        
    }, error = function(e){
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
write.csv(resdt_hc, file='metrics/metrics_hclust2.txt')


#######KMEDOIDS##########
resdt_kmed <- rbindlist(lapply(mPaths[1:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_kmed_metrics(E, geneset, dist(E))
        
        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'kmed')
        
        return(res)
        
        
        
    }, error = function(e){
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
write.csv(resdt_kmed, file='metrics/metrics_kmeds.txt')



#######WGCNA##########
for (i in 15:100) {
  resdt_wgcna <- rbindlist(lapply(mPaths[i:i], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
      load(mPath)
      
      E <- t(scale(t(E), scale = FALSE))
      
      km <- get_WGCNA_metrics(E, geneset, dist(E))
      
      res <- data.table(gse = GSE, ncol = ncol(E), 
                        intNC = km$intNC, 
                        intMetric = km$intMetric, 
                        extNC = km$extNC, 
                        extMetric = km$extMetric, 
                        method = 'wgcna')
      
      return(res)
      
      
      
    }, error = function(e){
      print(e)
      return(data.table(gse = character(), ncol = numeric(),
                        intNC = numeric(), intMetric = numeric(),
                        extNC = numeric(), extMetric = numeric(),
                        method = character()))
    })
    
  }))
  if (nrow(resdt_wgcna) == 1) {
    write.csv(resdt_wgcna, 
              file=paste('metrics/metrics_wgcna',i,'.txt', sep=''))
  }
}

enableWGCNAThreads()
resdt_wgcna <- rbindlist(lapply(mPaths[31:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_WGCNA_metrics(E, geneset, dist(E))
        
        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'wgcna')
        
        return(res)
        
        
        
    }, error = function(e){
        print(e)
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
resdt_wgcna
write.csv(resdt_wgcna, file='metrics/metrics_wgcna4.txt')


#######GMM##########
resdt_gmm <- rbindlist(lapply(mPaths[1:100], function(mPath){
    GSEGPL <- str_extract(mPath, "GSE\\d+-GPL\\d+")
    GSE <- str_extract(GSEGPL, "GSE\\d+")
    GPL <- str_extract(GSEGPL, "GPL\\d+")
    species <- listSpecies[listSpecies$geo_id == GSE]$organism[[1]]
    geneset <- listGseSignature[[GSE]][[1]]
    geneset <- convert_to_entrezid(geneset, species)
    tryCatch({
        load(mPath)
        
        E <- t(scale(t(E), scale = FALSE))
        
        km <- get_GMM_metrics(E, geneset, dist(E))
        
        res <- data.table(gse = GSE, ncol = ncol(E), 
                          intNC = km$intNC, 
                          intMetric = km$intMetric, 
                          extNC = km$extNC, 
                          extMetric = km$extMetric, 
                          method = 'gmm')
        
        return(res)
        
        
        
    }, error = function(e){
        return(data.table(gse = character(), ncol = numeric(),
                          intNC = numeric(), intMetric = numeric(),
                          extNC = numeric(), extMetric = numeric(),
                          method = character()))
    })
    
}))
write.csv(resdt_gmm, file='metrics/metrics_gmm.txt')



resdt_km <- read.csv('metrics/metrics_kmeans.txt', row.names = 1)
resdt_ica <- read.csv('metrics/metrics_ICA.txt', row.names = 1)
resdt_kmed <- read.csv('metrics/metrics_kmeds.txt', row.names = 1)

resdt_wgcna1 <- read.csv('metrics/metrics_wgcna2.txt', row.names = 1)
resdt_wgcna2 <- read.csv('metrics/metrics_wgcna3.txt', row.names = 1)
resdt_wgcna3 <- read.csv('metrics/metrics_wgcna4.txt', row.names = 1)
resdt_wgcna <- rbind(resdt_wgcna1, resdt_wgcna2, resdt_wgcna3)
rownames(resdt_wgcna) <- 1:nrow(resdt_wgcna)

resdt_hc <- read.csv('metrics/metrics_hclust2.txt', row.names = 1)
resdt_gmm <- read.csv('metrics/metrics_gmm.txt', row.names = 1)

pdt <- rbind(resdt_km, resdt_ica, resdt_kmed, resdt_wgcna, resdt_hc, resdt_gmm)
# ------------------------------------------------------------------------------
# Plot Results
# ------------------------------------------------------------------------------
ggplot(pdt, aes(x=method, y=intMetric)) +
    geom_boxplot(aes(fill=method))

ggplot(pdt, aes(x=method, y=extMetric)) +
    geom_boxplot(aes(fill=method))

ggplot(pdt, aes(x=method, y=extNC)) +
    geom_boxplot(aes(fill=method))

time <- data.table(method=c('km', 'ica', 'kmed', 'wgcna', 'hclust', 'gmm'),
                   datasets=c(nrow(resdt), nrow(resdt_ica), nrow(resdt_kmed), 
                              nrow(resdt_wgcna), nrow(resdt_hc), nrow(resdt_gmm)),
                   time=c(20, 26, 263, 840, 22, 117))

time$avg.time <- time$time / time$datasets

time








