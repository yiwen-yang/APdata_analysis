rm(list = ls())
library(dplyr)
library(tidyfst)
library(limma)
library(ggplot2)
library(openxlsx)

# simply read the corresponding files
toread_csv <- function(stage, stage_number){
  files <- data.frame()
  for (i in 1:stage_number){
    filedir <- paste0('lDNB_output//Max_dnb_score_module in ', stage, ' for sample ',i,'.csv')
    fileread <- read.csv(filedir)
    files <- rbind(files, fileread)
  }
  return(files)
}

# compute the DNB score in every stage
# Increase reliability, prevent interference of sample size, take average
tosumscore <- function(data, freqnum, stagenum, k){
  data <- na.omit(data)
  data %>% group_by(gene_name) %>% summarise(freq = n()) %>% filter(freq > freqnum*stagenum) -> data_tmp
  data_score <- data[which(data$gene_name %in% data_tmp$gene_name),]
  if(nrow(data_score) > 2){data_score <- avereps(data_score, ID = data_score$gene_name)}else {data_score = data_score}
  data_numeric <- apply(data_score[, -1], 2, as.numeric) %>% as.data.frame()
  if(nrow(data_numeric) >10){data_numeric <- data_numeric[order(data_numeric$score,decreasing = T),]}
  totalscore <- ifelse(nrow(data_numeric) > k, sum(data_numeric[1:k,]$score)/k, sum(data_numeric$score)/nrow(data_numeric))
  return(totalscore)
}

scores.sum <- data.frame()
for (i in seq(10,120,10)){
  k = i 
  Maxscore_health <- toread_csv('Health', 5)
  Maxscore_MAP <- toread_csv('MAP', 56)
  Maxscore_SAP <- toread_csv('SAP', 12)
  freqnum = 0.5
  health_score <- tosumscore(Maxscore_health, 0.1, 5 ,k)
  MAP_score <- tosumscore(Maxscore_MAP, freqnum, 56 ,k)
  SAP_score <- tosumscore(Maxscore_SAP, freqnum, 12 ,k)
  AP.score <- data.frame(k,health_score, MAP_score, SAP_score)
  scores.sum <- rbind(scores.sum, AP.score)
}

health_filter <- Maxscore_health %>%
  group_by(gene_name) %>%
  summarise(freq = n(), 
            score = sum(score) / freq) %>%
  filter(freq > 0.2 * 5) %>%
  na.omit()
map_filter <- Maxscore_MAP %>%
  group_by(gene_name) %>%
  summarise(freq = n(),
            score = sum(score) / freq) %>%
  filter(freq > 0.4 * 56) %>%
  na.omit()
sap_filter <- Maxscore_SAP %>%
  group_by(gene_name) %>%
  summarise(freq = n(),
            score = sum(score) / freq) %>%
  filter(freq > 0.5 * 12) %>%
  na.omit()

write.csv(health_filter, file = 'DNBproteinAndDNBscore_from_health.csv',quote = F)
write.csv(map_filter, file = 'DNBproteinAndDNBscore_from_map.csv',quote = F)
write.csv(sap_filter, file = 'DNBproteinAndDNBscore_from_sap.csv',quote = F)
 
# use DNB gene ----------------
tosumDNBscore <- function(data, freqnum, stagenum, k){
  data <- na.omit(data)
  data_score <- data[which(data$gene_name %in% map_filter$gene_name),]
  if(nrow(data_score) > 2){data_score <- avereps(data_score, ID = data_score$gene_name)}else {data_score = data_score}
  data_numeric <- apply(data_score[, -1], 2, as.numeric) %>% as.data.frame()
  if(nrow(data_numeric) >10){data_numeric <- data_numeric[order(data_numeric$score,decreasing = T),]}
  totalscore <- ifelse(nrow(data_numeric) > k, sum(data_numeric[1:k,]$score)/k, sum(data_numeric$score)/nrow(data_numeric))
  return(totalscore)
}
scoresDNB.sum <- data.frame()
for (i in seq(10,120,10)){
  k = i 
  freqnum = 0.5
  healthDNB_score <- tosumDNBscore(Maxscore_health, 0.1, 5 ,k)
  MAPDNB_score <- tosumDNBscore(Maxscore_MAP, freqnum, 56 ,k)
  SAPDNB_score <- tosumDNBscore(Maxscore_SAP, freqnum, 12 ,k)
  APDNB.score <- data.frame(k,healthDNB_score, MAPDNB_score, SAPDNB_score)
  scoresDNB.sum <- rbind(scoresDNB.sum, APDNB.score)
}

#read background network
allproteinlist <-  read.csv('raw_data/reference_sample.csv')
allproteinlist <- allproteinlist[,1:2]
allproteinlist$Health <- 0
allproteinlist$Healthscore <- 0
allproteinlist$map <- 0
allproteinlist$mapscore <- 0
allproteinlist$msap <- 0
allproteinlist$msapscore <- 0
allproteinlist$sap <- 0
allproteinlist$sapscore <- 0
allproteinlist <- allproteinlist[,-2]
healthallproteinscore <- na.omit(Maxscore_health) 
healthallproteinscore <- healthallproteinscore[healthallproteinscore$gene_name %in% health_filter$gene_name,] %>% avereps(., ID = .$gene_name) %>% as.data.frame() 
mapallproteinscore <- na.omit(Maxscore_MAP) 
mapallproteinscore <- mapallproteinscore[mapallproteinscore$gene_name %in% map_filter$gene_name,] %>% avereps(., ID = .$gene_name) %>% as.data.frame() 
sapallproteinscore <- na.omit(Maxscore_SAP)
sapallproteinscore <- sapallproteinscore[sapallproteinscore$gene_name %in% sap_filter$gene_name,] %>% avereps(., ID = .$gene_name) %>% as.data.frame() 
for (num in 1:nrow(allproteinlist)){
  if (allproteinlist[num,1] %in% health_filter$gene_name){
    allproteinlist[num,"Health"] <- 1
    allproteinlist[num, 'Healthscore'] <-  healthallproteinscore[healthallproteinscore$gene_name == allproteinlist[num, 'X'], 'score']}
}
  
for (num in 1:nrow(allproteinlist)){
  if (allproteinlist[num,1] %in% map_filter$gene_name){
    allproteinlist[num,"map"] <- 1
    allproteinlist[num, 'mapscore'] <-  mapallproteinscore[mapallproteinscore$gene_name == allproteinlist[num, 'X'], 'score']}
}

for (num in 1:nrow(allproteinlist)){
  if (allproteinlist[num,1] %in% sap_filter$gene_name){
    allproteinlist[num,"sap"] <- 1
    allproteinlist[num, 'sapscore'] <-  sapallproteinscore[sapallproteinscore$gene_name == allproteinlist[num, 'X'], 'score']}
}
sum(allproteinlist$Health)
sum(allproteinlist$map)
sum(allproteinlist$sap)

#------------edge-edge ---------------
ssnedgetoread_csv <- function(stage, stage_number){
  files <- data.frame()
  for (i in 1:stage_number){
    filedir <- paste0('lDNB_output//SSN for ', stage, ' in sample ',i,'.csv')
    fileread <- read.csv(filedir)
    fileread <- fileread[order(abs(fileread$deltaPCC),decreasing = T),]
    fileread <- fileread[abs(fileread$deltaPCC) > 0.1,]
    print(nrow(fileread))
    files <- rbind(files, fileread)
  }
  return(files)
}

tosumedge <- function(data, freqnum, stagenum){
  data <- na.omit(data)
  data$node <- paste0(data$node1,data$node2)
  data %>% group_by(node) %>% summarise(freq = n()) %>% filter(freq > freqnum*stagenum) -> data_tmp
  data_ssn <- data[which(data$node %in% data_tmp$node),]
  data_ssn <- data_ssn[!duplicated(data_ssn$node),]
    return(data_ssn)
}

health_ssn <- ssnedgetoread_csv('Health',5)
MAP_ssn <- ssnedgetoread_csv('MAP',56)
SAP_ssn <- ssnedgetoread_csv('SAP',12)
health_ssn.dnb <- tosumedge(health_ssn, 0.2, 5)
MAP_ssn.dnb <- tosumedge(MAP_ssn, 0.2, 56)
SAP_ssn.dnb <- tosumedge(SAP_ssn, 0.4, 12)
write.csv(health_ssn.dnb,file = 'health_ssn_dnb.csv',row.names = F, quote = F)
write.csv(MAP_ssn.dnb,file = 'MAP_ssn.dnb.csv',row.names = F, quote = F)
write.csv(SAP_ssn.dnb,file = 'SAP_ssn.dnb.csv',row.names = F, quote = F)

stringnetwork <- read.table('raw_data/string_interactions0.7AB.tsv') %>% .[,1:2]
names(stringnetwork) <- c('node1','node2')
stringnetwork$node12 <- paste0(stringnetwork$node1, stringnetwork$node2)
for (num in 1:nrow(stringnetwork)){
  if (stringnetwork[num,3] %in% health_ssn.dnb$node){
    stringnetwork[num,"health"] <- health_ssn.dnb[health_ssn.dnb$node == stringnetwork[num, 'node12'], 'deltaPCC']}
}
 
for (num in 1:nrow(stringnetwork)){
  if (stringnetwork[num,3] %in% MAP_ssn.dnb$node){
    stringnetwork[num,"MAP"] <- MAP_ssn.dnb[MAP_ssn.dnb$node == stringnetwork[num, 'node12'], 'deltaPCC']}
} 


for (num in 1:nrow(stringnetwork)){
  if (stringnetwork[num,3] %in% SAP_ssn.dnb$node){
    stringnetwork[num,"SAP"] <-  SAP_ssn.dnb[SAP_ssn.dnb$node == stringnetwork[num, 'node12'], 'deltaPCC']}
}
stringnetwork[is.na(stringnetwork)] <- 0
sum(stringnetwork$health != 0)
sum(stringnetwork$MAP != 0)
sum(stringnetwork$SAP != 0)
stringnetwork$health <- abs(stringnetwork$health)
stringnetwork$MAP <- abs(stringnetwork$MAP)
stringnetwork$SAP <- abs(stringnetwork$SAP)

write.csv(stringnetwork,file = 'AP_deltaPCC_raw.csv',quote = F,row.names = F)