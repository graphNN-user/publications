
library(dplyr)
library(enrichR)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- dbs$libraryName

run_enrichr <- function(weblive, dbs_in, csvfile_in) {
    genes <- read.csv(csvfile_in, header=T)
    genes <- genes$x
   if (weblive) {
    enrichr_res <- enrichr(genes, dbs_in )
    return(enrichr_res)
    }
   print("enrichr website not live!")
}

columns= c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes","database")
myRes = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(myRes) = columns

process_enrichr <- function(res_enrichr) {
  myRes = data.frame(matrix(nrow = 0, ncol = length(columns)))
  i = 1
  for (db in res_enrichr ) {
    dbnames <- names(res_enrichr[i])
    db$database <- dbnames
    i = i+1
    myRes <- rbind(myRes, db)
  }
  myRes2 <- select(myRes, database,Term,Overlap,P.value,Adjusted.P.value,Combined.Score,Genes)
  myRes3 <- filter(myRes2, Adjusted.P.value < 0.05) 
  myRes4 <- arrange(myRes3,Adjusted.P.value)
  myRes4
}

res1 <- run_enrichr(websiteLive, dbs, "subset4_C_neg_L_neg_M_pos_union_up_sorted_v6.csv" ) 
res1b <- process_enrichr(res1)
write.csv(res1b, file="subset4_C_neg_L_neg_M_pos_union_up_v6_enrichr.csv", row.names=F)

res2 <- run_enrichr(websiteLive, dbs, "subset4_C_neg_L_neg_M_pos_union_down_sorted_v6.csv" )
res2b <- process_enrichr(res2)
write.csv(res2b, file="subset4_C_neg_L_neg_M_pos_union_down_v6_enrichr.csv", row.names=F)








q()

if (websiteLive) dbs <- listEnrichrDbs()
dbs <- dbs$libraryName
genes <- read.csv("subset1_C_pos_L_pos_M_neg_union_down_sorted_v6.csv", header=T)
genes <- genes$x

if (websiteLive) {
    subset1_C_pos_L_pos_M_neg_union_down_v6 <- enrichr(genes, dbs)
}

genes <- read.csv("subset1_C_pos_L_pos_M_neg_union_up_sorted_v6.csv", header=T)
genes <- genes$x

if (websiteLive) {
    subset1_C_pos_L_pos_M_neg_union_up_v6 <- enrichr(genes, dbs)
}


genes <- read.csv("M_pos_C_neg_L_neg_doublecutoff_up_v2.csv", header=F)
genes <- genes$V1

if (websiteLive) {
    enrichr_M_pos_C_neg_L_neg_doublecutoff_up_v2 <- enrichr(genes, dbs)
}

genes <- read.csv("M_pos_C_neg_L_neg_doublecutoff_down_v2.csv", header=F)
genes <- genes$V1

if (websiteLive) {
    enrichr_M_pos_C_neg_L_neg_doublecutoff_down_v2 <- enrichr(genes, dbs)
}


columns= c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes","database") 
myRes = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(myRes) = columns

process_enrichr <- function(res_enrichr) {
  myRes = data.frame(matrix(nrow = 0, ncol = length(columns)))
  i = 1
  for (db in res_enrichr ) {
    dbnames <- names(res_enrichr[i])
    db$database <- dbnames
    i = i+1
    myRes <- rbind(myRes, db)
  }
  myRes2 <- select(myRes, database,Term,Overlap,P.value,Adjusted.P.value,Combined.Score,Genes)
  myRes3 <- filter(myRes2, Adjusted.P.value < 0.05) 
  myRes4 <- arrange(myRes3,Adjusted.P.value)
  myRes4
}

res1 <- process_enrichr(enrichr_M_pos_C_neg_L_neg_doublecutoff_down_v2)
write.table(res1, file="enrichr_M_pos_C_neg_L_neg_doublecutoff_down_v2.txt", sep="\t", row.names=F)

res2 <- process_enrichr(enrichr_M_pos_C_neg_L_neg_doublecutoff_up_v2)
write.table(res2, file="enrichr_M_pos_C_neg_L_neg_doublecutoff_up_v2.txt", sep="\t", row.names=F)

res3 <- process_enrichr(enrichr_M_pos_C_neg_L_neg_union_3list_up_v2)
write.table(res3, file="enrichr_M_pos_C_neg_L_neg_union_3list_up_v2.txt", sep="\t", row.names=F)

res4 <- process_enrichr(enrichr_M_pos_C_neg_L_neg_union_3list_down_v2)
write.table(res4, file="enrichr_M_pos_C_neg_L_neg_union_3list_down_v2.txt", sep="\t", row.names=F)

exit()



columns= c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value""Odds.Ratio","Combined.Score","Genes","database")                         
myRes = data.frame(matrix(nrow = 0, ncol = length(columns)))  
colnames(myRes) = columns

enrichr_M_pos_C_neg_L_neg_union_3list_up_v2_filtered2 <- lapply(enrichr_M_pos_C_neg_L_neg_union_3list_up_v2_filtered, function(x) filter(x, x$Adjusted.P.value < 0.05 ))



cutoff <- function(x) {
  x[1,"Adjusted.P.value"] < 0.05
}
enrichr_M_pos_C_neg_L_neg_union_3list_up_v2_filtered <- filter(lambda x: x[1,"Adjusted.P.value"] < 0.05, enrichr_M_pos_C_neg_L_neg_union_3list_up_v2) 





results <- c(enrichr_M_pos_C_neg_L_neg_union_3list_up_v2, enrichr_M_pos_C_neg_L_neg_union_3list_down_v2, enrichr_M_pos_C_neg_L_neg_doublecutoff_up_v2, enrichr_M_pos_C_neg_L_neg_doublecutoff_down_v2) 

saveRDS(results, file="enrichR_res.rds")






indices <- sapply(enrichr_M_pos_C_neg_L_neg_union_3list_up_v2, function(x) x[1,"Adjusted.P.value"] < 0.05 )
enrichr_M_pos_C_neg_L_neg_union_3list_up_v2_filtered <- enrichr_M_pos_C_neg_L_neg_union_3list_up_v2[indices]

