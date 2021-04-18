# We run original idaFast on our skeleton files
require(devtools)
require(graph)
load_all("../../pcalg/",
                          quiet=T,reset=T,export_all=T)
message(date()," >> Rversion: ")
print(version)

Ngenes <- c(100, 200, 500, 1000, 2000, 5000)
Ngenes <- c(200)
#############################################

# Read the design, so that we can separate Myc high from Myc Low
designFile <- paste0("01-005-00_Ergebnis_Uebersicht.csv")
descr.mat <- read.csv(designFile,header = TRUE, 
                      sep = "\t", dec = ",",row.names=1)
# get samples with Myc low and Myc high
myc_samples <- data.frame( High=as.numeric(rownames(descr.mat[descr.mat$Myc==1,])), Low =  as.numeric(rownames(descr.mat[descr.mat$Myc==0,])))


tiempoOPT <- c()
tiempo <- c()
i = 1
for (Ng in Ngenes){
   bigVarGenes <- readRDS(paste0("var_genes/Ng",Ng,"alpha0p5.rds"))
   for (myc in c("Low","High")){
     skel_file <- paste0("skeleton_files/Ng", 
                  Ng,".it001.Myc",myc,".alpha0p5.rds")

     MycHighLow <- bigVarGenes[myc_samples[,myc], ] 
     l <- readRDS(skel_file)
     subDat <- MycHighLow[l$ind,]
     Ns <- nrow(MycHighLow)   # number of samples
     subSize <- floor(Ns/2) 
     tiempoOPT <- c(tiempoOPT,Sys.time())
     message(Sys.time())
     mysol <- idaFastOpt(mcov=cov(subDat),graphEst=l$cpdag_graph) 
     message(Sys.time())
     tiempoOPT[i] <- Sys.time() - tiempoOPT[i]
     tiempo <- c(tiempo,Sys.time())
     oldsol <- list()
     for (j in 1:Ng){
         oldsol[[j]] <-  idaFast(j,1:Ng,cov(subDat),l$cpdag_graph)
     }
     message(Sys.time())
     tiempo[i] <- Sys.time() - tiempo[i]
     cat(Ng," Myc",myc,"idaFastOpt:", tiempoOPT[i],
         attr(tiempoOPT[i],"units"),"\n")
     cat(Ng," Myc",myc,"idaFast:", tiempo[i],
         attr(tiempo[i],"units"),"\n")
     i <- i + 1 
   }
}

