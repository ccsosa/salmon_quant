library(JASPAR2018)
library(rGADEM)
library(TFBSTools)
library(GNET2)
library(bears)


#https://compgenomr.github.io/book/motif-discovery.html


pwm_library = getMatrixSet(
  JASPAR2018,
  opts=list(
    collection = 'CORE',
    species    = 'Oryza sativa',
    matrixtype = 'PWM'
  ))



# find the most similar motif to our motif
pwm_sim = PWMSimilarity(
  
  # JASPAR library
  pwm_library, 
  
  # out motif
  unknown_pwm,
  
  # measure for comparison
  method = 'Pearson')

# library(BSgenome.Hsapiens.UCSC.hg19)
# path<- system.file("extdata/Test_100.fasta",package="rGADEM")
# FastaFile<-paste(path,sep="")
# Sequences <- readDNAStringSet(FastaFile, "fasta")
# gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)
# path<- system.file("extdata/jaspar2009.txt",package="rGADEM")
# seededPwm<-readPWMfile(path)
# grep("STAT1",names(seededPwm))
# STAT1.PWM=seededPwm[103]
# gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens,Spwm=STAT1.PWM)
