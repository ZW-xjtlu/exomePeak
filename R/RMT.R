

RMT <- function(
  INPUT_BAM,
  IP_BAM,
  INPUT2IP = NA, 
  GENE_ANNO_GTF = NA, 
  GENOME = NA, 
  UCSC_TABLE_NAME = "knownGene",
  TXDB = NA,
  EXOME_OUTPUT_DIR = NA,
  RMT_OUTPUT_DIR = NA,
  RMT_EXPERIMENT_NAME = "RMT_Result") {
  # Wrap parameters ##################################################
  PARAMETERS = list();
  PARAMETERS$INPUT_BAM = INPUT_BAM
  PARAMETERS$IP_BAM = IP_BAM
  PARAMETERS$INPUT2IP = INPUT2IP
  PARAMETERS$EXOME_OUTPUT_DIR = EXOME_OUTPUT_DIR
  PARAMETERS$GO_OUTPUT_DIR = RMT_OUTPUT_DIR
  PARAMETERS$GO_EXPERIMENT_NAME = RMT_EXPERIMENT_NAME
  PARAMETERS$GENE_ANNO_GTF=GENE_ANNO_GTF
  PARAMETERS$GENOME = GENOME
  PARAMETERS$UCSC_TABLE_NAME=UCSC_TABLE_NAME
  PARAMETERS$TXDB = TXDB
  
  
  # dependent variables
  if (is.na(PARAMETERS$EXOME_OUTPUT_DIR)) {
    PARAMETERS$EXOME_OUTPUT_DIR <- getwd()
  }
  if (is.na(PARAMETERS$GO_OUTPUT_DIR)) {
    PARAMETERS$GO_OUTPUT_DIR <- getwd()
  }
  
  # algrithm ##################################################
  #peak calling
  input_length <- length(PARAMETERS$INPUT_BAM)
  ip_length <- length(PARAMETERS$IP_BAM)
  all_filepath <- vector()
  
  # decide whether paired
  temp1 <- is.na(PARAMETERS$INPUT2IP)
  flag <- temp1[1]
  
  
  if(flag){
    for(i in 1:ip_length){
      experiment_name = paste(PARAMETERS$GO_EXPERIMENT_NAME,"/temp/IP_rep_",i,sep="")
      all_filepath[i] = PARAMETERS$IP_BAM[i]
      print(paste("Processing IP sample",i))
      temp <- exomepeak(GENE_ANNO_GTF = GENE_ANNO_GTF, 
                        GENOME = GENOME, UCSC_TABLE_NAME = UCSC_TABLE_NAME,
                        TXDB = TXDB, IP_BAM = c(PARAMETERS$IP_BAM[i]), 
                        INPUT_BAM = c(PARAMETERS$INPUT_BAM[i]), 
                        OUTPUT_DIR = PARAMETERS$EXOME_OUTPUT_DIR, 
                        POISSON_MEAN_RATIO = 1,
                        EXPERIMENT_NAME = experiment_name)
    }
  } else {
    for(i in 1:ip_length){
      experiment_name = paste(PARAMETERS$GO_EXPERIMENT_NAME,"/temp/IP_rep_",i,sep="")
      input_bam = PARAMETERS$INPUT_BAM[PARAMETERS$INPUT2IP[[i]][1]]
      if(length(PARAMETERS$INPUT2IP[[i]])>1){
        for(j in 2:length(PARAMETERS$INPUT2IP[[i]])){
          input_bam=c(input_bam,c(PARAMETERS$INPUT_BAM[PARAMETERS$INPUT2IP[[i]][j]]))
        }
      }      
      all_filepath[i] <- paste(PARAMETERS$EXOME_OUTPUT_DIR,experiment_name,sep = '/')
      print(paste("Processing IP sample",i))
      temp <- exomepeak(GENE_ANNO_GTF = GENE_ANNO_GTF, 
                        GENOME = GENOME, UCSC_TABLE_NAME = UCSC_TABLE_NAME,
                        TXDB = TXDB, IP_BAM = c(PARAMETERS$IP_BAM[i]), 
                        INPUT_BAM = input_bam, 
                        OUTPUT_DIR = PARAMETERS$EXOME_OUTPUT_DIR, 
                        POISSON_MEAN_RATIO = 1, 
                        EXPERIMENT_NAME = experiment_name)   
    }
  }
  
  #merge all peak
  for(i in 1:length(all_filepath)){
    path <- paste(all_filepath[i],"peak.xls",sep = '/')
    bam_file <- read.table(path,sep = "\t",header = FALSE)
    if(i == 1){
      all_peak <- bam_file
    } else {
      all_peak = rbind(all_peak,bam_file[2:nrow(bam_file), ])
    } 
  }
  dir.create(paste(PARAMETERS$GO_OUTPUT_DIR,
                   PARAMETERS$GO_EXPERIMENT_NAME,sep = '/'),
             recursive = TRUE,showWarnings = FALSE)
  dir = paste(PARAMETERS$GO_OUTPUT_DIR,
              PARAMETERS$GO_EXPERIMENT_NAME,sep = '/')
  write.table(all_peak,paste(dir,"all_peak.xls",sep = '/'),
              sep = "\t",quote = FALSE,col.names = FALSE,
              row.names = FALSE)
  
  #read all peaks and transform GRangesList
  filepath = paste(dir,"all_peak.xls",sep = '/')
  matrix_peak_read = read.table(filepath,header = TRUE,stringsAsFactors = FALSE)
  
  ori_peak <- .xls2Grangeslist(filepath)
  
  #remove overlaps
  overlaps <- findOverlaps(ori_peak,ori_peak)
  matrix_overlap <- matrix(0,length(overlaps),2)
  matrix_overlap[,1] <- queryHits(overlaps)
  matrix_overlap[,2] <- subjectHits(overlaps)
  
  num <- c(rep(0,length(ori_peak)))
  
  for(i in 1:nrow(matrix_overlap)){
    if(matrix_overlap[i, 1][1][[1]] != matrix_overlap[i, 2][1][[1]]){
      x <- ranges(ori_peak[matrix_overlap[i, 1]])[[1]]@width
      y <- ranges(ori_peak[matrix_overlap[i, 2]])[[1]]@width
      #num[i]=if(x<y) matrix_overlap[i,1][1][[1]] else matrix_overlap[i,2][1][[1]]
      if(x < y){
        num[matrix_overlap[i, 1][1][[1]]] <- matrix_overlap[i, 1][1][[1]] 
      } else if(x > y){
        num[matrix_overlap[i, 2][1][[1]]] <- matrix_overlap[i, 2][1][[1]]
      } else {
        num[matrix_overlap[i, 1][1][[1]]] <- matrix_overlap[i, 1][1][[1]]
        num[matrix_overlap[i, 2][1][[1]]] <- matrix_overlap[i, 2][1][[1]]
      }
    }
  }
  peak <- ori_peak[num == 0]
  tab <- which(num == 0)
  write.table(matrix_peak_read[tab, ],paste(dir,"merged_peak.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  write.table(matrix_peak_read[tab,1:12 ],paste(dir,"merged_peak.bed",sep = '/'),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  
  
  # read gene annotation from internet
  txdb <- .readTxDb2(PARAMETERS)
  exonRanges <- exonsBy(txdb, "gene")
  
  #rpkm in all inputbams and ipbams
  rpkm_input <- matrix(0,length(peak),length(PARAMETERS$INPUT_BAM))
  rpkm_ip <- matrix(0,length(peak),length(PARAMETERS$IP_BAM))
  count_input <- matrix(0,length(peak),length(PARAMETERS$INPUT_BAM))
  count_ip <- matrix(0,length(peak),length(PARAMETERS$IP_BAM))
  #rpkm in inutbams
  for(i in 1:length(PARAMETERS$INPUT_BAM)){
    filepath <- PARAMETERS$INPUT_BAM[i]
    aligns <- readGAlignments(filepath)
    para <- ScanBamParam(what="mapq")
    mapq <- scanBam(filepath, param=para)[[1]][[1]]
    # filter reads with mapq smaller than 30. 
    mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
    ID_keep <- (mapq >30)
    filtered <- aligns[ID_keep]
    id <- countOverlaps(filtered,exonRanges)
    transcriptome_filtered_aligns <- filtered[id>0]
    counts <- countOverlaps(peak, transcriptome_filtered_aligns)
    count_input[,i] <- counts
    numBases <- sum(width(peak))
    geneLengthsInKB <- numBases / 1000
    millionsMapped <- length(transcriptome_filtered_aligns) / 1000000
    rpm <- counts / millionsMapped
    rpkm_input[,i] <- rpm / geneLengthsInKB
    
  }
  write.table(count_input,paste(dir,"count_input.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  write.table(rpkm_input,paste(dir,"rpkm_input.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  #rpkm in ipbams
  for(i in 1:length(PARAMETERS$IP_BAM)){
    filepath <- PARAMETERS$IP_BAM[i]
    aligns <- readGAlignments(filepath)
    para <- ScanBamParam(what="mapq")
    mapq <- scanBam(filepath, param=para)[[1]][[1]]
    # filter reads with mapq smaller than 30. 
    mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
    ID_keep <- (mapq >30)
    filtered <- aligns[ID_keep]
    id <- countOverlaps(filtered,exonRanges)
    transcriptome_filtered_aligns <- filtered[id>0]
    counts <- countOverlaps(peak, transcriptome_filtered_aligns)
    count_ip[,i] <- counts
    numBases <- sum(width(peak))
    geneLengthsInKB <- numBases / 1000
    millionsMapped <- length(transcriptome_filtered_aligns) / 1000000
    rpm <- counts / millionsMapped
    rpkm_ip[,i] <- rpm / geneLengthsInKB 
  }
  write.table(count_ip,paste(dir,"count_ip.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  write.table(rpkm_ip,paste(dir,"rpkm_ip.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  #return to all conditions
  matrix_log2_fc = matrix(0,length(peak),ip_length)
  if(is.na(PARAMETERS$INPUT2IP)){
    for(i in 1:ip_length){
      matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(rpkm_input[,i]+0.01)
    }
  } else {
    for(i in 1:ip_length){
      if(length(PARAMETERS$INPUT2IP[[i]])==1){
        matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(rpkm_input[,PARAMETERS$INPUT2IP[[i]]]+0.01)
      } else {
        sum = matrix(0,length(peak),1)
        for(j in 1:length(PARAMETERS$INPUT2IP[i])){
          sum = sum + rpkm_input[,PARAMETERS$INPUT2IP[[i]][j]]
        }
        matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(sum + 0.01)
      }
    }
  }
  
  #feature selection
  v <- matrix_log2_fc[,1]
  for(i in 1:nrow(matrix_log2_fc)) {v[i] <- var(matrix_log2_fc[i,])}
  # hist(v)
  matrix_log2_fc <- cbind(tab,matrix_log2_fc[1:nrow(matrix_log2_fc),],matrix_peak_read[tab,])
  write.table(matrix_log2_fc,paste(dir,"all_information.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  # delete un-necessary files
  file.remove(paste(dir,"all_peak.xls",sep = '/'))
  file.remove(paste(dir,"all_information.xls",sep = '/')) 

  # notification
  print("***************************")
  print("***************************")
  print(paste("The combinatorial RNA methylome is generated under:",dir))
  print("Including:")
  print("1. Merged Peaks")
  print("2. Merged Peaks in BED format")
  print("3. Reads count for every peak in every bam file")
  print("4. RPKM for every peak in every bam file")
  }


# read gtf file
.readTxDb2 <- function(PARAMETERS){
  # download the annotation
  if ((!is.na(PARAMETERS$GENOME)) & (!is.na(PARAMETERS$UCSC_TABLE_NAME)) & is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB)) {
    op <- options(warn = (-1))
    txdb =makeTxDbFromUCSC(genome=PARAMETERS$GENOME,
                           tablename=PARAMETERS$UCSC_TABLE_NAME)
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB) ) {
    op <- options(warn = (-1))
    txdb=makeTxDbFromGFF(PARAMETERS$GENE_ANNO_GTF,format="gtf")
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$TXDB) ) {
    txdb=PARAMETERS$TXDB
  }
  
  
  # return data
  return(txdb)
}
.xls2Grangeslist <- function(filepath) {
    # read bed file
    a = read.table(filepath,sep="\t",header=TRUE,stringsAsFactors =FALSE)
    mcols_info =a[,13:length(a[1,])]
    a = a[,1:12]
    
    # get transcripts
    no_tx = length(a[,1])
    tx_id = 1:no_tx;
    tx_name = paste("exomepeak_",1:no_tx,sep="")
    tx_chrom = a[,1]
    tx_strand = a[,6]
    tx_start = a[,2]+1
    tx_end = a[,3]
    transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
    head(transcripts)
    
    # get splicings
    splicing = data.frame()
    for (i in 1:no_tx) {
      tx = a[i,]
      tx_id = i
      exon_rank=1:as.integer(tx[10])
      
      # get start
      temp = as.integer(strsplit(as.character(tx[12]), ",")[[1]]) + tx_start[i]
      exon_start=temp
      
      # get end
      temp = as.integer(strsplit(as.character(tx[11]), ",")[[1]])
      temp2 = temp + exon_start - 1
      exon_end=temp2
      
      # get CDS
      cds_start = exon_start
      cds_end = exon_end
      
      # get data frame
      splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
      
      # collect result
      splicing = rbind(splicing, splicing_tx)
    }
    
    # get genes
    tx_name = tx_name
    gene_id = as.character(a[,4])
    gene_id[is.na(gene_id)]="NA"
    gene=data.frame(tx_name,gene_id)
    
    # 
    
    peaks = makeTxDb(transcripts=transcripts, splicings=splicing,
                             genes=gene)
    tx <- exonsBy(peaks, "tx",use.names=TRUE)
    mcols(tx) <- mcols_info
    
    return(tx)
  }
