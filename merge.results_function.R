
###############################################################################
#This function takes an object of the type MEDIPS.meth and merge adjacent     #
#windows only if they belong to the same gene region. It gives results in a   #
#data.frame containing the following columns: chr, start, stop, windows,      #
#ensembl transcript name, HUGO name, and the rest of the columns from the     #
#MEDIPS.meth output                                                           #
###############################################################################

merge.regions<- function ( MEDIPS.object ){
  
  #Merge by promoter and take the mean of the values
  result.merge <- aggregate.data.frame( MEDIPS.object[,4:10], 
                                        by= list(MEDIPS.object$ROI), 
                                        FUN= mean)
  
  #Map those regions to their chromosome
  merge.chr <- aggregate.data.frame( MEDIPS.object$chr, 
                                     by= list(MEDIPS.object$ROI), 
                                     FUN= unique)
  
  #Take the start of the first window and the stop of the last window merged
  merge.start <- aggregate.data.frame( MEDIPS.object$start, 
                                       by= list(MEDIPS.object$ROI), 
                                       FUN= min)
  
  merge.stop <- aggregate.data.frame( MEDIPS.object$stop, 
                                      by= list(MEDIPS.object$ROI),
                                      FUN= max)
  
  merge.windows<- aggregate.data.frame( MEDIPS.object$ROI, 
                                        by= list(MEDIPS.object$ROI), 
                                        FUN= length)
  
  #sanity check:
  
  if ( identical(merge.chr$Group.1, merge.start$Group.1) &
       identical(merge.chr$Group.1, merge.stop$Group.1)  &
       identical(merge.chr$Group.1, merge.windows$Group.1) &
       identical(merge.chr$Group.1, result.merge$Group.1) == TRUE){
       result.merge$chr <- merge.chr$x
       result.merge$start <- merge.start$x
       result.merge$stop <- merge.stop$x
       result.merge$windows <- merge.windows$x
       
  } else {
    
    result.merge <- "error"
  }
  
  #############################################################################
  # In this part I make the table nicer and change the ENSEMBL names by       #
  # HUGO                                                                      #
  #############################################################################
  
  #Change the Ensembl number by the HUGO name
  
  library(biomaRt)
  
  transcript.names<-sub( "TSS_", "", result.merge$Group.1) #Find the transcript names
  result.merge$transcript<- transcript.names
  
  #Use biomaRt to obtain the HUGO names
  
  mart <- useDataset("hsapiens_gene_ensembl", 
                     useMart("ensembl"))
  
  hugo<- getBM ( filters= "ensembl_transcript_id", 
                 attributes= c("ensembl_transcript_id","hgnc_symbol"),
                 values = transcript.names, 
                 mart = mart)
  
  idx<- match(result.merge$transcript, hugo$ensembl_transcript_id)
  
  result.merge$gene<- hugo$hgnc_symbol[idx]
  
  
  #Sort them by logFC
  
  result.merge <- result.merge[order(result.merge$edgeR.logFC, 
                                     decreasing = TRUE),]
  
  return <- result.merge[, c(9:14,2:8)] #put columns in a logical order
  
}
  
