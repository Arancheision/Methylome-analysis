###############################################################################
#This function takes an object of the type MEDIPS.meth and merge adjacent     #
#windows only if they belong to the same gene. It gives results in a          #
#data.frame containing the following columns: chr, start, stop, windows,      #
#ensembl transcript name, HUGO name, and the rest of the columns from the     #
#MEDIPS.meth output                                                           #
###############################################################################

merge.by.gene<- function ( MEDIPS.object ){

  
  
  #Merge by gene and take the mean of the value
  result.merge <- aggregate.data.frame( MEDIPS.object[,4:14], 
                                        by= list(MEDIPS.object$gene), 
                                        FUN= mean)
  
  #Map those regions to their chromosome
  merge.chr <- aggregate.data.frame( MEDIPS.object$chr, 
                                     by= list(MEDIPS.object$gene), 
                                     FUN= unique)
  
  #Take the start of the first window and the stop of the last window merged
  merge.start <- aggregate.data.frame( MEDIPS.object$start, 
                                       by= list(MEDIPS.object$gene), 
                                       FUN= min)
  
  merge.stop <- aggregate.data.frame( MEDIPS.object$stop, 
                                      by= list(MEDIPS.object$gene),
                                      FUN= max)
  
  merge.windows<- aggregate.data.frame( MEDIPS.object$gene, 
                                        by= list(MEDIPS.object$gene), 
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
  
  
  #Sort them by logFC
  
  result.merge <- result.merge[order(result.merge$edgeR.logFC, 
                                     decreasing = TRUE),]
  #names(result.merge)[1]<- "gene"
  
  #return<- result.merge[, c(11:14,1:10)] #put columns in a logical order

  return<- result.merge
}

