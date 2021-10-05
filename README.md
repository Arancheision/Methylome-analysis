# Methylome-analysis#
Important!! it is work in progress!!
This is a set of scripts to analyse MDB-seq data using MEDIPs R package.
The script named MethylomeAnalysis starts out with the BAM files to generate the MEDIPs object "mr.edgeR" with the DMRs between control and treated.
The DMR_CPGs takes only DMRs that are inside CpG islands.
The CREs Analysis overlaps the DMRs of the mr.edgeR with the annotated chromatin regulatory regions from the ENDCODE project
