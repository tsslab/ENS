#!/usr/bin/Rscript
library(DiffBind)
trunk <- dba(sampleSheet="diffbind_all.csv")
trunk_count <- dba.count(trunk, minOverlap=3)
trunk_consensus <- dba.peakset(trunk, consensus=-DBA_REPLICATE, minOverlap=3)
#Example of contrast pairwise comparison between sample 
trunk_contrast <- dba.contrast(trunk_count, group1=trunk_count$masks$E2pos12,
                                 group2=trunk_count$masks$neg12,
                                 name1="E12", name2="neg")
trunk_analyze <- dba.analyze(trunk_contrast_2, method=DBA_DESEQ2)
#Then generate report
report <- dba.report(trunk_analyze, contrast=1, method=DBA_DESEQ2, th=1, bCounts=TRUE)
write.table(as.data.frame(report), file="report.txt", sep="\t", quote=F, row.names=F)
