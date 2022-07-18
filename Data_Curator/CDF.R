library ("tidyverse")
library ("ggplot2")
whitelist_raw<-read.table("/projectnb/bf528/users/hedgehog_2022/project_4/Data_curator/other_output_files/white_list.txt")

# Order by counts
whitelist_raw <- as_tibble(whitelist_raw[order(whitelist_raw$V1), , drop=FALSE])

#plot
#png("cdf.png")
#plot(ecdf(whitelist_raw$V1), main='CDF of Barcode Frequency',xlab='Barcodes', ylab = "Frequency", xlim=c(0,600))
#dev.off()

whitelist<-dplyr::filter(whitelist_raw,V1>100)
whitelist<-whitelist[,-1]

whitelist1<-whitelist_raw[100:nrow(whitelist_raw),]


# Writing  
#write_csv(whitelist,"whitelist_strict_filter.csv", col_names=F)
#writeLines(as.character(whitelist1, , drop=FALSE])), file('whitelist.txt'), sep="\n")

