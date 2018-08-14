wrong_sam <- strsplit(readLines("../../data/CMS_001_RNA_A_S1_aedes_aligned.sam"), "\t")
right_sam <- strsplit(readLines("../../data/CMS_001_RNA_A_S1_aligned.sam"), "\t")

wrong_sam_val <- sapply(wrong_sam, function (x) as.numeric(gsub("NM:i:", "", grep("NM:", x, value=TRUE))))
right_sam_val <- sapply(right_sam, function (x) as.numeric(gsub("NM:i:", "", grep("NM:", x, value=TRUE))))

summary(wrong_sam_val)
summary(right_sam_val)

t.test(wrong_sam_val, right_sam_val)

mean(wrong_sam_val==0)
mean(right_sam_val==0)
