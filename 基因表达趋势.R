BiocManager::install('Mfuzz')
library(Mfuzz)
library(marray)
x <- read.table("E:\\1_科研计划\\10_中国玫瑰孤儿基因OGs\\基因表达趋势分析\\salt.txt",header=T)
count <- data.matrix(x)
eset <- new("ExpressionSet",exprs = count)
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = 'mean')
eset <- filter.std(eset, min.std = 0)
eset <- standardise(eset)
c <- 6
m <- mestimate(eset)
cl <- mfuzz(eset, c = 6, m = m)
mfuzz.plot(eset, cl = cl, mfrow = c(2, 3), time.labels = c("CK","2h","24h","48h"))
mfuzzColorBar(main="Membership",cex.main=1)
#每个簇下基因数量
cl$size
#每个基因所属簇
head(cl$cluster)
#基因和 cluster 之间的 membership，用于判断基因所属簇，对应最大值的那个簇
head(cl$membership)
#整合关系输出
gene_cluster <- cbind(cl$cluster, cl$membership)
colnames(gene_cluster)[1] <- 'cluster'
write.table(gene_cluster, 'E:\\1_科研计划\\10_中国玫瑰孤儿基因OGs\\基因表达趋势分析\\salt_gene_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)
View(gene_cluster)
x <- read.table("E:\\1_科研计划\\10_中国玫瑰孤儿基因OGs\\基因表达趋势分析\\drought.txt",header=T)
count <- data.matrix(x)
eset <- new("ExpressionSet",exprs = count)
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = 'mean')
eset <- filter.std(eset, min.std = 0)
eset <- standardise(eset)
c <- 6
m <- mestimate(eset)
cl <- mfuzz(eset, c = 6, m = m)
mfuzz.plot(eset, cl = cl, mfrow = c(2, 3), time.labels = c("CK","30d","60d","90d"))
gene_cluster <- cbind(cl$cluster, cl$membership)
colnames(gene_cluster)[1] <- 'cluster'
write.table(gene_cluster, 'E:\\1_科研计划\\10_中国玫瑰孤儿基因OGs\\基因表达趋势分析\\drought_gene_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)