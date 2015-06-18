'Trabalho ECBD ' 



source("http://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.1.0.st.v1")
library(pd.hugene.1.0.st.v1)
#Diretoria dos ficheiros a analisar
setwd("~/GitHub/Data_Analysis/dataset")
########### 1- Carregamento dos Dados ############
library(oligo)
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)
#The Robust Multichip Average (RMA) 
eset <- rma(affyRaw)
# save the data to an output file (Data will be log2 transformed and normalized)
write.exprs(eset,file="data.txt")
my_frame <- data.frame(exprs(eset))
View(my_frame)
dim(my_frame)
class(my_frame)
featureNames(eset)[1:5]
sampleNames(eset)[1:5]
varMetadata(eset) # nao tem descrição das amostras
phenoData(eset) #nao tem imformação
annotation(eset) #"pd.hugene.1.0.st.v1"
experimentData(eset) #nao tem imformação
abstract(eset) #nao tem imformação


# The package contains an SQLite database. This database is accessible through a connection
conn <- db(pd.hugene.1.0.st.v1)
dbListTables(conn)
dbListFields(conn, 'featureSet')
dbListFields(conn, 'pmfeature')
sql <- 'SELECT * FROM pmfeature INNER JOIN featureSet USING(fsetid)'
probeInfo <- dbGetQuery(conn, sql)
probeInfo[1:10, 1:3]
head(probeInfo)

########### 2- Pre-processamento dos Dados ############
#Não temos valores omissos nos dados em analise.
sum(is.na(my_frame$GSM1446286_Can1.CEL))
sum(is.nan(my_frame$GSM1446286_Can1.CEL))
sum(is.na(my_frame$GSM1446287_Str1.CEL))
sum(is.nan(my_frame$GSM1446287_Str1.CEL))
sum(is.na(my_frame$GSM1446288_Tot1.CEL))
sum(is.nan(my_frame$GSM1446288_Tot1.CEL))
sum(is.na(my_frame$GSM1446289_Can2.CEL))
sum(is.nan(my_frame$GSM1446289_Can2.CEL))
sum(is.na(my_frame$GSM1446290_Str2.CEL))
sum(is.nan(my_frame$GSM1446290_Str2.CEL))
sum(is.na(my_frame$GSM1446291_Tot2.CEL))
sum(is.nan(my_frame$GSM1446291_Tot2.CEL))
sum(is.na(my_frame$GSM1446292_Can3.CEL))
sum(is.nan(my_frame$GSM1446292_Can3.CEL))
sum(is.na(my_frame$GSM1446293_Str3.CEL))
sum(is.nan(my_frame$GSM1446293_Str3.CEL))
sum(is.na(my_frame$GSM1446294_Tot3.CEL))
sum(is.nan(my_frame$GSM1446294_Tot3.CEL))

# Normalização, Background Corrections, Pm correction
library(genefilter)
sds=rowSds(my_frame)#calcula o desvio padrao por linha
sds[1:15]
m=median(sds)
m
mean(sds)
hist(sds, breaks=20, col="mistyrose")  						
abline(v=m, col="blue", lwd=4, lty=2)
abline(v=m*2, col="red", lwd=4, lty=2)
new_frame=my_frame[sds >= 2*m, ]


maximos=apply(my_frame,1,max)
maximos
#max value of gene expression 
max(maximos)
minimos=apply(my_frame,1,min)
#min value of gene expression 
min(minimos)
vl=maximos/minimos>2
new_frame2=my_frame[vl,]

## keep top 50 percent
filter=varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
frame_var_filter <- data.frame(exprs(filter))


## Expressão diferencial
require(Biobase)
object<-new("ExpressionSet", exprs=as.matrix(new_frame2))
object
tt = rowttests(object)
tt
#New dataframe ordered by column p value, ascending
pvalueorder = tt[order(tt$p.value),]


#cluster
eucD = dist(exprs(object[1:20])) 
cl.hier <- hclust(eucD)
plot(cl.hier) 

cl.hier <- hclust(eucD, method="single")
plot(cl.hier)

cl.hier <- hclust(eucD, method="average")
plot(cl.hier)


heatmap(exprs(object[1:20]), labCol = F)

km = kmeans(exprs(object[1:20]), 3) 
names(km)
km$cluster



