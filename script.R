source("http://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.1.0.st.v1")
library(pd.hugene.1.0.st.v1)
#biocLite("affy") # do not load this library, only use oligo for our case
#library(affy) # do not load this library, only use oligo for our case
#biocLite("oligo")

#Diretoria dos ficheiros a analisar
setwd("~/GitHub/Data_Analysis/dataset")

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


# pre processamento de dados
library(genefilter)
sds=rowSds(my_frame)#calcula o desvio padrao por linha
sds[1:15]
m=median(sds)
m
mean(sds)
hist(sds, breaks=20, col="mistyrose")
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
  						
abline(v=m, col="blue", lwd=4, lty=2)
abline(v=m*2, col="red", lwd=4, lty=2)
new_frame=my_frame[sds >= 2*m, ]


maximos=apply(my_frame,1,max)
maximos
max(maximos)
minimos=apply(my_frame,1,min)
min(minimos)


# #possivelmente lixo
# #media de valores
# mean(data$X) 
# 
# 
# 
# class(data) #classe do ficheiro
# typeof(data) #tipo do dataset
# nrow(data) #numero de linhas
# ncol(data) #numero de colunas
# sapply(data, class) #verificar a classe de cada coluna
# sapply(data, typeof) #verificar o tipo de valor que cada coluna tem
# sum(is.na(data)) #verificar se tem numeros omissos - Nao tem valores omissos
# summary(data$X) #resumo de todos os atributos
# 
# 
# # Distribuicao normal do valor  de pelo gene
# dens = density(amostra1)
# m = mean(data$X)
# hist(data$X,probability=T,col=gray(.9),main="Analysis of compartment-specific gene expression in breast cancer tumors",xlab="Value")
# lines(dens, col = "blue")
# abline(v=m, col = "green")
# curve(dnorm(x,mean(data$X),sd(data$X)),add=T,col="red")
# 
# 
# 
