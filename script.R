#source("http://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.1.0.st.v1")
#library(pd.hugene.1.0.st.v1)
#biocLite("affy")
#library(affy)
#data(pd.hugene.1.0.st.v1)
#biocLite("oligo")

#Diretoria dos ficheiros a analisar
setwd("~/GitHub/Data_Analysis/dataset")

library(oligo)
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)
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
annotation(eset) #nao tem imformação
experimentData(eset) #nao tem imformação
abstract(eset) #nao tem imformação

# pre processamento de dados
library(genefilter)
sds=rowSds(my_frame)#calcula o desvio padrao por linha
m=median(sds)
hist(sds, breaks=50, col="mistyrose")
sum(is.na(my_frame$GSM1446286_Can1.CEL))
sum(is.nan(myframe$GSM1446286_Can1.CEL))




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
