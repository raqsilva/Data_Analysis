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



#Cancer Cells 1
data1 = read.fwf(file="GSM1446286_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Stroma 1
data2 = read.fwf(file="GSM1446287_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Total 1  
data3 = read.fwf(file="GSM1446288_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Cancer Cells 2  
data4 = read.fwf(file="GSM1446289_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Stroma 2
data5 = read.fwf(file="GSM1446290_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Total 2
data6 = read.fwf(file="GSM1446291_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Cancer Cells 3
data7 = read.fwf(file="GSM1446292_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Stroma 3
data8 = read.fwf(file="GSM1446293_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#Total 3
data9 = read.fwf(file="GSM1446294_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
#amostra1
amostra1=data.frame(data3,data1$X,data2$X)
colnames(amostra1) <- c("ID_REF.VALUE  ", "Total 1 ","Cancer Cells 1","Stroma 1")

#amostra2
amostra2=data.frame(data6,data4$X,data5$X)
colnames(amostra2) <- c("ID_REF.VALUE  ", "Total 2 ","Cancer Cells 2","Stroma 2")
#amostra3
amostra3=data.frame(data9,data7$X,data8$X)
colnames(amostra3) <- c("ID_REF.VALUE  ", "Total 3 ","Cancer Cells 3","Stroma 3")
#media de valores
mean(data$X) 



class(data) #classe do ficheiro
typeof(data) #tipo do dataset
nrow(data) #numero de linhas
ncol(data) #numero de colunas
sapply(data, class) #verificar a classe de cada coluna
sapply(data, typeof) #verificar o tipo de valor que cada coluna tem
sum(is.na(data)) #verificar se tem numeros omissos - Nao tem valores omissos
summary(data$X) #resumo de todos os atributos


# Distribuicao normal do valor  de pelo gene
dens = density(amostra1)
m = mean(data$X)
hist(data$X,probability=T,col=gray(.9),main="Analysis of compartment-specific gene expression in breast cancer tumors",xlab="Value")
lines(dens, col = "blue")
abline(v=m, col = "green")
curve(dnorm(x,mean(data$X),sd(data$X)),add=T,col="red")



