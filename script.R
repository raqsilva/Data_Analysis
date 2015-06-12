# #source("http://bioconductor.org/biocLite.R")
# biocLite("pd.hugene.1.0.st.v1")
# library(pd.hugene.1.0.st.v1)
# # biocLite("affy")
# #library(affy)
# data(pd.hugene.1.0.st.v1)

#Diretoria dos ficheiros a analisar
m=setwd("~/GitHub/Data_Analysis/dataset_processed")
data1 = read.fwf(file="GSM1446286_sample_table.txt",sep=" ",header= TRUE,widths=c(8,10))
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
dens = density(data$X)
m = mean(data$X)
hist(data$X,probability=T,col=gray(.9),main="Analysis of compartment-specific gene expression in breast cancer tumors",xlab="Value")
lines(dens, col = "blue")
abline(v=m, col = "green")
curve(dnorm(x,mean(data$X),sd(data$X)),add=T,col="red")



