'Trabalho ECBD ' 

source("http://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.1.0.st.v1")
library(pd.hugene.1.0.st.v1)
#Diretoria dos ficheiros a analisar
setwd("~/GitHub/Data_Analysis/dataset")
### 1- Carregamento dos Dados ####
library(oligo)
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)
#The Robust Multichip Average (RMA) 
eset <- rma(affyRaw)#normalização dos dados
# Guarda os dados em txt,os dados são transformados em log2 e normalizados
write.exprs(eset,file="data.txt")
my_frame <- data.frame(exprs(eset))
#ver data framme
View(my_frame)
dim(my_frame)
class(my_frame)#33297     9
featureNames(eset)[1:5]# "7892501" "7892502" "7892503" "7892504" "7892505"
sampleNames(eset)[1:5]# "GSM1446286_Can1.CEL" "GSM1446287_Str1.CEL" "GSM1446288_Tot1.CEL" "GSM1446289_Can2.CEL" "GSM1446290_Str2.CEL"
varMetadata(eset) # nao tem descrição das amostras
phenoData(eset) #nao tem informação
annotation(eset) #"pd.hugene.1.0.st.v1"
experimentData(eset) #nao tem informação
abstract(eset) #nao tem informação


###2- Pre-processamento dos Dados ###
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
sds[1:15] #Ver os primeiros 15 desvios padroes
m=median(sds)#mediana de desvios de padroes 0.2639786 
m 
mean(sds)#media de desvios de padroes 0.2639786 
hist(sds, breaks=20, col="mistyrose") # histograma 						
abline(v=m, col="blue", lwd=4, lty=2)#mediana 
abline(v=m*2, col="red", lwd=4, lty=2)#limite superior 2 vezes a mediana
new_frame=my_frame[sds >= 2*m, ]


maximos=apply(my_frame,1,max)# maximos do valor de expressao dos genes
maximos
minimos=apply(my_frame,1,min)# minimo do valor de expressao dos genes
#maximo valor de gene expression
max(maximos)
#minimo valor de gene expression 
min(minimos)

vl=maximos/minimos>2
new_frame2=my_frame[vl,]#Data frame filtrado com  genes cujo rácio do máximo valor sobre o mínimo valor de expressão seja superior a 2


## keep top 50 percent
filter=varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
frame_var_filter <- data.frame(exprs(filter))


## Expressão diferencial
require(Biobase)
object<-new("ExpressionSet", exprs=as.matrix(new_frame2))# transformar o newframe2 em  expression set
object
tt = rowttests(object)#Realiza os t-tests e verificar os p-values
tt


## Obter nomes dos genes
ob=featureNames(object)
#biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)
unlist(mget(ob, hugene10sttranscriptclusterSYMBOL))
unlist(mget(ob, hugene10sttranscriptclusterGENENAME))
unlist(mget(ob, hugene10sttranscriptclusterENTREZID))



#Novo dataframe ordenado pela coluna de p value
pvalueorder= tt[order(tt$p.value),]
pvalueorder$p.value[1:20]# primeiros 20 resultados com menor p value
rr=rownames(pvalueorder)[1:20]
unlist(mget(rr, hugene10sttranscriptclusterENTREZID))



dim(object)
### Clustering ### 
eucD = dist(exprs(object[65:85])) 
cl.hier <- hclust(eucD)
plot(cl.hier) 
bb=rownames(exprs(object[65:85]))
unlist(mget(bb, hugene10sttranscriptclusterENTREZID))



cl.hier <- hclust(eucD, method="single")
plot(cl.hier)

cl.hier <- hclust(eucD, method="average")
plot(cl.hier)

#heatmap
heatmap(exprs(object[1:20]), labCol = F)
#kmeans
km = kmeans(exprs(object[1:20]), 3) 
names(km)
km$cluster


