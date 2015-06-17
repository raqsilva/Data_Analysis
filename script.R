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


# 
# exp_m2 = exprs(new_frame2)
# > exprs(ALLm2) = scale(exp_m2) 
# > design = model.matrix(~ALLm2$mol.biol)
# > fit = lmFit(ALLm2,design)
# > fit2 = eBayes(fit)
# > diff = topTable(fit2, coef=2, 100)
# > indexes = as.numeric(rownames(diff))
# > ALLcl = ALLm2[indexes,]
# > ALLcl
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 100 features, 79 samples 
# 
# > train = t(exprs(ALLcl[,1:60]))
# > test = t(exprs(ALLcl[,61:79]))
# 
# > library(class)
# > valores.previstos = knn(train,test,ALLcl$mol.biol[1:60])
# > valores.reais = ALLcl$mol.biol[61:79]
# > sum(valores.previstos == valores.reais)/length(ALLcl$mol.biol[61:79])
# [1] 0.9473684
# > table(valores.previstos, valores.reais)
# valores.reais
# valores.previstos BCR/ABL NEG
# BCR/ABL       8   0
# NEG           1  10
# > valores.previstos = knn(train,test,ALLcl$mol.biol[1:60], k=3)
# > table(valores.previstos, valores.reais)
# valores.reais
# valores.previstos BCR/ABL NEG
# BCR/ABL       8   2
# NEG           1   8
# 
# 
# > library(nnet)
# > ann = nnet(ALLcl$mol.biol[1:60]~.,data.frame(train),size=3)
# # weights:  307
# initial  value 42.759543 
# iter  10 value 11.086923
# .
# iter  40 value 0.001757
# final  value 0.000071 
# converged
# > valores.prev.ann = predict(ann, data.frame(test), type="class")
# > table(valores.prev.ann, valores.reais)
# valores.reais
# valores.prev.ann BCR/ABL NEG
# BCR/ABL       8   4
# NEG           1   6
# > sum(valores.prev.ann == valores.reais)/length(valores.reais)
# [1] 0.7368421
# 
# > library(rpart)
# > arv = rpart(ALLcl$mol.biol[1:60]~.,data.frame(train))
# > plot(arv, uniform=T, branch=0.4,margin=0.1,compress=T)
# > text(arv,use.n=T,cex=0.9)
# > classes.previstas.arv = predict(arv, data.frame(test), type="class")
# > sum(classes.previstas.arv == valores.reais)/length(valores.reais)
# [1] 0.7368421
# 
# 
# > library(MLInterfaces)
# > knnResult <- MLearn(mol.biol~., ALLcl, knnI(k=1), 1:60)
# > confuMat(knnResult)
# predicted
# given     BCR/ABL NEG
# BCR/ABL       8   1
# NEG           0  10
# 
# > nnetResultLOO <- MLearn(mol.biol~., ALLcl, nnetI, xvalSpec("LOO"), size=3, decay=0.01)
# > confuMat(nnetResultLOO)
# predicted
# given     BCR/ABL NEG
# BCR/ABL      32   5
# NEG           5  37
# 
# > treeResultCV <- MLearn(mol.biol~., ALLcl, rpartI, xvalSpec("LOG", 5, balKfold.xvspec(5)))
# > confuMat(treeResultCV)
# predicted
# given     BCR/ABL NEG
# BCR/ABL      29   8
# NEG          18  24
# 
# > filt = nsFilter(ALLs, require.entrez=T, remove.dupEntrez=T, var.func=IQR, var.cutoff=0.5, feature.exclude="^AFFX")
# > ALLf = filt$eset
# > affyUniverse = featureNames(ALLf)
# > entrezUniverse = unlist(mget(affyUniverse, hgu95av2ENTREZID))
# > ttests = rowttests(ALLfilt_bcrneg, "mol.biol")
# > smPV = ttests$p.value < 0.05
# > pvalFiltered = ALLf[smPV, ]
# > selectedEntrezIds = unlist(mget(featureNames(pvalFiltered),
#                                   hgu95av2ENTREZID))
# 
# > params = new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse, annotation="hgu95av2.db", ontology="MF", pvalueCutoff= 0.025, testDirection="over")
# > hgOver = hyperGTest(params)
# 
# > hgOver
# Gene to GO MF  test for over-representation 
# 1117 GO MF ids tested (34 have p < 0.025)
# Selected gene set size: 691 
# Gene universe size: 3985 
# Annotation package: hgu95av2 
# > summary(hgOver)
# GOMFID       Pvalue OddsRatio   ExpCount Count Size                                                                  Term
# 1  GO:0005509 2.156361e-06  2.537820 23.0622334    45  133                                                   calcium ion binding
# 2  GO:0003779 7.733430e-06  2.567301 19.7676286    39  114                                                         actin binding
# GO:0060589 3.999279e-04  1.965485 25.4898369    42  147                          nucleoside-triphosphatase regulator activity
# .
# 
# 
# 
# 
# # #
# # library(limma)
# # design = model.matrix(~new_frame)
# # 
# # fit = lmFit(new_frame$,design)
# # fit2 = eBayes(fit)
# # diff = topTable(fit2, coef=2, 10)
# # diff
# # 
# 
