## plateform
#unix
grep -v '^#' GPL91-30375.txt > PLATGPL91.txt

cut -f1,11 PLATGPL91.txt | cat > gpl91.txt

# R
gpl91<-read.table("gpl91.txt",h=T,sep="\t")

library(tidyr)

new<-separate(data = gpl91, col = Gene.Symbol, into = c("left", "right"), sep = "\\ ")
ok<-new[,1:2]
colnames(ok)<-c("probe","gene")
row.names(ok)<-ok$probe
gpl_91<-ok


## matrix
#unix
grep -v '^!' GSE96-GPL91_series_matrix.txt > GSE96.txt

#R
gse96<-read.table("GSE96.txt",h=T,sep="\t")
row.names(gse96)<-gse96$ID_REF

gse96<-gse96[,2:ncol(gse96)]


## annotations
grep -v '^!Sample_source_name_ch1' GSE96-GPL91_series_matrix.txt > GSE96annot.txt
cat GSE96-GPL91_series_matrix.txt | grep '^!Sample_source_name_ch1'

cat GSE96-GPL91_series_matrix.txt | grep '^!Sample_source_name_ch1' > samples.txt

cat GSE96-GPL91_series_matrix.txt | grep '^!Sample_geo_accession' > id.txt

cat id.txt samples.txt > annot.txt
#R
annot<-read.table("annot.txt",h=F,sep="\t")
annot<-t(annot)
colnames(annot)<-c("id","group")
annot<-annot[2:nrow(annot),]
annot<-data.frame(annot)
row.names(annot)<-annot$id

##
library(tidyverse)
annot$space <- annot$group
new.annot<-annot %>% unite(concat, id, group, sep="_")

## subset matrix with some sample types
ngse96<-gse96
colnames(ngse96)<-new.annot$concat
brain<-ngse96[,grep("whole.brain",colnames(ngse96))]
thalamus<-ngse96[,grep("thalamus",colnames(ngse96))]
df<-data.frame(brain,thalamus)


##analyse
input<-merge(ok,gse96,by="row.names")
input<-input[,3:ncol(input)]






#prepare affymetrix platform
affyplatform<- function(affy)
	{
   	if(!require(tidyr)){install.packages("tidyr")}
    	library(tidyr)
	input<-data.frame(affy)
	colnames(input)<-c("ID_REF","gene")
	new<-separate(data = input, col = gene, into = c("left", "right"), sep = " /// ")
	ok<-new[,1:2]
	colnames(ok)<-c("ID_REF","gene")
	row.names(ok)<-ok$probe
	gpl_new<-ok
	gpl_new
}


### annot expression matrix
library(dplyr)
newp%>%inner_join(matrix,by="ID_REF") %>% select(-ID_REF)%>% head()










##filtermatrix usage 1st column nammed "gene" with non uniq gene symbol
#final<-filtermatrix(data)


filtermatrix <- function(data)
	{
    #install require R packages if necessary and load them
    if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
        }

    nbc=ncol(data)
    data2<-data[,2:nbc]
    m2 <- apply(data2,1, mean)
    combined2<-data.frame(m2,data$gene,data)
    ord2<-combined2[with(combined2,order(-m2)),]
    ok<-ord2 %>% distinct(gene, .keep_all = T)
    nbcok=ncol(ok)
    row.names(ok)<-ok$gene
    final<-ok[,4:nbcok]
    return(final)
    }


##usage filter<-filtermean(data,0)
filtermean<-function(data,tb=0){
	meanr<-apply(data,1,mean)
	meanok<-meanr[!is.na(meanr)]
	#baseline<-mean(meanok)
	meantable<-data.frame(meanok)
	meantable
	combined<-merge(meantable,data,by="row.names")
	row.names(combined)<-combined$Row.names
	combord<-combined[with(combined,order(-meanok)),]
	colord<-ncol(combord)
	combord<-combord[,2:colord]
	sub<-subset(combord, combord$meanok > tb)
	colsub<-ncol(sub)
	sub<-sub[,2:colsub]
	return(sub)
        }



#usage variable matrix keep only variable genes from from dataframe with uniq gene symbol as row.names
#ok<-varsig(final)


varsig <- function(data,x=1)
        {
    	variable <- apply(data,1, var)
    	var<-variable[!is.na(variable)]	
    	tb<-mean(var)+(x*sd(var))
    	baseline<-mean(var)
    	vartable<-data.frame(var)
        vartable
    	combined<-merge(vartable,data,by="row.names")
	    row.names(combined)<-combined$Row.names
	    combord<-combined[with(combined,order(-var)),]
	    colord<-ncol(combord)
	    combord<-combord[,2:colord]
	    m <- apply(combord,1, mean)
    	plot(var,m)
    	abline(v=tb,col="blue",lty=3)
	    abline(v=baseline,col="red",lty=1)
    	sub<-subset(combord, combord$var > tb)
    	head(sub,n=30)[,1:3]
	    nrow(sub)
    	return(sub)
        }


#function to remove variance column from matrix with variable genes
#usage: final<-ppmatvar(ok)


ppmatvar<-function(data)
	{
	n=ncol(data)
	data2<-data[,2:n]
	return (data2)
	}




##transpose ok (varsig) matrix for PCA

data.pca<-transmat(ok)

transmat<-function(ok)
	{
	tmat<-t(ok)
	numberrow<-nrow(tmat)
	tmat<-tmat[2:numberrow,]
	tmat
	}

## perform limma 2 groups
## limmagroup(complete,groups=pheno$group,control="other")

limmagroup<-function(data,groups,control="HD")
	{
		{
    #install require R packages if necessary and load them
    	if(!require(dplyr)){
	if (!requireNamespace("BiocManager", quietly = TRUE))
   	 install.packages("BiocManager")
	 BiocManager::install("limma")
	    library(limma)
        		}
		}
	levels <- as.factor(groups)
	levels<-relevel(levels,ref=control)
	design <- model.matrix(~levels)
	rownames(design) = colnames(complete)
	
	# perform limma analysis
	tmp <- lmFit(complete,design=design)
	fit <- eBayes(tmp)
	res = topTable(fit,number = nrow(data),coef=2)
	res
}



merge2mat<- function(m1,m2)
		{
	dfm1<-data.frame(m1)
	dfm2<-data.frame(m2)
	df<-merge(dfm1,dfm2,by="row.names")
	row.names(df)<-df$Row.names
	df<-df[,2:ncol(df)]
	df
	}
