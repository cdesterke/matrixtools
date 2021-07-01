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

##transpose ok (varsig) matrix for PCA

data.pca<-transmat(ok)

transmat<-function(ok)
	{
	tmat<-t(ok)
	numberrow<-nrow(tmat)
	tmat<-tmat[2:numberrow,]
	tmat
	}


