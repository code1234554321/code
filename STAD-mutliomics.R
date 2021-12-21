#----------for STAD multi-omics:
library(iClusterPlus)
library(GenomicRanges)
library(lattice)
library(RTCGAToolbox)
library(bnstruct)
library(ggpubr)
library(NMF)
library(RColorBrewer)
library(data.table)

library(survival)
library(KMsurv)
library(RColorBrewer)
library(glmnet)
library(survivalROC)
library(survminer)
library(pheatmap)
library(survcomp)
cox_univariant_gene_regr<-function(myd,colums){
	Surv(as.numeric(myd$A1_OS),as.numeric(myd$Status))->myd.surv
	c()->univar_anova_p;
	c()->univar_coxph_HR;
	c()->univar_coxph_low95;
	c()->univar_coxph_high95;
	c()->univar_coxph_logtest;
	c()->univar_coxph_sctest;
	c()->univar_coxph_waldtest;
	c()->fpkm.mean;
	c()->fpkm.median;
	colnames(myd)[colums]->myd.names;
	print(myd.names)
	flush.console()
	for(i in myd.names){
		mean(myd[,i],na.rm=T)->tmp.mean;
		median(myd[,i],na.rm=T)->tmp.median;
		c(fpkm.mean,tmp.mean)->fpkm.mean;
		c(fpkm.median,tmp.median)->fpkm.median;
		as.formula(paste("myd.surv~",i))->tmp.formula;
		coxph(formula=tmp.formula,data=myd)->tmp.coxph;
		summary(tmp.coxph)->tmp.coxph.summary;
		c(univar_anova_p,tmp.coxph.summary$coefficients[,5])->univar_anova_p;
		c(univar_coxph_HR,tmp.coxph.summary$coefficients[,2])->univar_coxph_HR;
		c(univar_coxph_low95,tmp.coxph.summary$conf.int[,3])->univar_coxph_low95;
		c(univar_coxph_high95,tmp.coxph.summary$conf.int[,4])->univar_coxph_high95;
		c(univar_coxph_logtest,tmp.coxph.summary$logtest[3])->univar_coxph_logtest;
		c(univar_coxph_sctest,tmp.coxph.summary$sctest[3])->univar_coxph_sctest;
		c(univar_coxph_waldtest,tmp.coxph.summary$waldtest[3])->univar_coxph_waldtest;	
	}
	data.frame("gName"=myd.names,"Pvalue"=univar_anova_p,"HR"=univar_coxph_HR,"Low(95%CI)"=univar_coxph_low95,"High(95%CI)"=univar_coxph_high95,"Logrank"=univar_coxph_logtest,"Sctest"=univar_coxph_sctest,"Waldtest"=univar_coxph_waldtest,"fpkm_median"=fpkm.median,"fpkm_mean"=fpkm.mean)->myd.coxph.df;
	#myd.coxph.df[myd.coxph.df$fpkm_median>=1,]->myd.coxph.df;
	myd.coxph.df[order(myd.coxph.df$Logrank),]->myd.coxph.df;
	myd.coxph.df[!is.na(myd.coxph.df$Logrank),]->myd.coxph.df;
	return(myd.coxph.df);
}
cox_univariant_gene_regr(myd_exp_raw_filter,tmp_index)
draw_survial_curve_custom_status<-function(myd,column,bk,status_col){
	myd[myd[,column]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,column]),]->myd.rm;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd[,status_col]))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd[,status_col]))->myd.surv;
	}
	survfit(formula=myd.surv~myd.rm[,column])->myd.fit;
	survdiff(formula=myd.surv~myd.rm[,column],rho=0)->myd.diff;
	table(myd.rm[,column])->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=brewer.pal(length(myd.table),"Set1"),xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main="Method Kaplan Meier",xlim=c(0,max_xlim));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk));
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20));
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	#print(myd.diff)
	legend("topright",legend=paste(names(myd.table),paste("(n=",myd.table,")",sep="")),fill=brewer.pal(length(myd.table),"Set1"),bty="n");
	text(max_xlim*0.85,0.8,labels=paste("p=",round(pvalue,5),sep=""));
	return(c(pvalue,myd.table));
}
draw_survial_curve_custom<-function(myd,column,bk,myd_colors){
	myd[myd[,column]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,column]),]->myd.rm;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~myd.rm[,column])->myd.fit;
	survdiff(formula=myd.surv~myd.rm[,column],rho=0)->myd.diff;
	table(myd.rm[,column])->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main="Method Kaplan Meier",xlim=c(-100,max_xlim*1.1));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	#abline(h=seq(0.2,1,0.2),col=brewer.pal(9,"Greys")[3],lty=3)
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(N=",myd.table,")",sep="")),fill=myd_colors,bty="n",cex=1.2);
	if(pvalue<1e-5){
		text(x=bk/2,y=0.25,labels="p < 1e-5",bty="n",cex=1.2,pos=4,adj=0.5,font=3)
	}else{
		text(x=bk/2,y=0.25,labels=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.2,pos=4,adj=0.5,font=3)
	}
	return(c(pvalue,myd.table));
}

setwd("d:/stomach_Adenocarcinoma")
#--------------------read CNV data
read.table("CNV/nocnv/Copy Number Variation.nocnv.merge.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv
myd.cnv[myd.cnv$Chromosome!="X",]->myd.cnv#***********
as.numeric(myd.cnv$Chromosome)->myd.cnv$Chromosome;#************
CNregions(seg=myd.cnv,epsilon=0,adaptive=FALSE,frac.overlap=0.5,rmSmallseg=TRUE,nProbes=5)->myd.cnv_merged
data.frame("Samples"=rownames(myd.cnv_merged),myd.cnv_merged)->myd.cnv_merged
write.table(myd.cnv_merged,"multi-omics/CNV_merged.txt",quote=F,sep="\t",row.names=F)
readLines("multi-omics/common_samples")->myd_samples;
read.table("multi-omics/cnv_interval2genes.txt",header=T,sep="\t")->cnv_interval2genes
########read Methylation file
fread("multi-omics/STAD_450k_0.7cutoff.rmCrossMap.factors.res",header=T,sep="\t",stringsAsFactors=F)->myd.methy
as.data.frame(myd.methy)->myd.methy
which(gsub("\\.","-",colnames(myd.methy))%in%myd_samples)->methy_index;
myd.methy[,c(1:9,methy_index)]->myd.methy
#knn impute method to handle NA values
t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
knn.impute(myd.methy_t,k=10,cat.var=1:ncol(myd.methy_t),to.impute=1:nrow(myd.methy_t),using=1:nrow(myd.methy_t))->myd.methy_t
data.frame(myd.methy[,c(1:9)],t(myd.methy_t))->myd.methy;
write.table(myd.methy,"multi-omics/STAD_450k_07.methy_filter.txt",quote=F,sep="\t",row.names=F)
########read EXP file
read.table("multi-omics/HTSeq_FPKM_STAD.merge.SYMBOL.factors.txt",header=T,sep="\t",stringsAsFactors=F)->myd.exp
myd.exp$A0_Samples->rownames(myd.exp)
myd.exp[myd_samples,]->myd.exp;
#####################process EXP:
preprocess_expd<-function(expd,gStartColumn,aliveEvent){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove IRGP value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		length(expd[expd[,i]==0,i])->failed_samples;
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}
	expd[,c(1:(gStartColumn-1),filter_colums)]->expd.filter;
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
13->g_start_column
preprocess_expd(myd.exp,g_start_column,"Alive")->myd.exp_filter;
#--------------------------------------------------------
cox_univariant_gene_regr(myd.exp_filter,c(g_start_column:c(ncol(myd.exp_filter)-1)))->myd.exp_filter.df
myd.exp_filter.df[myd.exp_filter.df$Pvalue<0.05,]->myd.exp_filter.df
myd.exp_filter.df$gName->coxph_single_g.gName
#-----------------------------------------------------------------
as.character(myd.exp_filter$A0_Samples)->myd_samples;
myd.cnv_merged[myd_samples,]->myd.cnv_merged;
which(gsub("\\.","-",colnames(myd.methy))%in%myd_samples)->methy_index;
myd.methy[,c(1:9,methy_index)]->myd.methy
#----------do correlation calculate for CNV_EXP
#------------------------------------------------------------------for regions merged to one gene
fread("multi-omics/STAD.CNV_merged_regions2symbol.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv_region2symbol;
as.data.frame(myd.cnv_region2symbol)->myd.cnv_region2symbol;
myd.cnv_region2symbol$Sample->rownames(myd.cnv_region2symbol);
myd.cnv_region2symbol[myd_samples,]->myd.cnv_region2symbol;
cal_correlation_CNV_EXP_region2symbol<-function(cnvd,expd,cnvd_start_column){
	#-------calculate correlation
	cnvd$Sample->rownames(cnvd);
	expd$A0_Samples->rownames(expd);
	cnvd[myd_samples,]->cnvd;
	expd[myd_samples,]->expd;
	colnames(cnvd)->cnvd_names;
	c()->cnvd_expd_cor;
	c()->cnvd_expd_cor_pvalue;
	c()->cnvd_genes;
	c()->cnvd_regions;
	c()->expd_genes;
	for(i in cnvd_start_column:ncol(cnvd)){
		cnvd_names[i]->g_name;
		as.character(g_name)->g_name;
		if(length(is.na(expd[1,g_name]))==0){
			next;
		}
		cor.test(cnvd[,i],expd[,g_name])->tmp_cor;
		c(cnvd_regions,cnvd_names[i])->cnvd_regions;
		c(cnvd_genes,as.character(g_name))->cnvd_genes;
		c(expd_genes,as.character(g_name))->expd_genes;
		c(cnvd_expd_cor,tmp_cor$estimate)->cnvd_expd_cor;
		c(cnvd_expd_cor_pvalue,tmp_cor$p.value)->cnvd_expd_cor_pvalue;
	}
	data.frame("CnvRegion"=cnvd_regions,"CnvGene"=cnvd_genes,"ExpGene"=expd_genes,"Cor"=cnvd_expd_cor,"Cor.pvalue"=cnvd_expd_cor_pvalue)->res_df;
	return(res_df);
}
t(myd.cnv_region2symbol[,-1])->myd.cnv_region2symbol_t;
apply(myd.cnv_region2symbol_t,1,function(x){2^x})->test;
data.frame("Samples"=myd.cnv_region2symbol$Sample,test)->myd.cnv_region2symbol_logTrans;
colnames(myd.cnv_region2symbol)->colnames(myd.cnv_region2symbol_logTrans);

cal_correlation_CNV_EXP_region2symbol(myd.cnv_region2symbol_logTrans,myd.exp_filter,2)->CNV_region2symbol_EXP.cor;
unlist(lapply(CNV_region2symbol_EXP.cor[,4],function(x){log((1+x)/(1-x))*0.5}))->CNV_region2symbol_EXP.cor$Z_value;
CNV_region2symbol_EXP.cor[order(CNV_region2symbol_EXP.cor$Cor.pvalue),]->CNV_region2symbol_EXP.cor;
p.adjust(CNV_region2symbol_EXP.cor$Cor.pvalue)->CNV_region2symbol_EXP.cor$FDR;
CNV_region2symbol_EXP.cor[CNV_region2symbol_EXP.cor$FDR<0.05,]->CNV_region2symbol_EXP.cor_filter;
CNV_region2symbol_EXP.cor_filter[CNV_region2symbol_EXP.cor_filter$Z_value>0,]->CNV_region2symbol_EXP.cor_filter;
table(CNV_region2symbol_EXP.cor_filter$CnvGene)->CnvGene_table;
names(CnvGene_table[CnvGene_table!=0])->CNVCor_region2symbol_genes;
write.table(CNV_region2symbol_EXP.cor,"multi-omics/CNV_region2symbol_EXP.cor.txt",row.names=F,sep="\t",quote=F)
#-------do correlation calculate for MET_EXP
cal_correlation_MET_EXP<-function(metd,expd,metd_start_column,expd_start_column){
	#--------map gene to probes
	c()->g_probes;
	names(table(metd$Symbol))->col1_table;
	for(s in col1_table){
		which(metd$Symbol%in%s)->tmp_probes;
		c(g_probes,paste(tmp_probes,collapse="_"))->g_probes;
	}
	col1_table->names(g_probes);
	#-------prepare matrix;
	t(metd[,metd_start_column:ncol(metd)])->metd_t;
	gsub("\\.","-",rownames(metd_t),perl=T)->rownames(metd_t);
	metd_t[myd_samples,]->metd_t;
	metd$Probe->colnames(metd_t);
	#print(metd_t[1:10,1:10])
	#flush.console();
	#-------calculate correlation
	colnames(expd)->expd_names;
	c()->metd_expd_cor;
	c()->metd_expd_cor_pvalue;
	c()->metd_probes;
	c()->metd_genes;
	c()->expd_genes;
	for(i in expd_start_column:ncol(expd)){
		expd_names[i]->g_name;
		g_probes[g_name]->g_name_probes;
		if(is.na(g_name_probes)){
			next;
		}
		as.numeric(unlist(strsplit(g_name_probes,"_")))->g_name_probes;
		for(p in g_name_probes){	
			cor.test(expd[,i],metd_t[,p])->tmp_cor;
			c(metd_probes,as.character(colnames(metd_t)[p]))->metd_probes;
			c(metd_genes,g_name)->metd_genes;
			c(expd_genes,g_name)->expd_genes;
			c(metd_expd_cor,tmp_cor$estimate)->metd_expd_cor;
			c(metd_expd_cor_pvalue,tmp_cor$p.value)->metd_expd_cor_pvalue;
		}
	}
	data.frame("MetProbe"=metd_probes,"MetGene"=metd_genes,"ExpGene"=expd_genes,"Cor"=metd_expd_cor,"Cor.pvalue"=metd_expd_cor_pvalue)->res_df;
	return(res_df);
}
cal_correlation_MET_EXP(myd.methy,myd.exp_filter,10,g_start_column)->MET_EXP.cor
unlist(lapply(MET_EXP.cor[,4],function(x){log((1+x)/(1-x))*0.5}))->MET_EXP.cor$Z_value;
MET_EXP.cor[order(MET_EXP.cor$Cor.pvalue),]->MET_EXP.cor;
p.adjust(MET_EXP.cor$Cor.pvalue)->MET_EXP.cor$FDR;
MET_EXP.cor[MET_EXP.cor$FDR<0.05,]->MET_EXP.cor_filter
MET_EXP.cor_filter[MET_EXP.cor_filter$Z_value<0,]->MET_EXP.cor_filter
table(MET_EXP.cor_filter$MetGene)->MetGene_table;
names(MetGene_table[MetGene_table!=0])->METCor_genes;
write.table(MET_EXP.cor,"multi-omics/MET_EXP.cor.txt",row.names=F,sep="\t",quote=F)
#------------------------------------
intersect(coxph_single_g.gName,CNVCor_region2symbol_genes)->CNVCor_genes_coxph;
intersect(coxph_single_g.gName,METCor_genes)->METCor_genes_coxph;

CNV_region2symbol_EXP.cor_filter[which(CNV_region2symbol_EXP.cor_filter$CnvGene%in%CNVCor_genes_coxph),]->CNV_region2symbol_EXP.cor_filter1;
MET_EXP.cor_filter[which(MET_EXP.cor_filter$MetGene%in%METCor_genes_coxph),]->MET_EXP.cor_filter1;
########--------NMF cluster:
library(NMF);
library(IntNMF)
#---using CNVCor_genes_coxph------------------------------------------------
myd.exp_filter[,CNVCor_genes_coxph]->CNVCor_genes_EXP;
myd.exp_filter$A0_Samples->rownames(CNVCor_genes_EXP);
nmf(t(CNVCor_genes_EXP),2:10,nrun=50,seed=12345)->CNVCor_genes_nmf;
plot(CNVCor_genes_nmf);
consensusmap(CNVCor_genes_nmf,labCol=NA,labRow=NA,tracks=NA)
#---------------
myd.exp_filter[,METCor_genes_coxph]->METCor_genes_EXP;
myd.exp_filter$A0_Samples->rownames(METCor_genes_EXP)
nmf(t(METCor_genes_EXP),2:10,nrun=50,seed=12345)->METCor_genes_nmf;
plot(METCor_genes_nmf);
consensusmap(METCor_genes_nmf,labCol=NA,labRow=NA,tracks=NA)
################################################################################################
#------------------------------select the best rank: CNV->3,MET->3;using CNVCor_region2symbol_genes
retrive_cluster_names<-function(myd,myd_consensusmap,hvalue){
	rownames(myd)->sample_names;
	lapply(cut(myd_consensusmap$Colv,hvalue)$lower, function(l)rapply(l,function(i)i))->myd_cut_list;
	c()->cluster_sample_names;
	c()->tmp_cluster;
	1->c_index;
	for(i in myd_cut_list){
		c(cluster_sample_names,as.character(sample_names[unlist(i)]))->cluster_sample_names;
		c(tmp_cluster,rep(paste("C",c_index,sep=""),length(unlist(i))))->tmp_cluster;
		c_index+1->c_index;		
	}
	data.frame("Sample"=cluster_sample_names,"Cluster"=tmp_cluster)->cluster_df;
	return(cluster_df);
}
nmf(t(CNVCor_genes_EXP),3,nrun=50,seed=12345)->CNVCor_genes_nmf_3;#3 clusters
consensusmap(CNVCor_genes_nmf_3,labCol=NA,labRow=NA,tracks=NA)->CNVCor_genes_nmf_3_consensusmap
retrive_cluster_names(CNVCor_genes_EXP,CNVCor_genes_nmf_3_consensusmap,0.9)->CNVCor_genes_nmf_3_cluster;
#"C1"->CNVCor_genes_nmf_3_cluster[CNVCor_genes_nmf_3_cluster$Cluster=="C2",2]
table(CNVCor_genes_nmf_3_cluster$Cluster)
c("Sample","CNVCor_C")->colnames(CNVCor_genes_nmf_3_cluster)
paste("CNVCor",CNVCor_genes_nmf_3_cluster$CNVCor_C,sep="")->CNVCor_genes_nmf_3_cluster$CNVCor_C;

nmf(t(METCor_genes_EXP),3,nrun=50,seed=12345)->METCor_genes_nmf_3;#3 clusters
consensusmap(METCor_genes_nmf_3,labCol=NA,labRow=NA,tracks=NA)->METCor_genes_nmf_3_consensusmap;
retrive_cluster_names(METCor_genes_EXP,METCor_genes_nmf_3_consensusmap,0.9)->METCor_genes_nmf_3_cluster;
table(METCor_genes_nmf_3_cluster$Cluster)
c("Sample","METCor_C")->colnames(METCor_genes_nmf_3_cluster)
paste("METCor",METCor_genes_nmf_3_cluster$METCor_C,sep="")->METCor_genes_nmf_3_cluster$METCor_C;
#---for each cluster samples: draw survival curve;
merge(myd.exp_filter,CNVCor_genes_nmf_3_cluster,by.x="A0_Samples",by.y="Sample")->myd.exp_filter_CNVCor_cluster;
merge(myd.exp_filter,METCor_genes_nmf_3_cluster,by.x="A0_Samples",by.y="Sample")->myd.exp_filter_METCor_cluster;
brewer.pal(9,"Set1")->myd_colors
par(mfrow=c(2,1),mar=c(4,5,4,4))
draw_survial_curve_custom(myd.exp_filter_CNVCor_cluster,ncol(myd.exp_filter_CNVCor_cluster),500,myd_colors)#p<0.01
draw_survial_curve_custom(myd.exp_filter_METCor_cluster,ncol(myd.exp_filter_METCor_cluster),500,myd_colors)#p<0.01
##########for skewness test: Performs D'Agostino test for skewness in normally distributed data.
library(moments)
agostino.test(CNV_region2symbol_EXP.cor$Z_value, alternative = c("two.sided", "less", "greater"))#p-value < 2.2e-16,skew=1.2352
agostino.test(sample(MET_EXP.cor$Z_value,size=40000), alternative = c("two.sided", "less", "greater"))#p-value < 3.486e-14,skew=-0.37363
#################################################################
library(VennDiagram)
library(ggpubr);
rep("CNV",nrow(CNV_region2symbol_EXP.cor))->CNV_region2symbol_EXP.cor$Type;
rep("MET",nrow(MET_EXP.cor))->MET_EXP.cor$Type;
rbind(CNV_region2symbol_EXP.cor[,3:8],MET_EXP.cor[,3:8])->CNV_MET.cor_combined;
ggdensity(CNV_MET.cor_combined,x="Z_value",fill="Type",add="mean")
grid.newpage()
draw.pairwise.venn(area1=length(CNVCor_genes_coxph),area2=length(METCor_genes_coxph),cross.area=length(intersect(CNVCor_genes_coxph,METCor_genes_coxph)),scale=F,category=c("CNVCor(368 genes)","METCor(348 genes)"),fill=brewer.pal(11,"Set3")[1:2],cat.pos=0,lwd=5,lty="blank")
############################################################################################
####-----compare CNVCor_genes cluster with METCor_genes cluster:
merge(CNVCor_genes_nmf_3_cluster,METCor_genes_nmf_3_cluster,by.x="Sample",by.y="Sample")->CNVCor_METCor_nmf_cluster
table(CNVCor_METCor_nmf_cluster$METCor_C,CNVCor_METCor_nmf_cluster$CNVCor_C)->CNVCor_METCor_nmf_cluster_table;
layout(matrix(c(1,1,2,3,3,3),nrow=2,byrow=T))
barplot(apply(CNVCor_METCor_nmf_cluster_table,2,function(x){x/sum(x)}),col=brewer.pal(9,"Set1")[1:3],beside=F)
plot.new()
legend("topleft",legend=rownames(CNVCor_METCor_nmf_cluster_table),fill=brewer.pal(9,"Set1")[1:3])
library("gplots")
balloonplot(CNVCor_METCor_nmf_cluster_table,main="CNVCor genes clustering subset \n overlap with METCor genes clustering subset",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
chisq.test(CNVCor_METCor_nmf_cluster_table)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("chisq-p:","<1e-5"))
}else{
	legend("topleft",legend=paste("chisq-p:",round(table_pvalue,5)))
}
##################################################################################################################################
##################################################################################################################################
library(gplots)
library(iClusterPlus);
myd.cnv_region2symbol_logTrans[,as.character(CNVCor_genes_coxph)]->myd.cnv_merged_CNVCor;
myd.cnv_region2symbol_logTrans$Sample->rownames(myd.cnv_merged_CNVCor);
myd.cnv_merged_CNVCor[myd_samples,]->myd.cnv_merged_CNVCor;

t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
myd.methy$Probe->colnames(myd.methy_t);
gsub("\\.","-",rownames(myd.methy_t),perl=T)->rownames(myd.methy_t);
myd.methy_t[,as.character(MET_EXP.cor_filter1$MetProbe)]->myd.methy_t_METCor;
myd.methy_t_METCor[myd_samples,]->myd.methy_t_METCor;

unique(c(CNVCor_genes_coxph,METCor_genes_coxph))->CNVCor_METCor_genes;
myd.exp_filter[,CNVCor_METCor_genes]->CNVCor_METCor_genes_EXP;
myd.exp_filter$A0_Samples->rownames(CNVCor_METCor_genes_EXP)
CNVCor_METCor_genes_EXP[myd_samples,]->CNVCor_METCor_genes_EXP;
save.image("multi-omics/work-iclusterplus2.RData")
get_best_lambda<-function(myd_tune){
	myd_tune$fit->myd_tune_fit;
	myd_tune$lambda->myd_tune_lambda;
	which.min(unlist(lapply(myd_tune_fit,function(x){x$BIC})))->min_BIC_index;
	myd_tune_lambda[min_BIC_index,]->myd_tune_best_lambda;
	return(myd_tune_best_lambda);
}
#------------------------------------single lambda iclusterplusï¼5 clusters or 6 clusters
lapply(2:4,function(x){paste("icluster_tune_fit",x,sep="")->x_fit;get_best_lambda(eval(as.symbol(x_fit)))})#2: 0.651351351 0.002702703 0.856756757
##3:0.305405405 0.008108108 0.716216216
##4:0.95945946 0.01351351 0.57567568
c()->fit_list1;
for(i in 1:20){
	set.seed(i*100+i*3+i);
	iClusterPlus(dt1=as.matrix(myd.cnv_merged_CNVCor),dt2=as.matrix(myd.methy_t_METCor),dt3=as.matrix(CNVCor_METCor_genes_EXP),type=c("gaussian","gaussian","gaussian"),lambda=c(0.651351351,0.002702703,0.856756757),K=2,maxiter=20)->single_fit;
	c(fit_list1,list(single_fit))->fit_list1;
}
c()->fit_list2;
for(i in 1:20){
	set.seed(i*100+i*3+i);
	iClusterPlus(dt1=as.matrix(myd.cnv_merged_CNVCor),dt2=as.matrix(myd.methy_t_METCor),dt3=as.matrix(CNVCor_METCor_genes_EXP),type=c("gaussian","gaussian","gaussian"),lambda=c(0.305405405,0.008108108,0.716216216),K=3,maxiter=20)->single_fit;
	c(fit_list2,list(single_fit))->fit_list2;
}
generate_icluster_groups<-function(icluster_fit){
	icluster_fit$cluster->ic_cluster;
	table(ic_cluster)->ic_table;
	paste("iC",names(ic_table),sep="")->ic_table_groups;
	names(ic_table)[order(ic_table)]->ic_table_names;
	c()->new_clusters;
	for(j in 1:length(ic_cluster)){
		for(i in 1:length(ic_table_names)){
			as.numeric(ic_table_names[i])->ic;
			if(ic_cluster[j]==ic){
				c(new_clusters,ic_table_groups[i])->new_clusters;
			}
		}
	}
	return(new_clusters);
}
draw_icluster_repeat20_survival_curve<-function(fit_list,bk){
	c()->rank_repeat20_summary;
	par(mfrow=c(4,5))
	for(i in 1:20){
		fit_list[[i]]->single_fit;
		generate_icluster_groups(single_fit)->single_fit_ic;
		c(rank_repeat20_summary,table(single_fit_ic))->rank_repeat20_summary;
		data.frame("Sample"=rownames(myd.cnv_merged_CNVCor),single_fit_ic)->iCluster_df;
		merge(myd.exp_filter,iCluster_df,by.x="A0_Samples",by.y="Sample")->iCluster_merged;
		draw_survial_curve_custom(iCluster_merged,ncol(iCluster_merged),bk,brewer.pal(9,"Set1"))
	}
	matrix(rank_repeat20_summary,ncol=length(table(single_fit_ic)),byrow=T)->rank_repeat20_summary;
	names(table(single_fit_ic))->colnames(rank_repeat20_summary);
	data.frame("Repeat"=paste("Rep",1:20,sep=""),rank_repeat20_summary)->rank_repeat20_summary;
	return(rank_repeat20_summary);
}
draw_icluster_repeat20_survival_curve(fit_list1,500)->rank2_repeat20_summary
draw_icluster_repeat20_survival_curve(fit_list2,500)->rank3_repeat20_summary;
#k=2: 3 clusters, i=1/2;
fit_list1[[4]]->single_fit;
data.frame("Sample"=rownames(myd.cnv_merged_CNVCor),"iCluster"=generate_icluster_groups(single_fit))->iCluster_df;
change_group<-function(iclusters,oldGroups,newGroups){
	c()->new_clusters;
	for(i in 1:nrow(iclusters)){
		for(ic in 1:length(oldGroups)){
			if(iclusters[i,2]==oldGroups[ic]){
				c(new_clusters,newGroups[ic])->new_clusters;
			}
		}
	}
	new_clusters->iclusters$iCluster;
	return(iclusters);
}
merge(myd.exp_filter,iCluster_df,by.x="A0_Samples",by.y="Sample")->iCluster_merged;
par(mfrow=c(1,2))
draw_survial_curve_custom(iCluster_merged,ncol(iCluster_merged),500,brewer.pal(11,"Spectral")[c(4,10,11)])
subset_cluster<-function(myd,sub_clusters){
	c()->sub_clusters_index;
	for(i in sub_clusters){
		which(myd$iCluster==i)->tmp_index;
		c(sub_clusters_index,tmp_index)->sub_clusters_index;
	}
	data.frame(myd[sub_clusters_index,])->myd_sub;
	as.character(myd_sub$iCluster)->myd_sub$iCluster;
	return(myd_sub);
}
subset_cluster(iCluster_merged,c("iC1","iC2"))->test_;#log rank p-value:0.00037
draw_survial_curve_custom(test_,ncol(iCluster_merged),500,brewer.pal(11,"Spectral")[c(4,10)])
par(mfrow=c(1,2))
subset_cluster(iCluster_merged,c("iC1","iC3"))->test_;#log rank p-value: 1.722718e-05
draw_survial_curve_custom(test_,ncol(iCluster_merged),500,brewer.pal(11,"Spectral")[c(4,11)])
subset_cluster(iCluster_merged,c("iC2","iC3"))->test_;
draw_survial_curve_custom(test_,ncol(iCluster_merged),500,brewer.pal(11,"Spectral")[c(10,11)])
#--------------------compare CNVCor/METcor vs icluster----------------------------------------------------------------
merge(CNVCor_METCor_nmf_cluster,iCluster_df,by.x="Sample",by.y="Sample")->CNVCor_METCor_iC_cluster
write.table(CNVCor_METCor_iC_cluster,"multi-omics/CNVCor_METCor_iC_cluster.txt",quote=F,sep="\t",row.names=F)

table(CNVCor_METCor_iC_cluster$CNVCor_C,CNVCor_METCor_iC_cluster$iCluster)->table1###p-value < 1e-5
table(CNVCor_METCor_iC_cluster$METCor_C,CNVCor_METCor_iC_cluster$iCluster)->table2###p-value < 1e-5
par(mfrow=c(1,2))
balloonplot(table1,main="CNVCor genes clustering subset \n overlap with iCluster genes clustering subset",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
chisq.test(table1)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("chisq-p:","<1e-5"),cex=1.2)
}else{
	legend("topleft",legend=paste("chisq-p:",round(table_pvalue,5)),cex=1.2)
}
balloonplot(table2,main="METCor genes clustering subset \n overlap with iCluster genes clustering subset",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
chisq.test(table2)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("chisq-p:","<1e-5"),cex=1.2)
}else{
	legend("topleft",legend=paste("chisq-p:",round(table_pvalue,5)),cex=1.2)
}
#--------------------------------------------------------------------------------------------------
as.numeric(substr(iCluster_df$iCluster,3,3))->merged_iclusters;
preoder_samples<-function(myd,sample_clusters,cluster_type){
	c()->myd_order;
	for(i in cluster_type){
		which(sample_clusters==i)->tmp_index;
		hclust(dist(myd[tmp_index,]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(myd_order,tmp_index)->myd_order;	
	}
	return(myd[myd_order,]);
}
CNVCor_METCor_iC_cluster[,2:4]->col_CNVCor_METCor_iC_cluster
CNVCor_METCor_iC_cluster$Sample->rownames(col_CNVCor_METCor_iC_cluster);
generate_color<-function(x,myd_colors){
	names(table(x))->x_table;
	myd_colors[1:length(x_table)]->x_colors;
	x_table->names(x_colors);
	return(x_colors);
}
generate_gaps<-function(myd_clusters){
	table(myd_clusters)->cluster_count;
	c()->myd_gaps;
	0->tmp_gap;
	for(i in 1:(length(cluster_count)-1)){
		tmp_gap+cluster_count[i]->tmp_gap;
		c(myd_gaps,tmp_gap)->myd_gaps;
	}
	return(myd_gaps);
}
generate_gaps(merged_iclusters)->gaps_col;
list("CNVCor_C"=generate_color(CNVCor_METCor_iC_cluster[,2],myd_colors),"METCor_C"=generate_color(CNVCor_METCor_iC_cluster[,3],myd_colors),"iCluster"=generate_color(CNVCor_METCor_iC_cluster[,4],brewer.pal(11,"Spectral")[c(4,10,11)]))->col_colors;

preoder_samples(myd.cnv_merged_CNVCor,merged_iclusters,c(1,2,3))->test1;
#myd.cnv_merged_CNVCor
#c(colorRampPalette(c("#00FF00","#000000"))(50),colorRampPalette(c("#000000","#FF0000"))(200))->fill_colors;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(70),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
pheatmap(t(log2(test1+1)),cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,color=fill_colors,gaps_col=gaps_col,annotation_col=col_CNVCor_METCor_iC_cluster,annotation_colors=col_colors)
preoder_samples(myd.methy_t_METCor,merged_iclusters,c(1,2,3))->test2;
#myd.methy_t_METCor
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(250),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
pheatmap(t(test2),cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,color=fill_colors,gaps_col=gaps_col,annotation_col=col_CNVCor_METCor_iC_cluster,annotation_colors=col_colors)
#######################################################################################################
#-------------------------------------for CNV genes gain or loss, MET hyper or hypo
#myd.cnv_region2symbol[,CNVCor_genes_coxph]->myd.cnv_region2symbol_CNVCor;
apply(myd.cnv_region2symbol[,-1],1,function(x){c(length(x[x>=0.3]),length(x[x<=(-0.3)]))})->myd.cnv_merged_gainLoss;
as.data.frame(t(myd.cnv_merged_gainLoss))->myd.cnv_merged_gainLoss;
c("Gain","Loss")->colnames(myd.cnv_merged_gainLoss);
myd.cnv_region2symbol$Sample->myd.cnv_merged_gainLoss$Samples;
#for Methylation:
t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
myd.methy$Probe->colnames(myd.methy_t);
gsub("\\.","-",rownames(myd.methy_t),perl=T)->rownames(myd.methy_t);
apply(myd.methy_t,1,function(x){c(length(x[x>0.8]),length(x[x<0.2]))})->myd.methy_t_hyperHypo;
as.data.frame(t(myd.methy_t_hyperHypo))->myd.methy_t_hyperHypo;
c("MetHyper","MetHypo")->colnames(myd.methy_t_hyperHypo)
rownames(myd.methy_t)->myd.methy_t_hyperHypo$Samples;
merge(myd.cnv_merged_gainLoss,myd.methy_t_hyperHypo,by.x="Samples",by.y="Samples")->myd.CNV_MET_abnormal_frequency;
write.table(myd.CNV_MET_abnormal_frequency,"multi-omics/CNV_MET_abnormal_frequency.txt",quote=F,sep="\t",row.names=F)
#lm(Gain~Loss,data=myd.CNV_MET_abnormal_frequency)->g_l_lm;
#_-------draw points :
draw_genes_exp_point_v1<-function(expd,g1,g2){
	plot(expd[,g1],expd[,g2],xlab=paste(c("",g1),collapse=""),ylab=paste(c("",g2),collapse=""),pch=20,cex=0.5,bty="o");
	cor.test(expd[,g1],expd[,g2],method = c("spearman"))->g1_g2_cortest;
	g1_g2_cortest$estimate->cortest_cor;
	g1_g2_cortest$p.value->cortest_p;
	#-------------
	expd[,c(g1,g2)]->expd_g1_g2;
	#log2(expd_g1_g2+1)->expd_g1_g2;
	paste(g2,g1,sep="~")->g1_g2_formula;
	as.formula(g1_g2_formula)->g1_g2_formula;
	lm(g1_g2_formula,data=expd_g1_g2)->g1_g2_lmfit;
	abline(g1_g2_lmfit,lty=3,col="red",lwd=2)
	if(cortest_p < 0.001){
		"padj < 0.001"->cortest_p
	}else{
		paste("padj = ",round(cortest_p,3),sep="")->cortest_p;
	}
	paste("R2 = ",round(cortest_cor,3),sep="")->cortest_cor;
	paste(cortest_cor,cortest_p,sep=",")->leg;
	mtext(leg,cex=1.2,font=3);
	return(c(g1_g2_cortest$estimate,g1_g2_cortest$p.value));
}
par(mfrow=c(2,3))
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"Gain","Loss")->gainloss_cor;
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"Gain","MetHyper")->gainhyper_cor;
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"Gain","MetHypo")->gainhypo_cor;
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"MetHyper","Loss")->hyperloss_cor;
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"MetHypo","Loss")->hypoloss_cor;
draw_genes_exp_point_v1(myd.CNV_MET_abnormal_frequency,"MetHyper","MetHypo")->hyperhypo_cor;
rbind(gainloss_cor,gainhyper_cor,gainhypo_cor,hyperloss_cor,hypoloss_cor,hyperhypo_cor)->cnv_met_cor;
as.data.frame(cnv_met_cor)->cnv_met_cor;
c("R2","CorrP")->colnames(cnv_met_cor);
as.numeric(cnv_met_cor$CorrP)->cnv_met_cor$CorrP;
cnv_met_cor[order(cnv_met_cor$CorrP),]->cnv_met_cor;
p.adjust(cnv_met_cor$CorrP)->cnv_met_cor$FDR;
#######################################################################################################################
#-----------------------------immune landscape for iC3 and iC7: immune infiltrating,Leukocyte fraction,SNV-Neoantigen,Indel-Neoantigen
read.table("d:/script/data-set/TIMER-data/TCGA-immuneEstimation.txt",header=T,sep="\t",stringsAsFactors=F)->immuneEstimation;
merge(immuneEstimation,CNVCor_METCor_iC_cluster,by.x="barcode",by.y="Sample")->iC_cluster_immunescore;
write.table(iC_cluster_immunescore,"multi-omics/iC_cluster_immunescore.txt",quote=F,sep="\t",row.names=F)

calculate_immunescore_pvalue<-function(myd,columns,clusters,cluster_col){
	list()->c_groups;
	for(i in clusters){
		which(myd[,cluster_col]==i)->tmp_index;
		c(c_groups,list(tmp_name=tmp_index))->c_groups;
	}
	clusters->names(c_groups);
	c()->compare_clusters;
	c()->rank_test_pvalues;
	for(i in 1:(length(clusters)-1)){
		for(j in (i+1):length(clusters)){
			paste(clusters[i],clusters[j],sep="~")->tmp_name;
			c(compare_clusters,tmp_name)->compare_clusters;
			for(k in columns){
				as.numeric(unlist(c_groups[clusters[i]]))->c1_index;
				as.numeric(unlist(c_groups[clusters[j]]))->c2_index;
				if(length(clusters)>2){
					kruskal.test(list(myd[c1_index,k],myd[c2_index,k]))$p.value->tmp_pvalue;
				}else{
					wilcox.test(myd[c1_index,k],myd[c2_index,k])$p.value->tmp_pvalue;
				}
				c(rank_test_pvalues,tmp_pvalue)->rank_test_pvalues;
			}
		}
	}

	length(c_groups)->c_groups_len;
	matrix(rank_test_pvalues,c_groups_len*(c_groups_len-1)/2,length(columns),byrow=T)->res_df;	
	#print(c_groups_len*(c_groups_len-1)/2);
	#flush.console();
	compare_clusters->rownames(res_df);
	colnames(myd)[columns]->colnames(res_df);
	return(res_df);
}
calculate_immunescore_pvalue(iC_cluster_immunescore,2:7,c("iC1","iC2","iC3"),10)->iC_cluster_immunescore_ranktest;
colnames(iC_cluster_immunescore)[2:7]->immune_cells;
draw_boxplot_genes_by_factors<-function(expd,genes,f,myd_colors){
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	names(table(expd[,f]))->f_names;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	lapply(genes,function(x){expd[,x]->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",boxwex=0.5,pch=20,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),boxplot_max)*1.5)#->genes_f_values.boxplot;
	(boxplot_max-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	if(is.wholenumber(boxplot_max) && abs(boxplot_max) < 1){
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=seq(round(boxplot_min),boxplot_max*1.1,bk))
	}else{
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=round(seq(round(boxplot_min),boxplot_max*1.1,bk),3))
	}
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	print(bk);
	flush.console();
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.4,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	box();
	return(kruskal_test_res);
}
draw_boxplot_genes_by_factors_v2<-function(expd,genes,f,myd_colors){
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	names(table(expd[,f]))->f_names;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	lapply(genes,function(x){log2(expd[,x]+1)->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",boxwex=0.5,pch=20,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),round(boxplot_max)*1.5))#->genes_f_values.boxplot;
	(round(boxplot_max)-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	axis(side=2,at=seq(round(boxplot_min),round(boxplot_max)*1.1,bk),labels=seq(round(boxplot_min),round(boxplot_max)*1.1,bk))
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.2,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	#abline(v=seg_index,lty=3)
	box();
	return(kruskal_test_res);
}
#------------
iC_cluster_immunescore->expd;
par(mar=c(7,3,4,3))
"iCluster"->f;#change clinical features to draw boxplot 
draw_boxplot_genes_by_factors_v2(expd,immune_cells,f,brewer.pal(11,"Spectral")[c(4,10,11)])->test_;
#------------
merge(myd.exp_filter,iC_cluster_immunescore,by.x="A0_Samples",by.y="barcode")->myd.exp_filter_immunescore;
preorder_samples_iCluster_OS<-function(myd_merged,clusters){
	#myd_merged[,c(as.character(geneOrder),"RiskScore","RiskType")]->myd_order;
	#data.frame("A0_Samples"=myd_merged$A0_Samples,myd_order)->myd_order;
	#merge(myd_order,sample_clusters,by.x="A0_Samples",by.y="Sample")->myd_merged;
	c()->order_index;
	for(i in clusters){
		which(myd_merged$iCluster==i)->tmp_index;
		order(myd_merged[tmp_index,"A1_OS"])->tmp_index_hclust;
		tmp_index[tmp_index_hclust]->tmp_index;
		c(order_index,tmp_index)->order_index;
	}
	myd_merged[order_index,]->myd_merged;
	return(myd_merged)
}
preorder_samples_iCluster_OS(myd.exp_filter_immunescore,c("iC1","iC2","iC3"))->iC_cluster_immunescore_sort;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(70),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
generate_gaps(as.numeric(substr(iC_cluster_immunescore$iCluster,3,3)))->gaps_col;
list("CNVCor_C"=generate_color(CNVCor_METCor_iC_cluster[,2],myd_colors),"METCor_C"=generate_color(CNVCor_METCor_iC_cluster[,3],myd_colors),"iCluster"=generate_color(CNVCor_METCor_iC_cluster[,4],brewer.pal(11,"Spectral")[c(4,10,11)]),"A2_Event"=generate_color(iC_cluster_immunescore_sort$A2_Event,brewer.pal(9,"Greys")[c(3,9)]),"Age"=colorRampPalette(brewer.pal(9,"Purples"))(200),"A1_OS"=colorRampPalette(brewer.pal(9,"Greens"))(200))->col_colors;
pheatmap(t(iC_cluster_immunescore_sort[,immune_cells]),cluster_rows=F,color=fill_colors,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellheight=40,border_color=NA,annotation_col=iC_cluster_immunescore_sort[,c("CNVCor_C","METCor_C","iCluster","A1_OS","A2_Event")],annotation_colors=col_colors)
#---------------------
read.table("d:/script/data-set/immune-landscape.txt",header=T,sep="\t",stringsAsFactors=F)->immune_landscape;
retrive_TCGA_immuneValues<-function(samples,immuneSig){
	unlist(strsplit(as.character(samples[1]),split="-"))->s_split;
	if(length(s_split)>3){
		unlist(lapply(samples,function(x){as.character(x)->x;substr(x,1,nchar(x)-3)->s_res;s_res}))->tmp_samples;
	}
	tmp_samples->names(samples);
	unlist(lapply(tmp_samples,function(x){which(immune_landscape[,1]==x)->x_index;x_index}))->s_index;
	immune_landscape[s_index,immuneSig]->res;
	immune_landscape[s_index,1]->immune_landscape_s;
	data.frame(samples[immune_landscape_s],res)->res;
	c("Sample",immuneSig)->colnames(res);
	return(res);
}
#retrive_TCGA_immuneValues(CNVCor_METCor_iC_cluster$Sample,c("Leukocyte.Fraction","BCR.Shannon","TCR.Shannon"))->sample_immune_landscape;
retrive_TCGA_immuneValues(CNVCor_METCor_iC_cluster$Sample,c("SNV.Neoantigens","Silent.Mutation.Rate","Nonsilent.Mutation.Rate"))->sample_immune_landscape;
merge(sample_immune_landscape,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="Sample")->iC_cluster_immune_landscape;
write.table(iC_cluster_immune_landscape,"multi-omics/iC_cluster_immune_landscape.txt",quote=F,sep="\t",row.names=F)
calculate_immunescore_pvalue(iC_cluster_immune_landscape,2:4,c("iC1","iC2","iC3"),7)->iC_cluster_immune_landscape_ranktest;
par(mfrow=c(1,3));
colnames(iC_cluster_immune_landscape)->iC_names;
for(i in 2:4){
	as.formula(paste(iC_names[i],"iCluster",sep="~"))->boxplot_formula;
	max(iC_cluster_immune_landscape[,iC_names[i]],na.rm=T)->boxplot_res_max;
	min(iC_cluster_immune_landscape[,iC_names[i]],na.rm=T)->boxplot_res_min;
	boxplot(formula=boxplot_formula,data=iC_cluster_immune_landscape,col=brewer.pal(9,"Set1")[1:4],pch=20,ylim=c(boxplot_res_min*1.1,boxplot_res_max*1.5))->boxplot_res;
	paste(rownames(iC_cluster_immune_landscape_ranktest),round(iC_cluster_immune_landscape_ranktest[,iC_names[i]],5),sep=":")->legend_ranktest_p;
	legend("topleft",legend=c(paste(iC_names[i]," score",sep=""),legend_ranktest_p))
}
###################################################################################
#---compare clinical features in iC3 and iC4:-------------------------------------
library(dgof)
factors_barplot_clustering<-function(myd,sample_clusters,column){
	merge(myd,sample_clusters,by.x="A0_Samples",by.y="Sample")->myd.merged;
	myd.merged[!is.na(myd.merged[column]),]->myd.merged;
	table(myd.merged$iCluster,myd.merged[,column])->riskType_vs_cluster.table;
	dim(riskType_vs_cluster.table)->table_size;
	print(riskType_vs_cluster.table)
	#ks.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	chisq.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	max(riskType_vs_cluster.table)*0.9->max_ylim;
	barplot(riskType_vs_cluster.table,beside=T,border=F,col=brewer.pal(9,"Set1")[1:table_size[1]],ylab="Count")#for high and low risk boxplot 
	legend("topright",legend=rownames(riskType_vs_cluster.table),fill=brewer.pal(9,"Set1")[1:table_size[2]],horiz=F)
	text(x=1,y=max_ylim,labels=paste("chisq p:",round(riskType_vs_cluster.pvalue,5),sep=""),pos=4)
	mtext(colnames(myd)[column])
}

factors_barplot_ratio_clustering<-function(myd,sample_clusters,column){
	merge(myd,sample_clusters,by.x="A0_Samples",by.y="Sample")->myd.merged;
	myd.merged[!is.na(myd.merged[column]),]->myd.merged;
	table(myd.merged$iCluster,myd.merged[,column])->riskType_vs_cluster.table;
	dim(riskType_vs_cluster.table)->table_size;
	c()->table_ratio;
	for(i in 1:nrow(riskType_vs_cluster.table)){
		for(j in 1:ncol(riskType_vs_cluster.table)){
			riskType_vs_cluster.table[i,j]/sum(riskType_vs_cluster.table[i,])->res;
			c(table_ratio,res)->table_ratio;
		}
	}
	matrix(table_ratio,nrow=table_size[2],byrow=F)->table_ratio;
	rownames(riskType_vs_cluster.table)->colnames(table_ratio)
	colnames(riskType_vs_cluster.table)->rownames(table_ratio)
	print(table_ratio)
	#ks.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	chisq.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	barplot(table_ratio,beside=F,border=T,col=brewer.pal(9,"Set1")[1:table_size[2]],ylab="Ratio")#for high and low risk boxplot 
	legend(x="bottom",legend=colnames(riskType_vs_cluster.table),fill=brewer.pal(9,"Set1")[1:table_size[2]],horiz=T,bty="n",inset=c(0,-0.15),xpd=T)
	paste("Chisq-p:",round(riskType_vs_cluster.pvalue,5),sep="")->chisq_p;
	mtext(paste(colnames(myd)[column],chisq_p,sep="\n"))
}
#---change_values for clinical features;
#----------------change values-----------------------------------------------------
myd.exp_filter_immunescore->tmp.filter;
change_values<-function(myd,column,values){
	myd[!is.na(myd[,column]),]->myd
	colnames(myd)[column]->column_name;
	table(myd[,column])->N_stage.table;
	as.factor(names(N_stage.table))->N_stage.names;
	data.frame(column_name=N_stage.names,"Value"=values)->N_stage.df;
	c()->tmp.value;
	for(i in 1:nrow(myd)){
		for(j in 1:nrow(N_stage.df)){
			if(myd[i,column]==N_stage.df[j,1]){
				c(tmp.value,as.character(N_stage.df[j,2]))->tmp.value;
			}
		}
	}
	tmp.value->myd[,column];
	return(myd);
}
change_values(tmp.filter,4,c("T1","T1","T1","T2","T2","T2","T3","T4","T4","T4"))->tmp.filter;
change_values(tmp.filter,5,c("NX","N0","N1","N2","N3","N3","N3","NX"))->tmp.filter;
change_values(tmp.filter,6,c("M0","M1","MX"))->tmp.filter;
change_values(tmp.filter,7,c("X","I","I","I","II","II","II","III","III","III","III","IV"))->tmp.filter;
change_values(tmp.filter,8,c("G1","G2","G3","GX"))->tmp.filter;
#------------------------------------------------------------------
split_factor<-function(myd,column,values){
	c()->cut_values;
	c()->range_name;
	for(k in 1:(length(values)-1)){
		c(values[k],values[k+1])->row.value;
		rbind(cut_values,row.value)->cut_values;
		c(range_name,paste(values[k],values[k+1],sep="~"))->range_name;
	}
	c("Start","End")->colnames(cut_values);
	as.data.frame(cut_values)->cut_values;
	range_name->cut_values$Name;
	c()->test.values;
	for(j in 1:nrow(myd)){
		for(i in 1:nrow(cut_values)){
			if(myd[j,column]>=cut_values[i,1] && myd[j,column]<cut_values[i,2]){
				c(test.values,cut_values[i,3])->test.values;
			}
		}
	}
	test.values->myd[,column]
	return(myd);
}
0->tmp.filter[which(is.na(tmp.filter[,9])),9]
as.numeric(tmp.filter[,9])->tmp.filter[,9]
split_factor(tmp.filter,9,c(0,50,60,70,80,100))->tmp.filter;
#--------------------------------------------------
par(mfrow=c(2,3))
for(i in c(4:9)){
	factors_barplot_ratio_clustering(tmp.filter,CNVCor_METCor_iC_cluster,i)
}
#-----------------------------------------------------------------------------------------------------------
Sample_info_summary<-function(myd,columns){
	list()->summary.res;
	for(i in columns){
		as.data.frame(table(myd[i]))->tmp.df;
		c(colnames(myd)[i],"Values")->colnames(tmp.df);
		print(tmp.df);
	}
}
Sample_info_summary(tmp.filter,c(3:10));
Sample_info_summary(subset_cluster(tmp.filter,c("iC3")),c(3:10));
#------------------------------------------------------for CNVCor genes:
CNVCor_genes_count<-function(myd,intervals){
	names(table(as.character(myd$CnvGene)))->genes;
	lapply(genes,function(x){which(intervals$Gene%in%x)[1]})->cnv_intervals_index;
	intervals[unlist(cnv_intervals_index),]->res;
	lapply(res$CNV_region,function(x){unlist(strsplit(as.character(x),split="\\."))->res;res[1]})->cnv_genes_chr;
	unlist(cnv_genes_chr)->res$Chr;
	merge(myd,res,by.x="CnvGene",by.y="Gene")->res;
	return(res);
}
CNVCor_genes_count(CNV_region2symbol_EXP.cor_filter1,cnv_interval2genes)->CNVCor_genes_chr_count;
as.numeric(gsub("chr","",CNVCor_genes_chr_count$Chr))->CNVCor_genes_chr_count$ChrN
table(CNVCor_genes_chr_count$ChrN)->chr_table;
chr_table/sum(chr_table)->chr_table_ratio;
paste("Chr",names(chr_table),sep="")->names(chr_table_ratio);
layout(matrix(c(1,2,2),3,1,byrow=T))
barplot(chr_table_ratio,las=2,col=brewer.pal(12,"Set3")[12],cex.names=1.5);
boxplot(Cor~ChrN,data=CNVCor_genes_chr_count,col=brewer.pal(12,"Set3")[12],pch=20,cex.names=2.5,names=paste("Chr",names(chr_table),sep=""))
#ggviolin(CNVCor_genes_chr_count,x="Chr",y="Cor",add=c("mean_se"),fill=brewer.pal(11,"Set3")[2],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(0.7,0.8),linetype=2)

#--------------------------------------------------------------------------------------------------------------
METCor_genes_count<-function(myd,myd_cor){
	merge(myd_cor,myd[,1:9],by.x="MetProbe",by.y="Probe")->myd_merged;
	table(myd_merged$MetGene)->met_genes;
	met_genes[met_genes!=0]->met_genes;
	names(met_genes)->met_genes;
	c()->filter_index;
	for(g in met_genes){
		for(i in 1:nrow(myd_merged)){
			if(g==myd_merged[i,2]){
				c(filter_index,i)->filter_index;
				break;
			}
		}
	}
	return(myd_merged[filter_index,]);
}
METCor_genes_count(myd.methy,MET_EXP.cor_filter1)->METCor_genes_summary_count;
ggboxplot(METCor_genes_summary_count,x="Chr",y="Cor",fill=brewer.pal(12,"Set3")[12],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(-0.5,-0.6),linetype=2)
ggboxplot(METCor_genes_summary_count,x="Chr",y="TSSDist",fill=brewer.pal(12,"Set3")[12],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(-200,500,1000,1500),linetype=2)
#par(mfrow=c(2,1))
#boxplot(Cor~Chr,data=METCor_genes_summary_count,col=brewer.pal(11,"Set3")[1],pch=20)
#boxplot(TSSDist~Chr,data=METCor_genes_summary_count,col=brewer.pal(11,"Set3")[2],pch=20)

draw_METCor_genes_summary<-function(myd){
	layout(matrix(c(1,1,1,1,2,2,2,3),2,4,byrow=T));
	table(myd$Chr)->chr_table;
	c(paste("chr",seq(1,22,1),sep=""),"chrX")->chr_order;
	chr_table/sum(chr_table)->chr_table_ratio;
	names(chr_table)->names(chr_table_ratio);
	chr_table_ratio[chr_order]->chr_table_ratio;
	barplot(chr_table_ratio,col=brewer.pal(12,"Set3")[12],cex.names=1.5,las=2);
	table(myd$GeneType)->gtype_table;
	gtype_table[!is.na(names(gtype_table))]->gtype_table;
	gtype_table/sum(gtype_table)->gtype_table_ratio;
	names(gtype_table)->names(gtype_table_ratio);
	gtype_table_ratio[order(gtype_table_ratio,decreasing=T)[1:5]]->gtype_table_ratio;
	barplot(gtype_table_ratio,horiz=T,names.arg=F,axes=F,col=brewer.pal(12,"Set3")[12])->gtype_mp;
	axis(1);
	text(x=0.01,y=gtype_mp,labels=names(gtype_table_ratio),cex=1.5,pos=4);
	#axis(2,at=gtype_mp,mgp=c(0,-15,0),cex=2,tick=F,labels=names(gtype_table_ratio),las=2)
	table(myd$CGI)->cgi_table;
	which(names(cgi_table)%in%"")->null_index;
	c(null_index,which(names(cgi_table)%in%"NA"))->null_index;
	cgi_table[-null_index]->cgi_table;
	cgi_table/sum(cgi_table)->cgi_table_ratio;
	names(cgi_table)->names(cgi_table_ratio);
	barplot(cgi_table_ratio,col=brewer.pal(12,"Set3")[12],names.arg=F,axes=F,horiz=T)->cgitype_mp;
	axis(1);
	text(y=cgitype_mp,x=0.01,labels=names(cgi_table_ratio),cex=1.5,pos=4);

}
draw_METCor_genes_summary(METCor_genes_summary_count)
####################################################################################
############################################
#----DESeq use TCGA readcount:
library(DESeq2)
read.table("multi-omics/HTSeq-Counts_STAD.merge.SYMBOL.txt",header=T,sep="\t",stringsAsFactors=F)->myd_ht_count;
as.matrix(myd_ht_count[,-c(1,2)])->myd_readcount_matrix;
as.character(myd_ht_count[,2])->rownames(myd_readcount_matrix)
data.frame("SampleID"=gsub("-","\\.",CNVCor_METCor_iC_cluster$Sample),"Condition"=CNVCor_METCor_iC_cluster$iCluster)->iC_group_condition;
generate_group<-function(myd_clusters,target_clusters){
	c()->target_index;
	for(i in target_clusters){
		
		which(myd_clusters$Condition==i)->tmp_index;
		c(target_index,tmp_index)->target_index;
	}
	myd_clusters[target_index,]->myd_clusters_subset;
	data.frame("SampleID"=myd_clusters_subset[,1],"Condition"=as.character(myd_clusters_subset[,2]))->myd_clusters_subset;
	return(myd_clusters_subset);
}
#----------------------for iC1~iC2------------------
generate_group(iC_group_condition,c("iC1","iC2"))->iC_filter_group;
myd_readcount_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC2/iC1)
library(parallel)
detectCores()->no_cores;
makeCluster(no_cores-1)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-rowMeans(counts(dds))>=5
	dds[keep,]->dds.keep
	DESeq(dds.keep)->dds.keep.deseq
	results(dds.keep.deseq)->dds_res
	dds_res[order(dds_res$pvalue),]->dds_res_order;
	as.data.frame(dds_res_order)->dds_res_order_df;
	data.frame("TransID"=rownames(dds_res_order_df),dds_res_order_df)->dds_res_order_df;##log2((iC2 exp)/(iC1 exp))
	return(dds_res_order_df);
}
sample_distance_estimate<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-rowMeans(counts(dds))>=5
	dds[keep,]->dds.keep
	assay(vst(dds.keep),blind = FALSE)->dds_keep_vst;
	dist(t(dds_keep_vst))->sampleDists
	as.matrix(sampleDists)->sampleDistMatrix
	colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)->colors
	data.frame("Condition"=colData$Condition)->annot_col;
	colData$SampleID->rownames(annot_col)
	brewer.pal(9,"Set1")[c(1,2)]->annot_colors;
	names(table(annot_col))->names(annot_colors)
	list("Condition"=annot_colors)->annot_colors;
	pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,show_colnames=F,show_rownames=F,annotation_col=annot_col,annotation_row=annot_col,annotation_names_row=F,annotation_names_col=F,annotation_colors=annot_colors)->p;
	return(p);
}
do_DESeq(myd_readcount_matrix,iC_filter_group)->C2_C1_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
merge(C2_C1_deseq.res,myd_ht_count[,1:2],by.x="TransID",by.y="Tags")->C2_C1_deseq.res
C2_C1_deseq.res[order(C2_C1_deseq.res$padj),]->C2_C1_deseq.res;
write.table(C2_C1_deseq.res,"multi-omics/iC2_iC1_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
#----------------------for iC1~iC3------------------
generate_group(iC_group_condition,c("iC1","iC3"))->iC_filter_group;
as.matrix(myd_ht_count[,-c(1,2)])->myd_readcount_matrix;
as.character(myd_ht_count[,2])->rownames(myd_readcount_matrix)
myd_readcount_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC3/iC1)
library(parallel)
detectCores()->no_cores;
makeCluster(no_cores-1)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq(myd_readcount_matrix,iC_filter_group)->C3_C1_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
write.table(C3_C1_deseq.res,"multi-omics/iC3_iC1_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
#----------------------for iC2~iC3------------------
generate_group(iC_group_condition,c("iC2","iC3"))->iC_filter_group;
as.matrix(myd_ht_count[,-c(1,2)])->myd_readcount_matrix;
as.character(myd_ht_count[,2])->rownames(myd_readcount_matrix)
myd_readcount_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC3/iC2)
library(parallel)
detectCores()->no_cores;
makeCluster(no_cores-1)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq(myd_readcount_matrix,iC_filter_group)->C3_C2_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
write.table(C3_C2_deseq.res,"multi-omics/iC3_iC2_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
#########################################################
#--------------------draw volcano figures for C3_C2 and C4_C3:
draw_volcano_figure<-function(deseq,cut_p,cut_lfc,p_title){
	brewer.pal(9,"Set1")->myd_colors;
	unlist(lapply(1:nrow(deseq),function(x){if(deseq[x,3]<(-cut_lfc) && deseq[x,7]<cut_p){myd_colors[2]->res;}else if(deseq[x,3]>cut_lfc && deseq[x,7]<cut_p){myd_colors[1]->res;}else{myd_colors[9]->res};res}))->deseq_color;
	plot(x=deseq$log2FoldChange,y=-log(deseq$padj),pch=20,cex=0.5,col=deseq_color,xlim=c(-3,3),xlab="Log2 fold change",ylab="-Log10(FDR)")
	legend("top",legend=c("Down","Up"),pch=20,cex=1.5,col=myd_colors[c(2,1)],title=p_title,border=NA)
}
par(mfrow=c(2,2))
draw_volcano_figure(C2_C1_deseq.res,0.05,1,"C2/C1")
draw_volcano_figure(C3_C1_deseq.res,0.05,1,"C3/C1")
draw_volcano_figure(C3_C2_deseq.res,0.05,1,"C3/C2")
#-------------------------compare C3_C2 and C4_C3 DGEs: keep shared DGEs
library(limma)
generate_DEGs<-function(deseq_list,x){
	c()->tmp_res;
	for(y in deseq_list){
		"F"->res;
		data.frame("gName"=rownames(y),y)->y;
		for(j in 1:nrow(y)){	
			if(y[j,1]== as.character(x) && abs(y[j,3])>cut_lfc && y[j,7]<cut_p){
				"T"->res;
			}
		}
		c(tmp_res,res)->tmp_res;
	}
	return(tmp_res);
}
prepare_vennCount<-function(deseq_list,cut_p,cut_lfc){
	unlist(lapply(deseq_list,function(x){rownames(x)}))->total_transIDs;
	unique(total_transIDs)->total_transIDs;
	detectCores()->no_cores;
	makeCluster(no_cores-1)->c1;	
	clusterExport(c1,c("generate_DEGs","deseq_list","cut_p","cut_lfc"));#
	parLapply(c1,total_transIDs,function(x)generate_DEGs(deseq_list,x))->res_list;
	stopCluster(c1);
	unlist(res_list)->res_list;
	matrix(res_list,ncol=length(deseq_list),byrow=T)->res_list_matrix
	names(deseq_list)->colnames(res_list_matrix);
	total_transIDs->rownames(res_list_matrix);
	return(res_list_matrix);
}
list("C2_C1_deseq"=C2_C1_deseq.res,"C3_C1_deseq"=C3_C1_deseq.res,"C3_C2_deseq"=C3_C2_deseq.res)->deseq_list;
0.05->cut_p;
1->cut_lfc;#log2(1.5) not good; log2(2) also not good;
save.image("multi-omics/work-deseq-multiomics-STAD.RData")#use linux service to do this work;
#prepare_vennCount(deseq_list,cut_p,cut_lfc)->deseq_list.table;
read.table("TME/deseq_list.table.txt",header=T,sep="\t",stringsAsFactors=F,row.names=1)->deseq_list.table;
as.data.frame(deseq_list.table)->deseq_list.table
for(i in 1:ncol(deseq_list.table)){
	as.logical(deseq_list.table[,i])->deseq_list.table[,i];
}
vennCounts(deseq_list.table)->deseq_list.vennCount;
vennDiagram(deseq_list.vennCount,names=c("C2/C1","C3/C1","C3/C2"),cex=1.2,lwd=1.5,circle.col=brewer.pal(9,"Set1"));
#--------------
library(gplots)
C2_C1_deseq.res[abs(C2_C1_deseq.res$log2FoldChange)>cut_lfc,]->test1_
test1_[test1_$padj<cut_p,"TransID"]->test1_transID
C3_C1_deseq.res[abs(C3_C1_deseq.res$log2FoldChange)>cut_lfc,]->test2_;
test2_[test2_$padj<cut_p,"TransID"]->test2_transID;
C3_C2_deseq.res[abs(C3_C2_deseq.res$log2FoldChange)>cut_lfc,]->test3_;
test3_[test3_$padj<cut_p,"TransID"]->test3_transID;
list("C3/C2"=test1_transID,"C3/C1"=test2_transID,"C3/C2"=test3_transID)->deseq_TransID_list;#fold change > 1.5
venn(deseq_TransID_list)->deseq_TransID_list.table;
attr(deseq_TransID_list.table,"intersections")->deseq_TransID_list.table_instersections;
unlist(lapply(deseq_TransID_list.table_instersections,function(x){length(x)}))->venn_intersections;
deseq_TransID_list.table_instersections[["C3/C2:C3/C1:C3/C2"]]->CNVCor_METCor_df_genes;
myd_ht_count[which(myd_ht_count$Tags%in%CNVCor_METCor_df_genes),1]->CNVCor_METCor_df_genes
#----------------compare iC1 and iC2 CNV: CNVCor_METCor_iC_cluster+CNVCor_METCor_df_genes+myd.cnv_merged
which(colnames(myd.cnv_region2symbol)%in%CNVCor_METCor_df_genes)->cnv_intervals_index;
myd.cnv_region2symbol[,cnv_intervals_index]->myd.cnv_merged_changed;
for(i in 1:length(cnv_intervals_index)){
	for(j in 1:nrow(myd.cnv_merged_changed)){
			if(myd.cnv_merged_changed[j,i]>0.3){
				1->myd.cnv_merged_changed[j,i];
			}else if(myd.cnv_merged_changed[j,i]<(-0.3)){
				(-1)->myd.cnv_merged_changed[j,i];
			}else{
				0->myd.cnv_merged_changed[j,i];
			}
	}
}
data.frame("Sample"=myd.cnv_region2symbol$Sample,myd.cnv_merged_changed)->cnv_genes_changed;
merge(cnv_genes_changed,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="Sample")->cnv_genes_changed_merged;

compare_icluster_cnv_count<-function(myd,colums,clusters,clusters_col){
	data.frame("iC1_amp"=as.numeric(),"iC2_amp"=as.numeric(),"iC1_del"=as.numeric(),"iC2_del"=as.numeric(),"iC1_Normal"=as.numeric(),"iC2_Normal"=as.numeric(),"FisherP"=as.numeric())->res_df;
	for(i in colums){
		c()->cnv_del_count;
		c()->cnv_amp_count;
		c()->cnv_n_count;
		table(myd[,clusters_col])->total_count;
		total_count[clusters]->total_count;
		for(j in clusters){
			which(myd[,clusters_col]==j)->ic_index;
			length(which(myd[ic_index,i]==1))->amp_count;
			length(which(myd[ic_index,i]==-1))->del_count;
			length(which(myd[ic_index,i]==0))->n_count;
			c(cnv_del_count,del_count)->cnv_del_count;
			c(cnv_amp_count,amp_count)->cnv_amp_count;
			c(cnv_n_count,n_count)->cnv_n_count;
		}
		matrix(c(cnv_amp_count,cnv_del_count,cnv_n_count),2,3,byrow=F)->n_y_count;
		fisher.test(n_y_count)$p.value->n_y_count.pvalue;
#		print(n_y_count);
#		print(paste("fisher-test: ",n_y_count.pvalue,sep=""));
#		flush.console();
		rbind(res_df,c(cnv_amp_count,cnv_del_count,cnv_n_count,n_y_count.pvalue))->res_df;
	}
	c("iC1_amp","iC2_amp","iC1_del","iC2_del","iC1_Normal","iC2_Normal","FisherP")->colnames(res_df);
	data.frame("ID"=colnames(myd)[colums],res_df)->res_df;
	res_df[order(res_df$FisherP),]->res_df;
	p.adjust(res_df$FisherP)->res_df$FDR;
	return(res_df);
}
compare_icluster_cnv_count(cnv_genes_changed_merged,2:88,c("iC1","iC2"),92)->CNVCor_METCor_df_genes_cnvCount;
write.table(CNVCor_METCor_df_genes_cnvCount,"multi-omics/CNVCor_METCor_df_genes_cnvCount.txt",row.names=F,sep="\t",quote=F)
#----------------compare iC1 and iC2 MET: CNVCor_METCor_iC_cluster+CNVCor_METCor_df_genes+myd.methy_t
myd.methy[which(myd.methy$Symbol%in%CNVCor_METCor_df_genes),]$Probe->methy_probes;
myd.methy_t[,as.character(methy_probes)]->myd.methy_t_changed;
for(i in 1:ncol(myd.methy_t_changed)){
	for(j in 1:nrow(myd.methy_t_changed)){
			if(myd.methy_t_changed[j,i]>0.8){
				1->myd.methy_t_changed[j,i];##methylation hyper
			}else if(myd.methy_t_changed[j,i]<0.2){
				(-1)->myd.methy_t_changed[j,i];##methylation hypo
			}else{
				0->myd.methy_t_changed[j,i];##normal
			}
	}
}
data.frame("Sample"=rownames(myd.methy_t),myd.methy_t_changed)->methy_t_changed;
merge(methy_t_changed,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="Sample")->methy_t_changed_merged;
compare_icluster_cnv_count(methy_t_changed_merged,2:1269,c("iC1","iC2"),1273)->CNVCor_METCor_df_genes_metCount;

myd.methy[which(myd.methy$Probe%in%CNVCor_METCor_df_genes_metCount$ID),]->myd.methy_METCor_filter;
c()->test_gnames;
for(i in CNVCor_METCor_df_genes_metCount$ID){
	for(j in 1:nrow(myd.methy_METCor_filter)){
		if(i==myd.methy_METCor_filter[j,4]){
			c(test_gnames,as.character((myd.methy_METCor_filter[j,5])))->test_gnames;
			break;
		}
	}
}
test_gnames->CNVCor_METCor_df_genes_metCount$gName;
write.table(CNVCor_METCor_df_genes_metCount,"multi-omics/CNVCor_METCor_df_genes_metCount.txt",quote=F,sep="\t",row.names=F)
################--------------draw CNVCor_METCor_df_genes heatmap
preorder_cnvcor_metcor_df_genes_heatmap<-function(myd,sample_clusters,clusters,geneOrder){
	myd[,as.character(geneOrder)]->myd_order;
	data.frame("Sample"=myd$Sample,myd_order)->myd_order;
	merge(myd_order,sample_clusters,by.x="Sample",by.y="Sample")->myd_merged;
	c()->order_index;
	for(i in clusters){
		which(myd_merged$iCluster==i)->tmp_index;
		hclust(dist(myd_merged[tmp_index,as.character(geneOrder)]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(order_index,tmp_index)->order_index;
	}
	myd_merged[order_index,]->myd_merged;
	return(myd_merged)
}
brewer.pal(11,"RdYlBu")[c(11,7,1)]
preorder_cnvcor_metcor_df_genes_heatmap(cnv_genes_changed,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3"),CNVCor_METCor_df_genes_cnvCount$ID)->cnv_genes_changed_order;
generate_gaps(as.numeric(substr(cnv_genes_changed_order$iCluster,3,3)))->gaps_col;
list("CNVCor_C"=generate_color(cnv_genes_changed_order$CNVCor_C),"METCor_C"=generate_color(cnv_genes_changed_order$METCor_C),"iCluster"=generate_color(cnv_genes_changed_order$iCluster))->col_colors;
pheatmap(t(cnv_genes_changed_order[,2:88]),color=c("#313695","#F0F0F0","#A50026"),cluster_rows=F,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellheight=4,border_color=NA,annotation_col=cnv_genes_changed_order[,c("CNVCor_C","METCor_C","iCluster")],annotation_colors=col_colors,silent=F,labels_row=CNVCor_METCor_df_genes_cnvCount$ID)->cnv_heatmap_gtable;
cnv_heatmap_gtable#,border_color=brewer.pal(9,"Greys")[2]
#--------for MET: 
preorder_cnvcor_metcor_df_genes_heatmap(methy_t_changed,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3"),CNVCor_METCor_df_genes_metCount$ID)->met_genes_changed_order;
pheatmap(t(met_genes_changed_order[,2:101]),color=c("#313695","#F0F0F0","#A50026"),cluster_rows=F,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellheight=4,border_color=NA,annotation_col=met_genes_changed_order[,c("CNVCor_C","METCor_C","iCluster")],annotation_colors=col_colors,labels_row=CNVCor_METCor_df_genes_metCount$gName[1:100])#->met_heatmap_gtable;
#--------for EXP:
as.matrix(myd_ht_count[,-c(1,2)])->myd_readcount_matrix_icluster;
as.character(myd_ht_count[,1])->rownames(myd_readcount_matrix_icluster)

myd_readcount_matrix_icluster[CNVCor_METCor_df_genes,]->myd_readcount_matrix_iC1_iC2;
t(myd_readcount_matrix_iC1_iC2)->myd_readcount_matrix_iC1_iC2;
data.frame("Sample"=gsub("\\.","-",rownames(myd_readcount_matrix_iC1_iC2)),myd_readcount_matrix_iC1_iC2)->myd_readcount_matrix_iC1_iC2
preorder_cnvcor_metcor_df_genes_heatmap(myd_readcount_matrix_iC1_iC2,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3"),colnames(myd_readcount_matrix_iC1_iC2)[-1])->myd_readcount_matrix_iC1_iC2_order;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(130),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
pheatmap(t(log2(myd_readcount_matrix_iC1_iC2_order[,2:149]+1)),color=fill_colors,cluster_rows=T,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,show_rownames=F,cellheight=3,border_color=NA,annotation_col=myd_readcount_matrix_iC1_iC2_order[,c("CNVCor_C","METCor_C","iCluster")],annotation_colors=col_colors)
############################################################################################################
#----draw cnv count, met count, exp in iC1 ad iC2,KM in TCGA and KMplot: 
CNVCor_METCor_df_genes_cnvCount[CNVCor_METCor_df_genes_cnvCount$FDR<0.05,"ID"]->cnvCount_filter_genes;
CNVCor_METCor_df_genes_metCount[CNVCor_METCor_df_genes_metCount$FDR<0.05,"gName"]->metCount_filter_genes;
intersect(intersect(cnvCount_filter_genes,metCount_filter_genes),CNVCor_METCor_df_genes)->cnvCount_metCount_deseq_filter_genes;
#--------------------------------------------------------for survival analysis
########read EXP file
fread("multi-omics/HTSeq_FPKM_STAD.merge.SYMBOL.factors.txt",header=T,sep="\t",stringsAsFactors=F)->myd_exp_raw
as.data.frame(myd_exp_raw)->myd_exp_raw
#####################process EXP:
preprocess_expd<-function(expd,gStartColumn,aliveEvent){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove IRGP value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		length(expd[expd[,i]==0,i])->failed_samples;
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}
	expd[,c(1:(gStartColumn-1),filter_colums)]->expd.filter;
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
13->g_start_column
preprocess_expd(myd_exp_raw,g_start_column,"Alive")->myd_exp_raw_filter;
which(colnames(myd_exp_raw_filter)%in%cnvCount_metCount_deseq_filter_genes)->tmp_index;
cox_univariant_gene_regr(myd_exp_raw_filter,tmp_index)->shared_genes.coxph
draw_gene_survial_curve_custom<-function(myd,g,bk,myd_colors){
	myd[myd[,g]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,g]),]->myd.rm;
	median(myd[,g],rm.na=T)->g_median;
	unlist(lapply(myd[,g],function(x){if(x>g_median){"H"}else{"L"}}))->g_HL;
	g_HL->myd.rm$ExpLevel;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~ExpLevel,data=myd.rm)->myd.fit;
	survdiff(formula=myd.surv~ExpLevel,rho=0,data=myd.rm)->myd.diff;
	table(myd.rm$ExpLevel)->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("Method Kaplan Meier(",g,")",sep=""),xlim=c(0,max_xlim));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(n=",myd.table,")",sep="")),fill=myd_colors,bty="n");
	if(pvalue<1e-5){
		legend(x=bk/2,y=0.2,legend="p < 1e-5",bty="n",cex=1.2,pos=2,adj=0.5,font=2)
	}else{
		legend(x=bk/2,y=0.2,legend=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.2,pos=2,adj=0.5,font=2)
	}
	return(c(pvalue,myd.table));
}
draw_gene_survial_curve_custom_cut<-function(myd,g,bk,cuts,myd_colors){
	myd[myd[,g]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,g]),]->myd.rm;
	as.numeric(myd.rm[,g])->myd.rm[,g];
	sort(myd.rm[,g],na.last=NA)->g_x;
	c(min(g_x,na.rm=T))->g_cuts;
	c("L0")->g_cuts_level;
	for(i in 1:(cuts-1)){
		round(length(g_x)*(i/cuts))->i_index;
		c(g_cuts,g_x[i_index])->g_cuts;
		c(g_cuts_level,paste("L",i,sep=""))->g_cuts_level;
	}
	c(g_cuts,max(g_x,na.rm=T))->g_cuts;
	c(g_cuts_level,paste("L",cuts,sep=""))->g_cuts_level;
	c()->g_HL;
	for(i_g in myd.rm[,g]){
		for(i in 1:(length(g_cuts)-1)){
			if(i_g>=g_cuts[i] && i_g<=g_cuts[i+1]){
				c(g_HL,g_cuts_level[i+1])->g_HL;
				break;
			}
		}
	}
	g_HL->myd.rm$ExpLevel;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~ExpLevel,data=myd.rm)->myd.fit;
	survdiff(formula=myd.surv~ExpLevel,rho=0,data=myd.rm)->myd.diff;
	table(myd.rm$ExpLevel)->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("Method Kaplan Meier(",g,")",sep=""),xlim=c(0,max_xlim));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(n=",myd.table,")",sep="")),fill=myd_colors,bty="n");
	if(pvalue<1e-5){
		text(x=bk/2,y=0.1,labels="p < 1e-5",bty="n",cex=1.1,pos=4,adj=0.5,font=3)
	}else{
		text(x=bk/2,y=0.1,labels=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.1,pos=4,adj=0.5,font=3)
	}
	return(c(pvalue,myd.table));
}
#draw_gene_survial_curve_custom_cut(myd_exp_processed,"IDO1",2000,3,brewer.pal(11,"Spectral")[c(4,10,11)])->tmp_res;
do_gene_sa<-function(genes){
	c()->km_pvalue;
	for(g in genes){
		draw_gene_survial_curve_custom_cut(myd_exp_raw_filter,g,500,3,myd_colors)->res;
		c(km_pvalue,res[1])->km_pvalue;
	}
	data.frame("gName"=genes,"KMsurv"=km_pvalue)->res_df;
	return(res_df);
}
do_gene_sa(cnvCount_metCount_deseq_filter_genes)->CNVCor_METCor_df_genes_sa;
CNVCor_METCor_df_genes_sa[order(CNVCor_METCor_df_genes_sa$KMsurv),]->CNVCor_METCor_df_genes_sa
#---------------------------
CNVCor_METCor_df_genes_sa[CNVCor_METCor_df_genes_sa$KMsurv<0.05,1]->filtered_gnames;
draw_gene_survial_curve_custom_cut(myd_exp_raw_filter,c("KCNB1"),500,3,myd_colors)->tmp_res;
list()->res;
#par(mfrow=c(4,5))
for(g in filtered_gnames){
	draw_gene_survial_curve_custom_cut(myd_exp_raw,g,1000,3)->tmp_res;
	c(res,list(tmp_res))->res;
}
generate_group(iC_group_condition,c("iC1","iC2"))->iC_filter_group;
draw_cnvMET_count_exp<-function(cnv_count,met_count,exp_count,groups,g,legend_pos){
	brewer.pal(11,"Set3")->myd.colors;
	par(mfrow=c(1,4));
	which(cnv_count$ID==g)[1]->g_index;
	sum(cnv_count[g_index,c(2,4,6)])->ic1_sum;
	sum(cnv_count[g_index,c(3,5,7)])->ic2_sum;
	c(cnv_count[g_index,2]/ic1_sum,cnv_count[g_index,4]/ic1_sum)->cnv_ic1_count;
	c(cnv_count[g_index,3]/ic2_sum,cnv_count[g_index,5]/ic2_sum)->cnv_ic2_count;
	barplot(cnv_ic1_count,ylim=c(-1,1),col=myd.colors[1],axes=F,main=paste(g," CNV Level",sep=""),ylab="Ratio");
	barplot(-cnv_ic2_count,add=T,col=myd.colors[2],axes=F)->bar_mg;
	axis(side=2,at=seq(-1,1,0.2),labels=c(seq(1,0,-0.2),seq(0.2,1,0.2)))
	text(x=bar_mg,y=-0.9,labels=c("CNV Gain","CNV Loss"))
	cnv_count$FDR[g_index]->tmp_pvalue;
	if(tmp_pvalue<1e-5){
		legend("topleft",legend=paste("FDR: ","< 1e-5",sep=""))
	}else{
		legend("topleft",legend=paste("FDR: ",round(tmp_pvalue,5),sep=""))
	}
	legend("topright",legend=c("iC1","iC2"),fill=myd.colors[1:2])
	#--
	which(met_count$gName==g)[1]->g_index;
	sum(met_count[g_index,c(2,4,6)])->ic1_sum;
	sum(met_count[g_index,c(3,5,7)])->ic2_sum;
	c(met_count[g_index,2]/ic1_sum,cnv_count[g_index,4]/ic1_sum)->met_ic1_count;
	c(met_count[g_index,3]/ic2_sum,cnv_count[g_index,5]/ic2_sum)->met_ic2_count;
	barplot(met_ic1_count,ylim=c(-1,1),col=myd.colors[1],axes=F,main=paste(g," Methylation Level",sep=""),ylab="Ratio");
	barplot(-met_ic2_count,add=T,col=myd.colors[2],axes=F)->bar_mg;
	axis(side=2,at=seq(-1,1,0.2),labels=c(seq(1,0,-0.2),seq(0.2,1,0.2)))
	text(x=bar_mg,y=-0.9,labels=c("Hyper Methy","Hypo Methy"))
	met_count$FDR[g_index]->tmp_pvalue;
	if(tmp_pvalue<1e-5){
		legend("topleft",legend=paste("FDR: ","< 1e-5",sep=""))
	}else{
		legend("topleft",legend=paste("FDR: ",round(tmp_pvalue,5),sep=""))
	}
	legend("topright",legend=c("iC1","iC2"),fill=myd.colors[1:2])
	#-----
	which(exp_count$gName==g)->g_index;
	groups[which(groups$Condition=="iC1"),1]->group_normal;
	groups[which(groups$Condition=="iC2"),1]->group_tumor;
	as.numeric(exp_count[g_index,as.character(group_normal)])->group_normal_exp;
	as.numeric(exp_count[g_index,as.character(group_tumor)])->group_tumor_exp;
	C2_C1_deseq.res[which(C2_C1_deseq.res$gName==g),]->g_dds_res;
	data.frame("Expr"=c(group_normal_exp,group_tumor_exp),"Group"=rep(c("iC1","iC2"),c(length(group_normal),length(group_tumor))))->group_nt_exp;
	g_dds_res[,7]->group_nt_exp_ranktest_pvalue;
	g_dds_res[,3]->log2_fc;
	boxplot(formula=log2(Expr)~Group,data=group_nt_exp,col=brewer.pal(9,"Set3")[1:2],pch=20,main=paste(g," Expression Level",sep=""),ylab="Log2(Expr)");
	if(group_nt_exp_ranktest_pvalue<1e-5){
		c(paste("DESeq2 padj: ","< 1e-5",sep=""),paste("iC2/iC1 FC: ",round(2^log2_fc,5),sep=""))->group_legend;
	}else{
		c(paste("DESeq2 padj: ",round(group_nt_exp_ranktest_pvalue,5),sep=""),paste("iC2/iC1 FC: ",round(2^log2_fc,5),sep=""))->group_legend;
	}
	legend(legend_pos,legend=group_legend)
	#---os km
	draw_gene_survial_curve_custom_cut(myd_exp_raw_filter,g,500,3,brewer.pal(11,"Spectral")[c(4,10,11)])->res;
}
draw_cnvMET_count_exp(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,myd_ht_count,iC_filter_group,"KCNB1","bottomleft")#
draw_cnvMET_count_exp(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,myd_ht_count,iC_filter_group,"PLCXD3","topleft")#
draw_cnvMET_count_exp(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,myd_ht_count,iC_filter_group,"JPH3","bottomleft")#
draw_cnvMET_count_exp(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,myd_ht_count,iC_filter_group,"RIMS4","topleft")
draw_cnvMET_count_exp(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,myd_ht_count,iC_filter_group,"ARHGAP40","topleft")
#------------
retrive_gene_CNV_MET_EXP_pvalue<-function(cnv_count,met_count,exp_count,genes){
	data.frame("V1"=as.character(0),"V2"=as.character(0),"V3"=as.character(0),"V4"=as.character(0),"V5"=as.character(0),"V6"=as.character(0),"V7"=as.character(0),"V8"=as.character(0))->res;
	for(g in genes){
		cnv_count[which(cnv_count$ID==g),]->cnv_p;
		print(cnv_p);
		met_count[met_count$gName==g,]->met_p;
		print(met_p);
		exp_count[exp_count$gName==g,]->exp_p;
		print(exp_p);
		flush.console();
	}
	#return(res)
}
retrive_gene_CNV_MET_EXP_pvalue(CNVCor_METCor_df_genes_cnvCount,CNVCor_METCor_df_genes_metCount,C2_C1_deseq.res,cnvCount_metCount_deseq_filter_genes);
###################################################################################################################
#------------GEO data：
fread("multi-omics/GSE62254_exp_symbol.merged.factors.txt",header=T,sep="\t",stringsAsFactors=F)->GSE62254_exp;
as.data.frame(GSE62254_exp)->GSE62254_exp;
0->GSE62254_exp[which(GSE62254_exp$A2_Event==1),3]
preprocess_expd(GSE62254_exp,19,0)->GSE62254_exp_processed;
par(mfrow=c(1,2))
draw_gene_survial_curve_custom_cut(GSE62254_exp_processed,"KCNB1",500,2,brewer.pal(11,"Spectral")[c(4,10,11)])->res;#not good
draw_gene_survial_curve_custom_cut(GSE62254_exp_processed,"JPH3",500,2,brewer.pal(11,"Spectral")[c(4,10,11)])->res;#not good
#draw_gene_survial_curve_custom_cut(GSE62254_exp_processed,"PLCXD3",500,2,brewer.pal(11,"Spectral")[c(4,10,11)])->res;#not good
#draw_gene_survial_curve_custom_cut(GSE62254_exp_processed,"RIMS4",500,2,brewer.pal(11,"Spectral")[c(4,10,11)])->res;#not good
#draw_gene_survial_curve_custom_cut(GSE62254_exp_processed,"ARHGAP40",500,2,brewer.pal(11,"Spectral")[c(4,10,11)])->res;#not good
#------------
tmp.filter->expd;
c("Stage_T","Stage_N","Stage_M","Stage","Grade","Age","Gender")->clinical_features;
par(mar=c(8,3,4,3))
"Gender"->f;#change clinical features to draw boxplot
draw_boxplot_genes_by_factors_v2(expd,cnvCount_metCount_deseq_filter_genes,f,brewer.pal(9,"Set1"))->test_;
###################################################################################################################
read.table("multi-omics/STAD.mutect2.snv_vcf2matrix.res",header=T,sep="\t",stringsAsFactors=F)->myd_snv;
myd_snv[myd_snv$Variant_Classification!="Silent",]->myd_snv
myd_snv[myd_snv$Variant_Classification!="Intron",]->myd_snv

cal_snv_exp_cor<-function(snvd,expd,genes,groups){
	c()->gName;
	c()->snvGene;
	c()->snv_pos;
	c()->CorrP;
	c()->Corr;
	for(g in genes){
		#snvd[which(snvd$Hugo_Symbol==g),]->snv_values;
		snvd[,as.character(groups$SampleID)]->snv_state;
		gsub("\\.","-",groups$SampleID)->tmp_samples;
		expd[tmp_samples,as.character(g)]->exp_values;
		nrow(snv_state)->snv_count;
		for(i in 1:snv_count){
			if(length(snv_state[i,])==length(exp_values)){
				cor.test(as.numeric(snv_state[i,]),exp_values)->snv_exp.cor;
				c(gName,g)->gName;
				c(snvGene,snvd$Hugo_Symbol[i])->snvGene;
				c(snv_pos,snvd$HGVSc[i])->snv_pos;
				c(CorrP,snv_exp.cor$p.value)->CorrP;
				c(Corr,snv_exp.cor$estimate)->Corr;
			}
		}
	}
	data.frame("gName"=gName,"snvGene"=snvGene,"SNV"=snv_pos,"CorrP"=CorrP,"Corr"=Corr)->res;
	#res[!is.na(res$CorrP),]->res;
	#res[order(res$CorrP),]->res;
	return(res);
}
"JPH3"->snv_exp_cor_g;
cal_snv_exp_cor(myd_snv,CNVCor_genes_EXP,c(snv_exp_cor_g),iC_filter_group)->snv_exp_cor_g_cor;
snv_exp_cor_g_cor[!is.na(snv_exp_cor_g_cor$CorrP),]->snv_exp_cor_g_cor;
write.table(snv_exp_cor_g_cor,paste(c("multi-omics/",snv_exp_cor_g,"_snv_exp_cor.txt"),collapse=""),quote=F,sep="\t",row.names=F)

#-----------------------------------------------------------------------
cal_snv_group_ftest<-function(snvd,groups){
	c()->gName;
	c()->snv_pos;
	c()->group1_Mut;
	c()->group1_Nor;
	c()->group2_Mut;
	c()->group2_Nor;
	c()->FisherP;
	groups[groups$Condition=="iC1",1]->group1;
	groups[groups$Condition=="iC3",1]->group2;
	for(i in 1:nrow(snvd)){
		snvd[i,as.character(group1)]->group1_values;
		snvd[i,as.character(group2)]->group2_values;
		c(length(group1_values[group1_values!=0]),length(group1_values[group1_values==0]))->group1_count;
		c(length(group2_values[group2_values!=0]),length(group2_values[group2_values==0]))->group2_count;
		matrix(c(group1_count,group2_count),nrow=2,byrow=F)->group1_group2_table;
		#print(group1_group2_table);
		#flush.console();
		fisher.test(group1_group2_table)->group1_group2_table_test;
		#p.adjust(group1_group2_table_test$p.value,n=sum(group1_count,group2_count))->test_padj;
		c(gName,snvd$Hugo_Symbol[i])->gName;
		c(snv_pos,snvd$HGVSc[i])->snv_pos;
		c(FisherP,group1_group2_table_test$p.value)->FisherP;
		c(group1_Mut,group1_count[1])->group1_Mut;
		c(group1_Nor,group1_count[2])->group1_Nor;
		c(group2_Mut,group2_count[1])->group2_Mut;
		c(group2_Nor,group2_count[2])->group2_Nor;
		#c(padj,test_padj)->padj;
	}
	data.frame("gName"=gName,"SNV"=snv_pos,"FisherP"=FisherP,"Mut_iC1"=group1_Mut,"Nor_iC1"=group1_Nor,"Mut_iC3"=group2_Mut,"Nor_iC3"=group2_Nor)->res;
	p.adjust(res$FisherP,n=nrow(res))->res$padj;
	res[order(res$FisherP),]->res;
	return(res);
}
cal_snv_group_ftest(myd_snv,iC_filter_group)->myd_snv_group_ftest;
write.table(myd_snv_group_ftest,"multi-omics/myd_snv_group_ftest.txt",quote=F,sep="\t",row.names=F)
summary_gene_snv<-function(snvd,sample_cols){
	snvd[snvd$Variant_Classification!="Intron",]->snvd_filter;
	snvd_filter[snvd_filter$Variant_Classification!="Silent",]->snvd_filter
	names(table(as.character(snvd_filter$Hugo_Symbol)))->genes;
	c()->snvd_counts;
	for(g in genes){
		snvd_filter[snvd_filter$Hugo_Symbol==g,]->g_snvd_filter;
		for(i in sample_cols){
			g_snvd_filter[,i]->g_snvd_values;
			0->g_snvd_count;
			if(length(g_snvd_values[g_snvd_values!=0])!=0){
				length(g_snvd_values[g_snvd_values!=0])->g_snvd_count;
			}
			c(snvd_counts,g_snvd_count)->snvd_counts;
		}
	}
	matrix(snvd_counts,nrow=length(genes),ncol=length(sample_cols),byrow=T)->tmp_matrix;
	as.data.frame(tmp_matrix)->tmp_matrix;
	colnames(snvd)[sample_cols]->colnames(tmp_matrix);
	data.frame("gName"=genes,tmp_matrix)->tmp_matrix;
	return(tmp_matrix);
}
#summary_gene_snv(myd_snv,18:ncol(myd_snv))->gene_snv_merged;
#save.image("multi-omics/work-20190227.RData")
read.table("multi-omics/STAD.mutect2.gene_vcf2matrix.res",header=T,sep="\t",stringsAsFactors=F)->gene_snv_merged
cal_gene_snv_group_ftest<-function(snvd,groups){
	snvd->tmp_matrix;
	#--------------------------------
	c()->gName;
	c()->group1_Mut;
	c()->group1_Nor;
	c()->group2_Mut;
	c()->group2_Nor;
	c()->FisherP;
	#c()->padj;
	groups[groups$Condition=="iC1",1]->group1;
	groups[groups$Condition=="iC2",1]->group2;
	for(i in 1:nrow(tmp_matrix)){
		tmp_matrix[i,as.character(group1)]->group1_values;
		tmp_matrix[i,as.character(group2)]->group2_values;
		c(length(group1_values[group1_values!=0]),length(group1_values[group1_values==0]))->group1_count;
		c(length(group2_values[group2_values!=0]),length(group2_values[group2_values==0]))->group2_count;
		matrix(c(group1_count,group2_count),nrow=2,byrow=F)->group1_group2_table;
		#print(group1_group2_table);
		#flush.console();
		fisher.test(group1_group2_table)->group1_group2_table_test;
		#p.adjust(group1_group2_table_test$p.value,n=sum(group1_count,group2_count))->test_padj;
		c(FisherP,group1_group2_table_test$p.value)->FisherP;
		c(group1_Mut,group1_count[1])->group1_Mut;
		c(group1_Nor,group1_count[2])->group1_Nor;
		c(group2_Mut,group2_count[1])->group2_Mut;
		c(group2_Nor,group2_count[2])->group2_Nor;
		#c(padj,test_padj)->padj;
	}
	data.frame("gName"=tmp_matrix$gName,"FisherP"=FisherP,"Mut_iC1"=group1_Mut,"Nor_iC1"=group1_Nor,"Mut_iC2"=group2_Mut,"Nor_iC2"=group2_Nor)->res;
	res[order(res$FisherP),]->res;
	p.adjust(res$FisherP,n=nrow(res))->res$padj;
	return(res);
}
cal_gene_snv_group_ftest(gene_snv_merged,iC_filter_group)->gene_snv_count;
write.table(gene_snv_count,"multi-omics/gene_snv_count.txt",quote=F,sep="\t",row.names=F)

gene_snv_count[gene_snv_count$padj<0.1,1]->genes_snv_count_filter;
unlist(lapply(genes_snv_count_filter,function(x){which(gene_snv_merged$gName==x)->i;i}))->tmp_index;
gene_snv_merged[tmp_index,]->gene_snv_merged_filter;

preoder_snv_gene_samples<-function(snvd,groups,clusters){
	data.frame("gName"=snvd$gName)->res;
	for(g in groups){
		unlist(lapply(g,function(x){which(clusters$Condition==x)->i;i}))->tmp_index;
		snvd[,as.character(clusters[tmp_index,1])]->snvd_filter
		unlist(lapply(1:ncol(snvd_filter),function(x){which(snvd_filter[,x]!=0)->x_index;length(snvd_filter[x_index,x])}))->tmp_index_count;
		cbind(res,snvd_filter[,order(tmp_index_count)])->res;
	}
	res[,-1]->res;
	return(res);
}
preoder_snv_gene_samples(gene_snv_merged_filter,c("iC1","iC2","iC3"),iC_group_condition)->gene_snv_merged_filter_sort;
unlist(lapply(c("iC1","iC2","iC3"),function(x){which(iC_group_condition$Condition==x)->i;i}))->tmp_index;
CNVCor_METCor_iC_cluster[tmp_index,c(2,3,4)]->tmp_annotation_col;
gsub("-","\\.",as.character(CNVCor_METCor_iC_cluster[tmp_index,1]))->rownames(tmp_annotation_col)

pheatmap(gene_snv_merged_filter_sort,cluster_rows=F,cluster_cols=F,labels_row=genes_snv_count_filter,color=c(brewer.pal(9,"Greys")[2],colorRampPalette(brewer.pal(9,"Greys")[5:9])(30)),cellheight=15,gaps_col=gaps_col,show_colnames=F,annotation_col=tmp_annotation_col,annotation_colors=col_colors,border_color="white")




















