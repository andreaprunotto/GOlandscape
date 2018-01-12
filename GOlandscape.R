#################################################
## This is GOlandscape, the Gene Ontology tool ##
## conceived and developed by Andrea Prunotto  ##
##       Lausanne 2016 - Stockholm 2017        ## 
##   © Unil & Karolinska Institutet (2017)     ##      
#################################################

 ###############################################
 # This is the developemental version n. 7     #
 # See Documentation in the annexed folder     #
 # for detailed explanations about the methods #
 ###############################################

# INPUT: 
# (1) A (possibly tab separated) text file with at least 2 fields
# 1 field with gene ids, 1 field with a ranking variable
# Optionally, one field with a selection variable (FC)
# The file is stored in the same folder where GOlandscape.R is located

# OUTPUT: 
# (1) A list of tab separated text files goseq.*.tab which will be stored
# in the folder easygo_support, a subfolder of the working directory. 
# These are information that can be used later to save computational time.
# (2) A heatmap GOlandscape.*.GO_Landscape.pdf containing the stepwise thresholding 
# (3) A heatmap GOlandscape.*.Gene_Landscape.pdf connecting
# the most relevant DE genes with the most likley associated GO terms
# (4) A heatmap GOlandscape.*.complete.pdf containing both (2) and (3)
# (2,3,4) are saved in easygo_result, a subfolder of the working folder

# OPTIONS: See "Input Session"

# REMARKS: You need a connection to the Internet to work with this tool
# Gene IDs must be UNIQUE
# DE studies with too low significance (e.g. min. p-value = 0.01) won't work!
# Results' significance is given in 

 
################################################
################### Have Fun! ##################
################################################

# Loading essential libraries

library(goseq)
library(KEGG.db)
library(ComplexHeatmap)
library(circlize)

# General options

options(stringsAsFactors=FALSE, showWarnings=FALSE)
dir.create('easygo_results')
dir.create('easygo_support')


####################################################################
########################### INPUT SESSION ##########################
############## Here you can specify the input options ############## 
####################################################################

# Main cycle over the file names

for(filename in c('FILENAME.tab')) 
{
# Is there an header? [mandatory] default: TRUE
withhead=TRUE
# Field separator of $filename? (e.g. "\t", " ") (character) default: '\t'
separsym='\t'
# Is there any selection variable? e.g. logFC default: TRUE
addselvar=TRUE
# Which field contains the (unique) labels? (e.g. gene names) (character)  default: 1
# LABELS SHOULD BE UNIQUE, otherwise the software stop! 
labfield=1
# What field contains the ranking variable (e.g. p-values) (numeric) [mandatory] default: 2
ranfield=4
# How many significant digits you would like to trust (suggestion based on max. likelihood: 2) defaut: 2
sidigits=2
# What field contains the (signed) selection variable (if any)? (e.g. logFC) (numeric) [optional] default: 3
if(addselvar){selfield=2}
# Would you like to -log10 normalize the ranking variable? default: TRUE
log10ran=TRUE
# Would you like to print an overall distribution of n. of labels vs. rank? default: TRUE
plotdist=TRUE
# Would you like to print go size vs. rank? default: TRUE
plotgosrank=TRUE

# Genome and genes info
# Which genome? 
# Which input gene id?
# Which output gene id?
# Allowed IDs (so far): ENSEMBL SYMBOL REFSEQ ENTREZ

genomeid='Rattus.norvegicus'
ingeneid='ENSEMBL'
ougeneid='SYMBOL'

# Number of categories to show in the landscape
# BP: Biological Process CC: Cellular Component
# MF: Molecular Function KEGG: Pathways
# You can also add the relevance of the not annot. genes
# Number of genes appearing in the (gene) landscape

nbpcat=30
ncccat=30
nmfcat=30
nggcat=30
nuncat=0
ngenes=30

# Set the minimal value of the ranking variable
# The lower the number, the slower the calculation
# See Documentation!

miniran=2

# Automatic stepwise threshold sampling definition
# bygene: (keeping the n. of genes in each step constant) 
# byde: (keeping the amount of DE in each step constant)
# See Documentation before changing the default!

autostep=TRUE
automethod='byde' 
maxnstep=15
deincrement=5 

# Manual stepwise threshold sampling 
# To use only for compare with GSEA!
stepran=2

# Run the goseq analysis?
# Override previous geseq data?
# Just produce the data or also print the maps?

rungoseq=FALSE
forcegoseq=TRUE
printmap=TRUE

###########################################################################
########################### end of input session ##########################
###########################################################################




# Preparing the table which connects genes and GO
# For each genome (so far: Mus.musculus, Rattus.norvegicus, Homo.sapiens)

if(genomeid=='Mus.musculus'){     
    library(org.Mm.eg.db)
    library(Mus.musculus)
    ense=toTable(org.Mm.egENSEMBL)
    symb=toTable(org.Mm.egSYMBOL)
    refs=toTable(org.Mm.egREFSEQ)
    entr=ense; entr$ensembl_id=entr$gene_id
    entr=unique(entr);colnames(entr)[2]='entrez_id'
    goseqgenome=c('mm9')
    goca=toTable(org.Mm.egGO)
    path=toTable(org.Mm.egPATH)
}

if(genomeid=='Rattus.norvegicus'){
    library(org.Rn.eg.db)
    library(Rattus.norvegicus)
    ense=toTable(org.Rn.egENSEMBL)
    symb=toTable(org.Rn.egSYMBOL)
    refs=toTable(org.Rn.egREFSEQ)
    entr=ense; entr$ensembl_id=entr$gene_id
    entr=unique(entr);colnames(entr)[2]='entrez_id'
    goseqgenome=c('rn4')
    goca=toTable(org.Rn.egGO)
    path=toTable(org.Rn.egPATH)
}

if(genomeid=='Homo.sapiens'){ 
    library(org.Hs.eg.db)
    library(Homo.sapiens)
    ense=toTable(org.Hs.egENSEMBL)
    symb=toTable(org.Hs.egSYMBOL)
    refs=toTable(org.Hs.egREFSEQ)
    entr=ense; entr$ensembl_id=entr$gene_id
    entr=unique(entr);colnames(entr)[2]='entrez_id'
    goseqgenome=c('hg38')
    goca=toTable(org.Hs.egGO)
    path=toTable(org.Hs.egPATH)
}

idlist=list(ense,symb,refs,entr)

if(ingeneid=='ENSEMBL'){inputid=idlist[[1]];goseqid='ensGene'}
if(ingeneid=='SYMBOL'){ inputid=idlist[[2]];goseqid='geneSymbol'}
if(ingeneid=='REFSEQ'){ inputid=idlist[[3]];goseqid='refGene'}
if(ingeneid=='ENTREZ'){ inputid=idlist[[4]];goseqid='knownGene'}

if(ougeneid=='ENSEMBL'){shownid=idlist[[1]]}
if(ougeneid=='SYMBOL'){ shownid=idlist[[2]]}
if(ougeneid=='REFSEQ'){ shownid=idlist[[3]]}
if(ougeneid=='ENTREZ'){ shownid=idlist[[4]]}

# Read the file with at least (unique) labels and a numeric variable

file=read.table(filename,header=withhead,sep=separsym,check.names=FALSE,quote='@"')
file=file[!duplicated(file[,labfield]),]
rownames(file)=file[,labfield]
file0=file

# Build the raw data file
# Distinguishing up/down/both regulation
# datalist will contain all the files to be processed by goseq

if(addselvar)
{
	data=as.data.frame(cbind(rownames(file),file[,ranfield],file[,selfield]))
	colnames(data)=c('lab','ran','sel')
	data$lab=as.character(data$lab)
	data$ran=as.numeric(data$ran)
	data$sel=as.numeric(data$sel)
	both=data
	plus=data[data$sel>0,]
	minu=data[data$sel<0,]
	datalist=list(both,plus,minu)
}

if(!addselvar)
{
	data=as.data.frame(cbind(rownames(file),file[,ranfield]))
	colnames(data)=c('lab','ran')
	data$lab=as.character(data$lab)
	data$ran=as.numeric(data$ran)
	both=data
	datalist=list(both)
}

# Perform the analysis on datalist (both plus minus or only both)
# workdata contains each item of datalist

for(w in 1:length(datalist))
{

if(w==1){workdata=both;workname='up-and-down-regulated';filelabel='up- and down-regulated genes'}
if(w==2){workdata=plus;workname='up-regulated';filelabel='only up-regulated genes'}
if(w==3){workdata=minu;workname='down-regulated';filelabel='only down-regulated genes'}

workdata=workdata[,1:2]

# Normalize the p-value with log10 (if requested) 
# and approximate pvalues down to sidigiti significative digits. 
# It should not be more than 2, however (See Max. likelihood theorem).

if(log10ran){
	if(nrow(workdata[workdata$ran==0,])>0)
	{
		workdata[workdata$ran==0,]$ran=0.1*min(workdata[workdata$ran>0,]$ran)
	}
	workdata$ran=-log10(workdata$ran)
	workdata$ran=round(workdata$ran,sidigits)
	}else{
		workdata$ran=round(workdata$ran,sidigits)
	}

head(workdata,5)

showdegorig=workdata 

# Gene Id (in/out) conversion

if(ingeneid=='ENSEMBL'){showdegsymb=merge(showdegorig,ense,by.x='lab',by='ensembl_id',all.x=TRUE)}
if(ingeneid=='SYMBOL'){ showdegsymb=showdegorig}
if(ingeneid=='REFSEQ'){ showdegsymb=merge(showdegorig,refs,by.x='lab',by='accession',all.x=TRUE)}
if(ingeneid=='ENTREZ'){ showdegsymb=merge(showdegorig,entr,by.x='lab',by='entrez_id', all.x=TRUE)}

# Ontology database building  

labbackg=inputid[,2]
labtogo=merge(inputid,goca,by='gene_id')
labtopa=merge(inputid,path,by='gene_id')
labtogo=as.data.frame(cbind(labtogo[,2:3],labtogo[,5]))
colnames(labtogo)=c('gene','category','ontology')
labtopa=as.data.frame(cbind(labtopa[,2:3],'KEGG'))
colnames(labtopa)=c('gene','category','ontology')
genestogo=rbind(labtopa,labtogo)
k=as.data.frame(merge(genestogo,as.data.frame(labbackg),by.x='gene',by.y='labbackg',all.y=TRUE))
if(nrow(k[is.na(k$category),])!=0)
{
	k[is.na(k$category),]$category="unknown"
	k[k$category=="unknown",]$ontology="UNKN"
}

# Remove overall categories
k=k[k$category!="GO:0003674" & k$category!="GO:0005575" & k$category!="GO:0008150" & k$category!="01100",]

genestogo=k

##################################
# Defining the stepwise sampling #
##################################

# Manual stepwise sampling (to compare with GSEA. Not to be used in analysis context)

if(autostep==FALSE)
{

N=nrow(workdata[workdata$ran>=miniran,])
M=N
r=miniran
m=0

c=as.data.frame(cbind(r,M))

while(N!=M | m<max(workdata$ran))
{
  N=nrow(workdata[workdata$ran>=r,])
  M=nrow(workdata[workdata$ran>=r+stepran,])
  m=min(workdata[workdata$ran>=r,]$ran)
  r=r+stepran
  c0=cbind(r,M)
  c=rbind(c0,c)
}

steplist=data.frame()
prec=-1
for(i in 1:nrow(c))
{
  if(c$M[i]!=prec){s0=cbind(c$r[i],c$M[i]);steplist=rbind(s0,steplist)};
  prec=c$M[i]
}

steplist=as.numeric(steplist$V1)
steplist
}


#####################################
# Automatic step list (incremental) #
#####################################

if(autostep){

# Calculation of the amout of DE in the range of DE/ran variation

if(automethod=='byde'){
x=workdata
x=x[order(x$ran,decreasing=TRUE),]
a=rle(x$ran)
b=data.frame(sDE=a$values, nDEg=a$lengths)
b=b[order(b$sDE,decreasing=TRUE),]
b=b[b$sDE>=miniran,]
c=cbind(b,b$sDE*b$nDEg); colnames(c)=c('sDE','nDEg','DE')

ngps=miniran 
nstep=nrow(b)
s=0
steplist=as.data.frame(b$sDE);colnames(steplist)='V5'

while(nstep>maxnstep){
i=1
steplist=data.frame()
n=0
while(!is.na(c$nDEg[i+1])){
	s=0
	while(s<ngps & !is.na(c$nDEg[i+1])){
		s=s+c$DE[i]
		n=n+b$nDEg[i]
		i=i+1
	}
	coll=cbind(maxnstep,nstep,s,n,b$sDE[i])
	steplist=rbind(coll,steplist)
}
print(steplist)
nstep=nrow(steplist)
ngps=ngps+deincrement
}
steplist=as.numeric(steplist$V5)

}

#############################################
# Automatic step list (n. gene incremental) #
#############################################

if(automethod=='bygene'){

# Better not to mess up with workdata (memory)

x=workdata
x=x[order(x$ran,decreasing=TRUE),]
a=rle(x$ran)
b=data.frame(sDE=a$values, nDEg=a$lengths)
b=b[order(b$sDE,decreasing=TRUE),]
b=b[b$sDE>=miniran,]

# Now we run the autosampling...

ngps=1 ; nstep=nrow(b); s=0
steplist=as.data.frame(b$sDE);colnames(steplist)='V5'
while(nstep>maxnstep){
i=1; ss=0; steplist=data.frame()
while(!is.na(b$nDEg[i+1])){
	s=0
	while(s<ngps & !is.na(b$nDEg[i+1])){
		s=s+b$nDEg[i]
		i=i+1
	}
	ss=ss+s
	coll=cbind(maxnstep,nstep,s,ss,b$sDE[i])
	steplist=rbind(coll,steplist)
}
# It is nice to see the software looking for 
# the best sampling at each step!
print(steplist) 
nstep=nrow(steplist)
ngps=ngps+deincrement
}
steplist=as.numeric(steplist$V5)

}

# We got the steplist!
}

# Let's plot it: it is nice to see how many genes per step we have

if(plotdist)
{
x=workdata
x=x[order(x$ran,decreasing=TRUE),]
a=rle(x$ran)
b=data.frame(sDE=a$values, nDEg=a$lengths)
b=b[order(b$sDE,decreasing=TRUE),]
b=b[b$sDE>=miniran,]

d=data.frame()
for(s in c(steplist,max(x$ran))){d0=cbind(s,nrow(x[x$ran>=s,]));d=rbind(d0,d)}
d$V2=as.numeric(d$V2)
d$s=as.numeric(d$s)
d=d[order(d$s,decreasing=FALSE),]

xlim=c(miniran,max(x$ran))
ylim=c(1,max(d$V2)+10)
log='xy'
ptext=data.frame();
for(i in 1:(nrow(d))){ptext=rbind(ptext,d$s[i]+(d$s[i+1]-d$s[i])/2)}
ltext=d$V2

if(workname=='down-regulated'){plotcolor=c('cornflowerblue')}
if(workname=='up-regulated'){plotcolor=c('salmon')}
if(workname=='up-and-down-regulated'){plotcolor=c('green3')}

pdf(paste('./easygo_support/sDE_vs_nDEg.',filename,'.',workname,'.pdf',sep=''))
plot(b,xlim=xlim,ylim=ylim,log=log,col='grey',pch=19,axes=FALSE,xlab='',ylab='');
par(new=TRUE);
plot(d,xlim=xlim,ylim=ylim,log=log,col=plotcolor,type='s',lty=1,axes=FALSE,xlab='',ylab='');
par(new=TRUE);
plot(d,xlim=xlim,ylim=ylim,log=log,col=plotcolor,type='h',lty=1,
	xlab=expression(-log[10](p-value)),
	ylab='number of genes',
	main=paste(filename,'\nnumber of DE genes vs. DE p-value (',filelabel,')',sep='')
	);
text(cbind(ptext,ltext),labels=ltext,font=1,cex=0.8,pos=3,col=plotcolor)
dev.off()
}

#######################################################
####################### MAIN ##########################
#######################################################




if(printmap & length(steplist)>2){

# Calculate the GO p-values with goseq 
# Of course if not already present in easygo_support

if(rungoseq) 
{
w0=data.frame()
for(step in steplist)
{
	if(!file.exists(paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab',sep='')) | forcegoseq)
	{
    	selec=as.vector(as.integer(as.numeric(workdata$ran)>=step));
    	names(selec)=as.character(workdata$lab);
    	pwf=nullp(selec,goseqgenome,goseqid,plot.fit=FALSE);
    	goseqresult=goseq(pwf,goseqgenome,goseqid,use_genes_without_cat=TRUE,
        gene2cat=genestogo) 
    	if(nrow(goseqresult[goseqresult$category=='unknown',]!=0))
    	{
    		goseqresult[goseqresult$category=='unknown',]$ontology='UNKN'
    		goseqresult[goseqresult$category=='unknown',]$term='not annotated'
    	}
    	keggtable=toTable(KEGGPATHID2NAME);
    	tempkegg1=as.data.frame(merge(goseqresult,keggtable,by.y='path_id',by.x='category'));
    	tempkegg1$term=tempkegg1$path_name
    	tempkegg1$ontology='KEGG'
    	tempkegg1=tempkegg1[,1:(ncol(tempkegg1)-1)]
	    goseqresult=goseqresult[!is.na(goseqresult$ontology),]
    	goseqresult=rbind(goseqresult,tempkegg1)
    	w1=cbind(goseqresult,step);
    	w1$over_represented_pvalue=signif(w1$over_represented_pvalue,sidigits)
    	write.table(w1,paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab', sep=''),sep='\t');
	}else{
		w1=read.table(paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab', sep=''), sep='\t',quot='"',header=TRUE,
 		colClass=c('character','character','numeric','numeric','numeric','numeric','character','character','numeric'));
	}
 	w0=rbind(w0,w1)
} # end of cycle on writing goseq result step

}else # if no rungoseq, read the already saved goseq output files in easygo_support
{
w0=data.frame()
for(step in steplist)
{
	if(file.exists(paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab',sep='')) | forcegoseq)
	{
		w1=read.table(paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab', sep=''), sep='\t',quot='"',header=TRUE,
 		colClass=c('character','character','numeric','numeric','numeric','numeric','character','character','numeric'));
 	}else{
 		selec=as.vector(as.integer(as.numeric(workdata$ran)>=step));
    	names(selec)=as.character(workdata$lab);
    	pwf=nullp(selec,goseqgenome,goseqid,plot.fit=FALSE);
    	goseqresult=goseq(pwf,goseqgenome,goseqid,use_genes_without_cat=TRUE,
        gene2cat=genestogo) 
    	if(nrow(goseqresult[goseqresult$category=='unknown',]!=0))
    	{
    		goseqresult[goseqresult$category=='unknown',]$ontology='UNKN'
    		goseqresult[goseqresult$category=='unknown',]$term='not annotated'
    	}
	    keggtable=toTable(KEGGPATHID2NAME);
    	tempkegg1=as.data.frame(merge(goseqresult,keggtable,by.y='path_id',by.x='category'));
    	tempkegg1$term=tempkegg1$path_name
    	tempkegg1$ontology='KEGG'
    	tempkegg1=tempkegg1[,1:(ncol(tempkegg1)-1)]
    	goseqresult=goseqresult[!is.na(goseqresult$ontology),]
    	goseqresult=rbind(goseqresult,tempkegg1)
    	w1=cbind(goseqresult,step);
    	w1$over_represented_pvalue=signif(w1$over_represented_pvalue,sidigits)
    	write.table(w1,paste('./easygo_support/goseqresult.',filename,'.',workname,'.rankstep.',step,'.tab', sep=''),sep='\t');
	}
    w0=rbind(w0,w1)    
} # fi of cycle on steplist for no rungoseq

}#  fi of rungoseq no (reading goseqresult)

# Remove overall categories...

w0=w0[w0$category!="GO:0003674" & w0$category!="GO:0005575" & w0$category!="GO:0008150" & w0$category!="01100",]

# ...and Cleaning abusurd results

if(nrow(w0[w0$over_represented_pvalue=='Inf',])){w0[w0$over_represented_pvalue=='Inf',]=1}
if(nrow(w0[w0$over_represented_pvalue>1,])){w0[w0$over_represented_pvalue>1,]=1}
if(nrow(w0[w0$over_represented_pvalue<=0,])){w0[w0$over_represented_pvalue<=0,]$over_represented_pvalue=min(w0[w0$over_represented_pvalue>0,]$over_represented_pvalue)}

# We got the full association with p-values and DE genes in gomerge!

gomerge=w0

# Retriving the go sizes

gos=unique(as.data.frame(cbind(gomerge$category,gomerge$numInCat)))
colnames(gos)=c('category','gosize')
gos$gosize=as.numeric(gos$gosize)

# Establish the goarea (see Documentation), remove totally uninteresting terms

allgoarea=as.data.frame(aggregate(-log10(over_represented_pvalue)*step~category+term+ontology,data=gomerge,mean));
colnames(allgoarea)=c('name','description','ontology','area');
allgoarea=allgoarea[allgoarea$area>=0,];
allgoarea$area=sqrt(allgoarea$area)

# Selecting the most relevant terms according to type

gobparea=allgoarea[allgoarea$ontology=='BP',];  
gobparea=head(gobparea[order(gobparea$area,decreasing=TRUE),],nbpcat)
goccarea=allgoarea[allgoarea$ontology=='CC',];  
goccarea=head(goccarea[order(goccarea$area,decreasing=TRUE),],ncccat)
gomfarea=allgoarea[allgoarea$ontology=='MF',];  
gomfarea=head(gomfarea[order(gomfarea$area,decreasing=TRUE),],nmfcat)
keggarea=allgoarea[allgoarea$ontology=='KEGG',];
keggarea=head(keggarea[order(keggarea$area,decreasing=TRUE),],nggcat)
unotarea=head(allgoarea[allgoarea$ontology=='UNKN',],nuncat)

# Building the core of the landscape heatmaps

goarea=rbind(gobparea,goccarea,gomfarea,keggarea,unotarea)
goarea=goarea[order(goarea$area,decreasing=TRUE),];
goredu=merge(goarea,gos,by.x='name',by.y='category',sort=FALSE)
colnames(goredu)=c('name','description','ontology','area','gosize')
goredu=goredu[goredu$area!=0,]
gocore=xtabs(over_represented_pvalue~category+step,data=gomerge) #,sparse=TRUE)
gocore=gocore[goredu$name,]
if(length(gocore[gocore==0])){gocore[gocore==0]=0.1*min(gocore[gocore>0])}
gocore=-log10(gocore)

# Cleaning terms with no relevance

gocore=gocore[,colSums(gocore)>0]
steplist=colnames(gocore)

# Filling the GO landscape cells info (n. of genes in/out the GO cat.)

infogoland=data.frame()
infointern=data.frame()
for(step in steplist)
{
    selec=workdata[workdata$ran>=as.numeric(step),]
    ndgst=nrow(selec)
    n1=unique(merge(genestogo,selec,by.x='gene',by.y='lab'))
    n2=aggregate(gene~category,data=n1,length)
    n3=merge(goarea,n2,by.x='name',by.y='category')
    n4=as.data.frame(cbind(n3$name,step,paste(step,' (',ndgst,')',sep=''),paste(n3$gene,'/',ndgst-n3$gene,sep='')))
    n5=as.data.frame(cbind(n3$name,step,ndgst,n3$gene))
    colnames(n4)=c('category','step','step.ndeg','ndegin.ndegout')
    infogoland=as.data.frame(rbind(n4,infogoland))
    colnames(n5)=c('category','step','ndeg','ndegin')
    infointern=as.data.frame(rbind(n5,infointern))
}
infogoland$step=as.numeric(infogoland$step)
infointern$step=as.numeric(infointern$step) 
infointern$ndeg=as.numeric(infointern$ndeg) 
infointern$ndegin=as.numeric(infointern$ndegin) 

goland=gocore
plotgo=gocore
plotgocells=goland
for(c in rownames(goland))
{
    for(step in steplist)
    {
        if(goland[c,step]==0)
        {
            plotgocells[c,step]=paste('0\n0/',unique(infointern[infointern$step==step,]$ndeg),sep='')
        }
        else
        {
            plotgocells[c,step]=paste(round(goland[c,step],2),'\n',
                infogoland[infogoland$step==step & infogoland$category==c,]$ndegin.ndegout,sep='')
        }
    }
}
plotgocells=plotgocells[rownames(plotgo),colnames(plotgo)]
colnames(plotgo)=unique(infogoland[order(infogoland$step),]$step.ndeg)
rownames(plotgo)=paste(
    goredu[goredu$name %in% rownames(goland), ]$ontology,': ',
    goredu[goredu$name %in% rownames(goland), ]$description,' (',
    goredu[goredu$name %in% rownames(goland), ]$gosize,')',
    sep='')
abar=as.matrix(goredu$area)

# Coloring of the cell: Red up, Blue down, Green both

if(workname=='down-regulated'){coloring=c('grey','lightblue','cornflowerblue');capt='down-reg.\n'}
if(workname=='up-regulated'){coloring=c('grey','pink','salmon'); capt='up-reg.\n'}
if(workname=='up-and-down-regulated'){coloring=c('grey','lightgreen','green3');capt=''}

# Plot the GO rank vs. GO size. To compare with GSEA.

if(plotgosrank)
{
pdf(paste('./easygo_support/gorank_gosize','.',filename,'.',workname,'.pdf',sep=''))
g2=merge(as.data.frame(rownames(goland)),goredu,by.x='rownames(goland)',by.y='name',sort=FALSE)
g3=as.data.frame(cbind(rownames(g2),g2$area,g2$gosize,g2$ontology))
nggcat=nrow(g3[g3$V4=='KEGG',])
nbpcat=nrow(g3[g3$V4=='BP',])
ncccat=nrow(g3[g3$V4=='CC',])
nmfcat=nrow(g3[g3$V4=='MF',])
g4=data.frame()
if(nggcat>0){colgg=cbind(g3[g3$V4=='KEGG',],'purple');colnames(colgg)[5]='V5';g4=rbind(colgg,g4)}
if(nbpcat>0){colbp=cbind(g3[g3$V4=='BP',],'red');     colnames(colbp)[5]='V5';g4=rbind(colbp,g4)}
if(ncccat>0){colcc=cbind(g3[g3$V4=='CC',],'green4');  colnames(colcc)[5]='V5';g4=rbind(colcc,g4)}
if(nmfcat>0){colmf=cbind(g3[g3$V4=='MF',],'blue');    colnames(colmf)[5]='V5';g4=rbind(colmf,g4)}
if(nuncat==1){ngos=nrow(goland)-1}else{ngos=nrow(goland)}
g4$V2=as.numeric(g4$V2)
g4$V3=as.numeric(g4$V3)
correl=round(cor(g4$V2,g4$V3),3)
fitr=lm(log10(V2)~log10(V3),data=g4)

plot(x=g4$V3,y=g4$V2,col=g4$V5,pch=22,bg=g4$V5,
	ylab='S (global significance of the GO term)',
	xlab='GO size (number of genes in the category) [log]',
	log='xy',cex=1,
	main=paste(filename,'\nGO rank vs. GO size\n',ngos,' terms, ',filelabel,' corr: ',correl,sep='')
	);abline(fitr$coefficients[1],fitr$coefficients[2],col='black')
legend('topleft', c('BP','CC','MF','KEGG'),pch=22,col=c('red','green4','blue','purple','grey'),
	pt.bg=c('red','green4','blue','purple'))

dev.off()
}

###########################
# Print the GO landscape  #
###########################

# Here we use ComplexHeatMaps

selegenestogo=genestogo[genestogo$category %in% rownames(goland),]
pvalselegenestogo=unique(merge(selegenestogo,workdata,by.x='gene',by.y='lab'))
genelandfull=xtabs(ran~category+gene,data=pvalselegenestogo)
genestoshow=head(unique(pvalselegenestogo[order(pvalselegenestogo$ran,decreasing=TRUE),]$gene),ngenes)
geneland=genelandfull[rownames(goland),colnames(genelandfull) %in% genestoshow]
rownames(geneland)=rownames(plotgo)
geneland=sqrt(geneland*goredu$area)
geneland=geneland[,order(apply(geneland,2,max))]

# Convert the genes ID

toconv=as.data.frame(colnames(geneland));colnames(toconv)='gene'
topass=merge(toconv,ense,by.x='gene',by.y='ensembl_id',all.x=TRUE,sort=FALSE)
if(length(topass[is.na(topass$gene_id),])>1){topass[is.na(topass$gene_id),]$gene_id=topass[is.na(topass$gene_id),]$gene}
tosymb=merge(topass,symb,by='gene_id',sort=FALSE,all.x=TRUE)
if(length(tosymb[is.na(tosymb$symbol),])>1){tosymb[is.na(tosymb$symbol),]$symbol=tosymb[is.na(tosymb$symbol),]$gene}
duplic=tosymb[!duplicated(tosymb$gene),]
duplic=duplic[match(colnames(geneland),duplic$gene),]
colnames(geneland)=duplic$symbol

# Printing the pdf files. Trim the description. Default 60 chars.

trimchar=60
for(j in 1:nrow(goredu))
{
	if(nchar(goredu$description[j])>trimchar)
    {goredu$description[j]=paste(strtrim(goredu$description[j],trimchar-3),'...',sep='')}
}
ncolors=100

cf=0.01388889
generowfont=8
legeleng=max(nchar(goredu$description))
wht0=legeleng*generowfont*cf/2
ht0=Heatmap(cbind(goredu$area),
    cluster_rows=FALSE,
    show_heatmap_legend=FALSE,
    show_column_names=FALSE,
    cell_fun=function(j,i,x,y,w,h,col){grid.text(gp=gpar(fontsize=8),
        as.data.frame(goredu$description)[i,j], x, y)},
    col=colorRampPalette(coloring)(ncolors),
    name='ht0', width=unit(wht0,'inch'),column_title ='Category Description',column_title_gp = gpar(fontsize = 8)
    )

abar=cbind(goredu$area)
colnames(abar)='global sGO'
whta=0.25
hta=Heatmap(abar, col = colorRampPalette(coloring)(ncolors),
	show_column_names=FALSE,
    cell_fun=function(j,i,x,y,w,h,col){grid.text(gp=gpar(fontsize=5),round(cbind(goredu$area)[i,j],2),x,y)},
    name='hta', width=unit(whta,'inch'),
    column_title ='S',column_title_gp = gpar(fontsize = 8),
    show_heatmap_legend=FALSE)


annb0=as.data.frame(rev(unique(infointern$ndeg)));colnames(annb0)='n.genes'
annb1=as.data.frame(rev(unique(infointern$step)));colnames(annb1)='DE.step'
annb0=cbind(annb0,annb1)
annb=columnAnnotation(annb0,height=unit(0.5,"cm"),width=unit(0.5,"cm"),
	col=list(
		n.genes=colorRamp2(c(min(annb0$n.genes), (min(annb0$n.genes)+max(annb0$n.genes))/2, max(annb0$n.genes)), c(coloring)),
		DE.step=colorRamp2(c(min(annb0$DE.step), (min(annb0$DE.step)+max(annb0$DE.step))/2, max(annb0$DE.step)), c(coloring))
		),
	annotation_name_gp= gpar(fontsize = 6,font=2),gap = unit(c(1,1), "mm"),
	show_annotation_name=TRUE,
	annotation_name_side='left',
	show_legend=FALSE)

wht1=ncol(plotgo)*0.25
ht1=Heatmap(cbind(plotgo),
    show_heatmap_legend=FALSE,
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_column_names=TRUE,
    show_row_names=FALSE,
    row_names_side =c('left'),
    top_annotation=annb,
    col=colorRampPalette(coloring)(ncolors),
    cell_fun = function(j, i, x, y, w, h, col){grid.text(gp=gpar(fontsize = 3),plotgocells[i, j], x, y)},
    name='ht1', width=unit(wht1,'inch'),
    column_title ='GO Landscape',column_title_gp = gpar(fontsize = 8)
    )
if(addselvar)
{
anfc0=as.data.frame(data[rev(match(duplic$gene,data$lab)),]$sel);colnames(anfc0)='FC'
anfc0$FC=2^(anfc0$FC)
anna0=as.data.frame(workdata[rev(match(duplic$gene,workdata$lab)),]$ran)
colnames(anna0)='sDE'
anna0=cbind(anfc0,anna0)
anna=columnAnnotation(anna0,height=unit(0.5,"cm"),width=unit(0.5,"cm"),
	col=list(
		FC= colorRamp2(c(min(anna0$FC), (min(anna0$FC)+max(anna0$FC))/2, max(anna0$FC) ),c(coloring)),
		sDE=colorRamp2(c(min(anna0$sDE),(min(anna0$sDE)+max(anna0$sDE))/2,max(anna0$sDE)),c(coloring))
		),
	show_legend=TRUE, annotation_name_gp= gpar(fontsize = 6,font=2),gap = unit(c(1,1), "mm"),
	show_annotation_name=TRUE,
	annotation_name_side='right'
	)
}else{
anna0=as.data.frame(workdata[rev(match(duplic$gene,workdata$lab)),]$ran)
colnames(anna0)='sDE'
anna=columnAnnotation(anna0,height=unit(0.5,"cm"),width=unit(0.5,"cm"),
	col=list(
		sDE=colorRamp2(c(min(anna0$sDE),(min(anna0$sDE)+max(anna0$sDE))/2,max(anna0$sDE)),c(coloring))
		),
	show_legend=TRUE, annotation_name_gp= gpar(fontsize = 6,font=2),gap = unit(c(1,1), "mm"),
	show_annotation_name=TRUE,
	annotation_name_side='right'
	)

}

###########################
# Print the genelandscape #
###########################

wht2=ncol(geneland)*0.25
ht2=Heatmap(cbind(geneland[,rev(colnames(geneland))]),
    show_heatmap_legend=FALSE,
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_column_names=TRUE,
    col=colorRampPalette(coloring)(ncolors),
    cell_fun=function(j,i,x,y,w,h,col){grid.text(gp=gpar(fontsize=5),round(cbind(geneland[,rev(colnames(geneland))])[i, j],2), x, y)},
    name='ht2', width=unit(wht2,'inch'),
    show_row_names=FALSE,
    top_annotation=anna,
    column_title=paste('Gene Landscape (',filelabel,')',sep=''),column_title_gp = gpar(fontsize = 6)
    )


ann=as.data.frame(cbind(goredu$area,goredu$ontology,paste(goredu$name,sep='')))
hz=Heatmap(as.numeric(ann$V1),
    col=colorRampPalette(coloring)(ncolors),
    cell_fun=function(j,i,x,y,w,h,col){grid.text(gp=gpar(fontsize=5),as.data.frame(ann$V3)[i, j], x, y)},
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_column_names=FALSE,
    show_row_names=FALSE,
    show_heatmap_legend=FALSE,
    name='ontology',
    width=unit(1.5, "cm"),
    column_title ='GO ID',column_title_gp = gpar(fontsize = 8)
    )

ontology=as.data.frame(ann$V2)
colnames(ontology)='type'
ha=rowAnnotation(df=ontology,col=list(type=c("BP"="red","MF"="blue","CC"="green4","KEGG"="purple","UNKN"="grey")),
    width=unit(0.25,"cm"),
    show_legend=TRUE,
    annotation_legend_param=list(title='',fontsize = 8),
    name='ontology')

hs=Heatmap(cbind(goredu$gosize),col=colorRampPalette(c('grey','darkgrey'))(ncolors),
	cell_fun=function(j,i,x,y,w,h,col){grid.text(gp=gpar(fontsize=5),as.data.frame(goredu$gosize)[i, j], x, y)},
    width=unit(0.5,"cm"),
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_column_names=FALSE,
    show_row_names=FALSE,
    show_heatmap_legend=FALSE,
    column_title ='GO size',column_title_gp = gpar(fontsize = 8),
    name='gosize')

pdf(paste('./easygo_results/GOlandscape.',filename,'.',workname,'.complete.pdf', sep=''),
    width=unit(wht0+wht1+whta+wht2+2+1,'inch'));
draw(ht1+hta+ha+hs+hz+ht0+ht2,
    heatmap_legend_side = "left",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title=paste('GOlandscape results of "',filename,'" : ',filelabel,sep='')
    );
dev.off()

pdf(paste('./easygo_results/GOlandscape.',filename,'.',workname,'.Gene_Landscape.pdf', sep=''),
	width=unit(wht0+wht2+whta+2,'inch'));
    draw(ha+ht0+hta+ht2,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),heatmap_legend_side = "left",
    column_title=paste('GOlandscape results of "',filename,'"',sep='')
    );
dev.off()

pdf(paste('./easygo_results/GOlandscape.',filename,'.',workname,'.GO_Landscape.pdf', sep=''),
    width=unit(wht1+wht0+whta+3,'inch'));
draw(ht1+hta+hs+hz+ht0+ha,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title=paste('GOlandscape results of "',filename,'"',sep='')
    );
dev.off()

} # fi of printmap

if(length(steplist)<=2){print(paste('No clustering possible!'))}

} # cycle on workdata (up/down/both regulation data)

} # cycle on filenames

# End
