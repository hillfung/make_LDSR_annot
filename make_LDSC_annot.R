##########################################################################################################
##																										##
##			R function to make annotations for analysis using Stratified LD Score Regression			##
##			Version 2.4.6 (last modified on 05-04-2018)													##
##																										##
##########################################################################################################

library(data.table)

make.LDSR.annot<-function(input,out,annot.name,template.dir,input.type,sep=NULL,add.windows=T,GeneCode=NULL,ID.extension=T,gene.shoulder=NULL){
	begin.time<-Sys.time()

	log.file<-file(description=paste0(out,annot.name,".log"),open="wt")
	#error.file<-file(description=paste0(out,annot.name,".error"),open="wt")
	sink(file=log.file,append=F,type="output",split=T)
	sink(file=log.file,append=T,type="message")
	
	cat(rep("#",64),"\n##",rep(" ",60),"##\n##",rep(" ",19),"make.LDSR.annot v2.4.6",rep(" ",19),"##\n##",rep(" ",60),"##\n",rep("#",64),"\n\nAnalysis started at ",format(x=begin.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")
	
	if(missing(annot.name)){
		stop("Did not specify name for the new annotation\n")
	}
	cat("The new annotation will be called: ",annot.name,"\n",sep="")
	
	if(missing(out)){
		cat("WARNING: did not specify output directory. Files will be saved in the working directory\n")
		out<-paste0(getwd(),"/")
	}else{
		if(substr(x=out,start=nchar(out),stop=nchar(out))!="/"){
			out<-paste0(out,"/")
		}
		cat("Output will be saved at: ",out,"\n",sep="")
	}
	out<-paste0(out,annot.name)
	
	if(missing(template.dir)){
		cat("WARNING: Template directory was not specified. Attempting to look for 22 .annot[.gz]-files in the working directory...")
		template.dir<-paste0(getwd(),"/*")
	}else{
		cat("Looking for ",template.dir,".[1-22].annot[.gz]...",sep="")
	}
	template.list<-Sys.glob(paste0(template.dir,".",1:22,".annot"))
	if(length(template.list)==22){
		cat("done\n")
	}else{
		template.list<-Sys.glob(paste0(template.dir,".",1:22,".annot.gz"))
		if(length(template.list)==22){
			cat("done\n")
		}else{
			cat("ERROR\n")
			stop(paste0("Could not find ",template.dir,".[1-22].annot[.gz]\n",sep=""))
		}
	}
	
	if(is.numeric(add.windows) | (sum(is.na(suppressWarnings(as.numeric(add.windows))))==0 & is.logical(add.windows)==F)){
		add.windows<-as.numeric(add.windows)
		if(length(add.windows)==1){
			cat("A ",add.windows,"BP window was added to the new annotation\n",sep="")
		}else{
			cat("A ",add.windows[1],"BP",sep="")
			if(length(add.windows)>2){
				for(i in 2:(length(add.windows)-1)){
					cat(", ",add.windows[i],"BP",sep="")
				}
				cat(", and ",add.windows[length(add.windows)],"BP window were added to the new annotation\n",sep="")
			}else{
				cat(" and ",add.windows[2],"BP window were added to the new annotation\n",sep="")
			}
		}
	}else if(length(add.windows)>1){
		stop(cat("). Flag accepts: TRUE/T, FALSE/F, and numerical inputs\n",cat(add.windows,cat("Did not recognize add.windows=c(",sep=""),sep=", "),sep=""))
	}else if(add.windows==T){
		cat("A 100BP and 500BP window were added to the new annotation\n")
		add.windows<-c(100,500)
	}else if(add.windows==F){
		cat("No windows were added to new annotation (are you sure about this?)\n")
	}else{
		stop(cat("Did not recognize add.windows=",add.windows,". Flag accepts: TRUE/T, FALSE/F, and numerical inputs\n",sep=""))
	}
	
	if(missing(input.type)){
		stop("Did not specify the input type\n")
	}
	input.type<-match.arg(arg=input.type,choices=c("rs","id","name","bp","pos"))
	if(input.type=="rs"){
		rs.LDSR.function(input=input,out=out,annot.name=annot.name,template.list=template.list,add.windows=add.windows)
	}else if(input.type=="pos"){
		pos.LDSR.function(input=input,out=out,annot.name=annot.name,template.list=template.list,sep=sep,add.windows=add.windows)
	}else if(input.type%in%c("id","name")){
		gene.LDSR.function(input=input,out=out,annot.name=annot.name,template.list=template.list,input.type=input.type,add.windows=add.windows,GeneCode=GeneCode,ID.extension=ID.extension,gene.shoulder=gene.shoulder)
	}else{
		bp.LDSR.function(input=input,out=out,annot.name=annot.name,template.list=template.list,add.windows=add.windows,gene.shoulder=gene.shoulder)
	}
	
	end.time<-Sys.time()
	cat("Analysis ended at ",format(x=end.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")
	analysis.duration<-difftime(time1=end.time,time2=begin.time,units="secs")
	analysis.mins<-floor(x=floor(x=(analysis.duration/60)))
	analysis.secs<-analysis.duration-analysis.mins*60
	cat("Analysis took ",analysis.mins," minute(s) and ",analysis.secs," seconds\n",sep="")
	closeAllConnections()
}

rs.LDSR.function<-function(input,out,annot.name,template.list,add.windows){
	cat("User indicated a list of RSIDs was supplied as input\n")
	if(class(input)=="list"){
		cat("Input had class \"list\"\n")
		if(length(input)>1){
			cat("WARNING: input had length>1. The first element was used as input\n")
		}else if(length(input)==0){
			stop("Input had length 0\n")
		}
		input<-input[[1]]
	}
	
	if(class(input)%in%c("data.frame","matrix")){
		cat("Input had class \"data.frame\" or \"matrix\"\n")
		if(ncol(input)==1){
			cat("Input contained only one column. Proceeded with this as input\n")
			input<-input[,1]
		}else{
			cat("WARNING: Input contained multiple columns. Attempted to find an single column of RSIDs based on the column names\n")
			colnames(input)<-toupper(colnames(input))
			if(sum(colnames(input)%in%c("SNP","RS","RSID","ID","RS_ID","RSNUM","NUM","RS_NUM","RSNUMBER","NUMBER","RS_NUMBER"))==1){
				cat("The name of column number ",which(colnames(input)%in%c("SNP","RS","RSID","ID","RS_ID","RSNUM","NUM","RS_NUM","RSNUMBER","NUMBER","RS_NUMBER"))," uniquely matched the search criteria. Proceeded with this as input\n",sep="")
				input<-input[,colnames(input)%in%c("SNP","RS","RSID","ID","RS_ID","RSNUM","NUM","RS_NUM","RSNUMBER","NUMBER","RS_NUMBER")]
			}else{
				if(sum(colnames(input)%in%c("SNP","RS","RSID","ID","RS_ID","RSNUM","NUM","RS_NUM","RSNUMBER","NUMBER","RS_NUMBER"))>1){
					cat("Multiple column names matched the search criteria. ")
				}else{
					cat("None of the column names matched the search criteria. ")
				}
				cat("Scanning through each column, trying to find an \"rs\"-pattern\n")
				columns.found<-c()
				for(i in 1:ncol(input)){
					dum<-input[,i]
					if(is.character(dum)){
						dum<-tolower(dum)
					}
					if(length(grep(pattern="^rs",x=dum))>=1){
						columns.found<-c(columns.found,i)
					}
					rm(dum)
				}
				if(length(columns.found)==1){
					cat("Found one column with entries appearing to be RSIDs. Proceeded with this as input\n")
					input<-input[,c(columns.found)]
				}else{
					stop("Could not identify a single column of RSIDs. Please supply a vector as input to avoid problems\n")
				}
			}
		}
	}else if(is.vector(input)){
		cat("A vector was supplied as input\n")
	}
	input<-tolower(input)
	if(length(grep(pattern="^rs",x=input))<1){
		stop("Input did not appear to contain RSIDs\n")
	}
	input<-unique(input)
	cat("After removing duplicate entries, ",length(input)," entries remained\n",sep="")
	
	m<-0
	included<-NULL
	for(i in template.list){
		template<-prepareTemplate(i)
		chr<-template[1,]$CHR
		template$NEW.ANNOT<-0
		if(nrow(template[template$SNP%in%input,])==0){
			cat("WARNING: No SNPs were found on chromosome ",chr,"\n",sep="")
		}else{
			template[template$SNP%in%input,]$NEW.ANNOT<-1
			cat(sum(template$NEW.ANNOT)," SNPs were found on chromosome ",chr,"\n",sep="")
			m<-m+sum(template$NEW.ANNOT)
			included<-c(included,template[template$NEW.ANNOT==1,]$SNP)
		}
		if(is.numeric(add.windows)){
			for(j in add.windows){
				template<-make.windows(x=template,annot.name=annot.name,window.size=j)
			}
		}
		colnames(template)[colnames(template)=="NEW.ANNOT"]<-annot.name
		write.table(x=template,file=gzfile(paste0(out,".",chr,".annot.gz")),quote=F,sep="\t",row.names=F)
		rm(template,chr)
	}
	if(m==0){
		stop("No SNPs were included in the new annotation\n")
	}
	cat("In total, ",m,"/",length(input)," SNPs were included in the new annotation. These are saved in ",out,".included\n",sep="")
	write.table(x=included,file=paste0(out,".included"),quote=F,sep="\t",row.names=F,col.names=F)
	excluded<-input[(input%in%included)==F]
	if(length(excluded)!=0){
		cat(length(excluded)," entries were excluded from the new annotation. These are saved in ",out,".excluded\n",sep="")
		write.table(x=excluded,file=paste0(out,".excluded"),quote=F,sep="\t",row.names=F,col.names=F)
	}
}

pos.LDSR.function<-function(input,out,annot.name,template.list,sep,add.windows){
	cat("User indicated a list of chromosome/basepair-positions was supplied as list of input\n")
	if(class(input)=="list"){
		cat("Input had class \"list\"\n")
		if(length(input)>1){
			cat("WARNING: input had length>1. The first element was used as input\n")
		}else if(length(input)==0){
			stop("Input had length 0\n")
		}
		input<-input[[1]]
	}
	
	if(class(input)%in%c("data.frame","matrix")){
		cat("Input had class \"",class(input),"\"\n",sep="")
		if(ncol(input)==1){
			cat("Input contained only one column. Attempted to split the entries into two columns: one containing the chromosome and one with base-pair-positions\n")
			input<-pos.vectorToInput(x=input[,1],sep=sep)
		}else{
			cat("Input contained ",ncol(input)," columns. An attempt was made to find a column named \"CHR\" and one named \"BP\"\n",sep="")
			if(sum(c("CHR","BP")%in%colnames(input))>=2){
				input<-input[,c("CHR","BP")]
			}else{
				stop("Could not find the two necessary columns. It's easiest to provide a vector and the \"sep\" flag\n")
			}
		}
	}else if(is.vector(input)){
		cat("Input was a vector. Attempted to split the entries into two columns: one containing the chromosome and one with base-pair-positions\n")
		input<-pos.vectorToInput(x=input,sep=sep)
	}
	input.2<-input<-cbind(data.frame(ROW.NUMBER=1:nrow(input)),input,stringsAsFactors=F)
	input.2$CHR<-checkChromosomeColumn(input.2$CHR)
	input.2$BP<-sapply(X=input.2$BP,FUN=function(x){
		if(sum(is.na(suppressWarnings(as.numeric(x))))==0){
			return(as.numeric(x))
		}else{
			return(NA)
		}
	})
	input.excluded<-input.2[(input.2$CHR%in%c(1:22))==F | is.na(input$BP),]
	if(nrow(input.excluded)!=0){
		cat(nrow(input.excluded)," entries contained non-numerical or otherwise strange values. These were removed and saved in ",out,".excluded\n",sep="")
		input.excluded.2<-input[input$ROW.NUMBER%in%input.excluded$ROW.NUMBER,]
		input.excluded$CHR<-input.excluded.2[match(input.excluded$ROW.NUMBER,input.excluded.2$ROW.NUMBER),]$CHR
		input.excluded$BP<-input.excluded.2[match(input.excluded$ROW.NUMBER,input.excluded.2$ROW.NUMBER),]$BP
		write.table(x=input.excluded,file=paste0(out,".excluded"),quote=F,sep="\t",row.names=F)
		input.2<-input.2[input.2$CHR%in%c(1:22) & is.na(input.2$BP)==F,]
		rm(input.excluded,input.excluded.2)
		gc()
	}
	input.2<-input.2[,c("CHR","BP")]
	m<-0
	included<-NULL
	excluded<-NULL
	for(i in template.list){
		template<-prepareTemplate(i)
		chr<-template[1,]$CHR
		template$NEW.ANNOT<-0
		dum<-input.2[input.2$CHR==chr,]
		if(nrow(dum)==0){
			cat("WARNING: No SNPs were included from the input on chromosome ",chr,"\n",sep="")
		}else if(sum(template$BP%in%dum$BP)==0){
			cat("WARNING: ",nrow(dum)," SNPs were found for chromosome ",chr,", however, none of the positions matched with the template. Make sure the basepair positions in the input matches with the template build\n",sep="")
		}else{
			template[template$BP%in%dum$BP,]$NEW.ANNOT<-1
			cat(sum(template$NEW.ANNOT)," SNPs were found on chromosome ",chr,"\n",sep="")
			m<-m+sum(template$NEW.ANNOT)
			included<-c(included,template[template$NEW.ANNOT==1,]$SNP)
			excluded<-c(excluded,rownames(dum[(dum$BP%in%template$BP)==F,]))
		}
		if(is.numeric(add.windows)){
			for(j in add.windows){
				template<-make.windows(x=template,annot.name=annot.name,window.size=j)
			}
		}
		colnames(template)[colnames(template)=="NEW.ANNOT"]<-annot.name
		write.table(x=template,file=gzfile(paste0(out,".",chr,".annot.gz")),quote=F,sep="\t",row.names=F)
		rm(template,chr)
	}
	if(m==0){
		stop("No SNPs were included in the new annotation\n")
	}
	cat("In total, ",m,"/",nrow(input.2)," SNPs were included in the new annotation. These are saved in ",out,".included\n",sep="")
	write.table(x=included,file=paste0(out,".included"),quote=F,sep="\t",row.names=F,col.names=F)
	if(length(excluded)!=0){
		cat(length(excluded)," entries were excluded from the new annotation. These are saved in ",out,".excluded\n",sep="")
		write.table(x=excluded,file=paste0(out,".excluded"),quote=F,sep="\t",row.names=F,col.names=F)
	}
}

pos.vectorToInput<-function(x,sep){
	if(is.null(sep)){
		cat("WARNING: Did not specify the separator. An attempt was made to split the input using \"_\" or \":\" as separator\n")
		if(length(grep(pattern="_",x=x))==length(x)){
			sep<-"_"
		}else if(length(grep(pattern=":",x=x))==length(x)){
			sep<-":"
		}else{
			stop("Could not find \"CHR:BP\" or \"CHR_BP\" format in the input\n")
		}
	}
	cat("Using \"",sep,"\" as separator\n",sep="")
	x<-as.data.frame(t(sapply(X=x,FUN=splitInput,sep)))
	colnames(x)<-c("CHR","BP")
	if(sum(is.na(x$BP))==nrow(x)){
		stop("None of the entries were")
	}
	if(sum(is.na(x$BP))!=0){
		cat("WARNING: some entries were not split\n")
	}
	return(x)
}

splitInput<-function(x,sep){
	out<-unlist(strsplit(x=x,split=sep))
	if(length(out)>2){
		stop("(Some) entries were split into more than two elements. Please make sure you have specified the separator correctly\n")
	}
	if(sum(is.na(suppressWarnings(as.numeric(out))))==0){
		return(as.numeric(out))
	}else{
		return(c(x,NA))
	}
}

gene.LDSR.function<-function(input,out,annot.name,template.list,input.type,add.windows,GeneCode,ID.extension,gene.shoulder){
	if(input.type=="name"){
		merge.by<-"GENE_NAME"
		cat("User indicated a list of gene names was supplied as input\n")
		if(ID.extension==F){
			cat("WARNING: ID.extension==F does not work with input.type==\"name\" and will be ignored\n")
		}
	}else{
		merge.by<-"GENE_ID"
		cat("User indicated a list of gene IDs was supplied as input\n")
		if(ID.extension==F){
			cat("The gene IDs do not contain extensions. These will be removed from the GeneCode-file\n")
		}
	}
	input.type<-toupper(input.type)
	if(class(input)=="list"){
		cat("Input had class \"list\"\n")
		if(length(input)>1){
			cat("WARNING: input had length>1. The first element was used as input\n")
		}else if(length(input)==0){
			stop("Input had length 0\n")
		}
		input<-input[[1]]
	}
	if(class(input)%in%c("data.frame","matrix")){
		cat("Input had class \"",class(input),"\"\n",sep="")
		if(ncol(input)==1){
			cat("Input contained only one column. Proceeded with this as input\n")
			input<-input[,1]
		}else{
			cat("WARNING: Input contained multiple columns. Attempted to find a single column of gene ",input.type,"s based on the column names\n",sep="")
			colnames(input)<-toupper(colnames(input))
			if(sum(colnames(input)%in%c("GENE","GENEID","ID","GENE_ID","GENE.ID","GENENAME","NAME","GENE_NAME","GENE.NAME"))==1){
				cat("the name of column number ",which(colnames(input)%in%c("GENE","GENEID","ID","GENE_ID","GENE.ID","GENENAME","NAME","GENE_NAME","GENE.NAME"))," uniquely matched the search criteria. Proceeded with this as input\n")
				input<-input[,colnames(input)%in%c("GENE","GENEID","ID","GENE_ID","GENE.ID","GENENAME","NAME","GENE_NAME","GENE.NAME")]
			}else{
				stop(cat("Could not identify a single column of gene ",input.type,"s. Make sure you have specified the correct input type\n",sep=""))
			}
		}
	}else if(is.vector(input)){
		cat("A vector was supplied as input\n")
	}
	input<-unique(input)
	cat("Input contained ",length(input)," entries\n",sep="")

	cat("Checking the GeneCode-file\n")
	if(is.null(GeneCode)){
		stop("Did not provide GeneCode-file\n")
	}
	colnames(GeneCode)<-toupper(colnames(GeneCode))
	if(sum(c(merge.by,"CHROMOSOME","START","END")%in%colnames(GeneCode))!=4){
		stop("The GeneCode-file did not contain all necessary columns\n")
	}
	GeneCode<-GeneCode[,c(merge.by,"CHROMOSOME","START","END")]
	colnames(GeneCode)<-c(merge.by,"CHR","START","END")
	if(ID.extension==F & merge.by=="GENE_ID"){
		cat("Removing gene ID extensions from the input\n")
		input<-sapply(X=strsplit(x=input,split="\\."),FUN="[[",1)
		cat("Removing gene ID extensions from the GeneCode-file\n")
		GeneCode$GENE_ID<-sapply(X=strsplit(x=GeneCode$GENE_ID,split="\\."),FUN="[[",1)
	}
	GeneCode$CHR<-checkChromosomeColumn(x=GeneCode$CHR)
	GeneCode<-checkStartEnd(x=GeneCode,gene.shoulder=gene.shoulder)
	GeneCode<-GeneCode[GeneCode$CHR%in%c(1:22) & GeneCode$START!=-9999,]
	
	cat("Merging input and GeneCode-file\n")
	input.in<-input[input%in%GeneCode[,c(merge.by)]]
	if(length(input.in)==0){
		stop("None of the input genes were found in the GeneCode-file\n")
	}
	GeneCode.in<-GeneCode[GeneCode[,c(merge.by)]%in%input.in,]
	if(nrow(GeneCode.in)==0){
		stop("None of the input genes were found in the GeneCode file\n")
	}
	input.out<-input[(input%in%GeneCode.in[,c(merge.by)])==F]
	if(length(input.out)>0){
		cat(length(input.out),"/",length(input)," entries in the input-file could not be found in the GeneCode-file. These are saved in ",out,".excluded\n",sep="")
		write.table(x=input.out,file=paste0(out,".excluded"),quote=F,sep="\t",row.names=F,col.names=F)
	}
	cat(nrow(GeneCode.in)," entries were used to make the new annotation. These are saved in ",out,".included\n",sep="")
	write.table(x=GeneCode.in,file=paste0(out,".included"),quote=F,sep="\t",row.names=F)
	if(nrow(GeneCode.in)>length(input)){
		cat("WARNING: the number of entries used to make the new annotation was higher than the length of the input. This indicates that multiple entries in the ",merge.by," column matched the input\n",sep="")
	}
	
	templateToOut(x=GeneCode.in,out=out,annot.name=annot.name,template.list=template.list,add.windows=add.windows)
}

bp.LDSR.function<-function(input,out,annot.name,template.list,add.windows,gene.shoulder){
	cat("User indicated a list of base pair ranges was supplied as input\n")
	if(class(input)=="list"){
		cat("Input had class \"list\"\n")
		if(length(input)>1){
			cat("WARNING: list had length>1. The first element will be used as input\n")
		}else if(length(input)==0){
			stop("Input had length 0\n")
		}
		input<-input[[1]]
	}
	if((class(input)%in%c("data.frame","matrix"))==F){
		stop("Input must be a data frame or matrix\n")
	}
	colnames(input)<-toupper(colnames(input))
	if(sum(c("CHR","START","END")%in%colnames(input))!=3){
		stop("Input did not contain all necessary columns. Make sure they are named \"CHR\", \"START\", and \"END\"\n")
	}
	input<-input[,c("CHR","START","END")]
	input.2<-input<-cbind(data.frame(ROW.NUMBER=1:nrow(input)),input,stringsAsFactors=F)
	cat("Input contained ",nrow(input.2)," entries\n",sep="")
	input.2$CHR<-checkChromosomeColumn(x=input.2$CHR)
	input.2<-checkStartEnd(x=input.2,gene.shoulder=gene.shoulder)
	input.out<-input.2[(input.2$CHR%in%c(1:22))==F | input.2$START==-9999,]
	if(nrow(input.out)!=0){
		input.out[input.out$START==-9999,c("START","END")]<-input[input$ROW.NUMBER%in%input.out[input.out$START==-9999,]$ROW.NUMBER,c("START","END")]
		input.out.chr<-input[input$ROW.NUMBER%in%input.out$ROW.NUMBER,c("ROW.NUMBER","CHR")]
		input.out$CHR<-input.out.chr[match(input.out$ROW.NUMBER,input.out.chr$ROW.NUMBER),]$CHR
		input.out<-input.out[order(input.out$ROW.NUMBER),]
		cat(nrow(input.out)," entries contained non-numerical or otherwise strange data. These are saved in ",out, ".excluded\n",sep="")
		write.table(x=input.out,file=paste0(out,".excluded"),quote=F,sep="\t",row.names=F)
		input.2<-input.2[input.2$CHR%in%c(1:22) & input.2$START!=-9999,]
		rm(input.out,input.out.chr)
	}
	rm(input)
	input.2<-input.2[,c("CHR","START","END")]
	cat(nrow(input.2)," entries were used to make the new annotation. These are saved as ",out,".included\n",sep="")
	write.table(x=input.2,file=paste0(out,".included"),quote=F,sep="\t",row.names=F)
	
	templateToOut(x=input.2,out=out,annot.name=annot.name,template.list=template.list,add.windows=add.windows)
}

prepareTemplate<-function(x){
	if(substr(x=x,start=nchar(x)-2,stop=nchar(x))==".gz"){
		dum<-fread(input=paste("zcat",x),header=T,showProgress=F,data.table=F)
	}else if(substr(x=x,start=nchar(x)-5,stop=nchar(x))==".annot"){
		dum<-fread(input=x,header=T,showProgress=F,data.table=F)
	}else{
		stop("Did not recognize template file extension\n")
	}
	if(sum(c("CHR","BP","SNP","CM")%in%colnames(dum))!=4){
		stop("Template did not contain all necessary columns\n")
	}
	dum<-dum[,c("CHR","BP","SNP","CM")]
	return(dum)
}

make.windows<-function(x,annot.name,window.size){
	if(sum(x$NEW.ANNOT)==0){
		x$TEMP.WINDOW<-0
	}else{
		window.start<-window.end<-x[x$NEW.ANNOT==1,]$BP
		window.start<-ifelse(test=(window.start-window.size)<min(x$BP),yes=min(x$BP),no=window.start-window.size)
		window.end<-ifelse(test=(window.end+window.size)>max(x$BP),yes=max(x$BP),no=window.end+window.size)
		dum<-data.frame(START=window.start,END=window.end,stringsAsFactors=F)
		x$TEMP.WINDOW<-bpRangeToAnnot(x=x,bp.range=dum)
		rm(window.start,window.end,dum)
	}
	colnames(x)[colnames(x)=="TEMP.WINDOW"]<-paste0(annot.name,".extend.",window.size)
	return(x)
}

make.windows.bp<-function(x,bp.range,annot.name,window.size){
	x$dum<-0
	if(is.null(bp.range)==F){
		bp.range$START<-ifelse((bp.range$START-window.size)<min(x$BP),min(x$BP),bp.range$START-window.size)
		bp.range$END<-ifelse((bp.range$END+window.size)>max(x$BP),max(x$BP),bp.range$END+window.size)
		annot.range<-bp.range[1,c("START","END")]
		if(nrow(bp.range)>1){
			for(j in 2:nrow(bp.range)){
				if(bp.range[j,]$START<=annot.range[nrow(annot.range),]$END){
					annot.range[nrow(annot.range),]$END<-bp.range[j,]$END
				}else{
					annot.range<-rbind(annot.range,bp.range[j,c("START","END")])
				}
			}
		}
		annot.range$RANGE<-as.numeric(apply(X=annot.range,MARGIN=1,FUN=function(y){length(y[1]:y[2])}))
		annot.range.2<-rep(x=NA,times=sum(annot.range$RANGE))
		counter.end<-0
		for(i in 1:nrow(annot.range)){
			counter.start<-counter.end+1
			counter.end<-counter.start+annot.range[i,]$RANGE-1
			annot.range.2[counter.start:counter.end]<-annot.range[i,]$START:annot.range[i,]$END
		}
		if(nrow(x[x$BP%in%annot.range.2,])!=0){
			x[x$BP%in%annot.range.2,]$dum<-1	
		}
	}
	colnames(x)[colnames(x)=="dum"]<-paste0(annot.name,".extension.",window.size)
	return(x)
}

checkChromosomeColumn<-function(x){
	cat("Checking the chromosome column\n")
	if(length(table(x))>30){
		stop("The chromosome column contained more than 30 different inputs, which seems implausible. Please remove non-autosomal entries in the input and try again\n")
	}
	if(is.numeric(x)==F){
		cat("WARNING: the chromosome column contained non-numerical entries. Attempting to change the values...",sep="")
		x<-toupper(as.character(x))
		x<-gsub(pattern="CHR",replacement="",x=x)
		x<-sapply(x,checkChromosomeFunction)
		cat("done\n",sep="")
	}else if(is.integer(x)){
		x<-as.numeric(x)
	}
	return(x)
}

checkChromosomeFunction<-function(x){
	dum<-strsplit(x=x,split=NULL)[[1]]
	dum.out<-c()
	for(i in 1:length(dum)){
		if(is.na(suppressWarnings(as.numeric(dum[[i]])))==F){
			dum.out<-paste0(dum.out,dum[[i]])
		}
	}
	if(length(dum.out)!=0){
		return(as.numeric(dum.out))
	}else{
		return(x)
	}
}

checkStartEnd<-function(x,gene.shoulder){
	for(i in c("START","END")){
		cat("Checking the ",i," column...",sep="")
		x[,c(i)]<-sapply(X=x[,c(i)],FUN=function(x){
			if(sum(is.na(suppressWarnings(as.numeric(x))))==0){
				return(as.numeric(x))
			}else{
				return(-9999)
			}
		})
		cat("done\n")
	}
	if(nrow(x[x$START>=x$END,])==0){
		x$direction<-1
	}else{
		cat("WARNING: Some rows contained higher start positions than end positions. These values will be swapped\n")
		x$direction<-0
		x[x$START<=x$END,]$direction<-1
		x[x$START>x$END,]$direction<--1
		x[x$direction==-1,c("START","END")]<-t(apply(X=x[x$direction==-1,c("START","END")],MARGIN=1,FUN=function(x){return(x[c(2,1)])}))
	}
	x<-x[order(x$CHR,x$START,x$END),]
	if(is.null(gene.shoulder)==F){
		cat("Adding a ",format(gene.shoulder,scientific=F)," bp shoulder around the genes...",sep="")
		x$START<-x$START-gene.shoulder
		x$END<-x$END+gene.shoulder
		cat("done\n")
	}
	return(x[,colnames(x)!="direction"])
}

bpRangeToAnnot<-function(x,bp.range){
	x$dum<-0
	if(is.null(bp.range)==F){
		bp.range<-bp.range[,c("START","END")]
		bp.range$RANGE<-apply(X=bp.range,MARGIN=1,FUN=function(y){length(y[1]:y[2])})
		annot.range<-bp.range[1,]
		annot.range.2<-rep(NA,times=sum(bp.range$RANGE))
		annot.range.2[1:bp.range[1,]$RANGE]<-bp.range[1,]$START:bp.range[1,]$END
		if(nrow(bp.range)>1){
			counter.start<-1
			counter.end<-bp.range[1,]$RANGE
			for(i in 2:nrow(bp.range)){
				if(bp.range[i,]$START<=annot.range[nrow(annot.range),]$END){
					annot.range[nrow(annot.range),]$END<-bp.range[i,]$END
					dum<-annot.range[nrow(annot.range),]$START:annot.range[nrow(annot.range),]$END
					counter.end<-counter.start+length(dum)-1
					annot.range.2[counter.start:counter.end]<-dum
					rm(dum)
				}else{
					counter.start<-counter.end+1
					annot.range<-bp.range[i,]
					counter.end<-counter.start+annot.range[1,]$RANGE-1
					annot.range.2[counter.start:counter.end]<-annot.range[1,]$START:annot.range[1,]$END
				}
			}
		}
		annot.range.2<-annot.range.2[!is.na(annot.range.2)]
		if(length(annot.range.2)>1000000){
			seq.end<-seq.start<-seq(from=1,to=length(annot.range.2),by=1000000)
			seq.end[1:(length(seq.end)-1)]<-seq.end[1:(length(seq.end)-1)]+999999
			seq.end[length(seq.end)]<-length(annot.range.2)
			for(i in 1:length(seq.start)){
				if(nrow(x[x$BP%in%annot.range.2[seq.start[i]:seq.end[i]],])>0){
					x[x$BP%in%annot.range.2[seq.start[i]:seq.end[i]],]$dum<-1
				}
			}
		}else if(nrow(x[x$BP%in%annot.range.2,])!=0){
			x[x$BP%in%annot.range.2,]$dum<-1
		}
	}
	return(x$dum)
}

templateToOut<-function(x,out,annot.name,template.list,add.windows){
	for(template.dum in template.list){
		template<-prepareTemplate(x=template.dum)
		chr<-template[1,]$CHR
		if(nrow(x[x$CHR==chr,])==0){
			cat("No entries were found for chromosome ",chr,"\n",sep="")
			template$NEW.ANNOT<-0
			dum<-NULL
		}else{
			dum<-x[x$CHR==chr,]
			cat(nrow(dum)," entries were found for chromosome ",chr,"\n",sep="")
			template$NEW.ANNOT<-bpRangeToAnnot(x=template,bp.range=dum)
		}
		if(is.numeric(add.windows)){
			for(j in add.windows){
				template<-make.windows.bp(x=template,bp.range=dum,annot.name=annot.name,window.size=j)
			}
		}
		colnames(template)[colnames(template)=="NEW.ANNOT"]<-annot.name
		write.table(x=template,file=gzfile(paste0(out,".",chr,".annot.gz")),quote=F,sep="\t",row.names=F)
	}
}

collapse.annots<-function(annot.dir,out,annot.name,add.windows=T,collapse.type){
	begin.time<-Sys.time()
	
	log.file<-file(description=paste0(out,annot.name,".log"),open="wt")
	#error.file<-file(description=paste0(out,annot.name,".error"),open="wt")
	sink(file=log.file,append=F,type="output",split=T)
	sink(file=log.file,append=T,type="message")
	
	cat(rep("#",64),"\n##",rep(" ",60),"##\n##",rep(" ",19),"make.LDSR.annot v2.4.6",rep(" ",19),"##\n##",rep(" ",60),"##\n",rep("#",64),"\n\nAnalysis started at ",format(x=begin.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")
	
	if(missing(collapse.type)){
		stop("Did not specify how to collapse the annotations\nPlease enter \"union\" or \"intersection\" (written full and case-sensitive)")
	}
	collapse.type<-match.arg(arg=collapse.type,choices=c("union","intersection"))
	cat("Function was called to make the",collapse.type,"between the annotations\n")
	
	if(missing(annot.name)){
		stop("Did not specify name for the new annotation\n")
	}
	cat("The new annotation will be called:",annot.name,"\n")
	
	if(missing(annot.dir)){
		stop("No input was supplied via \"annot.dir\"\n",call.=F)
	}
	cat(length(annot.dir),"input(s) were supplied via \"annot.dir\"\n")
	inputs<-c()
	for(i in 1:length(annot.dir)){
		dum<-annot.dir[i]
		if(substr(x=dum,start=nchar(dum),stop=nchar(dum))=="/"){
			cat("WARNING: element",i,"of annot.dir appears to be a directory. Looking for .annot[.gz]-files in this directory...")
			dum.list<-paste0(dum,list.files(path=dum,pattern="*\\.*\\.anno*"))
			dum.list<-gsub("\\.gz$","",dum.list,perl=T)
			dum.list<-gsub("\\..\\.annot$","",dum.list,perl=T)
			dum.list<-gsub("\\...\\.annot$","",dum.list,perl=T)
			dum.list<-unique(dum.list)
			for(j in 1:length(dum.list)){
				if(length(Sys.glob(paste0(dum.list[j],".",1:22,".anno*")))!=22){
					stop(cat("Could not find ",dum.list[j],".[1-22].annot[.gz]\n",sep=""))
				}
			}
			inputs<-c(inputs,Sys.glob(paste0(dum.list,".1.anno*")))
			cat("done\n")
		}else{
			if(length(Sys.glob(paste0(dum,".",1:22,".anno*")))!=22){
				stop(cat("Could not find ",dum,".[1-22].annot[.gz]\n",sep=""))
			}
			inputs<-c(inputs,Sys.glob(paste0(dum,".1.anno*")))
		}
		rm(dum)
	}
	cat("Found",length(inputs),"annotations\n")
	inputs<-gsub(".gz$","",inputs,perl=T)
	inputs<-gsub(".1.annot$","",inputs,perl=T)
	if(length(inputs)!=length(unique(inputs))){
		inputs<-unique(inputs)
		cat("WARNING: the list of annotations contained duplicate entries. After removing them,",length(inputs),"annotations remain\n")
	}
	input.names<-sapply(strsplit(x=inputs,split="/"),tail,n=1)
	if(length(input.names)!=length(unique(input.names))){
		cat("The following inputs appeared to be the same annotation in different directories",inputs[input.names==input.names[duplicated(input.names)]],sep="\n")
		stop("Duplicated annotations were found. Please check \"annot.dir\" Function was aborted",call.=F)
	}
	if(length(inputs)==1){
		stop("Only one annotation was found",call.=F)
	}
	
	if(missing(out)){
		cat("WARNING: did not specify output directory. Files will be saved in the working directory\n")
		out<-paste0(getwd(),"/")
	}else{
		if(substr(x=out,start=nchar(out),stop=nchar(out))!="/"){
			out<-paste0(out,"/")
		}
		cat("Output will be saved at:",out,"\n")
	}
	out<-paste0(out,annot.name)
	
	if(is.numeric(add.windows) | (sum(is.na(suppressWarnings(as.numeric(add.windows))))==0 & is.logical(add.windows)==F)){
		add.windows<-as.numeric(add.windows)
		if(length(add.windows)==1){
			cat("A ",add.windows,"BP window was added to the new annotation\n",sep="")
		}else{
			cat("A ",add.windows[1],"BP",sep="")
			if(length(add.windows)>2){
				for(i in 2:(length(add.windows)-1)){
					cat(", ",add.windows[i],"BP",sep="")
				}
				cat(" and ",add.windows[length(add.windows)],"BP window were added to the new annotation\n",sep="")
			}else{
				cat(" and ",add.windows[2],"BP window were added to the new annotation",sep="")
			}
		}
	}else if(length(add.windows)>1){
		stop(cat("). Flag accepts: TRUE/T, FALSE/F, and numerical inputs\n",cat(add.windows,cat("Did not recognize add.windows=c(",sep=""),sep=", "),sep=""))
	}else if(add.windows==T){
		cat("A 100BP and 500BP window were added to the new annotation\n")
		add.windows<-c(100,500)
	}else if(add.windows==F){
		cat("No windows were added to the new annotation (are you sure about this?)\n")
	}else{
		stop(cat("Did not recognize add.windows=",add.windows,". Flag accepts: TRUE/T, FALSE/F, and numerical inputs\n",sep=""))
	}
	
	m<-0
	snp.list<-c()
	for(i in 1:22){
		dat<-readAnnot(Sys.glob(paste0(inputs[1],".",i,".anno*")))
		colnames(dat)[5]<-input.names[1]
		for(j in 2:length(inputs)){
			dum<-readAnnot(Sys.glob(paste0(inputs[j],".",i,".anno*")))
			if(identical(dat$SNP,dum$SNP)){
				dat<-cbind(dat,dum[,5])
				colnames(dat)[ncol(dat)]<-input.names[j]
			}else{
				#should never happen
				stop(paste(inputs[j],"did not contain the same SNPs as the previous input(s)\nIf this annotation was made using make.LDSR.annot, please contact Hill"))
			}
			rm(dum)
		}
		dat$NEW.ANNOT<-rowSums(dat[,5:ncol(dat)])
		if(sum(dat$NEW.ANNOT)==0){
			cat("WARNING: 0 SNPs were included on chromosome",i,"\n")
		}else{
			if(collapse.type=="union"){
				dat[dat$NEW.ANNOT>=1,]$NEW.ANNOT<-1
			}else{
				dat[dat$NEW.ANNOT!=length(inputs),]$NEW.ANNOT<-0
				dat[dat$NEW.ANNOT==length(inputs),]$NEW.ANNOT<-1
			}
			cat(nrow(dat[dat$NEW.ANNOT==1,]),"SNPs were included on chromosome",i,"\n")
			m<-m+nrow(dat[dat$NEW.ANNOT==1,])
			snp.list<-c(snp.list,dat[dat$NEW.ANNOT==1,]$SNP)
		}
		dat<-dat[,c(1:4,ncol(dat))]
		if(is.numeric(add.windows)){
			for(j in add.windows){
				dat<-make.windows(x=dat,annot.name=annot.name,window.size=j)
			}
		}
		colnames(dat)[5]<-annot.name
		write.table(x=dat,file=gzfile(paste0(out,".",i,".annot.gz")),quote=F,sep="\t",row.names=F)
	}
	cat(m,"SNPs were included in the annotation. A list can be found at",paste0(out,".included.gz\n"))
	write.table(x=snp.list,file=gzfile(paste0(out,".included.gz")),quote=F,sep="\t",row.names=F)
	end.time<-Sys.time()
	cat("Analysis ended at ",format(x=end.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")
	analysis.duration<-difftime(time1=end.time,time2=begin.time,units="secs")
	analysis.mins<-floor(x=floor(x=(analysis.duration/60)))
	analysis.secs<-analysis.duration-analysis.mins*60
	cat("Analysis took ",analysis.mins," minute(s) and ",analysis.secs," seconds\n",sep="")
	closeAllConnections()
}

readAnnot<-function(x){
	if(substr(x=x,start=nchar(x)-2,stop=nchar(x))==".gz"){
		dum<-fread(input=paste("zcat",x),header=T,showProgress=F,data.table=F)
	}else if(substr(x=x,start=nchar(x)-5,stop=nchar(x))==".annot"){
		dum<-fread(input=x,header=T,showProgress=F,data.table=F)
	}else{
		stop("Did not recognize template file extension\n")
	}
	if(sum(c("CHR","BP","SNP","CM")%in%colnames(dum))!=4){
		stop("Template did not contain all necessary columns\n")
	}
	if(!is.numeric(dum[,5])){
		stop(paste("The fifth column of",x,"must be numeric\nIf this annotation was made using make.LDSR.annot, please contact Hill"))
	}
	tmp<-unique(dum[,5])
	if(length(tmp)==1 & sum(tmp%in%c(0,1))!=1){
		stop(paste("The fifth column of",x,"must be a vector of 0's and 1's\nIf this annotation was made using make.LDSR.annot, please contact Hill"))
	}else if(length(tmp)==2 & sum(tmp%in%c(0,1))!=2){
		stop(paste("The fifth column of",x,"must be a vector of 0's and 1's\nIf this annotation was made using make.LDSR.annot, please contact Hill"))
	}else if(length(tmp)!=2 & length(tmp)!=1){
		stop(paste("The fifth column of",x,"must be a vector of 0's and 1's\nIf this annotation was made using make.LDSR.annot, please contact Hill"))
	}
	dum<-dum[,c("CHR","BP","SNP","CM",colnames(dum)[5])]
	return(dum)	
}