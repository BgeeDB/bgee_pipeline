#Marta Rosikiewicz 10/07/2012

###Invoking:
#R CMD BATCH --no-save --no-restore '--args celfile_path="celfile_path" chipType_path="chipType_path" output_path="output_path"' get_celfile_info.R Rout_path


###Arguments to provide:

#Rout_path="path to  .Rout file"

#celfile_path="path to a CEL file"
#chipType_info_path="path to file chipType_info.txt"
#output_file="path to an output file"


cmd_args = commandArgs(TRUE);

print (cmd_args)
if(length(cmd_args)==0){stop("no arguments provided")
} else {
  for(i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}


suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(methods))

###############################

get_qc_param<-function(data){

	qc_vector<-vector(length=2)

	#Affymetrix mas5
	m5call<-exprs(mas5calls(data,verbose=FALSE))
	qc_vector[1]<-as.character(round(table(m5call)["P"]/length(m5call)*100,digits=2))

	#arIQR
	nr_probesets<-length(unique(probeNames(data)))

	probe_rank<-split(rank(pm(data)),probeNames(data))
	probeset_average_rank<-sapply(probe_rank,mean)
	qc_vector[2]<-as.character(round(IQR(probeset_average_rank),digits=2))

	#names(qc_vector)<-c("percent_present","IQR_rank_tranform")
	return(qc_vector)
}

################################

#constructing matrix celfile_info that will contain all results:
celfile_info<-vector(length=11)
names(celfile_info)<-c("full_path","array_name","folder_name", "cdfname", "scan_date", "chipTypeId","percent_present","IQR_rank_tranform","unique_ID","error_type","error")


#obtaining the name of the directory from the full path to directory:
folder_name<-rev(unlist(strsplit(celfile_path,"/")))[2]
celfile_name<-rev(unlist(strsplit(celfile_path,"/")))[1]

celfile_info["full_path"]<-celfile_path
celfile_info["array_name"]<-celfile_name
celfile_info["folder_name"]<-folder_name

#reading in chipType
chipType_info<-as.matrix(read.table(chipType_info_path,sep="\t", header=T,colClasses="character",comment.char="",as.is=T,check.names=F,strip.white=T))
all_cdfname<-chipType_info[,"#cdfName"]
compatible_cdfname<-chipType_info[chipType_info[,"ChipTypeStatus"]=="compatible",1]

#reading in cel file inside try function (which prevent stoping execution of the script in case of error):
try_output<-try(cel_data<-ReadAffy(filenames=celfile_path))

#if the error occured during reading in cel file try_output contain error message
if(class(try_output)=="try-error"){
		celfile_error<-try_output[1]
		celfile_error<-as.character(celfile_error)
		celfile_error<-chartr("\n"," ",celfile_error)
        celfile_info["error_type"]<-"reading error"
		celfile_info["error"]<-celfile_error
        write.table(t(celfile_info),file=output_file,quote=F,row.names=F,sep="\t")
        quit("no")
}else{
	#obtaining the cdf name and scantime from succesfully read in cel file:

		celfile_info["cdfname"]<-cdfName(cel_data)
		celfile_info["scan_date"]<-protocolData(cel_data)$ScanDate

}

#checking if chipType is already known
if (!(celfile_info["cdfname"]%in%all_cdfname)){
        celfile_info["error_type"]<-"unknown cdfname"
        write.table(t(celfile_info),file=output_file,quote=F,row.names=F,sep="\t")
        quit("no")
}else{
		celfile_info["chipTypeId"]<-chipType_info[chipType_info[,"#cdfName"]==celfile_info["cdfname"],2]
}

#checking if chipType is compatible with the pipeline
if (!(celfile_info["cdfname"]%in%compatible_cdfname)){
        celfile_info["error_type"]<-"incompatible chipType"
        write.table(t(celfile_info),file=output_file,quote=F,row.names=F,sep="\t")
        quit("no")
}

#computing quality scores
try(celfile_info[c("percent_present","IQR_rank_tranform")]<-get_qc_param(cel_data))

#catching errors occured during computing quality scores
if(class(try_output)=="try-error"){
		celfile_error<-try_output[1]
		celfile_error<-as.character(celfile_error)
		celfile_error<-chartr("\n"," ",celfile_error)
        celfile_info["error_type"]<-"processing error"
		celfile_info["error"]<-celfile_error
        write.table(t(celfile_info),file=output_file,quote=F,row.names=F,sep="\t")
        quit("no")

}

#computing unique id
celfile_info["unique_ID"]<-paste(celfile_info["percent_present"],celfile_info["IQR_rank_tranform"],sep="_")

write.table(t(celfile_info),file=output_file,quote=F,row.names=F,sep="\t")
quit(save="no")

