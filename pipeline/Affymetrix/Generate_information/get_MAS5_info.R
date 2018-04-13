#Marta Rosikiewicz 10/07/2012

###Invoking:
#R CMD BATCH --no-save --no-restore '--args MAS5_path="MAS5_path" output_path="output_path"' get_MAS5_info.R Rout_path


###Arguments to provide:

#Rout_path="path to  .Rout file"

#MAS5_path="path to a MAS5 file"
#output_file="path to an output file"

#reading in command arguments
cmd_args = commandArgs(TRUE);

print (cmd_args)

if(length(cmd_args)==0){stop("no arguments provided")
} else {
  for(i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

#constructing vector MAS5_info that will contain all results:
MAS5_info<-vector(length=7)
names(MAS5_info)<-c("full_path","array_name","folder_name", "percent_present","unique_ID","error_type","error")

#obtaining the name of the directory from the full path to directory:
folder_name<-rev(unlist(strsplit(MAS5_path,"/")))[2]
MAS5_name<-rev(unlist(strsplit(MAS5_path,"/")))[1]

MAS5_info["full_path"]<-MAS5_path
MAS5_info["array_name"]<-MAS5_name
MAS5_info["folder_name"]<-folder_name

#correspondance between different call names type
corresp_call<-c(
        'present'='present',
        'absent'='absent',
        'marginal'='marginal',
        'undefined'='undefined',
        'p'='present',
        'a'='absent',
        'm'='marginal',
        'u'='undefined',
        'Present'='present',
        'Absent'='absent',
        'Marginal'='marginal',
        'Undefined'='undefined',
        'P'='present',
        'A'='absent',
        'M'='marginal',
        'U'='undefined',
        'RP'='undefined'
        )

#the script assume that the MAS5 calls are in the 2 column of file
MAS5_call_col_nr=2

try_output<-try( MAS5_file<-as.matrix(read.delim(MAS5_path,sep="\t",header=T,comment.char="",colClasses="character")))


#if the error occured during reading in cel file try_output contain error message
if(class(try_output)=="try-error"){
		MAS5_error<-try_output[1]
		MAS5_error<-as.character(MAS5_error)
		MAS5_error<-chartr("\n"," ",MAS5_error)
		MAS5_info["error_type"]<-"reading error"
		MAS5_info["error"]<-MAS5_error
		write.table(t(MAS5_info),file=output_file,quote=F,row.names=F,sep="\t")
       quit("no")

}

nr_genes<-nrow(MAS5_file)

#computing the number of occurence for each call type
MAS5_table<-table(MAS5_file[,MAS5_call_col_nr])
MAS5_table_names<-names(MAS5_table)

#changing the call type names from file to the standard one
names(MAS5_table)<-corresp_call[MAS5_table_names]

#computing the MAS5 percent present
MAS5_info["percent_present"]<-as.character(round(MAS5_table["present"]/nr_genes*100,digits=2))

#computing the unique id
calls_names<-c("present",'absent','marginal','undefined')
MAS5_info["unique_ID"]<-paste(MAS5_table[calls_names],collapse="_")


write.table(t(MAS5_info),file=output_file,quote=F,row.names=F,sep="\t")
quit(save="no")





