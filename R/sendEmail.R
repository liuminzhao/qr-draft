sendEmail <- function(subject, text, address, file = FALSE, filename = "MyRFile"){
  if(file == FALSE){
    sys.arg <- paste("echo '", text, "' | mail -s '", subject,  "' ", address, sep = "")
  } else{
    ##		filename <- cat(filename, "." , unlist(strsplit(text, split="\\."))[2], sep = "")
    ##		sys.arg <- paste("uuencode ", text, " ", filename, "| mail -s ", subject, " ", address, sep = "")
    sys.arg <- paste("mail -s ", subject, " ", address, " < ", text, sep = "")
  }

  system(sys.arg)
}
