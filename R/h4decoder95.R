#' generate_barcode
#'
#' It generates a list of barcodes of 9 bases, based on quaternary Hamming(9,5).
#'
#'@examples 
#'print("Generating barcodes...")
#'barcodes=generate_barcode()
#'print("done.")
#'print(paste("We just generated", length(barcodes), "barcodes."))
#'
#'@export

generate_barcode = structure(function ## It generates a list of barcodes. 
### It generates a list of barcodes of 9 bases, based on quaternary Hamming(9,5).
(){
barcodes = matrix(nrow=(4^5),ncol=9)
colnames(barcodes) = c("p1","p2","d1","p3","d2","d3","d4","p4","d5")
c = 1 
for (d1 in 0:3) {
  for (d2 in 0:3) {
    for (d3 in 0:3) {
      for (d4 in 0:3) {
        for (d5 in 0:3) {
          p1 = (4-((d1+d2+d4+d5)%%4))%%4
          p2 = (4-((d1+d3+d4)%%4))%%4
          p3 = (4-((d2+d3+d4)%%4))%%4
          p4 = (4-((d5)%%4))%%4
          barcodes[c,c(1:9)] = c(p1,p2,d1,p3,d2,d3,d4,p4,d5)
          c = c + 1
        }
      }
    }
  }
}
barcodes[barcodes == 0] = "A"
barcodes[barcodes == 1] = "C"
barcodes[barcodes == 2] = "G"
barcodes[barcodes == 3] = "T"
barcodeslist = ""
for (i in 1:9) {
  barcodeslist = paste(barcodeslist,barcodes[,i],sep="")
}
barcodeslist
# Returns a list of 9 bp barcodes.
}, ex=function(){
  print("Generating barcodes...")
  barcodes=generate_barcode()
  print("done.")
  print(paste("We just generated", length(barcodes), "barcodes."))
})
 
##############

#' filter_barcode
#'
#' It filters barcodes list according to GC contents. GC content calculation is based on Ringo package from bioconductor.
#'
#'@param barcodeslist A list of barcodes to be filtered
#'
#'@examples 
#'print("Filtering barcodes...")
#'barcodes = generate_barcode()
#'filtered_barcodes = filter_barcode(barcodes)
#'print("done.")
#'print(paste("Only",length(filtered_barcodes$barcodeslist), "barcodes remain. "))
#'
#'@importFrom Ringo "compute.gc"
#'
#'@export


filter_barcode = structure(function ## It filters barcodes list according to GC contents.
### GC content calculation is based on Ringo package from bioconductor.
(barcodeslist ##<< A list of barcodes to be filtered.
){
library(Ringo)
GC<-compute.gc(barcodeslist)
dfbarcodes<-data.frame(cbind(barcodeslist, GC), stringsAsFactors=FALSE)
dfbarcodes<-subset(dfbarcodes, GC > 0.4 & GC < 0.6)
return(dfbarcodes)

# It returns a dataframe containing a filtered barcode list with according GC contents (between 40% and 60% of GCs).
}, ex=function(){
  print("Filtering barcodes...")
  barcodes = generate_barcode()
  filtered_barcodes = filter_barcode(barcodes)
  print("done.")
  print(paste("Only",length(filtered_barcodes$barcodeslist), "barcodes remain. "))
})

############

#' compile_primer
#'
#' It generates full sequence of sequencing primers according to your design. It fuses adaptor, barcode and primer sequence of your choice.
#'
#'@param up_seq Primer 5' sequence up to the barcode.
#'@param down_seq Primer 3' sequence down to the barcode.
#'@param barcodes Dataframe of filtered barcodes.
#'
#'@examples 
#'print("Compiling barcodes...")
#'barcodes = generate_barcode()
#'filtered_barcodes = filter_barcode(barcodes)
#'compiled_primer = compile_primer("AAA", "CCC", filtered_barcodes[,1])
#'print("done.")
#'
#'@export

compile_primer = structure(function ## It generates full sequence of sequencing primers according to your design.
### It fuses adaptor, barcode and primer sequence of your choice.
(
up_seq, ##<< Primer 5' sequence up to the barcode. 
down_seq, ##<< Primer 3' sequence down to the barcode.
barcodes ##<< Dataframe of filtered barcodes. 
){
paste(up_seq, barcodes,down_seq,sep="")
# It returns a list of sequencing primers. 
}, ex=function(){
  print("Compiling barcodes...")
  barcodes = generate_barcode()
  filtered_barcodes = filter_barcode(barcodes)
  compiled_primer = compile_primer("AAA", "CCC", filtered_barcodes[,1])
  print("done.")
})

##############

#' correct_barcode
#'
#' It corrects barcodes obtained by sequencing. The capacity of correction of the algorithm is one mismatch maximum.
#'
#'@param index Vector of barcodes obtained from sequencing.
#'
#'@examples 
#'print("Correcting index...")
#'index=c("GACATTGGG","CCTCTACAA","CATATGTGG","ACTCTCGCT","GTTTTAGGA") 
#'# correct barcodes are GATATTGGG, CCTCTTCAA, CATATGTGG, ACTCTCTCT, GTTTTAGGG 
#'corrected_index = correct_barcode(index)
#'print("done.")
#'print(paste(corrected_index[[3]],"errors have been corrected within the",corrected_index[[2]], "index list provided"))
#'
#'@export


correct_barcode = structure(function ## It corrects barcodes obtained by sequencing.
### The capacity of correction of the algorithm is one mismatch maximum. 
(
index ##<< Vector of barcodes obtained from sequencing.
){
bin2dec <- function(x) {
x <- as.character(as.numeric(x))
b <- as.numeric(unlist(strsplit(x, "")))
pow <- 2 ^ ((length(b) - 1):0)
sum(pow[b == 1])
} 
cor_index=apply(t(index),2,function(index_cur){
	x = as.matrix(rbind(unlist(strsplit(index_cur,""))))
	x[x == "A"] <- 0
	x[x == "C"] <- 1
	x[x == "G"] <- 2
	x[x == "T"] <- 3
	p1 = (as.numeric(x[1])+as.numeric(x[3])+as.numeric(x[5])+as.numeric(x[7])+as.numeric(x[9]))%%4
	p2 = (as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[6])+as.numeric(x[7]))%%4
	p3 = (as.numeric(x[4])+as.numeric(x[5])+as.numeric(x[6])+as.numeric(x[7]))%%4
	p4 = (as.numeric(x[8])+as.numeric(x[9]))%%4
	p_max = max(p1,p2,p3,p4)
	if (p1>0) {p1=1};
	if (p2>0) {p2=1};
	if (p3>0) {p3=1};
	if (p4>0) {p4=1};
	bin_pos = paste(p4,p3,p2,p1, sep="")
	err_pos = bin2dec(bin_pos)
	corr_value = (as.numeric(x[err_pos])-p_max)%%4
	corr_value[corr_value == 0] = "A"
	corr_value[corr_value == 1] = "C"
	corr_value[corr_value == 2] = "G"
	corr_value[corr_value == 3] = "T"
	if (err_pos > 0) { 
		substr(index_cur,(err_pos),(err_pos))<-(corr_value)
	} 
	return(index_cur)
})
error_count=sum(!(cor_index == index))
return(list(cor_index,length(index),error_count))
# It returns a list containing corrected barcodes, total number of index and number of index containing 1 mismatch. 
}, ex=function(){
  print("Correcting index...")
  index=c("GACATTGGG","CCTCTACAA","CATATGTGG","ACTCTCGCT","GTTTTAGGA") # correct barcodes are GATATTGGG, CCTCTTCAA, CATATGTGG, ACTCTCTCT, GTTTTAGGG 
  corrected_index = correct_barcode(index)
  print("done.")
  print(paste(corrected_index[[3]],"errors have been corrected within the",corrected_index[[2]], "index list provided"))
})

