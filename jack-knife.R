#Script for jack-knife

getFullTable <- function(path) {
	name_txt <- data.frame()
	z<-data.frame()
	temp = list.files(path = path, pattern="*.breakage")   
	for (i in 1:length(temp)) {
		a<-read.table(paste0(path,temp[i],sep=""), sep="\t", row.names=1, stringsAsFactors=F)
		if (length(name_txt) == 0) {
			name_txt <- rownames(a)
		} else {
			name_txt <- intersect(name_txt, rownames(a))
		}
	}
	for (i in 1:length(temp)) {
		a<-read.table(paste0(path,temp[i],sep=""), sep="\t", row.names=1, stringsAsFactors=F)
		if (length(z) == 0) {
			z <- a[name_txt,]
		} else {
			z<-cbind(z, a[name_txt,])
		}
	}
	rownames(z) <- name_txt
	return(z)
}
#First Argument is directory with control breakages
#Second Argument is directory with tumor breakages
args <- commandArgs()
con <- getFullTable(args[6])
exp <- getFullTable(args[7])
intintint <- intersect(rownames(con), rownames(exp))   
print(intintint)
full_table <- cbind(con[intintint,], exp[intintint,])  
#For script speed up we use only significant differ CpG islands
ff<-unlist(lapply(intintint, function(x) { return(wilcox.test(exp[x,], con[x,])$p.value < 0.05) }))

full_table <- full_table[ff,]
tp <- 0
fp <- 0
tn <- 0
fn <- 0
res_or<-c(rep.int(0, length(colnames(con))), rep.int(1, length(colnames(exp))))
library(e1071)
for (i in seq(1, length(res_or), 1)) {
	res_ss <- res_or[-c(i)]
	a<-svm(res_ss ~ ., t(full_table[,-c(i)]), kernel="linear")
	res <- predict(a, t(full_table[,i]))
	print(res)
	print(res_or[i])
	if ((res < 0.5)&&(res_or[i] == 0)) {
		tn <- tn + 1
	}
	if ((res < 0.5)&&(res_or[i] == 1)) {
		fn <- fn + 1
	}
	if ((res > 0.5)&&(res_or[i] == 1)) {
		tp <- tp + 1
	}
	if ((res > 0.5)&&(res_or[i] == 0)) {
		fp <- fp + 1
	}
}
write.table(file='file.txt', x=full_table)
print(tp)
print(fp)
print(tn)
print(fn)
