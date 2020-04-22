#R code for calculating concordance between duplicates at each snp 

args = commandArgs(trailingOnly=TRUE)
duplicates_filename = args[1]
duplicates_to_drop = args[2]
threshold = args[3]
d = read.table(duplicates_filename, header=TRUE)

n_pairs = (dim(d)[2]-6)/2
C = c()
var_pairs = list()  
for (i in 1:n_pairs) {
	index_1 = 2*i-1+6
	index_2=2*i+6
	var1 = d[,index_1]
	var2 = d[,index_2]
	
	var_pairs[[i]] = c(names(d)[index_1],names(d)[index_2])
	
	if ((sum(is.na(var1))<1500 && mean(var1, na.rm=TRUE)>0.01) && (sum(is.na(var2))<1500 && mean(var2, na.rm=TRUE)>0.01)) {
		C[i] = cor(var1,var2, use="complete.obs")
	} else{
		C[i] = NA
	}
 	
}

#var_pairs[which(C<0.99)]
#C[which(C<0.99)]
#1-length(C[which(C<0.99)])/length(C)
#head(d[,unlist(var_pairs[which(C<0.99)])])

remove_list = unlist(var_pairs[which(C<threshold)])
write(remove_list, file=duplicates_to_drop)