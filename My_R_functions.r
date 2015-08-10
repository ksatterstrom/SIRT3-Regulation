# Functions in this file:
#
# norm_sarika_array
# count_GO_terms
# get_GO_enrichment
# filter_probes
# combine_probes
# combine_expression
# filter_by_GO
# filter_by_GSEA
# match_refseqs
# get_Lever_auc
# get_Lever_pvals
# plot_Lever_output

##################################################
norm_sarika_array <- function(sarika_table) {
##################################################

 # norm_sarika_array function:
 # -This assumes that you've got the M17 samples in columns 5 to 8
 # -This assumes that you've got the wt samples in columns 13 to 16
 # -This assumes that you the gene column is labeled "gene"
 # -This normalizes by the average of the wt samples and returns the
 # M17 and wt samples (wt samples first, then M17 samples)

 temp_table = data.frame(sarika_table[,13:16], sarika_table[,5:8], row.names = sarika_table$gene)

 for (i in 1:nrow(temp_table)) {
   current_mean = mean(c(temp_table[i,1], temp_table[i,2], temp_table[i,3], temp_table[i,4]))
   for (j in 1:8) {
      temp_table[i,j] = temp_table[i,j] / current_mean
   }
 }

 return(temp_table)

}

##################################################
count_GO_terms <- function(corr_table, GO_terms) {
##################################################

 # count_GO_terms function:
 # -On the input matrix, the GO column *must* be labeled "GO"
 # -If there is more than one GO number there, they must be
 # separated by " | "

 current_length = length(GO_terms)
 running_counts = rep(0, current_length)

 # Only keep the columns that have one of the GO terms that we want
 for (i in 1:nrow(corr_table)) {
   current_GO = corr_table$GO[i]
   current_GO = as.numeric(unlist(strsplit(current_GO, split = " | ", fixed = T)))
   for (j in 1:current_length) {
      GO_term_J = GO_terms[j]
      if (GO_term_J %in% current_GO) running_counts[j] = running_counts[j]+1 
   }
 }

 return(running_counts)

}


################################################################
get_GO_enrichment <- function(corr_table, GO_terms, GO_counts) {
################################################################

 # get_GO_enrichment function:
 # -On the input matrix, the GO column *must* be labeled "GO"
 # -If there is more than one GO number there, they must be
 # separated by " | "

 neg_scores = GO_counts/nrow(corr_table)
 pos_scores = rep(1, length(neg_scores)) - neg_scores

 current_length = length(GO_terms)
 E_scores = rep(0, current_length)
 Cum_E_scores = E_scores

 corr_table = corr_table[order(corr_table$distance), ]  

 for (i in 1:nrow(corr_table)) {
   current_GO = corr_table$GO[i]
   current_GO = as.numeric(unlist(strsplit(current_GO, split = " | ", fixed = T)))
   for (j in 1:current_length) {
      GO_term_J = GO_terms[j]
      if (GO_term_J %in% current_GO) {
	    E_scores[j] = E_scores[j] + pos_scores[j] } else {
	    E_scores[j] = E_scores[j] - neg_scores[j] } 
   }
   Cum_E_scores = rbind(Cum_E_scores, E_scores)
 }

 return(Cum_E_scores)

}


#######################################
filter_probes <- function(corr_table) {
#######################################

 # filter_probes function:
 # -On the input matrix, the gene column *must* be labeled "gene"
 # -Likewise, the ID column *must* be labeled "id"
 # -This function is only capable of removing _s_ and _x_ probes
 # -It leaves _s_ and/or _x_ probes if there is no normal probe for a gene
 # -The returned matrix is ordered by gene

 # Remove un-annotated cells
 corr_table = subset(corr_table, gene != "")

 # Order by gene
 corr_table = corr_table[order(corr_table$gene), ] 

 # Start going through, checking for ones to remove
 num_deleted = 1
 num_cycles = 0

 while (num_deleted != 0) {   # iterate as long as I'm deleting
 num_cycles = num_cycles + 1  # I want to know how many iterations in takes
 j = 1				# this will be the index for the rows to delete
 rows_to_delete = 0		# this will store the rows to delete
 num_deleted = 0      		# store the number deleted in the current cycle
	for (i in 2:nrow(corr_table)) {
		# If gene is the same as the previous one
 		if (corr_table$gene[i] == corr_table$gene[i-1]) {
 			# If the previous ID does not have _x_ or _s_
			if (grepl("_x_",corr_table$id[i-1]) == FALSE & grepl("_s_",corr_table$id[i-1]) == FALSE) {
				# If the current ID does have _x_ or _s_
				if (grepl("_x_",corr_table$id[i]) == TRUE | grepl("_s_",corr_table$id[i]) == TRUE) {
					# Store the current row as one to be deleted
					rows_to_delete[j] = i
					num_deleted = num_deleted + 1
					j = j+1
				}
			}
		}
 	}
	# If something needs to be deleted, store the new Correlation Table as itself, minus the rows to delete
	if(num_deleted == 0) break else corr_table <- corr_table[-c(rows_to_delete),] 
 }

 # Reverse the row order
 corr_table <- corr_table[nrow(corr_table):1,]

 # Go through again, with things in reverse order
 num_deleted = 1
 num_cycles = 0

 while (num_deleted != 0) {   # iterate as long as I'm deleting
 num_cycles = num_cycles + 1  # I want to know how many iterations in takes
 j = 1				# this will be the index for the rows to delete
 rows_to_delete = 0		# this will store the rows to delete
 num_deleted = 0      		# store the number deleted in the current cycle
	for (i in 2:nrow(corr_table)) {
		# If gene is the same as the previous one
 		if (corr_table$gene[i] == corr_table$gene[i-1]) {
 			# If the previous ID does not have _x_ or _s_
			if (grepl("_x_",corr_table$id[i-1]) == FALSE & grepl("_s_",corr_table$id[i-1]) == FALSE) {
				# If the current ID does have _x_ or _s_
				if (grepl("_x_",corr_table$id[i]) == TRUE | grepl("_s_",corr_table$id[i]) == TRUE) {
					# Store the current row as one to be deleted
					rows_to_delete[j] = i
					num_deleted = num_deleted + 1
					j = j+1
				}
			}
		}
 	}
	# If something needs to be deleted, store the new Correlation Table as itself, minus the rows to delete
	if(num_deleted == 0) break else corr_table <- corr_table[-c(rows_to_delete),] 
 }

 # Put the row order back right
 corr_table <- corr_table[nrow(corr_table):1,]

 return(corr_table)

}


########################################
combine_probes <- function(corr_table) {
########################################

 # combine_probes function:
 # -On the input matrix, the gene column *must* be labeled "gene"
 # -Likewise, the distance column *must* be labeled "distance"
 # -The returned matrix is ordered by gene
 # -I have not yet implemented an option for combination strategies
 # other than averaging

 # Filter out un-annotated genes 
 corr_table = corr_table[complete.cases(corr_table$gene),]

 # Order by gene
 corr_table = corr_table[order(corr_table$gene), ]

 # I'll cycle through the list
 counter = 0			# I'll need a running counter
 j = 1				# this will be the index for the rows to delete
 rows_to_delete = 0		# this will store the rows to delete
 num_deleted = 0

 for (i in 2:nrow(corr_table)) {
	# If gene is different from the previous one
	if (corr_table$gene[i] != corr_table$gene[i-1]) {
		counter = 1	# reset the counter
	}
	# If gene is the same as the previous one
	if (corr_table$gene[i] == corr_table$gene[i-1]) {
		# Mark the previous row for deletion
		rows_to_delete[j] = i-1
		num_deleted = num_deleted+1
		j = j+1
		# Change current correlation value to weighted average of correlation value
		corr_table$distance[i] = ((corr_table$distance[i] + counter*corr_table$distance[i-1])/(counter + 1))
		# Increment running counter
		counter = counter+1
 	}
 }
 # If something needs to be deleted, store the new Correlation Table as itself, minus the rows to delete
 if(num_deleted > 0) corr_table <- corr_table[-c(rows_to_delete),] 

 return(corr_table)

}


########################################
combine_expression <- function(corr_table, first_col, second_col) {
########################################

 # combine_probes function:
 # -On the input matrix, the gene column *must* be labeled "gene"
 # -The expression values should probably be natural scale
 # -The returned matrix is ordered by gene
 # -I have not yet implemented an option for combination strategies
 # other than geometric mean

 # Filter out un-annotated genes 
 corr_table = corr_table[complete.cases(corr_table$gene),]

 # Order by gene
 corr_table = corr_table[order(corr_table$gene), ]

 # I'll cycle through the list
 counter = 0			# I'll need a running counter
 rows_to_delete = 0		# this will store the rows to delete 
 k = 1				# this will be the index for the rows to delete
 num_deleted = 0

 # Temporarily add extra rows to corr_table
 corr_table[nrow(corr_table)+1,] = corr_table[1,] # add 1 extra row
 corr_table[nrow(corr_table)+1,] = corr_table[1,] # add second extra row
 corr_table[2:(nrow(corr_table)-1),] = corr_table[1:(nrow(corr_table)-2),] # re-center so there's a row buffer either side
 corr_table[1,] = corr_table[(nrow(corr_table)-1),] # the extra first row is the actual last; the extra last row is the actual first
 
 # Mark first row for deletion later 
 rows_to_delete[k] = 1
 k = 2 

 # Now go ahead
 for (i in 2:(nrow(corr_table)-1)) {

	# If gene is different from the previous one
	if (corr_table$gene[i] != corr_table$gene[i-1]) {
		counter = 1 # reset counter
		row_holder = corr_table[i,] # store the row's expression value
	}

	# If gene is the same as the previous one
	if (corr_table$gene[i] == corr_table$gene[i-1]) {
		# Mark the previous row for deletion
		rows_to_delete[k] = i-1
		num_deleted = num_deleted+1
		k = k+1
		# Increment running counter
		counter = counter+1		
		# Keep storing the row
		row_holder[counter,] = corr_table[i,] 
 	}

	# If the next gene is different
	if (corr_table$gene[i] != corr_table$gene[i+1]) {
		if (counter > 1) { # if counter is greater than 1
			print(i)
			# store the geometric means of the columns in row_holder as the new expression values
			for (j in first_col:second_col) {
				corr_table[i,j] = geometric.mean(row_holder[,j])
			}
		}
	}
 }
 # If something needs to be deleted, store the new Correlation Table, minus the rows to delete
  
 rows_to_delete[k] = nrow(corr_table) # get rid of that last row that I added
 corr_table <- corr_table[-c(rows_to_delete),]

 return(corr_table)

}


########################################
combine_expression_avg <- function(corr_table, first_col, second_col) {
########################################

 # combine_probes function:
 # -On the input matrix, the gene column *must* be labeled "gene"
 # -The expression values should probably be natural scale
 # -The returned matrix is ordered by gene
 # -This combines expression values using an average

 # Filter out un-annotated genes 
 corr_table = corr_table[complete.cases(corr_table$gene),]

 # Order by gene
 corr_table = corr_table[order(corr_table$gene), ]

 # I'll cycle through the list
 counter = 0			# I'll need a running counter
 rows_to_delete = 0		# this will store the rows to delete 
 k = 1				# this will be the index for the rows to delete
 num_deleted = 0

 # Temporarily add extra rows to corr_table
 corr_table[nrow(corr_table)+1,] = corr_table[1,] # add 1 extra row
 corr_table[nrow(corr_table)+1,] = corr_table[1,] # add second extra row
 corr_table[2:(nrow(corr_table)-1),] = corr_table[1:(nrow(corr_table)-2),] # re-center so there's a row buffer either side
 corr_table[1,] = corr_table[(nrow(corr_table)-1),] # the extra first row is the actual last; the extra last row is the actual first
 
 # Mark first row for deletion later 
 rows_to_delete[k] = 1
 k = 2 

 # Now go ahead
 for (i in 2:(nrow(corr_table)-1)) {

	# If gene is different from the previous one
	if (corr_table$gene[i] != corr_table$gene[i-1]) {
		counter = 1 # reset counter
		row_holder = corr_table[i,] # store the row's expression value
	}

	# If gene is the same as the previous one
	if (corr_table$gene[i] == corr_table$gene[i-1]) {
		# Mark the previous row for deletion
		rows_to_delete[k] = i-1
		num_deleted = num_deleted+1
		k = k+1
		# Increment running counter
		counter = counter+1		
		# Keep storing the row
		row_holder[counter,] = corr_table[i,] 
 	}

	# If the next gene is different
	if (corr_table$gene[i] != corr_table$gene[i+1]) {
		if (counter > 1) { # if counter is greater than 1
			print(i)
			# store the geometric means of the columns in row_holder as the new expression values
			for (j in first_col:second_col) {
				corr_table[i,j] = mean(row_holder[,j])
			}
		}
	}
 }
 # If something needs to be deleted, store the new Correlation Table, minus the rows to delete
  
 rows_to_delete[k] = nrow(corr_table) # get rid of that last row that I added
 corr_table <- corr_table[-c(rows_to_delete),]

 return(corr_table)

}


################################################
filter_by_GO <- function(corr_table, GO_terms) {
################################################

 # filter_by_GO function:
 # -On the input matrix, the GO column *must* be labeled "GO"
 # -If there is more than one GO number there, they must be
 # separated by " | "

 keep_rows = 0

 # Only keep the columns that have one of the GO terms that we want
 for (i in 1:nrow(corr_table)) {
   current_GO = corr_table[i,]$GO
   current_GO = as.character(unlist(strsplit(current_GO, split = " | ", fixed = T)))
   for (j in 1:length(current_GO)) {
      current_GO_J = current_GO[j]
      cond = current_GO_J %in% GO_terms
      if (cond) {
	   keep_rows = c(keep_rows, i)
         break }  # this stops looking for GO terms for a gene after it finds the first one successfully
   }
}

 keep_rows = keep_rows[-1] # this gets rid of the first zero
 corr_table = corr_table[keep_rows,] 

 return(corr_table)

}


################################################
filter_by_GSEA <- function(corr_table, GSEA_terms) {
################################################

 # filter_by_GSEA function:
 # -On the input matrix, the GSEA column *must* be labeled "GSEA"
 # -If there is more than one GSEA number there, they must be
 # separated by " | "

 keep_rows = 0

 # Only keep the columns that have one of the GSEA terms that we want
 for (i in 1:nrow(corr_table)) {
   current_GSEA = corr_table[i,]$GSEA
   current_GSEA = as.character(unlist(strsplit(current_GSEA, split = " | ", fixed = T)))
   for (j in 1:length(current_GSEA)) {
      current_GSEA_J = current_GSEA[j]
      cond = current_GSEA_J %in% GSEA_terms
      if (cond) {
	   keep_rows = c(keep_rows, i)
         break }  # this stops looking for terms for a gene after it finds the first one successfully
   }
}

 keep_rows = keep_rows[-1] # this gets rid of the first zero
 corr_table = corr_table[keep_rows,] 

 return(corr_table)

}


#################################################
match_refseqs <- function(corr_table, linetags) {
#################################################

 # match_refseqs function:
 # -On the input matrix, the refseqs column *must* be labeled "refseqs"
 # -If there is more than one refseq for a gene, they must be separated by " | "
 # -You must also pass it a list of refseq IDs from the sequence file
 # -It returns corr_rows, which is a list of the lines in the input matrix
 # which have a matching refseq in the sequence file, and sequence_rows, which
 # is the row number in the sequence file of the match

 take = corr_table$refseqs
 corr_rows = 0 
 sequence_rows = 0

 # for each id of interest (in take), find corresponding lines from the sequence file
 for (i in 1:length(take)) {
   thisone = take[i]
   thisone = as.character(unlist(strsplit(thisone, split = " | ", fixed = T)))
   for (j in 1:length(thisone)) {
     thisoneJ = thisone[j]
     if (thisoneJ %in% linetags) {
	  corr_rows = c(corr_rows, i)
        sequence_rows = c(sequence_rows, match(thisoneJ, linetags))
        break  }
   } 
 } # this stops looking for accessions for a gene after it finds the first one successfully
 corr_rows = corr_rows[-1]
 sequence_rows = sequence_rows[-1] # this gets rid of the first zero
 
 output <- list(corr_rows, sequence_rows)
 return(output) 

}

#################################################
get_Lever_auc <- function(input_file) {
#################################################

auc_vals = read.delim(input_file, header = F, stringsAsFactors = FALSE, skip = 1)
auc_vals = auc_vals[,6]
auc_vals = as.numeric(gsub("stat=", "", auc_vals))

return(auc_vals)

}


#################################################
get_Lever_pvals <- function(input_file) {
#################################################

p_vals = read.delim(input_file, header = F, stringsAsFactors = FALSE, skip = 1)
p_vals = p_vals[,5]
p_vals = as.numeric(gsub("pval=", "", p_vals))

return(p_vals)

}


##################################################
plot_Lever_output <- function(in_file, out_file, motif_file) {
##################################################

XX = read.delim(in_file, header = F, stringsAsFactors = FALSE, skip = 1)

set = XX[,2]
set = as.numeric(gsub("SET=", "", set))
index = XX[,4]
index = gsub("COMBO INDEX=(", "", index, fixed = T)
index = gsub(")", "", index, fixed = T)
index = as.numeric(index)
nullmean = XX[,7]
nullmean = as.numeric(gsub("mean=", "", nullmean))
nullsd = XX[,8]
nullsd = as.numeric(gsub("sd=", "", nullsd))
nulllow = nullmean - nullsd
nullhigh = nullmean + nullsd
auc = XX[,6]
auc = as.numeric(gsub("stat=", "", auc))
qvalue = XX[,5]
qvalue = as.numeric(gsub("pval=", "", qvalue))

auclimits = XX[,9]
auclimits = gsub("error=(", "", auclimits, fixed = T)
auclimits = gsub(")", "", auclimits, fixed = T)
auclow = 0
auchigh = 0
for (i in 1:length(auclimits)) {
  thisonee = auclimits[i]
  thisonee = as.numeric(unlist(strsplit(thisonee, split = ",")))
  auclow[i] = thisonee[1]
  auchigh[i] = thisonee[2]  }

motifs = read.delim(motif_file, header = F, stringsAsFactors = FALSE, flush = T)
labels = 0
options("warn" = -1)
for (i in 1:nrow(motifs)) {
 thisonee = motifs[i,1]
 cond = is.na(as.numeric(thisonee))
 if (cond) labels = c(labels, thisonee) }
labels = labels[-1]
options("warn" = 0) # here, labels are as given in the motif file
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# add consensus sites:
consensus = 0
for (i in 1:length(labels)) {
 thisonee = labels[i]
 value = match(thisonee, motifs[,1])
 theA = motifs[value+1,]
 theC = motifs[value+2,]
 theG = motifs[value+3,]
 theT = motifs[value+4,]
 # ---------------------
 theA = theA[,complete.cases(as.numeric(theA))]
 theC = theC[,complete.cases(as.numeric(theC))]
 theG = theG[,complete.cases(as.numeric(theG))]
 theT = theT[,complete.cases(as.numeric(theT))]
 # ----------------------
 mat = rbind(theA, theC, theG, theT)
 mat = data.frame(mat, stringsAsFactors = F)
 row.names(mat) = c("A", "C", "G", "T")
 letter = 0
 for (j in 1:ncol(mat)) {
   thiscol = as.numeric(mat[,j])
   maxposition = which.max(thiscol)
   letter[j] = row.names(mat)[maxposition] }
 letter = paste(letter, collapse = "") 
 consensus[i] = letter }
stopifnot(length(consensus) == length(labels))
# --------------------------------------------
# --------------------------------------------
# make new labels with consensus sites:
labels2 = 0
for (i in 1:length(labels)) {
  thelab = labels[i]
  thesite = consensus[i]
  thelab = as.character(unlist(strsplit(thelab, split = "_")))
  newlab = thelab[1]
  newlab = paste(newlab, "(", sep = " ")
  newlab = paste(newlab, thesite, sep = "")
  newlab = paste(newlab, ")", sep = "")
  # ----------------------------------------
  cond = length(grep("_MA", labels[i])) == 1
  if (cond) {
  add = thelab[grep("MA", thelab)]
  newlab = paste(newlab, "[", sep = " ")
  newlab = paste(newlab, add, sep = "")
  newlab = paste(newlab, "]", sep = "") } 
labels2[i] = newlab  }

XX = data.frame(set, index, labels2, nullmean, nulllow, nullhigh, auc, auclow, auchigh, qvalue, stringsAsFactors = FALSE)

# ========================================================================================

take = 25
set = subset(XX, set == 0)
set = set[rev(order(set$auc)),]
set = set[1:take,]
set = set[nrow(set):1,]

pdf(out_file, width = 8.50, height = 5.50, onefile = T)

par(mai = c(0.80, 3.00, 0.20, 1.40))
xmin = min(min(set$nulllow), min(set$auclow))
xmax = max(max(set$nullhigh), max(set$auchigh))
plot.new()
plot.window(xlim = c(xmin, xmax), ylim = c(1, take))
rect(xmin - 100, -100, xmax + 100, take + 100, col = "yellow")
abline(v = 0.50, col = 2, lwd = 2)

vsize = 0.50
hsize = 0.20
for (i in 1:take) {
 null1 = set[i,]$nulllow
 null2 = set[i,]$nullhigh
 stat = set[i,]$auc
 stat1 = set[i,]$auclow
 stat2 = set[i,]$auchigh
 theq = set[i,]$qvalue
 theq = round(theq, 1)
 if (theq == 0) theq = "0.00"
 if (theq == 1) theq = "1.00"
 theq = paste("(P = ", theq, sep = "")
 theq = paste(theq, ")", sep = "")
 theq = paste(round(stat,3), theq, sep = " ")
 thelab = set[i,]$labels
 abline(h = i, lty = 2) 
 rect(null1, i - vsize, null2, i + vsize, col = "grey")
 points(stat, i, pch = 16, col = 1, cex = 1.40)
 segments(stat1, i, stat2, i, col = 1, lwd = 2)
 segments(stat1, i - hsize, stat1, i + hsize, col = 1, lwd = 2)
 segments(stat2, i - hsize, stat2, i + hsize, col = 1, lwd = 2)  
 mtext(theq, at = i, side = 4, line = 0.20, font = 2, las = 2, adj = 0) 
 mtext(thelab, at = i, side = 2, line = 0.20, font = 2, las = 2, adj = 1, cex = 0.70) }

box(lwd = 2)
axis(1, line = -0.5, lwd = 0, tick = F, col = "white", cex.axis = 0.80, las=2, font=2)
axis(1, labels = FALSE, tick = T, tcl = -0.3, font=2)
mtext("AUC Statistic (Enrichment)", side = 1, line = 2.40, font = 2, cex = 1.30)

graphics.off()

}