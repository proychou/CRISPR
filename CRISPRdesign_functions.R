#Scripts for CRISPR guide design
#Pavitra Roychoudhury
#Sep 2017

#Return a list of guides matching a pam sequence and specified length
#ref_seq is a DNAString object
#pattern_fwd is character, e.g. 'NGG' default
find_guides<-function(ref_seq,pattern_fwd='NGG',guide_length=20){
	require('Biostrings'); 
	require('dplyr');
	require('stringr')
	
	ref_seq;
	pattern_fwd<-DNAString(pattern_fwd)
	
	#Create reverse complement
	ref_seq_rev<-reverseComplement(ref_seq);
	ref_seq_rev;
	
	#Look for PAM sequence matches: Search reference sequence for PAM sequence 
	# and return hits in forward and reverse direction.
	matches_fwd<-matchPattern(pattern_fwd,ref_seq,fixed=F);
	matches_rev<-matchPattern(pattern_fwd,ref_seq_rev,fixed=F);
	matches_fwd;
	matches_rev;
	print('Initial hits:')
	print(length(matches_fwd)+length(matches_rev))
	
	#Find candidate targets: Get target seqs (fwd and rev) for each PAM sequence match from above.
	#Forward
	upstream_inds_fwd<-IRanges(start(matches_fwd)-guide_length,
														 width=guide_length+nchar(pattern_fwd));
	upstream_inds_fwd<-upstream_inds_fwd[start(upstream_inds_fwd)>0];
	guide_seqs_fwd<-sapply(upstream_inds_fwd,function(x)return(as.character(ref_seq[x])));
	rm(upstream_inds_fwd);
	
	#Reverse
	upstream_inds_rev<-IRanges(start(matches_rev)-guide_length,
														 width=guide_length+nchar(pattern_fwd));
	upstream_inds_rev<-upstream_inds_rev[start(upstream_inds_rev)>0];
	guide_seqs_rev<-sapply(upstream_inds_rev,function(x)return(as.character(ref_seq_rev[x])));
	rm(upstream_inds_rev);
	
	#Preview candidate guides (with the pam sequence)
	head(guide_seqs_fwd,10);
	head(guide_seqs_rev,10);
	
	#Remove duplicates
	guide_seqs_fwd<-unique(guide_seqs_fwd);
	guide_seqs_rev<-unique(guide_seqs_rev);
	print('Unique hits:')
	print(length(guide_seqs_fwd)+length(guide_seqs_rev))
	
	#Analyze guide frequency: For each candidate guide, record direction it was found, length
	# (they're all the same length),and whether there's a G in position 1, 2 and 3. 
	# Then look for the number of occurrences of each guide in the forward and reverse 
	# directions within the list of sequences for all group M.
	
	#Analyze guides
	candidate_guides<-data.frame(
		gseq=c(guide_seqs_fwd,guide_seqs_rev),
		dir=c(rep('F',length(guide_seqs_fwd)),
					rep('R',length(guide_seqs_rev))),
		stringsAsFactors=F)%>%
		mutate(length=str_length(gseq),
					 guide=substr(gseq,1,guide_length),
					 pam=substr(gseq,guide_length+1,guide_length+nchar(pattern_fwd)));
	head(candidate_guides);
	return(candidate_guides)
	
}
