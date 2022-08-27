#' Prepare data for GUST classification
#' return a list of data frames, containing fitch rates of protein positions.
#' @examples initialize.gust();
initialize.gust <- function() {
	cat('initializing ...\n'); flush.console();
	fitch.block <- list();
	for(k in 1:10) {
		fitch.block.name <- paste('fitch.block.', k, sep='');
		fitch.block[[k]] <- get(fitch.block.name)
	}
	assign('fitch.block', fitch.block, envir=.GlobalEnv);
#	cat(sum(unlist(lapply(fitch.block, nrow))), '\n');
}

#' Find hypo- and hyper-mutated samples
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write the output files
#' @param output.prefix: prefix used to name the output file as prefix.outlier.txt .
#' @return NULL
#' @examples find.outliers(input.file.name='./examples/TCGA.ACC.mutect.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.ACC');
find.outliers <- function(input.file.name, output.folder, output.prefix) {
	cat('finding outliers ...\n', '   ', input.file.name, '...'); flush.console();
	outlier.file.name <- paste(output.folder, '/', output.prefix, '.outlier.txt', sep='');
	if(nchar(output.prefix) == 0) {
		outlier.file.name <- paste(working.dir, 'outlier.txt', sep='');
	}

	mut.all <- read.table(input.file.name, sep='\t', header=T, quote='', stringsAsFactors=F);
	cat(nrow(mut.all), 'lines\n'); flush.console();
	mut.filtered <- mut.all[, c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'dbSNP_RS', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'One_Consequence', 'Hugo_Symbol', 'Gene', 'Feature', 'ENSP', 'HGVSc', 'HGVSp_Short', 'Amino_acids', 'Codons', 'ENSP', 'RefSeq', 'Entrez_Gene_Id')];
	colnames(mut.filtered) <- c('sample', 'chr', 'pos', 'rsid', 'ref', 'alt', 'filter', 'type', 'symbol', 'ensembl.gene', 'cds.id', 'protein.id', 'cds.comp', 'prot.comp', 'aa', 'codon', 'prot.id', 'refseq', 'gene.id');
	mut.filtered <- mut.filtered[which(mut.filtered$type %in% c('synonymous_variant', 'missense_variant', 'stop_gained', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')), ];
	agg <- aggregate(chr ~ sample, data=mut.filtered, 'length')
	agg <- agg[order(agg$chr), ]
	# lower.cutoff <- m - 3*d;
	# upper.cutoff <- m + 3*d;	
	x <- log(agg$chr, 2)
	m <- mean(x)
	d <- sd(x)
	q1 <- quantile(x, 0.25);
	q3 <- quantile(x, 0.75);
	iqr <- IQR(x)
	lower.cutoff <- 2^(q1 - iqr*1.5);
	upper.cutoff <- 2^(q3 + iqr*1.5);
	lower.outlier <- which(agg$chr < lower.cutoff)
	upper.outlier <- which(agg$chr > upper.cutoff)
	cat('     mean (', round(2^m), '), std (', round(2^d), '), cutoffs (', lower.cutoff, ' - ', upper.cutoff, ')\n'); flush.console();
	cat('     hypo- (', length(lower.outlier), ') ... hyper (', length(upper.outlier), ')\n'); flush.console();
	
	agg$outlier <- 'no';
	agg[lower.outlier, 'outlier'] <- 'hypo';
	agg[upper.outlier, 'outlier'] <- 'hyper';

	write.table(agg, outlier.file.name, sep='\t', quote=F, row.names=F);
}

#' Parse the VCF file to aggregate variants by genes and mutation types
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write output files
#' @param output.prefix: prefix used to name the output files. Output files include prefix.error.txt (mutations cannot be mapped), prefix.mut.filtered.txt (somatic variants with matching types, prefix.mut.cnt.txt (somatic variants aggregated by sample, gene and mutational type), prefix.mut.summary.txt (somatic variants aggregated by gene and mutational type across all samples), prefix.symbol_2_cds_id.txt (mapping gene symbols to ensembl transcript ids), and prefix.log.txt (log information).
#' @return NULL
#' @examples parse.aggregated.mut(input.file.name='./examples/TCGA.ACC.mutect.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.ACC')
parse.aggregated.mut <- function(input.file.name, output.folder='', output.prefix='') {
	cat('parsing aggregated calls ...\n');
	cat('   ', input.file.name, '...'); flush.console();

	error.file.name <- paste(output.folder, '/', output.prefix, '.error.txt', sep='');
	outlier.file.name <- paste(output.folder, '/', output.prefix, '.outlier.txt', sep='');
	log.file <- file(paste(output.folder, '/', output.prefix, '.log.txt', sep=''), open='a');
	mut.filtered.file.name <- paste(output.folder, '/', output.prefix, '.mut.filtered.txt', sep='');
	mut.cnt.file.name <- paste(output.folder, '/', output.prefix, '.mut.cnt.txt', sep='');
	mut.summary.file.name <- paste(output.folder, '/', output.prefix, '.mut.summary.txt', sep='');
	map.file.name <- paste(output.folder, '/', output.prefix, '.symbol_2_cds_id.txt', sep='');
	if(nchar(output.prefix) == 0) {
		error.file.name <- paste(working.dir, 'error.txt', sep='');
		outlier.file.name <- paste(working.dir, 'outlier.txt', sep='');
		log.file.name <- paste(working.dir, 'log.txt', sep='');
		mut.filtered..file.name <- paste(working.dir, 'mut.filtered..txt', sep='');
		mut.cnt.file.name <- paste(working.dir, 'mut.cnt.txt', sep='');
		mut.summary.file.name <- paste(working.dir, 'mut.summary.txt', sep='');
		map.file.name <- paste(working.dir, 'symbol_2_cds_id.txt', sep='');
	}

	mut.all <- read.table(input.file.name, sep='\t', header=T, quote='', stringsAsFactors=F);
	cat(nrow(mut.all), 'lines\n'); flush.console();
	mut.filtered <- mut.all[, c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'dbSNP_RS', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'One_Consequence', 'Hugo_Symbol', 'Gene', 'Feature', 'ENSP', 'HGVSc', 'HGVSp_Short', 'Amino_acids', 'Codons', 'ENSP', 'RefSeq', 'Entrez_Gene_Id')];
	colnames(mut.filtered) <- c('sample', 'chr', 'pos', 'rsid', 'ref', 'alt', 'filter', 'type', 'symbol', 'ensembl.gene', 'cds.id', 'protein.id', 'cds.comp', 'prot.comp', 'aa', 'codon', 'prot.id', 'refseq', 'gene.id');
	if(nchar(outlier.file.name) > 0) {
		outliers <- read.table(outlier.file.name, sep='\t', header=T, stringsAsFactors=F);
		mut.filtered <- mut.filtered[which(mut.filtered$sample %in% outliers[which(outliers$outlier == 'no'), 'sample']), ];
	}
	cat('    mutations after removing outlier samples', nrow(mut.filtered), '\n'); flush.console();	
	
	mut.filtered <- mut.filtered[which(mut.filtered$type %in% c('synonymous_variant', 'missense_variant', 'stop_gained', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')), ];
	mut.filtered$type <- ifelse(mut.filtered$type == 'synonymous_variant', 'syn', ifelse(mut.filtered$type == 'missense_variant', 'missense', ifelse(mut.filtered$type == 'stop_gained', 'nonsense', ifelse(mut.filtered$type =='frameshift_variant', 'fs', ifelse(mut.filtered$type =='inframe_insertion', 'inf_ins', 'inf_del')))));
	fs.idx <- which(mut.filtered$type == 'fs');
	fs.idx.ins <- intersect(fs.idx, grep('ins', mut.filtered$cds.comp));
	fs.idx.del <- intersect(fs.idx, grep('del', mut.filtered$cds.comp));
	fs.idx.dup <- intersect(fs.idx, grep('dup', mut.filtered$cds.comp));
	if(sum(length(fs.idx.ins), length(fs.idx.del), length(fs.idx.dup)) != length(fs.idx)) {
		cat('!!!!!!!!!!!!!!!!!!! Error !!!!!!!!!!!!! CHECK FS INDELS!!!!!!!!!\n');  flush.console();	
	}
	mut.filtered[fs.idx.ins, 'type'] <- 'fs_ins'
	mut.filtered[fs.idx.del, 'type'] <- 'fs_del'
	mut.filtered[fs.idx.dup, 'type'] <- 'fs_ins'	
	cat('    coding mutations', nrow(mut.filtered ), '; ');  flush.console();
	
	pattern.cds.1 <- str_match(mut.filtered$cds.comp, 'c\\.([0-9]+)([A-Z]+)>([A-Z]+)')
	pattern.cds.2 <- str_match(mut.filtered$cds.comp, 'c\\.([0-9]+)([ins|del|dup|\\+])(.+)')
	pattern.cds.3 <- str_match(mut.filtered$cds.comp, 'c\\.([0-9]+)_([0-9]+)([ins|del|dup|\\+])(.+)')
	pattern.prot.1 <- str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)([A-Z*])$')
	pattern.prot.2 <- str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)([A-Z])fs(.+)')
	pattern.prot.3 <- str_match(mut.filtered$prot.comp, 'p\\.([A-Z*])([0-9]+)([A-Z*a-z_0-9]+)')
	pattern.codon <- str_match(mut.filtered$codon, '([A-Za-z\\-]+)/([A-Za-z]+)')
	pattern.aa <- str_match(mut.filtered$aa, '([A-Z\\-]+)/([A-Z*]+)')
	mut.filtered$cds.pos <- pattern.cds.1[, 2]
	mut.filtered$cds.ref <- pattern.cds.1[, 3]
	mut.filtered$cds.alt <- pattern.cds.1[, 4]
	idx <- which(is.na(mut.filtered$cds.alt));
	mut.filtered[idx, c('cds.pos', 'cds.ref', 'cds.alt')] <- pattern.cds.2[idx, 2:4]
	idx <- which(is.na(mut.filtered$cds.alt));
	mut.filtered[idx, c('cds.pos', 'cds.ref', 'cds.alt')] <- pattern.cds.3[idx, c(2, 4, 5)]
	mut.filtered$prot.pos <- pattern.prot.1[, 3]
	mut.filtered$prot.ref <- pattern.prot.1[, 2]
	mut.filtered$prot.alt <- pattern.prot.1[, 4]
	idx <- which(is.na(mut.filtered$prot.alt) & mut.filtered$type %in% c('fs_ins', 'fs_del'));
	mut.filtered[idx, c('prot.pos', 'prot.ref', 'prot.alt')] <- pattern.prot.2[idx, c(3, 2, 4)]
	idx <- which(is.na(mut.filtered$prot.alt) & mut.filtered$type %in% c('inf_ins', 'inf_del'));
	mut.filtered[idx, c('prot.pos', 'prot.ref', 'prot.alt')] <- pattern.prot.3[idx, c(3, 2, 4)]
				
	e <- t(apply(mut.filtered[, c('cds.id', 'cds.pos')], 1, function(x) { i<-which(seq_bp$cds.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c('', '', ''); } else { p<-as.numeric(x[2]); s=seq_bp[i, 'cds.seq']; prev=substr(s, p-1, p-1); current=substr(s, p, p); post=substr(s, p+1, p+1); r<-c(prev, current, post);} }))
	f <- unlist(apply(mut.filtered[, c('prot.id', 'prot.pos')], 1, function(x) { i<-which(seq_aa$pep.id.2==x[1]); r <- c(); if(length(i)==0) { r <- c(''); } else { p<-as.numeric(x[2]); s=seq_aa[i, 'pep.seq']; current=substr(s, p, p); r<-c(current); } }))
	mut.filtered$seq.prev <- e[, 1];
	mut.filtered$seq.current <- e[, 2];
	mut.filtered$seq.post <- e[, 3];
	mut.filtered$prot.current <- f;
	index.error.1 <- which(!mut.filtered$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & mut.filtered$seq.current != mut.filtered$cds.ref);
	index.error.2 <- which(!mut.filtered$type %in% c('fs_ins', 'fs_del', 'inf_ins', 'inf_del') & mut.filtered$prot.current != mut.filtered$prot.ref);
	index.error <- unique(c(index.error.1, index.error.2));
	errors <- c();
	if(length(index.error) > 0) {
		errors <- mut.filtered[index.error, ]
		mut.filtered <- mut.filtered[-index.error, ]
	}
	cat(' with matching seq', nrow(mut.filtered), '\n'); flush.console();
	if(nrow(mut.filtered) > 0) {
		mut.filtered$mut_type <- paste(mut.filtered$cds.ref, mut.filtered$cds.alt, sep='-');
		mut.filtered$mut_type <- ifelse(mut.filtered$mut_type == 'C-T' & mut.filtered$seq.post == 'G', 'C-T-h', ifelse(mut.filtered$mut_type == 'G-A' & mut.filtered$seq.prev == 'C', 'G-A-h', mut.filtered$mut_type));
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'fs_ins', 'fs_ins', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'fs_del', 'fs_del', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'inf_ins', 'inf_ins', mut.filtered$mut_type);
		mut.filtered$mut_type <- ifelse(mut.filtered$type == 'inf_del', 'inf_del', mut.filtered$mut_type);
		mut.filtered <- merge(mut.filtered, annotation, by.x='ensembl.gene', by.y='ensembl.gene.id', all.x=T, all.y=F);
		mut.filtered$gene.id <- as.integer(as.character(mut.filtered$gene.id.y));
		mut.filtered$symbol <- mut.filtered$symbol.y;
		write.table(mut.filtered, mut.filtered.file.name, sep='\t', row.names=F, quote=F);
	}
	write.table(errors, error.file.name, sep='\t', row.names=F, quote=F);
		
	agg.error <- aggregate(chr ~ symbol + pos + ensembl.gene + cds.id + prot.id, data=errors, 'length');
	agg.error <- agg.error[order(-agg.error$chr), ]
	agg.error <- aggregate(chr ~ symbol +  ensembl.gene + cds.id + prot.id, data=errors, 'length');
	agg.error <- agg.error[order(-agg.error$chr), ]
	cat('    unmatched mutations:', length(index.error), 'in', length(unique(agg.error$symbol)), 'genes, max cnt in unmatched genes (', agg.error[1, 'symbol'], ':', agg.error[1, 'chr'], ')\n'); flush.console();

	mut.list <- split(mut.filtered, mut.filtered$sample)
	cat('    total samples after removing outliers', length(mut.list), '\n'); flush.console();
	list.genes <- c();
	summary.mut <- c();
	symbol.2.cds.id <- c();
	for(ml in 1:length(mut.list)) {
		mut <- mut.list[[ml]];
		sample.name <- unique(mut$sample);
		# subject.name <- substr(sample.name, 1, 20);
		subject.name <- sample.name;
		cat('        ', ml, '...', subject.name, '...', nrow(mut), '\n', file=log.file); flush.console();
		if(!is.na(subject.name) & nrow(mut) > 0) {
			x <- aggregate(chr ~ symbol + type + mut_type, data=mut, 'length');
			m <- aggregate(chr ~ gene.id + symbol + cds.pos + cds.ref + cds.alt + mut_type + type + prot.pos + prot.ref + prot.alt, data=mut, 'length');
			summary.mut <- rbind(summary.mut, m);
			symbol.2.cds.id <- rbind(unique(symbol.2.cds.id), unique(mut[, c('cds.id', 'symbol')]));
			colnames(x) <- c('symbol', 'type', 'mut_type', subject.name);
			list.genes[[ml]] <- x;
		}
	}
	close(log.file);
	genes.samples <- Reduce(function(...) merge(..., all=T, by=c('symbol', 'type', 'mut_type')), list.genes)
	write.table(genes.samples, mut.cnt.file.name, sep='\t', row.names=F, quote=F);
	cat('    mutations by sample', nrow(genes.samples), 'x', ncol(genes.samples)-3, ';');  
		
	summary.mut <- aggregate(chr ~ gene.id + symbol + cds.pos + cds.ref + cds.alt + mut_type + prot.pos + prot.ref + prot.alt, data=summary.mut, 'sum');
	summary.mut <- summary.mut[order(-summary.mut$chr), ]
	write.table(summary.mut, mut.summary.file.name, sep='\t', row.names=F, quote=F);
	cat('unique mutations', nrow(summary.mut), ';');  

	map <- c();
	for(s in unique(symbol.2.cds.id$symbol)) {
		cds.ids <- symbol.2.cds.id[which(symbol.2.cds.id$symbol == s), 'cds.id'];
		cds.length <- seq.length[which(seq.length$cds.id %in% cds.ids), ];
		longest.id <- cds.length[which.max(cds.length$length), 'cds.id'][1]
		map <- rbind(map, c(s, longest.id));
	}
	colnames(map) <- c('symbol', 'cds.id');
	write.table(map, map.file.name, sep='\t', row.names=F, quote=F);	
	cat('unique genes', nrow(map), '\n');  
}

#' Organize a collection of variants in a gene into counts 
#' @param x: a data frame with four columns: name, type, mut_type, cnt.
#' @return: a data frame with counts of variants grouped by types
organize.obs <- function(x) {
	colnames(x) <- c('name', 'type', 'mut_type', 'cnt');
	mt <- unique(x$mut_type);
	syn.sum <- data.frame(mut_type=mt, sample=rep(0, length(mt)));
	missense.sum <- syn.sum;
	nonsense.sum <- syn.sum;
	fs.sum <- syn.sum;
	inframe.sum <- syn.sum;
	s <- which(x$type == 'syn');
	m <- which(x$type == 'missense');
	n <- which(x$type == 'nonsense');
	fi <- which(x$type == 'fs_ins');
	fd <- which(x$type == 'fs_del');
	ii <- which(x$type == 'inf_ins');
	id <- which(x$type == 'inf_del');
	if(length(s) > 0) {
		syn.sum <- x[s, c('mut_type', 'cnt')];
	}
	if(length(m) > 0) {
		missense.sum <- x[m, c('mut_type', 'cnt')];
	}
	if(length(n) > 0) {
		nonsense.sum <- x[n, c('mut_type', 'cnt')];
	}
	if(length(fi) + length(fd) > 0) {
		fs.sum <- x[c(fi, fd), c('mut_type', 'cnt')];
	}
	if(length(ii) + length(id) > 0) {
		inframe.sum <- x[c(ii, id), c('mut_type', 'cnt')];
	}
	colnames(syn.sum)[2] <- 'cnt.syn';
	colnames(missense.sum)[2] <- 'cnt.missense';
	colnames(nonsense.sum)[2] <- 'cnt.nonsense';
	colnames(fs.sum)[2] <- 'cnt.fs';
	colnames(inframe.sum)[2] <- 'cnt.inframe';
	sum.obs <- Reduce(function(...) merge(..., all=T, by='mut_type'), list(syn.sum, missense.sum, nonsense.sum, fs.sum, inframe.sum))
	colnames(sum.obs) <- c('mut_type', 'syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs');
	sum.obs[is.na(sum.obs)] <- 0;
	sum.obs$mut_type <- as.character(sum.obs$mut_type);
	return(sum.obs);	
}

#' optimization function to estimate independent probabilities of a specific mutation type 
prob <- function(theta, obs_exp, compare) {
	mut_type <- unique(obs_exp$mut_type);
	if(compare == 'missense') {
		mut_type <- mut_type[which(mut_type != 'fs')];
	}
	p <- 0;
	for(t in mut_type) {
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_compare <- 0;
		exp_compare <- 0;
		obs_total <- 0;
		if(compare == 'missense') {
			obs_compare <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
			exp_compare <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		} else if(compare == 'nonsense') {
			obs_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('nonsense.obs', 'fs.obs')]);  ## mut_type=='fs', nonsense.obs=0; otherwise, fs.obs=0;
			exp_compare <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		} else {
			obs_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('missense.obs', 'nonsense.obs', 'fs.obs')]);
			exp_compare <- sum(obs_exp[which(obs_exp$mut_type == t), c('missense.exp', 'nonsense.exp')]);
		}
		obs_total <- obs_syn + obs_compare;
		if(exp_syn == 0 | exp_compare == 0) {
			pp <- 0;
		} else {
			pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_compare) + obs_syn*log(exp_syn) + obs_compare*log(theta*exp_compare) - obs_total*log(exp_syn + theta * exp_compare);
		}
		p <- p + pp;
#		cat(t, '\t', pp, '\t', p, '\n');
	}
	return(-p);
}

#' optimization function to estimate joint probabilities of different types of variants. Indels are treated independently from substitutions.
prob.joint <- function(par, obs_exp) {
#	cat(par);
	phi <- par[1];
	psi <- par[2];
	if(phi < -5 | phi > 5 | psi < -5 | psi > 5) {
		return(NA);
	}
	mut_type <- unique(obs_exp$mut_type);
	p <- 0;
	for(t in mut_type) {
#		cat(t, ':' );
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
		exp_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		obs_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.obs'];
		exp_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		exp_syn <- ifelse(exp_syn == 0, 1, exp_syn);
		exp_missense <- ifelse(exp_missense == 0, 1, exp_missense);
		exp_nonsense <- ifelse(exp_nonsense == 0, 1, exp_nonsense);
		
		pp <- 0;
		if(t == 'inframe') {
		} else if(t == 'fs') {	
			obs_fs <- obs_exp[which(obs_exp$mut_type == t), 'fs.obs'];
			obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inframe'), 'inframe.obs'];
			if(obs_fs > 0 | obs_inframe > 0) {
				gene.length <- sum(obs_exp[which(obs_exp$mut_type == t), c('syn.exp', 'missense.exp', 'nonsense.exp')])/3;
				exp_fs <- gene.length*2/3;  ## expected fs is the length of the protein
	#			exp_fs <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp']; 
				exp_inframe <- gene.length*1/3;  ## expected fs is the length of the protein
				obs_total <- obs_inframe + obs_fs;
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_fs) + obs_syn*log(exp_syn) + obs_fs*log(psi*exp_fs) - obs_total*log(exp_syn + psi*exp_fs);
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_fs) + obs_syn*log(exp_syn) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_syn + exp(psi)*exp_fs);
				pp <- lfactorial(obs_total) - lfactorial(obs_inframe) - lfactorial(obs_fs) + obs_inframe*log(exp_inframe) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_inframe + exp(psi)*exp_fs);
			}
#			pp=0;
		} else {
			obs_total <- obs_syn + obs_missense + obs_nonsense;
			if(obs_total > 0) {
#				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*log(phi*exp_missense) + obs_nonsense*log(psi*exp_nonsense) - obs_total*log(exp_syn + phi*exp_missense + psi*exp_nonsense);
				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*phi + obs_missense*log(exp_missense) + obs_nonsense*psi + obs_nonsense*log(exp_nonsense) - obs_total*log(exp_syn + exp(phi)*exp_missense + exp(psi)*exp_nonsense);
			}
		}
		p <- p + pp;
#		cat(t, '\t', pp, '\t', p, '\n');
	}
#	cat(' ', -p, '\n');
	return(-p);
}

#' optimization function to estimate joint probabilities of different types of variants including indels.
prob.joint.indel <- function(par, obs_exp) {
#	cat(par);
	phi <- par[1];
	psi <- par[2];
	if(phi < -5 | phi > 5 | psi < -5 | psi > 5) {
		return(NA);
	}
	mut_type <- unique(obs_exp$mut_type);
	mut_type <- mut_type[which(mut_type != 'NA-NA')]
	p <- 0;
	for(t in mut_type) {
#		cat(t, ':' );
		obs_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.obs'];
		exp_syn <- obs_exp[which(obs_exp$mut_type == t), 'syn.exp'];
		obs_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.obs'];
		exp_missense <- obs_exp[which(obs_exp$mut_type == t), 'missense.exp'];
		obs_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.obs'];
		exp_nonsense <- obs_exp[which(obs_exp$mut_type == t), 'nonsense.exp'];
		exp_syn <- ifelse(exp_syn == 0, 1, exp_syn);
		exp_missense <- ifelse(exp_missense == 0, 1, exp_missense);
		exp_nonsense <- ifelse(exp_nonsense == 0, 1, exp_nonsense);
		
		pp <- 0;
		if(t %in% c('inf_ins', 'inf_del')) {
		} else if(t %in% c('fs_ins', 'fs_del')) {	
			gene.length <- sum(exp_syn, exp_missense, exp_nonsense)/3;
			exp_fs <- gene.length * 2/3;
			exp_inframe <- gene.length * 1/3;
			obs_fs <- obs_exp[which(obs_exp$mut_type == t), 'fs.obs'];
			obs_inframe <- 0; 
			if(t == 'fs_ins') {
				obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inf_ins'), 'inframe.obs'];
				exp_fs <- exp_fs * 0.34;
				exp_inframe <- exp_inframe * 0.34;
			} else {
				obs_inframe <- obs_exp[which(obs_exp$mut_type == 'inf_del'), 'inframe.obs'];
				exp_fs <- exp_fs * 0.66;
				exp_inframe <- exp_inframe * 0.66;
			}
			obs_total <- obs_inframe + obs_fs;
			if(obs_total > 0) {
				pp <- lfactorial(obs_total) - lfactorial(obs_inframe) - lfactorial(obs_fs) + obs_inframe*log(exp_inframe) + obs_fs*psi + obs_fs*log(exp_fs) - obs_total*log(exp_inframe + exp(psi)*exp_fs);
			}
		} else {
			obs_total <- obs_syn + obs_missense + obs_nonsense;
			if(obs_total > 0) {
				pp <- lfactorial(obs_total) - lfactorial(obs_syn) - lfactorial(obs_missense) - lfactorial(obs_nonsense) + obs_syn*log(exp_syn) + obs_missense*phi + obs_missense*log(exp_missense) + obs_nonsense*psi + obs_nonsense*log(exp_nonsense) - obs_total*log(exp_syn + exp(phi)*exp_missense + exp(psi)*exp_nonsense);
			}
		}
		p <- p + pp;
#		cat(t, '\t', pp, '\t', p, '\n');
	}
#	cat(' ', -p, '\n');
	return(-p);
}

#' Estimate selection coefficients of missense mutations and truncating mutations
#' @param obs_exp: A data frame with observed and expected counts
#' @param joint: A boolean flag to estimate joint probabilities
#' @param indel: A boolean flag to consider indel jointly with or independently from substitutions 
#' @return A vector with observed counts, end values and estimated probabilities
compute.selection <- function(obs_exp, joint=FALSE, indel=FALSE) {
	obs_exp[is.na(obs_exp)] <- 0;
	sums <- colSums(obs_exp[-1]);
	sum.obs <- sum(sums[c('syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs')]);
	sum.obs.disabling <- sum(sums[c('nonsense.obs', 'fs.obs')]);
	sum.exp <- sum(sums[c('syn.exp', 'missense.exp', 'nonsense.exp')]);
	opt.missense <- c(NA, NA);
	opt.nonsense <- c(NA, NA);
	opt.prot <- c(NA, NA);
	if(sum.obs >= 3 & sum.exp > 0) {
		if(joint) {
			opt.joint <- c();
			if(indel) {
				opt.joint <- unlist(optim(par=c(0, 0), fn=prob.joint.indel, obs_exp=obs_exp)); 
			} else {
				opt.joint <- unlist(optim(par=c(0, 0), fn=prob.joint, obs_exp=obs_exp)); 
			}
			opt.missense <- as.numeric(opt.joint[c(1,3)]);
			opt.nonsense <- as.numeric(opt.joint[c(2,3)]);
		} else {
			a <- sums[1]; b <- sums[2]; c <- sums[3]; d <- sums[4]; e <- sums[5];
			if(b > 1 & (a + b) >= 3) {
				opt.missense <- unlist(optimize(prob, interval=c(0, 100), obs_exp, 'missense')); 
			}
			if((c + d) > 1 & (a + c + d) >= 3) {
				opt.nonsense <- unlist(optimize(prob, interval=c(0, 100), obs_exp, 'nonsense')); 
			}
			if((b + c + d) > 1 & (a + b + c + d) >= 3) {
				opt.prot <- unlist(optimize(prob, interval=c(0, 100), obs_exp, 'prot')); 
			}
		}
	}
	opt <- c(sums[1:5], opt.missense, opt.nonsense, opt.prot);
	return(opt);
}

#' Get expected counts of variants for a vector of gene symbols
#' @param symbols.focus: A vector of gene symbols
#' @return A data frame with expected counts of variants
get.exp <- function(symbols.focus) {
	result <- c();
	for(s in symbols.focus) {
		id.focus <- map[which(map$symbol == s), 'cds.id']
		exp.focus <- results.seq_cat[which(results.seq_cat$cds.id == id.focus), ]
		exp.focus.cnt <- colSums(exp.focus[, 2:4])
		exp.focus.ratio <- c(exp.focus.cnt[2]/exp.focus.cnt[1], exp.focus.cnt[3]/exp.focus.cnt[1])
		result <- rbind(result, exp.focus.ratio);
	}
	rownames(result) <- seq(1, nrow(result));
	result <- as.data.frame(result);
	result$symbol <- symbols.focus;
	return(result);
}

#' Combine reciprocal mutation types
#' @param obs_exp: A data frame with observed and expected counts of variants
#' @return A data frame with observed and expected counts of variants after combining reciprocal mutation types
condense.type <- function(obs_exp) {
	AC.TG <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-C', 'T-G')), -1], na.rm=T);
	AG.TC <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-G', 'T-C')), -1], na.rm=T);
	AT.TA <- colSums(obs_exp[which(obs_exp$mut_type %in% c('A-T', 'T-A')), -1], na.rm=T);
	CA.GT <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-A', 'G-T')), -1], na.rm=T);
	CG.GC <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-G', 'G-C')), -1], na.rm=T);
	GA.CT <- colSums(obs_exp[which(obs_exp$mut_type %in% c('G-A', 'C-T')), -1], na.rm=T);
	CTh.GAh <- colSums(obs_exp[which(obs_exp$mut_type %in% c('C-T-h', 'G-A-h')), -1], na.rm=T);
	fs_ins <- colSums(obs_exp[which(obs_exp$mut_type %in% c('fs_ins')), -1]);
	fs_del <- colSums(obs_exp[which(obs_exp$mut_type %in% c('fs_del')), -1]);
	inf_ins <- colSums(obs_exp[which(obs_exp$mut_type %in% c('inf_ins')), -1]);
	inf_del <- colSums(obs_exp[which(obs_exp$mut_type %in% c('inf_del')), -1]);
	
	result <- rbind(AC.TG, AG.TC, AT.TA, CA.GT, CG.GC, GA.CT, CTh.GAh, fs_ins, fs_del, inf_ins, inf_del);
	result <- data.frame(mut_type=rownames(result), result, stringsAsFactors=F);
	result[is.na(result)] <- 0;
	return(result);
}

#' Estimate selection coefficients
#' @param data.all: A data frame with observed counts of variants
#' @param genes.opt.file.name: A file to write the computed selection coefficients.
#' @param map: A data frame containing symbol_2_cds_id mapping
#' @return A data frame with selection coefficients
get.opt.genes <- function(data.all, genes.opt.file.name='', indel=FALSE, map) {
	genes.all.opt <- c();
	counter <- 0;
	symbols.all <- unique(data.all$symbol);
	for(s in symbols.all) {
		counter <- counter + 1;
		if(counter %% 1000 == 0) { cat(counter, '...'); flush.console(); }
		data.s <- data.all[which(data.all$symbol == s), ]
		obs <- organize.obs(data.s);
		cds.id <- unique(map[which(map$symbol == s), 'cds.id'])[1];
		expected <- results.seq_cat[which(results.seq_cat$cds.id == cds.id), 1:4]
		total.cnt <- colSums(expected[, 2:4]);
		expected <- rbind(expected, data.frame(mut_type='fs_ins', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='fs_del', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='inf_ins', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		expected <- rbind(expected, data.frame(mut_type='inf_del', cnt_syn=total.cnt[1], cnt_nonsyn=total.cnt[2], cnt_nonsense=total.cnt[3]));
		obs_exp <- merge(obs, expected, by='mut_type', all=T)
		obs_exp <- condense.type(obs_exp);
		colnames(obs_exp) <- c('mut_type', 'syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs', 'syn.exp', 'missense.exp', 'nonsense.exp');

#		opt <- compute.selection(obs_exp)[c(1:4, 5, 7, 9)];	
		opt <- compute.selection(obs_exp, joint=T, indel=indel)[c(1:6, 8, 10)];	
		genes.all.opt <- rbind(genes.all.opt, opt);	
	}
	genes.all.opt <- as.data.frame(genes.all.opt);
	colnames(genes.all.opt) <- c('syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs', 'max.missense', 'max.nonsense', 'max.prot');
	rownames(genes.all.opt) <- symbols.all;
	genes.all.opt$symbol <- symbols.all;
	genes.all.opt <- merge(genes.all.opt, map, by='symbol', all.x=T, all.y=F);
	cat('\n'); flush.console();

	if(!is.null(genes.opt.file.name) && nchar(genes.opt.file.name) > 0) {
		write.table(genes.all.opt, genes.opt.file.name, sep='\t', row.names=F, quote=F);
	}
	return(genes.all.opt);
}

#' Compute statistics related to somatic selection
#' @param output.folder: A folder containing aggregated mutations produced by the parse.aggregated.mut() function.
#' @param output.prefix: Prefix used in the parse.aggregated.mut() function
#' @return NULL
#' @examples: compute.selection.stats(output.folder='./examples/', output.prefix='TCGA.ACC');
compute.selection.stats <- function(output.folder='', output.prefix='') {
	cat('computing somatic selection ...\n');
	working.dir <- paste(output.folder, '/', sep='');
	cat('    folder', output.folder, '...');
	mut.cnt.file.name <- paste(working.dir, output.prefix, '.mut.cnt.txt', sep='');
	mut.summary.file.name <- paste(working.dir, output.prefix, '.mut.summary.txt', sep='');
	genes.opt.file.name <- paste(working.dir, output.prefix, '.genes.opt.txt', sep='');
	map.file.name <- paste(working.dir, output.prefix, '.symbol_2_cds_id.txt', sep='');
	if(nchar(output.prefix) == 0) {
		mut.cnt.file.name <- paste(working.dir, 'mut.cnt.txt', sep='');
		mut.summary.file.name <- paste(working.dir, 'mut.summary.txt', sep='');
		genes.opt.file.name <- paste(working.dir, 'genes.opt.txt', sep='');
		map.file.name <- paste(working.dir, 'symbol_2_cds_id.txt', sep='');
	}
	genes.samples <- read.table(mut.cnt.file.name, header=T, sep='\t', stringsAsFactors=F);
	summary.mut <- read.table(mut.summary.file.name, header=T, sep='\t', stringsAsFactors=F);
	map <- read.table(map.file.name, header=T, sep='\t', stringsAsFactors=F);
	#map <- merge(map, gtex.skin.exp[, c('Description', 'avg.reads')], by.x='symbol', by.y='Description', all.x=T, all.y=F);
	#dim(map)  ## 17648     3  ## 19166     3
	map <- merge(map, seq.length, by='cds.id', all.x=T, all.y=F);
	data.all <- cbind(genes.samples[, 1:3], cnt=apply(genes.samples[, 4:ncol(genes.samples)], 1, sum, na.rm=T));
	cat(' genes.samples:', nrow(genes.samples), '; summary.mut:', nrow(summary.mut), '; map:', nrow(map), '; data.all:', nrow(data.all), '\n'); flush.console();

	cat('    processing ...'); flush.console();
	opts.result <- get.opt.genes(data.all, genes.opt.file.name, indel=T, map=map);
}

#' Detect mutational hotspots
#' @param x: a vector of mutated positions
#' @return A data frame with peak coordinates
find_peaks <- function (x, m=3, bw=bw){
	shape <- diff(sign(diff(x, na.pad = FALSE)))
	pks <- sapply(which(shape < 0), FUN = function(i){
		z <- i - m + 1
		z <- ifelse(z > 0, z, 1)
		w <- i + m + 1
		w <- ifelse(w < length(x), w, length(x))
		if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
	})
	pks <- unlist(pks)
	pks.end <- sapply(pks, function(pk) { 
		pk <- ifelse(pk - bw < 1, bw, ifelse(pk + bw > length(x), length(x) - bw, pk));
		a <- x[(pk-m):(pk+m)];
		if(min(a, na.rm=T) > 0) {
			x[(pk-m):(pk+m)] <- min(a, na.rm=T);
		}
		x.after <- x[pk:length(x)] - x[pk]; e <- min(which(x.after!=0)) + pk - m + 1; return(e)
	});
	pks <- ifelse(pks<=bw, 1, pks-bw);
	pks.end <- ifelse(pks.end + bw > length(x), length(x), pks.end+bw);
	return(data.frame(pk.start=pks, pk.end=pks.end));
}

#' Compute statistics related to mutational spatial distributions
#' @param output.folder: A folder containing a "mut.summary.txt" file produced by the parse.aggregated.mut() function.
#' @param output.prefix: A prefix string used label the various files. Spatial statistics is written to a file "prefix.spatial.txt".
#' @return NULL
#' @examples: compute.spatial.stat(folder='./examples/', prefix='TCGA.ACC');
compute.spatial.stats <- function(output.folder, output.prefix, bw=5, m=3) {
	cat('computing spatial distribution statistics ...\n'); flush.console();
	input.file.name <- paste(output.folder, output.prefix, '.mut.summary.txt', sep='');
	output.file.name <- paste(output.folder, output.prefix, '.spatial.txt', sep='');
	if(nchar(output.prefix) == 0) {
		input.file.name <- paste(output.folder, 'mut.summary.txt', sep='');
		output.file.name <- paste(output.folder, 'spatial.txt', sep='');
	}
	sample.mut <- read.table(input.file.name, sep='\t', header=T, stringsAsFactors=F);
	symbols <- unique(sample.mut$symbol);
	scores <- c();
	for(s in symbols) {
		subset.mut <- sample.mut[which(sample.mut$symbol == s), ]
		wt.stop <- 0;
		mut <- subset.mut[which(subset.mut$prot.alt == '*' | subset.mut$mut_type %in% c('fs_del', 'fs_ins')), ];
		if(nrow(mut) > 0) {
			wt.stop <- sum(mut$prot.pos * mut$chr)/sum(mut$chr)
		}

		mut <- subset.mut[which(subset.mut$prot.ref != subset.mut$prot.alt & nchar(subset.mut$mut_type) < 6 & subset.mut$prot.alt != '*'), ];
		scr.multiple.pct <- scr.max.pct <- NA;
		max.cnt <- max.peak.cnt <- multiple.cnt <- summit.start <- summit.end <- 0;
		if(nrow(mut) > 0) {
			agg.mut <- aggregate(chr ~ prot.pos, data=mut, 'sum');
			total.cnt <- sum(agg.mut$chr);
			max.cnt <- max.peak.cnt <- max(agg.mut$chr);
			summit.start <- summit.end <- agg.mut[which.max(agg.mut$chr), 'prot.pos'];
			single.cnt <- length(which(agg.mut$chr == 1));
			multiple.cnt <- sum(agg.mut[which(agg.mut$chr > 1), 'chr']);
			scr.multiple.pct <- multiple.cnt/total.cnt;
			scr.max.pct <- max.cnt/total.cnt;
			mut.distr <- unlist(apply(agg.mut, 1, function(x) { return(rep(x[1], x[2])) }));
			if(length(mut.distr) > 1 & bw > 0) {
				mut.density <- density(mut.distr, bw=bw, kernel='rectangular', n=max(mut.distr) + bw);
				xx <- round(mut.density$x);
				yy <- mut.density$y;
				yy <- round(yy/sum(yy), 4);
				xy <- data.frame(xx, yy);
				xy <- aggregate(yy ~ xx, data=xy, 'max');
				xx <- xy$xx;
				yy <- xy$yy;
				pks <- find_peaks(yy, m=m, bw=bw);
				pks$pk.start <- xx[pks$pk.start]
				pks$pk.end <- xx[pks$pk.end]
				pks.mut <- apply(pks, 1, function(x) { 
					x1 <- x[1];
					x2 <- x[2];
					aa <- ifelse(x1-2 < 1, 1, x1-2);
					bb <- ifelse(x2+2 > max(xx), max(xx), x2+2);
					i <- which(mut.distr >= aa & mut.distr <= bb); 
					if (length(i) > 1) return(i); 
				})
				pks.mut.cnt <- unlist(lapply(pks.mut, length));
				if(length(which(pks.mut.cnt > 0)) > 0) {
					max.peak.cnt <- max(pks.mut.cnt);
					if(max.peak.cnt > max.cnt) {
						scr.max.pct <- max.peak.cnt/total.cnt;
					}
					if(max.peak.cnt > 1) {
						coord <- mut.distr[pks.mut[[which.max(pks.mut.cnt)]]];
#						coord <- as.numeric(gsub('.prot.pos', '', coord))
						summit.start <- min(coord);
						summit.end <- max(coord);
					}
				}
				pks.mut.uniq <- unique(unlist(pks.mut))
				pks.multiple.cnt <- length(pks.mut.uniq);
				if(pks.multiple.cnt > multiple.cnt) {
					multiple.cnt <- pks.multiple.cnt;
					scr.multiple.pct <- multiple.cnt/total.cnt;
				}
			}
		}
		scores <- rbind(scores, c(max.cnt, scr.multiple.pct, scr.max.pct, max.peak.cnt, multiple.cnt, summit.start, summit.end, wt.stop));
	}
	scores <- data.frame(symbol=symbols, score=scores);
	write.table(scores, output.file.name, sep='\t', row.names=F, quote=F);	
}

#' Retreve fitch rate for each position of a protein.
#' @param id: An ID indicating the gene.
#' @param type: A flag indicating what type the "id" is. Possible values are "symbol" (gene symbol), "gene.id" (Entrez gene ID), "ensembl.cds" (Ensemble transcript ID), "ensembl.protein" (Ensemble protein ID), and "ensembl.gene" (Ensembl gene ID).
#' @return A vector containing fitch rate using 100 vertebrate alignments.
check.gene.consv.all <- function(id, type='symbol') {
	refseq.id <- id;
	if (type == 'symbol') {
		refseq.id <- unique(refseq_map[which(refseq_map$symbol == id), 'refseq.2'])
	} else if (type == 'gene.id') {
		refseq.id <- unique(refseq_map[which(refseq_map$gene.id == id), 'refseq.2'])
	} else if (type == 'ensembl.cds') {
		refseq.id <- unique(refseq_map[which(refseq_map$cds.id == id), 'refseq.2'])
	} else if (type == 'ensembl.protein') {
		refseq.id <- unique(refseq_map[which(refseq_map$protein.id == id), 'refseq.2'])
	} else if (type == 'ensembl.gene') {
		refseq.id <- unique(refseq_map[which(refseq_map$ensembl.gene == id), 'refseq.2'])
	}
	refseq.id <- refseq.id[1];
	k <- ref.block[which(ref.block$refseq == refseq.id), 'k']
	if(length(k) > 0) {
		fitch.k <- fitch.block[[k]]
		gene.fitch <- fitch.k[which(fitch.k$refseq == refseq.id), ]
		return(gene.fitch);
	}
}

#' Compute statistics related to evolutionary conservation
#' @param output.folder: A folder containing a "mut.summary.txt" file and a 'spatial.txt' file.
#' @param output.prefix: A prefix string used label the various files. Conservation statistics is written to a file "prefix.fitch.txt".
#' @return NULL
#' @examples: compute.fitch.stats(output.folder='./examples/', output.prefix='TCGA.ACC');
compute.fitch.stats <- function(output.folder, output.prefix) {
	cat('computing conservational statistics ...\n'); flush.console();
	mut.summary.file.name <- paste(output.folder, output.prefix, '.mut.summary.txt', sep='');
	spatial.file.name <- paste(output.folder, output.prefix, '.spatial.txt', sep='');
	output.file.name <- paste(output.folder, output.prefix, '.fitch.txt', sep='');
	if(nchar(output.prefix) == 0) {
		mut.summary.file.name <- paste(output.folder, 'mut.summary.txt', sep='');
		spatial.file.name <- paste(output.folder, 'spatial.txt', sep='');
		output.file.name <- paste(output.folder, 'fitch.txt', sep='');
	}

	mut.summary <- read.table(mut.summary.file.name, sep='\t', header=T, stringsAsFactors=F);
	mut.in <- mut.summary[which(mut.summary$prot.ref != mut.summary$prot.alt & nchar(mut.summary$mut_type) < 6 & mut.summary$prot.alt != '*'), ]
	agg.in <- aggregate(chr ~ symbol + prot.pos, data=mut.in, 'sum')
	agg.in.list <- split(agg.in, agg.in$symbol);
	consv <- do.call(rbind, lapply(agg.in.list, function(x) {
		consv.all <- check.gene.consv.all(id=unique(x$symbol), type='symbol');
		return(consv.all[which(consv.all$pos %in% x$prot.pos), c('refseq', 'pos', 'fitchratevertebrate')]);
	}));
	consv$symbol <- rownames(consv);	
	consv$symbol <- gsub('\\.[0-9]+', '', consv$symbol);
	colnames(consv) <- c('refseq', 'prot.pos', 'fitchratevertebrate', 'symbol');
	
	agg.anno <- merge(agg.in, consv, by=c('symbol', 'prot.pos'), all=T);
	bb <- unique(agg.anno);
	bb$fitch <- bb$fitchratevertebrate * bb$chr
	agg.bb <- aggregate(fitch ~ symbol, data=bb, 'sum', na.rm=T)
	agg.ss <- aggregate(chr ~ symbol, data=bb, 'sum', na.rm=T)
	agg.bs <- unique(merge(agg.bb, agg.ss, by='symbol'))
	agg.ff <- data.frame(symbol=agg.bs$symbol, fitch=agg.bs$fitch/agg.bs$chr)

	spatial <- read.table(spatial.file.name, sep='\t', header=T, stringsAsFactors=F);
	spatial <- spatial[which(spatial$score.6 > 0), ]
	spatial.bb <- merge(bb, spatial, by='symbol', all=F);
	spatial.bb <- spatial.bb[which(spatial.bb$prot.pos >= spatial.bb$score.6 & spatial.bb$prot.pos <= spatial.bb$score.7), ]
	agg.bb <- aggregate(fitch ~ symbol + score.6, data=spatial.bb, 'sum', na.rm=T)
	agg.ss <- aggregate(chr ~ symbol + score.6, data=spatial.bb, 'sum', na.rm=T)
	agg.bs <- unique(merge(agg.bb, agg.ss, by='symbol'))
	agg.spatial <- data.frame(symbol=agg.bs$symbol, fitch.summit=agg.bs$fitch/agg.bs$chr)
	
	agg <- merge(agg.ff, agg.spatial, by='symbol', all=T);
	write.table(agg, output.file.name, sep='\t', row.names=F, quote=F);
}

#' Merge selection coefficients, spatial concentration, conservation scores and gene expression data
#' @param output.folder: A folder containing a "genes.opt.txt" file and a 'spatial.txt' file, a 'fitch.txt' file.
#' @param output.prefix: A prefix string used label the various files. Combined data frame is written to a file "prefix.combined.txt".
#' @return NULL
#' @examples: append.info(output.folder='./examples/', output.prefix='TCGA.ACC');
append.info <- function(output.folder, output.prefix) {
	cat('appending somatic selection, spatial distribution and conservation stats ...\n'); flush.console();
	genes.opt.file.name <- paste(output.folder, output.prefix, '.genes.opt.txt', sep='');
	spatial.file.name <- paste(output.folder, output.prefix, '.spatial.txt', sep='');
	fitch.file.name <- paste(output.folder, output.prefix, '.fitch.txt', sep='');
	output.file.name <- paste(output.folder, output.prefix, '.combined.txt', sep='');
	if(nchar(output.prefix) == 0) {
		genes.opt.file.name <- paste(output.folder, 'genes.opt.txt', sep='');
		spatial.file.name <- paste(output.folder, 'spatial.txt', sep='');
		fitch.file.name <- paste(output.folder, 'fitch.txt', sep='');
		output.file.name <- paste(output.folder, 'combined.txt', sep='');
	}
	genes.opt <- read.table(genes.opt.file.name, sep='\t', header=T);
	cat('    selection data (', nrow(genes.opt), ')\n'); flush.console();
	genes.spatial <- read.table(spatial.file.name, sep='\t', header=T, stringsAsFactors=F);
	cat('    spatial data (', nrow(genes.spatial), ')\n'); flush.console();
	genes.fitch <- read.table(fitch.file.name, sep='\t', header=T, stringsAsFactors=F);
	cat('    conservation data (', nrow(genes.fitch), ')\n'); flush.console();

	genes.opt$total.mut.cnt <- rowSums(genes.opt[, c('syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs')]);
	genes.opt$prot.mut.cnt <- rowSums(genes.opt[, c('missense.obs', 'nonsense.obs', 'fs.obs')]);
	genes.opt$total.mut.pct <- genes.opt$total.mut.cnt/genes.opt$length*100;
	genes.opt$prot.mut.pct <- genes.opt$prot.mut.cnt/genes.opt$length*100;
	
	colnames(genes.spatial) <- c('symbol', 'score.max_cnt', 'score.mult_pct', 'score.max_pct', 'max.peak.cnt', 'multiple.cnt', 'summit.start', 'summit.end', 'wt.stop');

	genes.opt.fitch <- Reduce(function(...) merge(..., by='symbol', all=T, stringsAsFactors=F), list(genes.opt, genes.spatial, genes.fitch));
	genes.opt.fitch <- unique(genes.opt.fitch);
	genes.opt.fitch <- genes.opt.fitch[which(genes.opt.fitch$total.mut.cnt > 0), ];
	genes.opt.fitch$fitch.summit <- ifelse(genes.opt.fitch$max.peak.cnt > 1, genes.opt.fitch$fitch.summit, genes.opt.fitch$fitch);
	genes.opt.fitch$max.peak.cnt <- ifelse(genes.opt.fitch$max.peak.cnt == 1, 0, genes.opt.fitch$max.peak.cnt);
	write.table(genes.opt.fitch, output.file.name, sep='\t', row.names=F, quote=F);
}

#' Correct biases for genes with too few somatic variants.
#' @param opts.scs: A data frame containing all features.
#' @return A data frame containing all features after correction
correct.small.sample.bias <- function(opts.scs) {
	# max.nonsense check
	check.idx <- which(opts.scs$syn.obs == 0 & opts.scs$inframe.obs == 0 & opts.scs$missense.obs > opts.scs$cnt.disabling*10 & opts.scs$max.nonsense > 2);
	opts.scs[check.idx, 'max.nonsense'] <- -5;
#	cat(' nonsense... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$fs.obs >= 10 & (opts.scs$inframe.obs == 0 | opts.scs$fs.obs/opts.scs$inframe.obs >= 10 & opts.scs$max.nonsense < 4));
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 2;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$nonsense.obs >= 10 & (opts.scs$syn.obs == 0 | opts.scs$nonsense.obs/opts.scs$syn.obs >= 10 & opts.scs$max.nonsense < 4));
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 2;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$nonsense.obs >= 10 & opts.scs$nonsense.obs > opts.scs$syn.obs/2 & opts.scs$max.nonsense < 4);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$syn.obs > 0 & opts.scs$nonsense.obs/opts.scs$syn.obs >= 5 & opts.scs$max.nonsense < 4);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 2;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$inframe.obs > 0 & opts.scs$fs.obs <= opts.scs$inframe.obs*3 & opts.scs$nonsense.obs < opts.scs$syn.obs*2 & opts.scs$max.nonsense > 2);
	opts.scs[check.idx, 'max.nonsense'] <- 0;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$nonsense.obs >= (opts.scs$syn.obs + opts.scs$missense.obs)/3 & opts.scs$max.nonsense < -2);
	opts.scs[check.idx, 'max.nonsense'] <- 0;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.disabling >= opts.scs$total.mut.cnt/5 & opts.scs$max.nonsense < -2);
	opts.scs[check.idx, 'max.nonsense'] <- 0;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.disabling == 1 & opts.scs$total.mut.cnt >= 9 & opts.scs$max.nonsense > 0);
	opts.scs[check.idx, 'max.nonsense'] <- -4;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.disabling == 1 & opts.scs$max.nonsense > 0);
	opts.scs[check.idx, 'max.nonsense'] <- 0;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.disabling <= 2 & opts.scs$total.mut.cnt > 25 & opts.scs$max.nonsense > -1 & opts.scs$max.nonsense < 0);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] - 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$fs.obs == 0 & opts.scs$inframe.obs == 1 & opts.scs$max.nonsense > 1 & opts.scs$max.nonsense < 2);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$fs.obs >= 10 & opts.scs$inframe.obs == 0 & opts.scs$max.nonsense > 2.5 & opts.scs$max.nonsense < 3);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$fs.obs/opts.scs$inframe.obs >= 20 & opts.scs$max.nonsense > 2.5 & opts.scs$max.nonsense < 3);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$nonsense.obs > opts.scs$missense.obs & opts.scs$max.nonsense > 1 & opts.scs$max.nonsense < 2);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$disabling.cnt == 1 & opts.scs$total.mut.cnt > 5 & opts.scs$max.nonsense > 2);
	opts.scs[check.idx, 'max.nonsense'] <- 0;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.neutral >= opts.scs$cnt.disabling & opts.scs$max.nonsense > 2);
	opts.scs[check.idx, 'max.nonsense'] <- -5;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.neutral == 1 & opts.scs$cnt.disabling >= 10 & opts.scs$max.nonsense < 4);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 2;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.neutral >= 1 & opts.scs$cnt.disabling/opts.scs$cnt.neutral >= 5 & opts.scs$max.nonsense < 3);
	opts.scs[check.idx, 'max.nonsense'] <- opts.scs[check.idx, 'max.nonsense'] + 1.5;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$cnt.neutral >= 1 & opts.scs$cnt.disabling/opts.scs$cnt.neutral >= 5 & opts.scs$max.nonsense <= 0);
	opts.scs[check.idx, 'max.nonsense'] <- 3;
#	cat(' ... convert', length(check.idx), '\n'); flush.console();
	opts.scs[which(opts.scs$max.nonsense > 5), 'max.nonsense'] <- 5;
	
	# max.missense check	
	check.idx <- which(opts.scs$syn.obs >= 1 & opts.scs$missense.obs > opts.scs$syn.obs*10 & (opts.scs$nonsense.obs + opts.scs$fs.obs) <= 1 & opts.scs$max.missense < 2);
	opts.scs[check.idx, 'max.missense'] <- 4;	
#	cat(' missense ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$syn.obs >= 1 & opts.scs$missense.obs >= opts.scs$syn.obs*9 & opts.scs$max.missense < 2);
	opts.scs[check.idx, 'max.missense'] <- opts.scs[check.idx, 'max.missense'] + 1;	
#	cat(' missense ... convert', length(check.idx), '\n'); flush.console();
	check.idx <- which(opts.scs$max.peak.cnt >= 10 & opts.scs$score.summit > 0.4 & opts.scs$max.missense < 3);
	opts.scs[check.idx, 'max.missense'] <- 3.5;	
#	cat(' missense ... convert', length(check.idx), '\n'); flush.console();
	
	return(opts.scs);
}

#' Reformat a data frame to be suitable for modeling.
#' @param output.folder: A folder containing a "combined.txt" file.
#' @param output.prefix: A prefix string used to label the files. Reformatted data frame is written to a file "prefix.reformatted.txt".
#' @return A data frame containing all features after reformat
#' @examples: reformat(output.folder='./examples/', output.prefix='TCGA.ACC');
reformat <- function(output.folder, output.prefix) {
	cat('reformatting ...\n');
	input.file.name <- paste(output.folder, output.prefix, '.combined.txt', sep='');
	output.file.name <- paste(output.folder, output.prefix, '.reformatted.txt', sep='');
	if(nchar(output.prefix) == 0) {
		input.file.name <- paste(output.folder, 'combined.txt', sep='');
		output.file.name <- paste(output.folder, 'reformatted.txt', sep='');
	}
	
	opts.scs <- read.table(input.file.name, sep='\t', header=T, stringsAsFactors=F);
	
	opts.scs <- opts.scs[which(!is.na(opts.scs$syn.obs)), ];
	opts.scs <- opts.scs[which(!is.na(opts.scs$length)), ];
	opts.scs <- correct.small.sample.bias(opts.scs);
	opts.scs$length <- opts.scs$length/9;
	opts.scs$cnt.disabling <- rowSums(opts.scs[, c('nonsense.obs', 'fs.obs')]);
	opts.scs$fitch <- ifelse(is.na(opts.scs$fitch), mean(opts.scs$fitch, na.rm=T), opts.scs$fitch);
	opts.scs$fitch <- ifelse(opts.scs$fitch < 0.01, 0, opts.scs$fitch);
	opts.scs$fitch.summit <- ifelse(is.na(opts.scs$fitch.summit), opts.scs$fitch, opts.scs$fitch.summit);
	opts.scs$fitch <- round(opts.scs$fitch, 3);
	opts.scs$fitch.summit <- round(opts.scs$fitch.summit, 3);
	opts.scs$wt.stop.pct <- opts.scs$wt.stop/opts.scs$length;
	opts.scs$wt.stop.pct <- ifelse(opts.scs$wt.stop.pct > 1, 1, opts.scs$wt.stop.pct);
	opts.scs$max.missense <- ifelse(is.na(opts.scs$max.missense), 0, opts.scs$max.missense);
	opts.scs$max.nonsense <- ifelse(is.na(opts.scs$max.nonsense), 0, opts.scs$max.nonsense);
	opts.scs$max.missense <- round(opts.scs$max.missense, 2);
	opts.scs$max.nonsense <- round(opts.scs$max.nonsense, 2);
	
	colnames(opts.scs) <- ifelse(colnames(opts.scs)=='score.mult_pct', 'score.peak', ifelse(colnames(opts.scs)=='score.max_pct', 'score.summit', colnames(opts.scs)))
	opts.scs$score.max_cnt <- ifelse(is.na(opts.scs$score.max_cnt), 0, opts.scs$score.max_cnt);
	
	opts.scs[which(is.na(opts.scs$score.peak)), 'score.peak'] <- 0;
	opts.scs[which(is.na(opts.scs$score.summit)), 'score.summit'] <- 0;
	opts.scs$max.prot <- rowSums(opts.scs[, c('max.missense', 'max.nonsense')], na.rm=T);
	opts.scs <- opts.scs[which(opts.scs$total.mut.cnt >= 5 & !is.na(opts.scs$max.missense) & !is.na(opts.scs$max.missense) & !is.na(opts.scs$max.prot)), ]

	opts.scs$syn.obs.pct <- opts.scs$syn.obs/opts.scs$total.mut.cnt;
	opts.scs$missense.obs.pct <- opts.scs$missense.obs/opts.scs$total.mut.cnt;
	opts.scs$truncating.obs.pct <- (opts.scs$nonsense.obs + opts.scs$fs.obs)/opts.scs$total.mut.cnt;
	opts.scs$wt.stop.pct <- opts.scs$wt.stop/opts.scs$length;
	
	opts.scs <- merge(opts.scs, anno.length[, c('ensembl_transcript_id', 'ensembl.gene.id')], by.x='cds.id', by.y='ensembl_transcript_id', all.x=T, all.y=F);

	write.table(opts.scs, output.file.name, sep='\t', row.names=F, quote=F);
}

#' Predict oncogenes, tumor suppressor genes and passenger genes.
#' @param output.folder: A folder containing a "reformatted.txt" file.
#' @param output.prefix: A prefix string used to label the files. Predictions are written to a file "prefix.predictions.txt".
#' @return NULL
#' @examples: make.prediction(output.folder='./examples/', output.prefix='TCGA.ACC');
make.prediction <- function(output.folder, output.prefix) {
	cat('making predictions ...\n');
	input.file.name <- paste(output.folder, output.prefix, '.reformatted.txt', sep='');
	output.file.name <- paste(output.folder, output.prefix, '.predictions.txt', sep='');
	if(nchar(output.prefix) == 0) {
		input.file.name <- paste(output.folder, 'reformatted.txt', sep='');
		output.file.name <- paste(output.folder, 'predictions.txt', sep='');
	}
	all.data <- read.table(input.file.name, sep='\t', header=T, stringsAsFactors=F);
	col.sub <- c('missense.obs.pct', 'truncating.obs.pct', 'max.missense', 'max.nonsense', 'score.peak', 'score.summit', 'max.peak.cnt', 'wt.stop.pct', 'fitch', 'fitch.summit')

	predictors.all <- all.data
	na.row = which(is.na(predictors.all), arr.ind=T)
	if(length(na.row) > 0) predictors.all = predictors.all[-unique(na.row[, 1]), ];
	classes.pred <- predict(md, predictors.all[, col.sub])
	classes.pred.prob <- predict(md, predictors.all[, col.sub], 'prob')
	classes.pred.prob.max <- apply(classes.pred.prob, 1, max)
	colnames(classes.pred.prob) <- c('prob.CIG', 'prob.OG', 'prob.TSG');
	predictors.all$pred <- classes.pred;
	predictors.all$pred <- as.character(predictors.all$pred);
	predictors.all$pred <- ifelse(predictors.all$pred == 'CIG', 'PG', predictors.all$pred);
	predictors.all <- cbind(predictors.all, classes.pred.prob, 'prob.pred'=classes.pred.prob.max)
	predictors.all <- merge(predictors.all, cgc.drivers, by='symbol', all.x=T, all.y=F);
		
	predictors.all.orig <- predictors.all;
	dim(predictors.all.orig);  ## 91889
	predictors.all.orig$cgc.cat <- ifelse(is.na(predictors.all.orig$cgc.cat), 'passenger', as.character(predictors.all.orig$cgc.cat));
	predictors.all.orig$tumor <- gsub('TCGA.', '', output.prefix);
	predictors.all.orig <- merge(unique(annotation[, c('symbol', 'gene.id')]), predictors.all.orig, by='symbol', all.x=F, all.y=T);

	predictors.output <- predictors.all.orig[, c('symbol', 'ensembl.gene.id', 'gene.id', 'length', 'tumor', 'syn.obs', 'missense.obs', 'nonsense.obs', 'fs.obs', 'inframe.obs', 'max.missense', 'max.nonsense', 'missense.obs.pct',   'truncating.obs.pct', 'score.peak', 'score.summit', 'max.peak.cnt', 'wt.stop.pct', 'fitch', 'fitch.summit', 'cgc.cat', 'prob.CIG', 'prob.OG', 'prob.TSG', 'pred', 'prob.pred')]
	colnames(predictors.output) <- c('Symbol', 'Ensembl.gene', 'Entrez.gene', 'Gene.length', 'Tumor', 'Silent', 'Missense', 'Nonsense', 'Indel.fs', 'Indel.inframe', 'Sel.missense', 'Sel.truncating', 'R.missense', 'R.truncating', 'R.peak', 'R.summit', 'C.summit', 'R.length', 'E.gene', 'E.summit', 'CGC.cat', 'PG.prob', 'OG.prob', 'TSG.prob', 'GUST.pred', 'GUST.prob')
	predictors.output$R.missense <- round(predictors.output$R.missense, 2)
	predictors.output$R.truncating <- round(predictors.output$R.truncating, 2)
	predictors.output$R.peak <- round(predictors.output$R.peak, 2)
	predictors.output$R.summit <- round(predictors.output$R.summit, 2)
	predictors.output$R.length <- round(predictors.output$R.length, 2)
	write.table(predictors.output, output.file.name, sep='\t', row.names=F, quote=F);
}

#' Retreve somatic mutations of a gene from one or more tumor types.
#' @param gene.symbol: A gene symbol.
#' @param folder: A folder containing a "prefix.mut.summary.txt" file and a "prefix.symbol_2_cds_id.txt" file.
#' @param prefix: A prefix string used to label the files. Usually, it contains TCGA.tumor_abbreviation.
#' @return A list with three elements: the first element is a data frame with mutations counts aggregated on protein level; the second element is a data frame with mutations counts aggregated on nucleotide level; the third element is the Ensembl transcript ID.
check.gene.mut <- function(gene.symbol, folder, prefix) {
	mut.summary <- read.table(paste(folder, '/', prefix, '.mut.summary.txt', sep=''), header=T, sep='\t', stringsAsFactors=F);
	mut.summary.s <- mut.summary[which(mut.summary$symbol == gene.symbol), ];
	mut.summary.s.neutral <- mut.summary.s[which(mut.summary.s$mut_type %in% c('inf_del', 'inf_ins') | mut.summary.s$prot.ref == mut.summary.s$prot.alt), ]
	mut.summary.s.nn <- setdiff(mut.summary.s, mut.summary.s.neutral);
	
	if(nrow(mut.summary.s.neutral) > 0) {
		mut.summary.s.neutral$type <- 'neutral'
	}
	if(nrow(mut.summary.s.nn) > 0) {
		mut.summary.s.nn$type <- ifelse(mut.summary.s.nn$mut_type %in% c('fs_del', 'fs_ins') | mut.summary.s.nn$prot.alt == '*', 'nonsense', 'missense'); ## account for diff prot.alt at a pos
	}
	
	mut.summary.s <- rbind(mut.summary.s.neutral, mut.summary.s.nn);
	mut.summary.s <- mut.summary.s[order(mut.summary.s$prot.pos), ]
	result <- c();
	if(nrow(mut.summary.s) > 0) {
		agg <- aggregate(chr ~ prot.pos + prot.ref + prot.alt + type, data=mut.summary.s, 'sum')
		agg <- agg[order(-agg$chr), ];
		agg$tumor <- gsub('TCGA.', '', prefix);
		result <- rbind(result, agg)
	}
	
	colnames(result) <- c('prot.pos', 'prot.ref', 'prot.alt', 'type', 'cnt', 'tumor');
	result <- result[order(result$tumor, result$prot.pos), ]

	aa.length <- max(anno.length[which(anno.length$symbol == gene.symbol), 'length']);
	result$length <- aa.length;
	
	map <- read.table(paste(folder, '/', prefix, '.symbol_2_cds_id.txt', sep=''), header=T, sep='\t', stringsAsFactors=F);
	cds.id <- map[map$symbol==gene.symbol, 'cds.id'][1];

	return(list(result, mut.summary.s, cds.id));
}

#' Perform GUST classification of oncogens, tumor suppressor genes and passenger genes from exome sequencing data.
#' @param input.file.name: A VCF-formatted file to read SnpEff annotated somatic variants. This input file shall contain these fields: Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, Entrez_Gene_Id.
#' @param output.folder: the folder to write the temporary and final prediction files. 
#' @param output.prefix: prefix used to name the temporary and final prediction files. Temporary files include prefix.outlier.txt, prefix.mut.cnt.txt, prefix.mut.filtered.txt, prefix.mut.summary.txt, prefix.error.txt, prefix.log.txt, prefix.symbol_2_cds_id.txt, prefix.genes.opt.txt, prefix.spatial.txt, prefix.fitch.txt, prefix.combined.txt, prefix.reformatted.txt, prefix.gust.txt
#' @param steps: A vector of integers indicating which functions to execute. 1-find.outliers(), 2-parse.aggregated.mut(), 3-compute.selection.stats, 4-compute.spatial.stats(), 5-compute.fitch.stats(), 6-append.info(), 7-reformat(), 8-make.prediction(). Default to 1:8.
#' @return NULL
#' @examples gust(input.file.name='./examples/TCGA.ACC.mutect.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.ACC');
gust <- function(input.file.name, output.folder, output.prefix, steps=1:8) {
	output.folder <- paste(output.folder, '/', sep='');
	initialize.gust();
	if(1 %in% steps) {
		find.outliers(input.file.name=input.file.name, output.folder=output.folder, output.prefix=output.prefix);
	}
	if(2 %in% steps) {
		parse.aggregated.mut(input.file.name=input.file.name, output.folder=output.folder, output.prefix=output.prefix);
	}
	if(3 %in% steps) {
		compute.selection.stats(output.folder=output.folder, output.prefix=output.prefix);
	}
	if(4 %in% steps) {
		compute.spatial.stats(output.folder=output.folder, output.prefix=output.prefix);
	}
	if(5 %in% steps) {
		compute.fitch.stats(output.folder=output.folder, output.prefix=output.prefix);
	}
	if(6 %in% steps) {
		append.info(output.folder=output.folder, output.prefix=output.prefix)
	}
	if(7 %in% steps) {
		reformat(output.folder=output.folder, output.prefix=output.prefix)
	}
	if(8 %in% steps) {
		make.prediction(output.folder=output.folder, output.prefix=output.prefix)
	}
}
