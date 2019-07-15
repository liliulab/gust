#' Plot distribution of somatic mutations.
#' @param mut.in.tumor: A list returned from the check.gene.mut() function.
#' @param par: A boolean indicating if using current grid arrangement. Default is TRUE.
#' @param bw: bandwidth for density estimate. Default is 5.
#' @param main: title of the plot.
#' @param side: side to plot the protein position.
#' @return NULL
plot.mut <- function(mut.in.tumor, par=T, bw=5, main='', side=1) {
	cds.id <- mut.in.tumor[[3]];
	mut <- mut.in.tumor[[1]]
	max.pos <- max(mut$length);
	if(is.na(max.pos)) {
		max.pos <- max(mut$prot.pos) + 10;
	}
	tumors <- unique(mut$tumor);
	if(length(tumors) > 1) {
		main <- tumors;
	}
	par(mar=c(5,5,5,5));
	if(par) {
		par(mfrow=c(2, ceiling(length(tumors)/2)));
	}

	for(t in 1:length(tumors)) {
		this <- mut[which(mut$tumor == tumors[t]), ];
		y <- max(this$cnt) + 5;
		plot(0, type='n', xlim=c(1, max.pos), ylim=c(1, y), xlab='', ylab='Frequency', main=main[t], cex.main=2, cex.axis=1.5, cex.lab=1.5, las=1, xaxt='n');
		axis(side, cex.axis=1.5);
		mtext(side=side, line=2.5, 'Protein Position', cex=1.5); ## , cex=1.5
		
		missense <-  this[which(this$type == 'missense'), ];
		if(nrow(missense) > 0) {
			missense <- aggregate(cnt ~ prot.pos, data=missense, 'sum');
			points(x=missense$prot.pos, y=missense$cnt, type='h', col='blue', lwd=2);
			points(x=missense$prot.pos, y=missense$cnt, pch=19, col='blue');
		}
		nonsense <- this[which(this$type == 'nonsense'), ];
		if(nrow(nonsense) > 0) {
			nonsense <- aggregate(cnt ~ prot.pos, data=nonsense, 'sum');
			points(x=nonsense$prot.pos, y=nonsense$cnt, type='h', col='red', lwd=2);
			points(x=nonsense$prot.pos, y=nonsense$cnt, pch=4, col='red');
		}
		neutral <- this[which(this$type == 'neutral'), ];
		if(nrow(neutral) > 0) {
			neutral <- aggregate(cnt ~ prot.pos, data=neutral, 'sum');
			points(x=neutral$prot.pos, y=neutral$cnt, type='h', col='green', lwd=2);
			points(x=neutral$prot.pos, y=neutral$cnt, pch=3, col='green');
		}
		
		mut.distr <- as.numeric(unlist(apply(this[which(this$type != 'neutral'), ], 1, function(x) { return(rep(x[1], x[5])) })));
		mut.density <- density(mut.distr, bw=bw, kernel='rectangular', n=max(mut.distr) + bw);
		y = mut.density$y;
		y = round(y/sum(y), 3);
		x = mut.density$x;
		x = c(x, max.pos);
		y = c(y, 0);
		par(new=T);
		plot(x, y*100, type='l', col=alpha('gray30', 0.4), xlim=c(1, max.pos), axes=F, xlab=NA, ylab=NA, lwd=2);
		axis(4, col.axis='gray30', col.ticks='gray30', cex.axis=1.3, las=1)
		mtext(side=4, line=3, 'Density', col='gray30', cex=1.5)
	}	
}

#' Plot distribution of somatic mutations of a gene.
#' @param symbol: A gene symbol.
#' @param folder: A folder containing a "prefix.mut.summary.txt" file and a "prefix.symbol_2_cds_id.txt" file.
#' @param prefix: A prefix string used to label the files. Usually, it contains TCGA.tumor_abbreviation.
#' @param par: A boolean indicating if using current grid arrangement. Default is TRUE.
#' @param bw: bandwidth for density estimate. Default is 5.
#' @param main: title of the plot.
#' @param side: side to plot the protein position.
#' @return A data frame from the check.gene.mut() function
#' @examples m=plot.this.one('TP53', './examples/', 'TCGA.ACC')
plot.this.one <- function(symbol, folder, prefix, par=F, bw=5, main='', side=1) {
	mut.in.tumor <- check.gene.mut(symbol, folder, prefix)
	if(nchar(main)==0) {
		main <- paste(symbol, gsub('TCGA.', '', prefix), sep=':');
	}
	plot.mut(mut.in.tumor, par=F, bw=bw, main=main, side=side)
	return(mut.in.tumor);
}

#' Plot distribution of somatic mutations of multiple genes.
#' @param multiple: A data frame. Each row contains a gene symbol, a folder, a prefix and an optional class label.
#' @param par: A boolean indicating if using current grid arrangement. Default is TRUE.
#' @param bw: bandwidth for density estimate. Default is 5.
#' @param main: title of the plot.
#' @param side: side to plot the protein position.
#' @return A data frame from the check.gene.mut() function
#' @examples df=data.frame(rbind(c('TP53', './examples', 'TCGA.ACC', 'TSG'), c('CTNNB1', './examples/', 'TCGA.ACC', 'OG'), c('TTN', './examples/', 'TCGA.ACC', 'PG'))); m=plot.multiple(df)
plot.multiple <- function(multiple, par=F, bw=5, side=1) {
	n <- nrow(multiple);
	muts <- list();
	par(mfrow=c(3, ceiling(n/3)), mar=rep(1,4));
#	par(mfrow=c(3, 5), mar=rep(1,4));
	for(i in 1:n) {
		s <- multiple[i, 1];
		folder <- multiple[i, 2];
		prefix <- multiple[i, 3];
		label <- multiple[i, 4];
		mut.in.tumor <- check.gene.mut(s, folder, prefix)
		plot.mut(mut.in.tumor, par=F, bw=bw, main=paste(s, gsub('TCGA.', '', prefix), label, sep=':'), side=side)
		muts[[i]] <- mut.in.tumor
	}
	return(muts);
}
