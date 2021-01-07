%% Identifying Differentially Expressed Genes from RNA-Seq Data
% This example shows how to test RNA-Seq data for differentially expressed
% genes using a negative binomial model.

% Copyright 2010-2016 The MathWorks, Inc.

%% Introduction
% A typical differential expression analysis of RNA-Seq data consists of
% normalizing the raw counts and performing statistical tests to reject or
% accept the null hypothesis that two groups of samples show no significant
% difference in gene expression. This example shows how to inspect the
% basic statistics of raw count data, how to determine size factors for
% count normalization and how to infer the most differentially expressed
% genes using a negative binomial model.
%
% The dataset for this example comprises of RNA-Seq data obtained in the
% experiment described by Brooks et al. [1]. The authors investigated the
% effect of siRNA knock-down of pasilla, a gene known to play an important
% role in the regulation of splicing in _Drosophila melanogaster_. The
% dataset consists of 2 biological replicates of the control (untreated)
% samples and 2 biological replicates of the knock-down (treated) samples.

%% Inspecting Read Count Tables for Genomic Features
% The starting point for this analysis of RNA-Seq data is a count matrix,
% where the rows correspond to genomic features of interest, the columns
% correspond to the given samples and the values represent the number of
% reads mapped to each feature in a given sample. 
% 
% The included file |pasilla_count_noMM.mat| contains two tables with the
% count matrices at the gene level and at the exon level for each of the
% considered samples. You can obtain similar matrices using the function
% |featurecount|.

load pasilla_count_noMM.mat

%%

% preview the table of read counts for genes
geneCountTable(1:10,:)

%% 
% Note that when counting is performed without summarization, the
% individual features (exons in this case) are reported with their
% metafeature assignment (genes in this case) followed by the start and
% stop positions.

% preview the table of read counts for exons
exonCountTable(1:10,:)

%%
% You can annotate and group the samples by creating a logical vector as
% follows:

samples = geneCountTable(:,3:end).Properties.VariableNames;
untreated = strncmp(samples,'untreated',length('untreated'))
treated = strncmp(samples,'treated',length('treated'))

%% Plotting the Feature Assignments 
% The included file also contains a table |geneSummaryTable| with the
% summary of assigned and unassigned SAM entries. You can plot the basic
% distribution of the counting results by considering the number of reads
% that are assigned to the given genomic features (exons or genes for this
% example), as well as the number of reads that are unassigned (i.e. not
% overlapping any feature) or ambiguous (i.e. overlapping multiple
% features).

st = geneSummaryTable({'Assigned','Unassigned_ambiguous','Unassigned_noFeature'},:)
bar(table2array(st)','stacked');
legend(st.Properties.RowNames','Interpreter','none','Location','southeast');
xlabel('Sample')
ylabel('Number of reads')

%%
% Note that a small fraction of the alignment records in the SAM files is
% not reported in the summary table. You can notice this in the difference
% between the total number of records in a SAM file and the total number of
% records processed during the counting procedure for that same SAM file.
% These unreported records correspond to the records mapped to reference
% sequences that are not annotated in the GTF file and therefore are not
% processed in the counting procedure. If the gene models account for all
% the reference sequences used during the read mapping step, then all
% records are reported in one of the categories of the summary table.

geneSummaryTable{'TotalEntries',:} - sum(geneSummaryTable{2:end,:})


%% Plotting Read Coverage Across a Given Chromosome
% When read counting is performed without summarization using the function
% |featurecount|, the default IDs are composed by the attribute or metafeature
% (by default, gene_id) followed by the start and the stop positions of the
% feature (by default, exon). You can use the exon start positions to plot
% the read coverage across any chromosome in consideration, for example
% chromosome arm 2L.

% consider chromosome arm 2L
chr2L = strcmp(exonCountTable.Reference,'2L');
exonCount = exonCountTable{:,3:end};

% retrieve exon start positions
exonStart = regexp(exonCountTable{chr2L,1},'_(\d+)_','tokens');
exonStart = [exonStart{:}];
exonStart = cellfun(@str2num, [exonStart{:}]'); 

% sort exon by start positions
[~,idx] = sort(exonStart);

% plot read coverage along the genomic coordinates
figure;
plot(exonStart(idx),exonCount(idx,treated),'.-r',...
exonStart(idx),exonCount(idx,untreated),'.-b');
xlabel('Genomic position');
ylabel('Read count (exon level)');
title('Read count on Chromosome arm 2L');

% plot read coverage for each group separately 
figure;
subplot(2,1,1); 
plot(exonStart(idx),exonCount(idx,untreated),'.-r');
ylabel('Read count (exon level)');
title('Chromosome arm 2L (treated samples)');
subplot(2,1,2); 
plot(exonStart(idx),exonCount(idx,treated),'.-b');
ylabel('Read count (exon level)');
xlabel('Genomic position');
title('Chromosome arm 2L (untreated samples)');

%%
% Alternatively, you can plot the read coverage considering the starting
% position of each gene in a given chromosome. The file
% |pasilla_geneLength.mat| contains a table with the start and stop
% position of each gene in the corresponding gene annotation file.

% load gene start and stop position information
load pasilla_geneLength
geneLength(1:10,:)

%%

% consider chromosome 3 ('Reference' is a categorical variable)
chr3 = (geneLength.Reference == '3L') | (geneLength.Reference == '3R');
sum(chr3)

% consider the counts for genes in chromosome 3 
counts = geneCountTable{:,3:end};
[~,j,k] = intersect(geneCountTable{:,'ID'},geneLength{chr3,'ID'});
gstart = geneLength{k,'Start'};
gcounts = counts(j,:);

% sort according to ascending start position
[~,idx] = sort(gstart);

% plot read coverage by genomic position
figure;
plot(gstart(idx), gcounts(idx,treated),'.-r',...
	gstart(idx), gcounts(idx,untreated),'.-b');
xlabel('Genomic position')
ylabel('Read count');
title('Read count on Chromosome 3');

%% Normalizing Read Counts
% The read count in RNA-Seq data has been found to be linearly related to
% the abundance of transcripts [2]. However, the read count for a given
% gene depends not only on the expression level of the gene, but also on
% the total number of reads sequenced and the length of the gene
% transcript. Therefore, in order to infer the expression level of a gene
% from the read count, we need to account for the sequencing depth and the
% gene transcript length. One common technique to normalize the read count
% is to use the RPKM (Read Per Kilobase Mapped) values, where the read
% count is normalized by the total number of reads yielded (in millions)
% and the length of each transcript (in kilobases). This normalization
% technique, however, is not always effective since few, very highly
% expressed genes can dominate the total lane count and skew the expression
% analysis.

%%
% A better normalization technique consists of computing the effective
% library size by considering a size factor for each sample. By dividing
% each sample's counts by the corresponding size factors, we bring all the
% count values to a common scale, making them comparable. Intuitively, if
% sample A is sequenced N times deeper than sample B, the read counts of
% non-differentially expressed genes are expected to be on average N times
% higher in sample A than in sample B, even if there is no difference in
% expression.
% 
% To estimate the size factors, take the median of the ratios of observed
% counts to those of a pseudo-reference sample, whose counts can be
% obtained by considering the geometric mean of each gene across all
% samples [3]. Then, to transform the observed counts to a common scale,
% divide the observed counts in each sample by the corresponding size
% factor.

% estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(counts,2); 
nz = pseudoRefSample > 0; 
ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
sizeFactors = median(ratios,1)

%%

% transform to common scale
normCounts = bsxfun(@rdivide,counts,sizeFactors);
normCounts(1:10,:)

%%
% You can appreciate the effect of this normalization by using the function
% |boxplot| to represent statistical measures such as median, quartiles,
% minimum and maximum.

figure;

subplot(2,1,1)
maboxplot(log2(counts),'title','Raw Read Count','orientation','horizontal')
ylabel('sample')
xlabel('log2(counts)')

subplot(2,1,2)
maboxplot(log2(normCounts),'title','Normalized Read Count','orientation','horizontal')
ylabel('sample')
xlabel('log2(counts)')


%% Computing Mean, Dispersion and Fold Change 
% In order to better characterize the data, we consider the mean and the
% dispersion of the normalized counts. The variance of read counts is given
% by the sum of two terms: the variation across samples (raw variance) and
% the uncertainty of measuring the expression by counting reads (shot noise
% or Poisson). The raw variance term dominates for highly expressed genes,
% whereas the shot noise dominates for lowly expressed genes. You can plot
% the empirical dispersion values against the mean of the normalized counts
% in a log scale as shown below.

% consider the mean
meanTreated = mean(normCounts(:,treated),2); 
meanUntreated = mean(normCounts(:,untreated),2);

% consider the dispersion
dispTreated = std(normCounts(:,treated),0,2) ./ meanTreated; 
dispUntreated = std(normCounts(:,untreated),0,2) ./ meanUntreated;

% plot on a log-log scale
figure;
loglog(meanTreated,dispTreated,'r.');
hold on;
loglog(meanUntreated,dispUntreated,'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest');

%%
% Given the small number of replicates, it is not surprising to expect that
% the dispersion values scatter with some variance around the true value.
% Some of this variance reflects sampling variance and some reflects the
% true variability among the gene expressions of the samples. 

%% 
% You can look at the difference of the gene expression among
% two conditions, by calculating the fold change (FC) for each gene, i.e. the
% ratio between the counts in the treated group over the counts in the
% untreated group. Generally these ratios are considered in the log2 scale,
% so that any change is symmetric with respect to zero (e.g. a ratio of 1/2
% or 2/1 corresponds to -1 or +1 in the log scale).

% compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated) / 2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% plot mean vs. fold change (MA plot)
mairplot(meanTreated, meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')

%% 
% It is possible to annotate the values in the plot with the corresponding
% gene names, interactively select genes, and export gene lists to the
% workspace by calling the |mairplot| function as illustrated below:

mairplot(meanTreated,meanUntreated,'Labels',geneCountTable.ID,'Type','MA');

%%
% It is convenient to store the information about the mean value and fold
% change for each gene in a table. You can then access information about a
% given gene or a group of genes satisfying specific criteria by indexing
% the table by gene names.

% create table with statistics about each gene
geneTable = table(meanBase,meanTreated,meanUntreated,foldChange,log2FC);
geneTable.Properties.RowNames = geneCountTable.ID;

%%

% summary 
summary(geneTable)

%%

% preview 
geneTable(1:10,:)

%%

% access information about a specific gene
myGene = 'FBgn0261570';
geneTable(myGene,:)
geneTable(myGene,{'meanBase','log2FC'})

% access information about a given gene list
myGeneSet = {'FBgn0261570','FBgn0261573','FBgn0261575','FBgn0261560'};
geneTable(myGeneSet,:)

%% Inferring Differential Expression with a Negative Binomial Model
% Determining whether the gene expressions in two conditions are
% statistically different consists of rejecting the null hypothesis that
% the two data samples come from distributions with equal means. This
% analysis assumes the read counts are modeled according to a negative
% binomial distribution (as proposed in [3]). The function |nbintest|
% performs this type of hypothesis testing with three possible options to
% specify the type of linkage between the variance and the mean.

%% 
% By specifying the link between variance and mean as an identity, we
% assume the variance is equal to the mean, and the counts are modeled by
% the Poisson distribution [4].

tIdentity = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','Identity');
h = plotVarianceLink(tIdentity);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';
 
%% 
% Alternatively, by specifying the variance as the sum of the shot noise
% term (i.e. mean) and a constant multiplied by the squared mean, the
% counts are modeled according to a distribution described in [5]. The
% constant term is estimated using all the rows in the data.

tConstant = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','Constant');
h = plotVarianceLink(tConstant);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

%% 
% Finally, by considering the variance as the sum of the shot noise term
% (i.e. mean) and a locally regressed non-parametric smooth function of the
% mean, the counts are modeled according to the distribution proposed in
% [3]. 

tLocal = nbintest(counts(:,treated),counts(:,untreated),'VarianceLink','LocalRegression');
h = plotVarianceLink(tLocal);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

%%
% In order to evaluate which fit is the best for the data in consideration,
% you can compare the fitting curves in a single plot, as shown below.

h = plotVarianceLink(tLocal,'compare',true);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';
  
%%
% The output of |nbintest| includes a vector of P-values. A P-value
% indicates the probability that a change in expression as strong as the
% one observed (or even stronger) would occur under the null hypothesis,
% i.e. the conditions have no effect on gene expression. In the histogram
% of the P-values we observe an enrichment of low values (due to
% differentially expressed genes), whereas other values are uniformly
% spread (due to non-differentially expressed genes). The enrichment of
% values equal to 1 are due to genes with very low counts.

figure;
histogram(tLocal.pValue,100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment')

%%
% Filter out those genes with relatively low count to observe a more
% uniform spread of non-significant P-values across the range (0,1]. Note
% that this does not affect the distribution of significant P-values.

lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
histogram(tLocal.pValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')

%% Multiple Testing and Adjusted P-values
% Thresholding P-values to determine what fold changes are more significant
% than others is not appropriate for this type of data analysis, due to the
% multiple testing problem. While performing a large number of
% simultaneous tests, the probability of getting a significant result
% simply due to chance increases with the number of tests. In order to
% account for multiple testing, perform a correction (or adjustment) of
% the P-values so that the probability of observing at least one
% significant result due to chance remains below the desired significance
% level. 
%
% The Benjamini-Hochberg (BH) adjustment [6] is a statistical method that
% provides an adjusted  P-value answering the following question: what
% would be the fraction of false positives if all the genes with adjusted
% P-values below a given threshold were considered significant? Set a
% threshold of 0.1 for the adjusted P-values, equivalent to consider a 10%
% false positives as acceptable, and identify the genes that are
% significantly expressed by considering all the genes with adjusted
% P-values below this threshold.

% compute the adjusted P-values (BH correction)
padj = mafdr(tLocal.pValue,'BHFDR',true);

% add to the existing table
geneTable.pvalue = tLocal.pValue;
geneTable.padj = padj;

% create a table with significant genes
sig = geneTable.padj < 0.1;
geneTableSig = geneTable(sig,:);
geneTableSig = sortrows(geneTableSig,'padj');
numberSigGenes = size(geneTableSig,1)

%% Identifying the Most Up-regulated and Down-regulated Genes
% You can now identify the most up-regulated or down-regulated genes by
% considering an absolute fold change above a chosen cutoff. For example, a
% cutoff of 1 in log2 scale yields the list of genes that are up-regulated
% with a 2 fold change.

% find up-regulated genes
up = geneTableSig.log2FC > 1;
upGenes = sortrows(geneTableSig(up,:),'log2FC','descend');
numberSigGenesUp = sum(up)

% display the top 10 up-regulated genes
top10GenesUp = upGenes(1:10,:)

% find down-regulated genes
down = geneTableSig.log2FC < -1;
downGenes = sortrows(geneTableSig(down,:),'log2FC','ascend');
numberSigGenesDown = sum(down)

% find top 10 down-regulated genes
top10GenesDown = downGenes(1:10,:)

%%
% A good visualization of the gene expressions and their significance is
% given by plotting the fold change versus the mean in log scale and
% coloring the data points according to the adjusted P-values.

figure
scatter(log2(geneTable.meanBase),geneTable.log2FC,3,geneTable.padj,'o')
colormap(flipud(cool(256)))
colorbar;
ylabel('log2(Fold Change)')
xlabel('log2(Mean of normalized counts)')
title('Fold change by FDR')

%%
% You can see here that for weakly expressed genes (i.e. those with low
% means), the FDR is generally high because low read counts are dominated
% by Poisson noise and consequently any biological variability is drowned
% in the uncertainties from the read counting.

%% References
% [1] Brooks et al. Conservation of an RNA regulatory map between
%     Drosophila and mammals. Genome Research 2011. 21:193-202.
%
% [2] Mortazavi et al. Mapping and quantifying mammalian transcriptomes by
%     RNA-Seq. Nature Methods 2008. 5:621-628.
%
% [3] Anders et al. Differential expression analysis for sequence count
%     data. Genome Biology 2010. 11:R106.
%
% [4] Marioni et al. RNA-Seq: An assessment of technical reproducibility
%     and comparison with gene expression arrays. Genome Research 2008.
%     18:1509-1517.
%
% [5] Robinson et al. Moderated statistical test for assessing differences
%     in tag abundance. Bioinformatics 2007. 23(21):2881-2887.
%
% [6] Benjamini et al. Controlling the false discovery rate: a practical
%     and powerful approach to multiple testing. 1995. Journal of the Royal
%     Statistical Society, Series B57 (1):289-300.
