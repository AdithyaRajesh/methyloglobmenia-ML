methemoglobmenia_gbk = getgenbank('NM_001171660');
load mito
mito = methemoglobmenia_gbk.Sequence;
mito_length = length(mito)
first_300_bases = seqdisp(mito(1:300))
figure
ntdensity(mito)
bases = basecount(mito)
figure
basecount(mito,'chart','pie');
title('Distribution of Nucleotide Bases for Methemoglobmenia');
figure
dimers = dimercount(mito,'chart','bar')
title('Methemoglobmenia Genome Dimer Histogram');
[h,l] = featureview(methemoglobmenia_gbk,{'CDS','tRNA','rRNA','D_loop'},...
                                      [2 1 2 2 2],'Fontsize',16);
legend(h,l,'interpreter','none');
title('Homo sapiens methemoglobmenia, complete genome')
features = featureparse(methemoglobmenia_gbk,'Sequence',true)
coding_sequences = features.CDS;
coding_sequences_id = sprintf('%s ',coding_sequences.gene)
numCDS = numel(coding_sequences);
CDS = cell(numCDS,1);
for i = 1:numCDS
     seq = coding_sequences(i).Sequence;
     CDS{i} = seq(1:3*floor(length(seq)/3));
end
allCDS = [CDS{:}];
codoncount(allCDS)
figure
count = codoncount(allCDS,'figure',true,'geneticcode','Methemoglobmenia');
title('Human methemoglobmenia Genome Codon Frequency')
