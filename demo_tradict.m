% For this demo, we will use the first 90% of submissions to the SRA for A.
% thaliana to train a Tradict model. We will then use this trained model to
% predict the expression of the samples in the remaining 10% of submissions
% using only expression values of the selected markers obtained during 
% training.

%%% 1. SET UP THE ENVIRONMENT
% Load the quality filtered training transcriptome collection for A. thaliana.
% This will load three variables:
% - Y is a 21277 genes x 2597 samples matrix of transcripts per million (TPM)
%   expression values.
% - qt is a dataset that is an NCBI query table.
% - tids is a 21277 genes cell array of transcript ids.
% - sids is a 2597 samples cell array of sample ids.
%
% Change the below file path as needed.
% The file can be obtained from: 
% https://drive.google.com/file/d/0B252lj6tAx8XV0VCVFNfZXI5WVk/view?usp=sharing
load /Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered_public.mat

% Load the gene sets defining the transcriptional programs.
% This will load two variables of interest:
% - sets is a 150 transcriptional programs x 1 cell array where sets{i} is
%   another cell array of all the transcript IDs that belong to
%   transcriptional program i.
% - setnames is a 150 transcriptional programs x 2 cell array of names for
%   the transcriptional programs. setnames{i,1} and setnames{i,2} is the
%   short and long form names for transcriptional program i.
%
% Change the below file path as needed.
% The file can be obtained from: 
% https://drive.google.com/file/d/0B252lj6tAx8XV0VCVFNfZXI5WVk/view?usp=sharing
load ~/GitHub/transcriptome_compression/analysis/gene_ontology/Athaliana_representative_gene_set_02-Apr-2016.mat

% Add latent_log and this repository to your path. 
% Change the below file paths as needed.
% Code for the latent log can be cloned from GitHub
% e.g. git clone https://github.com/surgebiswas/latent_log.git
path(genpath('~/GitHub/latent_log'), path)
path(genpath('~/GitHub/tradict'), path)


% Partition the data into training and test sets.
% The below function will split the full collection of training
% transcriptomes by submission into a training set which contains the first
% (historically speaking) ~90% of training trainscriptomes and a test set
% which contains the latest (historically speaking) ~10% of training
% transcriptomes.
[ytrain, ytest, ktrain] = partition_data(Y', qt, 0.1);

%%% 2. TRAINING
% Set up the data for training. 
% ytrain is currently in units of TPMs; however, given samples were
% sequenced to different depths, the confidence in each TPM measurement
% varies from sample to sample. To factor in different sequencing depths,
% we first multiply each sample (row) of our training data by the
% sequencing depth (in millions of reads). This recasts the units of our
% training data into total measured transcripts. 
o = qt.spots(ktrain)/1000000; % sequencing depth of training samples (in millions of reads)
T = ytrain.*repmat(o, 1, size(ytrain,2) ); % multiply each row by depth

% Now train tradict. Type 'help tradict_train.m' for usage information.
model = tradict_train(T, o, tids, sets);

%%% 3. PREDICTION
% Make test set predictions.
% To do this we first multiply the TPM measurements for the marker genes by the
% sequencing depth to obtain total measured transcripts. 
% 
% model.S contains the indices of the selected markers.
o_test = qt.spots(~ktrain)/1000000;
t_test = ytest.*repmat(o_test, 1, size(ytest,2) );

% We then input only the expression values (in total measured transcripts)
% of the marker genes to construct the prediction.
% 
% - s_hat is a test_samples x transcriptional programs matrix of predicted 
%   expression values for the transcriptional programs.
% - z_hat is a test_samples x genes matrix of predicted log-TPM expression 
%   values for all genes.
%
% s_hat and z_hat can now be used for downstream analyses. 
[ s_hat, ~, z_hat ] = tradict_predict( t_test(:,model.S), o_test, model );


%%% 4. EVALUATE PERFORMANCE
% To evaluate the accuracy of these predictions, we must first calculate
% the actual expression values. 
%
% To compute log TPMs of all genes, we could naively take the usual 
% logarithm. However, due to zeros in the dataset this requires us to 
% choose a psuedocount, which we don't have a good idea of how to do. For
% this reason we compute the latent logarithm. Note that this behaves very
% similarly to log, however for low measured TPMs especially at low
% sequencing depths, latent log (lag) considers both the gene's a priori 
% abundance and the sequencing depth in order to better estimate its true 
% log TPM. 
% 
% z = log(t_test + 0.1) can be used as an alternative for the below in
% order to compare lag vs log. 
z = lag_dataset(t_test, o_test, 'priors', model.lag_priors);

% Given log TPMs we get the actual expression values of the transcriptional
% programs as follows:
zs = standardize(z, 'mu', model.train_mu, 'std', model.train_sig);
s = zs*model.geneset.coef;

% The test set contains mulitple submissions, which in aggregate, contain
% more biological signal than we would likely encounter in practice.
% Therefore we compute intra-submission performance by first performing a
% submission adjustment, which subtracts off the submission specific means
% from the expression values of all genes and transcriptional programs. We
% then compare z-score transformed expression so that different programs
% and genes are comparable. 

% Perform intra-submission adjustment 
% 
% tta and pta are test-set true and predicted, respectively, intra-submission 
% adjusted expression values for the transcriptional programs. These are
% test_samples x transcriptional programs matrices.
tsa = standardize(subadjust(s, qt.Submission(~ktrain)));
psa = standardize(subadjust(s_hat, qt.Submission(~ktrain)));

% tza and pza are test-set true and predicted, respectively, intra-submission 
% adjusted log TPM expression values for the genes. These are test_samples x genes 
% matrices.
tza = standardize(subadjust(z, qt.Submission(~ktrain)));
pza = standardize(subadjust(z_hat, qt.Submission(~ktrain)));

% Avg. PCC for transcriptional programs
disp('Avg. PCC for transcriptional programs:');
disp(corr(tsa(:), psa(:)));

% Avg. PCC for genes
disp('Avg. PCC for genes:');
disp(corr(tza(:), pza(:)));






