%% rgp custom workflow
% Use mLASSO-StARS to build a TRN from gene expression and prior at
% dynamic.
%% Modifications:
%1. removed all the zscore inside functions
%2. add options for using dYdt or PtXt
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
%% Guohua Yuan 190730

clear all
close all
restoredefaultpath

matlabDir = '../..';

addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% 1. Import gene expression data, list of regulators, list of target genes
% into a Matlab .mat object.
outdir = './outputs/test1/rgp_meta_dYdtYlogtpm1zscorelogistic_acctfazscorelogistic_exp20cel/peakinter';
geneExprTFAdir = fullfile(outdir, 'processedGeneExpTFA');
mkdir(geneExprTFAdir)
geneExprMat = fullfile(geneExprTFAdir,['coaccpeakgenergpfilhic_hocop1e_5_b_' 'geneExprGeneLists.mat']);
normGeneExprFile = './inputs/geneExpression/acctfazscorelogistic_p01rna_peakinter_coaccpeakgenergpfilhic_hocop1e_5_b.txt';
%normGeneExprFile = './inputs/geneExpression/atac_meta_chromvar_order_zscore.txt';
%normGeneExprFile = './inputs/geneExpression/atac_chromvarzscore_add_rna_logtpm1zscore_meta.txt';

dYdtYFile = './inputs/geneExpression/dYdtY_meta_logtpm1zscorelogistic_order.txt';
%dYdtYFile = './inputs/geneExpression/rna_meta_logtpm1zscorelogistic_order.txt';
targGeneFile = './inputs/targRegLists/targene_exp20c.txt';
potRegFile = './inputs/targRegLists/regulator_exp20c_hocop1e_5.txt';
tfaGeneFile = './inputs/targRegLists/targene_exp20c.txt';
%PtXtFile = './inputs/geneExpression/PtXt_pb_zscore.190911.txt';
PtXtopt = '';% options are 'PtXt' or ''

logistic_opt =''; %'' or 'logistic'
Yopt = '_dYdtY';% options are '' or '_dYdtY'; '' for use Y, '_dYdtY' use V

%V_grislimat = './inputs/geneExpression/V_grisli.mat';

try
    ls(geneExprMat)
catch
disp('1. importGeneExpGeneLists_custom.m')
% importGeneExpGeneLists_custom(normGeneExprFile,dYdtYFile,V_grislimat,Yopt,targGeneFile,potRegFile,...
%         tfaGeneFile,geneExprMat)
importGeneExpGeneLists_custom(normGeneExprFile,dYdtYFile,Yopt,targGeneFile,potRegFile,...
        tfaGeneFile,geneExprMat)
end


if PtXtopt
% input PtXt file
load(geneExprMat);
PtXt = importdata(PtXtFile);
PtXt_matrix = reshape(PtXt(:,4),[tottg,tottf,totSamps]);
clear PtXt
end

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
priorNamelist = {'prior_coaccpeakgenergpfilhic_hocop1e_5_b_1',
             'prior_coaccpeakgenergpfilhic_hocop1e_5_b_Trans',
             'prior_coaccpeakgenergpfilhic_hocop1e_5_b_TransVar'
             };
         
for pind = 1:1%length(priorNamelist)          
priorName = priorNamelist{pind};
priorFile = ['./inputs/priors/' priorName '.txt']; % ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

close all

try 
    ls(tfaMat)
catch

    disp('2. integratePrior_estTFA.m')
    integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
         minTargets, tfaMat)
end


%% Evaluate models for a variety of leave-out gene expression datasets, 
%% lambda bias values, and both TFA methods

%loInfo = {'','';'./inputs/leaveOutLists/lotest1.txt', '_lo1';};

loInfo = {'','';};

lambdaBiases = [1 .5 .25 .1]; % correspond to "no,moderate, and strong prior reinforcement
tfaOpts = {'_TFmRNA'}; % the two TFA options '','_TFmRNA'

totLos = size(loInfo,1);
totTfas = length(tfaOpts);

%% parameters for Step 3 (calculating network instabilities w/ bStARS)
totSS = 50;
targetInstability = .05;
lambdaMin = .01;
lambdaMax = 1;
extensionLimit = 1;
totLogLambdaSteps = 25; % will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5;
subsampleFrac = .63;
%% parameters for Step 4 (ranking TF-gene interactions)
meanEdgesPerGene = 15;
instabSource = 'Network';


loopt ='lo';% options are 'lo' or ''

for lind = 1:totLos  % you can use parfor loop here, if you have parallel computing toolbox  
    leaveOutSampleList = '';
    leaveOutInf = '';
    
    if loopt
    leaveOutSampleList = loInfo{lind,1};
    leaveOutInf = loInfo{lind,2};
    end

for tind = 1:totTfas
    tfaOpt = tfaOpts{tind};

    
for lambdaBias = lambdaBiases
    
%% 3. Calculate network instabilities using bStARS

instabilitiesDir = fullfile(outdir,strrep(['instabilities_targ' ...
    num2str(targetInstability) '_SS' num2str(totSS) leaveOutInf '_bS' num2str(bStarsTotSS)],'.','p'));
mkdir(instabilitiesDir)

netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
instabOutMat = fullfile(instabilitiesDir,netSummary);
if  PtXtopt
    disp('3. estimateInstabilitiesTRNbStARS_custom.m')
    estimateInstabilitiesTRNbStARS_custom(geneExprMat,tfaMat,PtXt_matrix,lambdaBias,tfaOpt,...
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,...
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
else
    disp('3. estimateInstabilitiesTRNbStARS_custom_noPtXt.m')
    estimateInstabilitiesTRNbStARS_custom_noPtXt(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,...
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
end

%% 4. For a given instability cutoff and model size, rank TF-gene
% interactions, calculate stabilities and network file for jp_gene_viz
% visualizations
priorMergedTfsFile = ['./inputs/priors/' priorName '_mergedTfs.txt'];
try % not all priors have merged TFs and merged TF files
    ls(priorMergedTfsFile) 
catch
    priorMergedTfsFile = '';
end
targetInstability = .05;
networkDir = strrep(instabilitiesDir,'instabilities','networks');
mkdir(networkDir);
networkSubDir = fullfile(networkDir,[instabSource ...
    strrep(num2str(targetInstability),'.','p') '_' ...
    num2str(meanEdgesPerGene) 'tfsPerGene']);
mkdir(networkSubDir)
trnOutMat = fullfile(networkSubDir,netSummary);
outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
networkHistDir = fullfile(networkSubDir,'Histograms');
mkdir(networkHistDir)

subsampHistPdf = fullfile(networkHistDir,[netSummary '_ssHist']);

if PtXtopt
disp('4. buildTRNs_mLassoStARS_custom.m')
buildTRNs_mLassoStARS_custom(instabOutMat,tfaMat,priorMergedTfsFile,...
    meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
    outNetFileSparse)
else
disp('4. buildTRNs_mLassoStARS_custom_noPtXt.m')
buildTRNs_mLassoStARS_custom_noPtXt(instabOutMat,tfaMat,priorMergedTfsFile,...
    meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
    outNetFileSparse)
end

%% 5. Calculate R2 for cross validation
if leaveOutInf
modSizes = 1:meanEdgesPerGene;
r2OutMat = [instabOutMat '_r2pred'];

%disp('5. calcR2predFromStabilities_custom_logistic')
%calcR2predFromStabilities_custom_logistic(instabOutMat,trnOutMat,r2OutMat,modSizes)
if PtXtopt
    disp('5. calcR2predFromStabilities_custom_PtXt')
calcR2predFromStabilities_custom_PtXt(instabOutMat,trnOutMat,r2OutMat,modSizes)
else
    if logistic_opt
disp('5. calcR2predFromStabilities_custom')
calcR2predFromStabilities_custom(instabOutMat,trnOutMat,r2OutMat,modSizes)
    else
disp('5. calcR2predFromStabilities_custom_logistic')
calcR2predFromStabilities_custom_logistic(instabOutMat,trnOutMat,r2OutMat,modSizes)
    end
end
end

end
end
end
end