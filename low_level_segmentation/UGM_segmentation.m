clearvars;

clc;

im_name='2_1_s.bmp';

% TODO: Update library path
% Add  library paths
basedir='UGM';
addpath(genpath(basedir));


%Set model parameters
%cluster color
K=4; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[10 1]; % Potts Model

%Load images
im = imread(im_name);


[rows, cols, channels] = size(im);
%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%
im = RGB2Lab(im);


%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term
% nodePot = P( color at pixel 'x' | Cluster color 'c' )  

%Preparing data for GMM fiting
im=double(im);
x=reshape(im,[rows*cols, channels]);
% gmm_color = fitgmdist(x,K);
gmm_color = gmdistribution.fit(x,K);
mu_color=gmm_color.mu;

% Estimate Unary potentials
nodePot = gmm_color.posterior(x);
data_term = nodePot;


%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
[edgePot,edgeStruct] = CreateGridUGMModel(rows, cols, K ,smooth_term);


if ~isempty(edgePot)

    % color clustering
    [~,c] = max(reshape(data_term,[rows*cols K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different UGM inference algorithms
    disp('Loopy Belief Propagation'); tic;
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot, edgePot, edgeStruct); toc;
    [~,c] = max(nodeBelLBP,[],2);
    im_lbp = reshape(mu_color(c,:),size(im));

    % Max-sum
    disp('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot, edgePot, edgeStruct);
    im_bp = reshape(mu_color(decodeLBP,:),size(im));
    toc;

    disp('ICMrestart'); tic;
    ICMrestartDecoding = UGM_Decode_ICMrestart(nodePot,edgePot,edgeStruct,50);
    im_icm = reshape(mu_color(ICMrestartDecoding,:),size(im));
    toc;

    disp('Max of Marginals'); tic;
    maxOfMarginalsMFdecode = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_MeanField);
    im_mM = reshape(mu_color(maxOfMarginalsMFdecode,:),size(im));
    toc;

    % Graph Cut
%     disp('Graph Cut'); tic;
%     decodeGC = UGM_Decode_GraphCut(nodePot, edgePot, edgeStruct);
% %     im_gc = reshape(mu_color(decodeGC,:),size(im));
%     im_gc = reshape(decodeGC, rows, cols);
%     toc;


    
    % TODO: apply other inference algorithms and compare their performance
    %
    % - Graph Cut
    % - Linear Programing Relaxation
    
    figure
    subplot(2,3,1),imshow(Lab2RGB(im));xlabel('Original');
    subplot(2,3,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
    subplot(2,3,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
    subplot(2,3,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');
    subplot(2,3,5),imshow(Lab2RGB(im_icm),[]);xlabel('ICMrestart');
    subplot(2,3,6),imshow(Lab2RGB(im_mM),[]);xlabel('Max of Marginals');
    
else
   
    error('You have to implement the CreateGridUGMModel.m function');

end