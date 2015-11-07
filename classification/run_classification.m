% Calculate one-shot classification error rate
function run_classification

    classdir = 'model_refits';
    load('items_classification','nrun','ntrain','ntest','cell_Y');

    trainset = num2cell(vec(1:ntrain));
    testset = num2cell(vec(1:ntest));

    fprintf(1,'one-shot classification results\n');
    perror = zeros(nrun,1);
    for r=1:nrun
        Y = cell_Y{r};
        fscore = @(itrain,itest) fclassify(itrain,itest,r,classdir);
        perror(r) = myclassify(trainset,testset,fscore,Y,'score');
        fprintf(1,' run %d (error %s%%)\n',r,num2str(perror(r),3));
    end
    fprintf(1,'average error: %s%%\n',num2str(mean(perror),3));

end

%
% Bayesian classification score
%
%  log P(I_T|I_C) + log(I_C|I_T) - log(I_C)
%   for training image I_C
%   and test image I_T
%
% Input
%  itrain: train index
%  itest: test index
%  irun: run index
%  classdir: file directory for "crossFit" files
%
% Output
%  log_score: [scalar] score
%
function log_score = fclassify(itrain,itest,irun,classdir)

    srun = num2str(irun);
    strain = num2str(itrain);
    stest = num2str(itest);

    fn_test_to_train = ['run',srun,'_fit_test' ,stest,'_to_image_train',strain];
    fn_train_to_test = ['run',srun,'_fit_train',strain,'_to_image_test',stest];

    % load file optimizing P(I_T|I_C)
    if ~exist(fullfile(classdir,[fn_train_to_test,'.mat']),'file')
       fprintf(1,'Please download pre-computed model results to use this feature. Program quiting...\n');
       assert false;
    end
    load(fullfile(classdir,fn_train_to_test),'fit_score','prior_score');
    pair.fit_score = fit_score;
    pair.prior_score = prior_score;
    pair_fit_train_to_test = pair;
    clear pair
    
    % load file optimizing P(I_C|I_T)
    load(fullfile(classdir,fn_test_to_train),'fit_score','prior_score');
    pair.fit_score = fit_score;
    pair.prior_score = prior_score;
    pair_fit_test_to_train = pair;
    clear pair
   
    % compute score
    [log_P_IT_given_IC,prior_scores] = log_post_pred(pair_fit_train_to_test);
    log_P_IC_given_IT = log_post_pred(pair_fit_test_to_train);    
    log_P_IC = logsumexp(prior_scores(:)); 
    log_score = log_P_IC_given_IT + log_P_IT_given_IC - log_P_IC;
    
end

%
% Log-posterior predictive score, log P(I_2 | I_1), using discrete approximation
%
% Input
%  pair: structure
%  
% Output
%  logscore: log conditional probability
%  prior_scores:  [k x 1] log joint probability (unnormalized) weights for
%       each parse of I_1
%
function [logscore,prior_scores] = log_post_pred(pair)

    % extra info
    logfit = pair.fit_score;
    prior_scores = pair.prior_score;
    
    % normalize weights
    logwt = prior_scores - logsumexp(prior_scores(:));
    
    % combine weights with fit term
    logv = logfit(:) + logwt;
    logscore = logsumexp(logv);
    
end