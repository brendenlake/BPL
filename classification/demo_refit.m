% Demo for refitting a set of motor programs, learned from one image, to a
% new image. Corresponds to Figure 4 in the paper.

% Parameters
use_precomputed = true; % (yes/no) use pre-computed results?
fast_mode = true; % (yes/no) only fit affine parameters not strokes?
irun = 9; % run index (1...20)
itrain = 9; % training item index (1...20)
itest = 9; % test item index (1...20)

if fast_mode && ~use_precomputed
    fprintf(1,'Fast mode refits a model using only an affine transformation, rather than fitting the motor proram.\n');
    warning_mode('Fast mode is for demo purposes only and was not used in paper results.')
end

load('items_classification','cell_train','cell_test');
lib = loadlib;

img_train = cell_train{irun}{itrain}.img;
img_test = cell_test{irun}{itest}.img;

% load previous model fit
fn_fit_train = fullfile('model_fits',makestr('run',irun,'_train',itrain,'_G'));
load(fn_fit_train,'G');
K = length(G.models); % add image to model class (removed to save memory)
for i=1:K
   G.models{i}.I = G.img; 
end

if use_precomputed
    clear pair
    fn_refit_test = fullfile('model_refits',makestr('run',irun,'_fit_train' ,itrain,'_to_image_test',itest));
    if ~exist([fn_refit_test,'.mat'],'file')
       fprintf(1,'Please download pre-computed model results to use this feature. Program quiting...\n');
       return
    end
    load(fn_refit_test,'pair');
    K = length(pair.Mbest);
    for i=1:K
       pair.Mbest{i}.I = img_test; 
    end
else    
    % refit models to new image
    Mbest = cell(K,1);
    fit_score = nan(K,1);
    prior_score = G.scores;
    for i=1:K
        fprintf(1,'re-fitting parse %d of %d\n',i,K);
        [Mbest{i},fscore] = FitNewExemplar(img_test,G.samples_type{i},lib,true,fast_mode);
        fit_score(i) = fscore;
    end
    
    % save output structure
    pair = struct;
    pair.Mbest = Mbest;
    pair.fit_score = fit_score;
    pair.prior_score = prior_score;
end

% visualize
figure;
sz = [307 1016];
pos = get(gcf,'Position');
pos(3:4) = sz;
set(gcf,'Position',pos);
subplot(6,2,1);
plot_image_only(img_train);
title('training image');
subplot(6,2,2);
plot_image_only(img_test);
title('test image');
count = 3;
for i=1:K
    subplot(6,2,count);
    count = count + 1;
    vizMP(G.models{i},'motor');
    title(makestr('parse ',i,'; fit score: ',pair.prior_score(i)));
    subplot(6,2,count);
    count = count + 1;
    vizMP(pair.Mbest{i},'motor');
    title(makestr('refit score: ',pair.fit_score(i)));
end