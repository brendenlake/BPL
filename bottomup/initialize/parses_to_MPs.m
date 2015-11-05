%
% Convert a set of coarse parses to MotorPrograms
%
% Input
%  I: image
%  S_walks: [nwalk x 1 cell] random walks on a skeleton graph
%  nint: number of parses to choose
%  lib: library
%  verbose: display details?
%
% Output
%  bestM_pass2: [ninit x 1 cell] list of best candidate motor programs
%  score_sorted: [n x 1] score of all parses considered in sorted order
%
function [bestM_pass2,score_sorted] = parses_to_MPs(I,PP,ninit,lib,verbose)
    
    if ~exist('verbose','var')
       verbose = false; 
    end
    ps = defaultps_bottomup;

    % Process the random walk
    PP.freeze;
    S_walks = PP.get_S; % return new random walks
    nwalk = numel(S_walks);

    % Transform sub-parsed random walks into motor programs
    fclassify = @(traj,scale) classify_traj_as_subid(traj,scale,lib);
    MP = cell(nwalk,1);
    if verbose, fprintf(1,'\nconverting to MotorPrograms...'); end
    parfor i=1:nwalk    
       MP{i} = sequence_to_MP(S_walks{i},I,fclassify,lib);       
    end
    if verbose, fprintf(1,'done.\n'); end
    
    % Remove all but the most promising parses
    K_order = ninit*ps.nwalk_order_mult;
    [bestM_pass1,score_sorted,score_indx] = pick_K_best(MP,lib,K_order);
%     if verbose
%         viz_parses(MP(score_indx),score_sorted);
%     end
    
    % Optimize the order of the strokes for each parse
    Norder = numel(bestM_pass1);
    if verbose, fprintf(1,'\noptimizing stroke order (for %d parses)\n',Norder); end  
    for i=1:Norder
        if verbose
            fprintf(1,'%d,',i);
            if mod(i,5)==0
                fprintf(1,'\n');
            end
        end
        optimize_order_MP(bestM_pass1{i},lib);
    end
    if verbose, fprintf(1,'done.\n'); end
    
    % Pick the best "ninit" samples
    [bestM_pass2,score_sorted,score_indx] = pick_K_best(bestM_pass1,lib,ninit);
    
end

 % Visualize the bag of parses
function viz_parses(list_M,scores)
   n = numel(list_M);   
   nrow = ceil(sqrt(n));
   figure(199);
   clf
   for i=1:n
      subplot(nrow,nrow,i);
      vizMP(list_M{i},'motor','start_size',8);
      title(makestr(scores(i)));
   end           
end

%
% Pick the K best hypotheses in a set
%
% Input
%  list_M: [n x 1 cell] of hypotheses
%  lib
%  K: how many to choose?
%
% Output
%  bestM: [k x 1 cell] of best hypotheses, from best to worst
%  score_sorted: [n x 1 scalar] score of all hypotheses in order
function [bestM,score_sorted,score_indx] = pick_K_best(list_M,lib,K)
    n = numel(list_M);
    scores = score_set(list_M,lib);
    npick = min(n,K);
    [score_sorted,score_indx] = randsort(scores,'descend');
    pick = score_indx(1:npick);
    bestM = list_M(pick);
end

%
% Score the set of hypotheses in list_M
%
function scores = score_set(list_M,lib)
    n = numel(list_M);
    scores = -inf(n,1);
    for i=1:n
        scores(i) = scoreMP_NoRel(list_M{i},lib,...
                    'type',true,'token',true,'image',false);
    end
end