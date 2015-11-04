classdef MCMC < BetterHandle
    %%MCMC Simple MH and Gibbs moves for the MotorProgram Class
    %%
    %%   This class provides some simple moves for sampling from
    %%   some variables, with fixed dimension. It does not provide
    %%   a sampler for the entire joint space.
    %%
    properties
        
       debug_full_score = false;
        
       %% book-keeping to count the number of times
       %% a move is accepted
       
       %% parameters that don't call the image evaluation
       accept_shape_type = 0;
       n_shape_type = 0;
       accept_scale_type = 0;
       n_scale_type = 0;
       accept_gpos = 0;
       n_gpos = 0;
       accept_seval = 0;
       n_seval = 0;
       accept_seval_token = 0;
       n_seval_token = 0;
       accept_sid = 0;
       n_sid = 0;
       accept_r = 0;
       n_r = 0;
       
       %% parameters that call the image evaluation function
       accept_shape_token = 0;
       n_shape_token = 0;
       accept_scale_token = 0;
       n_scale_token = 0;
       accept_pos_token = 0;
       n_pos_token = 0;
       
    end
    
    methods   
        
        function this = MCMC(debug_full_score)
            if exist('debug_full_score','var')
                this.debug_full_score = debug_full_score;
            end
        end
        
        %% print the percent of times various moves are accepted
        function acceptance_ratios(this)
           fprintf(1,'Acceptance ratios\n');
           fprintf(1,' shape type: %s%% (of %d)\n',num2str(this.accept_shape_type/this.n_shape_type*100,3),this.n_shape_type);
           fprintf(1,' scale type: %s%% (of %d)\n',num2str(this.accept_scale_type/this.n_shape_type*100,3),this.n_scale_type);
           fprintf(1,' global pos: %s%% (of %d)\n',num2str(this.accept_gpos/this.n_gpos*100,3),this.n_gpos);
           if this.n_seval > 0 
               fprintf(1,' seval type: %s%% (of %d)\n',num2str(this.accept_seval/this.n_seval*100,3),this.n_seval);
               fprintf(1,' seval tokn: %s%% (of %d)\n',num2str(this.accept_seval_token/this.n_seval_token*100,3),this.n_seval_token);
           end
           fprintf(1,' sub-ids   : %s%% (of %d)\n',num2str(this.accept_sid/this.n_sid*100,3),this.n_sid);
           fprintf(1,' R         : %s%% (of %d)\n',num2str(this.accept_r/this.n_r*100,3),this.n_r);
           
           fprintf(1,' scale tokn: %s%% (of %d)\n',num2str(this.accept_scale_token/this.n_scale_token*100,3),this.n_scale_token);
           fprintf(1,' shape tokn: %s%% (of %d)\n',num2str(this.accept_shape_token/this.n_shape_token*100,3),this.n_shape_token);
           fprintf(1,' pos   tokn: %s%% (of %d)\n',num2str(this.accept_pos_token/this.n_pos_token*100,3),this.n_pos_token);
        end
        
        %% SHAPE TYPE MOVE
        function [val_new,accept] = mh_shape_type(this,sid,bid,M,lib) 
            %% sid: stroke index
            %% bid: sub-stroke index
            
            %% make proposal
            curr_shape_type = M.S{sid}.shapes_type(:,:,bid);
            prop_shape_type = propose_shape_type(curr_shape_type,M.parameters);
            
            %% score proposal
            curr_score = score_shape_type(M,sid,bid,lib,M.parameters);
            if this.debug_full_score
               curr_score2 = score_shape_type(M,sid,bid,lib,M.parameters,true);
            end
            M.S{sid}.shapes_type(:,:,bid) = prop_shape_type;
            prop_score = score_shape_type(M,sid,bid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_shape_type(M,sid,bid,lib,M.parameters,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH Shape type passed score test\n');
                end
            end
                
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.shapes_type(:,:,bid) = curr_shape_type;
            end
            val_new = M.S{sid}.shapes_type(:,:,bid);
            
            this.accept_shape_type = this.accept_shape_type + accept;
            this.n_shape_type = this.n_shape_type + 1;
            
            if isinf(scoreMP(M,lib))
               error('score should not be negative'); 
            end
        end
        
        %% SHAPE TOKEN MOVE
        function [val_new,accept] = mh_shape_token(this,sid,bid,M,lib) 
            %% sid: stroke index
            %% bid: sub-stroke index
            
            %% make proposal
            curr_shape_token = M.S{sid}.shapes_token(:,:,bid);
            prop_shape_token = propose_shape_token(curr_shape_token,M.parameters);
            
            %% score proposal
            curr_score = score_shape_token(M,sid,bid,lib,M.parameters);
            if this.debug_full_score
               curr_score2 = score_shape_token(M,sid,bid,lib,M.parameters,true);
            end
            M.S{sid}.shapes_token(:,:,bid) = prop_shape_token;
            prop_score = score_shape_token(M,sid,bid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_shape_token(M,sid,bid,lib,M.parameters,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH Shape token passed score test\n');
                end
            end
                
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.shapes_token(:,:,bid) = curr_shape_token;
            end
            val_new = M.S{sid}.shapes_token(:,:,bid);
            
            this.accept_shape_token = this.accept_shape_token + accept;
            this.n_shape_token = this.n_shape_token + 1;
            
            if isinf(scoreMP(M,lib))
               error('score should not be negative'); 
            end
        end
        
        
        %% SCALE TYPE MOVE
        function [val_new,accept] = mh_scale_type(this,sid,bid,M,lib)
            %% sid: stroke index
            %% bid: sub-stroke index
           
            %% make proposal
            curr_invscale_type = M.S{sid}.invscales_type(bid);
            prop_invscale_type = propose_scale_type(curr_invscale_type,M.parameters);
            
            %% score proposal
            curr_score = score_scale_type(M,sid,bid,lib,M.parameters);
            if this.debug_full_score
                curr_score2 = score_scale_type(M,sid,bid,lib,M.parameters,true);
            end
            M.S{sid}.invscales_type(bid) = prop_invscale_type;
            prop_score = score_scale_type(M,sid,bid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_scale_type(M,sid,bid,lib,M.parameters,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH Scale type passed score test\n');
                end
            end
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
               M.S{sid}.invscales_type(bid) = curr_invscale_type;
            end            
            val_new = M.S{sid}.invscales_type(bid);
            
            this.accept_scale_type = this.accept_scale_type + accept;
            this.n_scale_type = this.n_scale_type + 1;
            
        end
        
        %% SCALE TOKEN MOVE
        function [val_new,accept] = mh_scale_token(this,sid,bid,M,lib)
            %% sid: stroke index
            %% bid: sub-stroke index
           
            %% make proposal
            curr_invscale_token = M.S{sid}.invscales_token(bid);
            prop_invscale_token = propose_scale_token(curr_invscale_token,M.parameters);
            
            %% score proposal
            curr_score = score_scale_token(M,sid,bid,lib,M.parameters);
            if this.debug_full_score
                curr_score2 = score_scale_token(M,sid,bid,lib,M.parameters,true);
            end
            M.S{sid}.invscales_token(bid) = prop_invscale_token;
            prop_score = score_scale_token(M,sid,bid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_scale_token(M,sid,bid,lib,M.parameters,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH Scale token passed score test\n');
                end
            end
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
               M.S{sid}.invscales_token(bid) = curr_invscale_token;
            end            
            val_new = M.S{sid}.invscales_token(bid);
            
            this.accept_scale_token = this.accept_scale_token + accept;
            this.n_scale_token = this.n_scale_token + 1;
            
        end

        %% GLOBAL POSITION FOR UNIHIST RELATION
        function [val_new,accept] = mh_gobal_position(this,sid,M,lib)
            
            %% sid: stroke index           
            assert(strcmp(M.S{sid}.R.type,'unihist'));
            
            %% make proposal
            curr_gpos = M.S{sid}.R.gpos;
            prop_gpos = propose_global_position(curr_gpos,M.parameters);
            
            %% score proposal
            curr_score = score_global_position(M,sid,lib);
            if this.debug_full_score
               curr_score2 = score_global_position(M,sid,lib,true); 
            end
            M.S{sid}.R.gpos = prop_gpos;
            prop_score = score_global_position(M,sid,lib);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_global_position(M,sid,lib,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH position passed score test\n');
                end
            end
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.R.gpos = curr_gpos;
            end            
            val_new = M.S{sid}.R.gpos;
            
            this.accept_gpos = this.accept_gpos + accept;
            this.n_gpos = this.n_gpos + 1;
            
        end
        
        %% POSITION FOR A STROKE
        function [val_new,accept] = mh_token_position(this,sid,M,lib)
            %% sid: stroke index
                       
            %% make proposal
            curr_pos_token = M.S{sid}.pos_token;
            prop_pos_token = propose_local_position(curr_pos_token,M.parameters);
            
            %% score proposal
            curr_score = score_local_position(M,sid,lib);
            if this.debug_full_score
               curr_score2 = score_local_position(M,sid,lib,true); 
            end
            M.S{sid}.pos_token = prop_pos_token;
            prop_score = score_local_position(M,sid,lib);
            
            if this.debug_full_score %% debug to make sure full-score matches
                prop_score2 = score_local_position(M,sid,lib,true);
                if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH token position passed score test\n');
                end
            end
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.pos_token = curr_pos_token;
            end            
            val_new = M.S{sid}.pos_token;
            
            this.accept_pos_token = this.accept_pos_token + accept;
            this.n_pos_token = this.n_pos_token + 1;
        end
        
        %% EVAL-SPOT TYPE
        function [val_new,accept] = mh_eval_spot_type(this,sid,M,lib)
            %% sid: stroke index
            assert(strcmp(M.S{sid}.R.type,'mid'));
            
            %% make proposal
            curr_eval_spot_type = M.S{sid}.R.eval_spot_type;
            assert(~isempty(curr_eval_spot_type));
            prop_eval_spot_type = propose_eval_spot(curr_eval_spot_type,M.parameters);
            
            %% score proposal
            curr_score = score_eval_spot_type(M,sid,lib,M.parameters);
            if this.debug_full_score
                curr_score2 = score_eval_spot_type(M,sid,lib,M.parameters,true);
            end
            M.S{sid}.R.eval_spot_type = prop_eval_spot_type;
            prop_score = score_eval_spot_type(M,sid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
               prop_score2 = score_eval_spot_type(M,sid,lib,M.parameters,true);
               if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH eval spot type passed score test\n');
                end
            end            
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.R.eval_spot_type = curr_eval_spot_type;
            end            
            val_new = M.S{sid}.R.eval_spot_type;
            
            this.accept_seval = this.accept_seval + accept;
            this.n_seval = this.n_seval + 1;
            
        end
        
        function [val_new,accept] = mh_eval_spot_token(this,sid,M,lib)
            %% sid: stroke index
            assert(strcmp(M.S{sid}.R.type,'mid'));
            
            %% make proposal
            curr_eval_spot_token = M.S{sid}.R.eval_spot_token;
            prop_eval_spot_token = propose_eval_spot(curr_eval_spot_token,M.parameters);
            
            %% score proposal
            curr_score = score_eval_spot_token(M,sid,lib,M.parameters);
            if this.debug_full_score
                curr_score2 = score_eval_spot_token(M,sid,lib,M.parameters,true);
            end
            M.S{sid}.R.eval_spot_token = prop_eval_spot_token;
            prop_score = score_eval_spot_token(M,sid,lib,M.parameters);
            
            if this.debug_full_score %% debug to make sure full-score matches
               prop_score2 = score_eval_spot_token(M,sid,lib,M.parameters,true);
               if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH eval spot token passed score test\n');
                end
            end            
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score);
            if ~accept
                M.S{sid}.R.eval_spot_token = curr_eval_spot_token;
            end            
            val_new = M.S{sid}.R.eval_spot_token;
            
            this.accept_seval_token = this.accept_seval_token + accept;
            this.n_seval_token = this.n_seval_token + 1;
            
        end
        
        %% SUB-STROKE IDS
        function [val_new,accept] = gibbs_substroke_id(this,sid,bid,M,lib)
            % sid: stroke index
            % bid: sub-stroke index
            curr_id = M.S{sid}.ids(bid);
            nclust = lib.N;
            score = score_all_subid(M,sid,bid,lib);            
           
            if this.debug_full_score
                score2 = zeros(nclust,1);
                for c=1:nclust
                    M.S{sid}.ids(bid) = c;
                    score2(c) = scoreMP(M,lib);
                end
                v1 = score-logsumexp(score(:));
                v2 = score2-logsumexp(score2(:));
                v1(isinf(v1)) = [];
                v2(isinf(v2)) = [];
                assert( aeq(v1,v2) );
                fprintf(1,'DEBUG: Gibbs sub-stroke ids passed score test\n');
            end
            
            %% normalize probability distribution
            score = score - logsumexp(score);
            pvec = exp(score);
            assert(aeq(sum(pvec),1));
            pvec = pvec ./ sum(pvec);
            
            %% sample the new id
            vr = mnrnd(1,pvec);
            assert(~any(isnan(vr)));
            new_id = find(vr);            
            M.S{sid}.ids(bid) = new_id;
            accept = new_id ~= curr_id;
            
            val_new = M.S{sid}.ids(bid);
            
            this.accept_sid  = this.accept_sid + accept;
            this.n_sid = this.n_sid + 1;
            
        end
        
        %% RELATION (WITH TOKEN PARAMETER FOR ATTACHMENT)
        function [val_new,accept] = mh_relation(this,sid,M,lib)
            
            %% make proposal
            PM = defaultps;
            curr_R = M.S{sid}.R;
            prop_R = propose_relation(M.S{sid}.pos_token,M.S(1:sid-1),lib,PM);
            
            %% score proposal
            curr_score = score_relation(M,sid,lib);
            g_prop_to_curr = pp_relation(M,sid,lib,PM);
            if this.debug_full_score
                curr_score2 = score_relation(M,sid,lib,true);
            end
            M.S{sid}.R = prop_R;
            prop_score = score_relation(M,sid,lib);
            
            g_curr_to_prop = pp_relation(M,sid,lib,PM);
            if this.debug_full_score
               prop_score2 = score_relation(M,sid,lib,true);
               if ~isinf(prop_score)
                    assert(aeq(prop_score-curr_score,prop_score2-curr_score2));
                    fprintf(1,'DEBUG: MH relation passed score test\n');
                end
            end
            
            %% accept or reject
            accept = mh_accept(prop_score,curr_score,g_curr_to_prop,g_prop_to_curr);
            if ~accept
                M.S{sid}.R = curr_R;
            end            
            val_new = M.S{sid}.R;
            
            %% record keeping
            this.accept_r = this.accept_r + accept;
            this.n_r = this.n_r + 1;
                 
        end
        
        
    end
end

%% general MH acceptance function
%% 
%% Input
%%  prop_score: log-probability of model with proposal
%%  curr_score: log-probability of current model
%%  g_curr_to_prop: log-prob. of transition from current to proposal model
%%  g_prop_to_curr: reverse log-prob
function accept = mh_accept(prop_score,curr_score,g_curr_to_prop,g_prop_to_curr)
    if exist('g_curr_to_prop','var')
        lr = prop_score - curr_score + g_prop_to_curr - g_curr_to_prop;
    else
        lr = prop_score - curr_score; 
    end
    
    if isinf(prop_score)
       accept = false;
       return
    end
    
    r = min(1,exp(lr));
    assert(~isinf(curr_score));
    accept = rand < r;
end

function R = propose_relation(pos_token,prev_strokes,lib,PM)
    %% Propose new relation by sampling from prior, except the global
    %% position is sampled from near the token position
    R = CPD.sample_relation_type(lib,prev_strokes);
    if strcmp(R.type,'mid')
        R.eval_spot_token = CPD.sample_relation_token(lib,R.eval_spot_type);
    elseif strcmp(R.type,'unihist')
        R.gpos(1) = pos_token(1) + PM.mcmc_prop_relpos_mlty*lib.rel.sigma_x*randn;
        R.gpos(2) = pos_token(2) + PM.mcmc_prop_relpos_mlty*lib.rel.sigma_y*randn;
    end
end

function ll = pp_relation(M,sid,lib,PM)
    %% The probability of proposing relation M.S{sid}.R
    ll = CPD.score_relation_type(lib,M.S{sid}.R);
    if strcmp(M.S{sid}.R.type,'mid')
       ll = ll + CPD.score_relation_token(lib,M.S{sid}.R.eval_spot_token,...
                                              M.S{sid}.R.eval_spot_type); 
    elseif strcmp(M.S{sid}.R.type,'unihist')
       gpos = M.S{sid}.R.gpos;
       pos_token = M.S{sid}.pos_token;
       ll = ll - lib.Spatial.score(gpos,sid); %% since this sample was not used
       ll = ll + mvnormpdfln(gpos(1),pos_token(1),PM.mcmc_prop_relpos_mlty*lib.rel.sigma_x);
       ll = ll + mvnormpdfln(gpos(2),pos_token(2),PM.mcmc_prop_relpos_mlty*lib.rel.sigma_y);
    end
end

%% RELATIONS

function ll = score_relation(M,sid,lib,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        ll = scoreMP(M,lib);
        return 
    end
    
    ll = CPD.score_relation_type(lib,M.S{sid}.R);
    if strcmp(M.S{sid}.R.type,'mid')
       ll = ll + CPD.score_relation_token(lib,M.S{sid}.R.eval_spot_token,...
                                              M.S{sid}.R.eval_spot_type); 
    end
    pos = M.S{sid}.pos_token;
    ll = ll + CPD.score_position(lib,pos,M.S{sid}.R,M.S(1:sid-1));
end

%% SHAPE

function shape_type = propose_shape_type(shape_type,PM)
    [~,dim,n] = size(shape_type);
    assert(dim==2 && n==1);
    shape_type = shape_type + PM.mcmc.prop_shape_sd .* randn(size(shape_type));
end

function shape_token = propose_shape_token(shape_token,PM)
    shape_token = propose_shape_type(shape_token,PM); 
end

function score = score_shape_type(M,sid,bid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    shape_type = M.S{sid}.shapes_type(:,:,bid);
    shape_token = M.S{sid}.shapes_token(:,:,bid);
    id = M.S{sid}.ids(bid);
    lp = CPD.score_shape_type(lib,shape_type,id); 
    ll = CPD.score_shape_token(lib,shape_token,shape_type);
    score = lp + ll;
end

function [score,out] = score_shape_token(M,sid,bid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        [score,out] = scoreMP(M,lib);        
        return 
    end
    shape_type = M.S{sid}.shapes_type(:,:,bid);
    shape_token = M.S{sid}.shapes_token(:,:,bid);
    
    %% account for future positions, which this parameter may influence
    ll_next_pos = 0;
    for i=sid+1:M.ns
       nextS = M.S{i};
       ll_next_pos = ll_next_pos + CPD.score_position(lib,nextS.pos_token,nextS.R,M.S(1:i-1));
    end
    
    %% sum up scores
    score = ll_next_pos + CPD.score_image(M.I,M.pimg) + CPD.score_shape_token(lib,shape_token,shape_type);    
end

%% SCALE

function invscales_type = propose_scale_type(invscales_type,PM)
    assert(isscalar(invscales_type));
    invscales_type = invscales_type + PM.mcmc.prop_scale_sd .* randn(size(invscales_type));
end

function invscales_token = propose_scale_token(invscales_token,PM)
    invscales_token = propose_scale_type(invscales_token,PM);
end

function score = score_scale_type(M,sid,bid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return
    end
    scales_type = M.S{sid}.invscales_type(bid);
    scales_token = M.S{sid}.invscales_token(bid);
    id = M.S{sid}.ids(bid);
    lp = CPD.score_invscale_type(lib,scales_type,id); 
    ll = CPD.score_invscale_token(lib,scales_token,scales_type);
    score = lp + ll;
end

function score = score_scale_token(M,sid,bid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    scales_type = M.S{sid}.invscales_type(bid);
    scales_token = M.S{sid}.invscales_token(bid);
    
    %% account for future positions, which this parameter may influence
    ll_next_pos = 0;
    for i=sid+1:M.ns
       nextS = M.S{i};
       ll_next_pos = ll_next_pos + CPD.score_position(lib,nextS.pos_token,nextS.R,M.S(1:i-1));
    end
    
    score = ll_next_pos + CPD.score_image(M.I,M.pimg) + CPD.score_invscale_token(lib,scales_token,scales_type);
end

%% POSITION

function pos = propose_local_position(pos,PM)
    pos = propose_global_position(pos,PM);
end

function gpos = propose_global_position(gpos,PM)
    gpos = gpos + PM.mcmc.prop_gpos_sd .* randn(size(gpos));
end

function score = score_global_position(M,sid,lib,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    R = M.S{sid}.R;    
    lp = CPD.score_relation_type(lib,R);
    ll = CPD.score_position(lib,M.S{sid}.pos_token,R,M.S(1:sid-1));
    score = lp + ll;
end

function score = score_local_position(M,sid,lib,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    R = M.S{sid}.R;
    ll = CPD.score_position(lib,M.S{sid}.pos_token,R,M.S(1:sid-1));
    
    %% account for future positions, which this parameter may influence
    ll_next_pos = 0;
    for i=sid+1:M.ns
       nextS = M.S{i};
       ll_next_pos = ll_next_pos + CPD.score_position(lib,nextS.pos_token,nextS.R,M.S(1:i-1));
    end
    
    score = CPD.score_image(M.I,M.pimg) + ll + ll_next_pos;
end

%% EVALUATION SPOT RELATION

function eval_spot_type = propose_eval_spot(eval_spot_type,PM)
    eval_spot_type = eval_spot_type + PM.mcmc.prop_relmid_sd .* randn;
end

function score = score_eval_spot_type(M,sid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    R = M.S{sid}.R;
    lp = CPD.score_relation_type(lib,R);
    ll = CPD.score_relation_token(lib,R.eval_spot_token,R.eval_spot_type);
    score = lp + ll;
end    

function score = score_eval_spot_token(M,sid,lib,~,fullscore)
    if ~exist('fullscore','var')
       fullscore = false; 
    end
    if fullscore
        score = scoreMP(M,lib);
        return 
    end
    R = M.S{sid}.R;
    lp = CPD.score_relation_token(lib,R.eval_spot_token,R.eval_spot_type);
    ll = CPD.score_position(lib,M.S{sid}.pos_token,R,M.S(1:sid-1));
    score = lp + ll;
end  