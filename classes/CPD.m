classdef CPD
    % StaticCPD 
    %  Defines the conditional probability distributions that make up the
    %  BPL model
    
    methods (Static)
        
        % -- 
        % Model of global statistics
        % --        
        function [logwt,x_ink,x_cov] = score_stat(libclass,M)
            % Compute the log-score of a motor program M, in terms
            % of the simple ink and coverage statistics.            
            x_ink = get_ink_overlap(M,'stroke');
            x_cov = get_ratio_coverage_type(M);                        
            log_ink = libclass.stat.q_ink_emp.pdf(x_ink) - libclass.stat.q_ink_prior.pdf(x_ink);
            log_canvas = libclass.stat.q_canvas_emp.pdf(x_cov) - libclass.stat.q_canvas_prior.pdf(x_cov);
            logwt = log_ink + log_canvas;
        end        
        
        % ---
        % Number of strokes model (Kappa)
        % ---
        
        function samps = sample_number(libclass,nsamp)
            if ~exist('nsamp','var')
               nsamp = 1; 
            end
            samps = zeros(nsamp,1);
            for i = 1:nsamp
                samps(i) = find(mnrnd(1,libclass.pkappa));
            end
        end
       
        function ll = score_number(libclass,data)
        %   data: vector [n x 1] of number of strokes
        % returns
        %   ll : vector [n x 1] of log-likelihood scores
            assert(isvector(data));
            if any(data>length(libclass.pkappa))
               ll = -inf;
               return
            end
            ll = log(libclass.pkappa(data));
        end 
        
        % ---
        % Sequence of sub-strokes model (z)
        % ---
        
        function nsub = sample_nsub(libclass,ns)  
        % ns    : number of strokes in the character    
            pvec = libclass.pmat_nsub(ns,:);
            nsub = find(mnrnd(1,pvec));
            assert(isscalar(nsub));
            assert(~isnan(nsub));
        end
        
        
        function samps = sample_sequence(libclass,ns,nsub,nsamp)
        % ns    : number of strokes in the character
        % nsub  : number of sub-strokes in the character (optional)      
        % nsamp : number of samples        
            if ~exist('nsamp','var')
               nsamp = 1; 
            end
            samps = cell(nsamp,1);

            % Produce each of the samples
            for i=1:nsamp
                
                %% number of sub-strokes
                if (~exist('nsub','var') || isempty(nsub))
                    nsub = CPD.sample_nsub(libclass,ns);
                end
                
                %% sample the sequence
                sq = zeros(nsub,1);                
                pStart = exp(libclass.logStart);
                sq(1) = find(mnrnd(1,pStart));
                for bid=2:nsub
                    prev = sq(bid-1);
                    sq(bid) = find(mnrnd(1,libclass.pT(prev)));
                end      
                samps{i} = sq;
            end

            if nsamp==1
               samps = samps{1}; 
            end
            
        end        
        
        function [ll,lv] = score_sequence(libclass,ns,ids)
        %  ids: [k x 1] data sequence
        %  ns: [scalar] how many strokes are in the character. If ns<0,
        %    this means to leave out scoring the number of sub-strokes
        %  output is a scalar log-likelihood value
            assert(isvector(ids));
            nsub = length(ids);
            nstate = libclass.N;
            assert(all(ids<=nstate));
            
            %% score the number of sub-strokes
            sz = size(libclass.pmat_nsub);
            if ns < 0 % do not include this part of the score...
               lp = 0;
            elseif ns > sz(1) || nsub > sz(2)
               lp = -inf; 
            else
               lp = log(libclass.pmat_nsub(ns,nsub));
            end                
            
            %% score the sub-stroke ids
            lv = zeros(nsub,1);
            lv(1) = libclass.logStart(ids(1));            
            for i=2:nsub
               lv(i) = libclass.logT(ids(i-1),ids(i)); 
            end                        
            ll = lp + sum(lv);
            assert(~any(isnan(ll)));
            
        end
        
        % ---
        % Shape model (x)
        % ---        
        function bspline_stack = sample_shape_type(libclass,subid)
        % Given a vector of cluster ids (id), sample
        % from the vanilla gaussian model
        % associated with it
        %
        % Input
        %  subid [k x 1]: the id (index) of the sub-stroke
        %
        % Output
        %  bspline_stack: [ncpt x 2 x k] sampled spline
             assert(isvector(subid));
             k = length(subid);
             
             if isunif(libclass)
               bspline_stack = CPDUnif.sample_shape_type(libclass,subid);
               return
             end            
             
             
             Cov = libclass.shape.Sigma(:,:,subid);
             rows_bspline = mvnrnd(libclass.shape.mu(subid,:),Cov);
             ncpt = size(rows_bspline,2)./2;
             bspline_stack = zeros(ncpt,2,k);
             for i=1:k
                 bspline_stack(:,:,i) = reshape(rows_bspline(i,:),[ncpt 2]);
             end
        end
        
        function ll = score_shape_type(libclass,bspline_stack,subid)
        % 
        % Input
        %  subid [k x 1]: the id (index) of the sub-stroke
        %  bspline_stack: (ncpt x 2 x k) shapes of bsplines
        %
        % Output
        %  ll : [k x 1]: vector of scores
            assert(isvector(subid));
            k = length(subid);
            assert(size(bspline_stack,3)==k);
            
            if isunif(libclass)
               ll = CPDUnif.score_shape_type(libclass,bspline_stack,subid);
               return
            end
            
            ll = score_shape_type_helper(libclass,bspline_stack,subid);
        end
        
        
        function bspline_stack = sample_shape_token(libclass,bspline_stack)
        %  bspline_stack: (ncpt x 2 x k) shapes of bsplines    
             sz = size(bspline_stack);
             bspline_stack = bspline_stack + libclass.tokenvar.sigma_shape .* randn(sz);
        end
        
        function ll = score_shape_token(libclass,bspline_stack_token,bspline_stack_type)
        % 
        % Input
        %  bspline_stack_type: (ncpt x 2 x k) shapes of bsplines
        %  bspline_stack_token: (ncpt x 2 x k) shapes of bsplines
        %
        % Output
        %  ll : [k x 1]: vector of scores
            k = size(bspline_stack_token,3);
            assert(k==size(bspline_stack_type,3));
            sz = size(bspline_stack_token);            
            X = bspline_stack_token(:)';
            mu = bspline_stack_type(:)';
            llvec = mvnormpdfln(X,mu,libclass.tokenvar.sigma_shape*ones(size(X)));
            llvec = reshape(llvec,sz);            
            ll = sum(llvec,1);
            ll = sum(ll,2);
            ll = squeeze(ll);                 
        end
        
        function ll = score_shape_marginalize(libclass,bspline_stack_token,subid)
        % 
        % Input
        %  subid [k x 1]: the id (index) of the sub-stroke
        %  bspline_stack: (ncpt x 2 x k) shapes of bsplines
        %
        % Output
        %  ll : [k x 1]: vector of scores
            assert(isvector(subid));
            k = length(subid);
            assert(size(bspline_stack_token,3)==k);
            
            if isunif(libclass)
               ll = CPDUnif.score_shape_marginalize(libclass,bspline_stack_token,subid);
               return
            end
            
            ll = score_shape_type_helper(libclass,bspline_stack_token,subid,libclass.tokenvar.sigma_shape^2);
        end
                
        % --
        % Scale model (y)
        % --        
        
        function invscales = sample_invscale_type(libclass,subid)
            % subid: [k x 1] vector of sub-stroke ids
            %
            assert(isvector(subid));
            
            if isunif(libclass)
               invscales = CPDUnif.sample_invscale_type(libclass,subid);
               return
            end
            
            theta = libclass.scale.theta(subid,:);    
            invscales = gamrnd(theta(:,1),theta(:,2));            
        end
        
        function lprob = score_invscale_type(libclass,invscales_type,subid)
            % subid: [k x 1] vector of sub-stroke ids
            assert(isvector(invscales_type));
            assert(isvector(subid));
            k = length(subid);
            assert(numel(invscales_type)==k);    
            
            if isunif(libclass)
                lprob = CPDUnif.score_invscale_type(libclass,invscales_type,subid);
                return
            end
            
            theta = libclass.scale.theta(subid,:);    
            prob = gampdf(invscales_type(:),theta(:,1),theta(:,2));
            lprob = log(prob);
            assert(numel(lprob)==k);
        end        
        
        function invscales_token = sample_invscale_token(libclass,invscales_type)
        % Gaussian noise, but don't allow negative scales.
        % Sampling is done by rejection sampling
            sz = size(invscales_type);
            invscales_token = invscales_type + libclass.tokenvar.sigma_invscale.*randn(sz);
            while any(isinf(CPD.score_invscale_token(libclass,invscales_token,invscales_type)))
                invscales_token = invscales_type + libclass.tokenvar.sigma_invscale.*randn(sz);
            end
        end
        
        function ll = score_invscale_token(libclass,invscales_token,invscales_type)
            % Score the token-level inverse scales.
            % Computes the right normalization constant,
            % since negative values are not allowed.
            %
            % Input
            %  invscales_token: [n x 1] vector
            %  invscales_type: [n x 1] vector
            %
            % Output
            %  ll: [n x 1] vector of log-likelihood        
            ll = mvnormpdfln(invscales_token(:)',invscales_type(:)',libclass.tokenvar.sigma_invscale);
            ll = ll(:);
            ll(invscales_token <= 0) = -inf; % don't allow invscales that are negative
            
            % correction for positive only invscales
            p_below = normcdf(0,invscales_type,libclass.tokenvar.sigma_invscale);
            p_above = 1-p_below; 
            ll = ll - log(p_above);
        end
        
        % --
        % Relation model (R)
        % ---
        function R = sample_relation_type(libclass,previous_strokes)
            nprev = length(previous_strokes);
            stroke_num = nprev+1;
            types = {'unihist','start','end','mid'};
            ncpt = libclass.ncpt;
            if nprev == 0
                indx = 1;
            else
                out = mnrnd(1,libclass.rel.mixprob);
                if any(isnan(out))
                   error('improper relation mixing prob.'); 
                end
                indx = find(out);                
            end
            type = types{indx};
            switch (type)
                case 'unihist'
                    gpos = libclass.Spatial.sample(stroke_num);
                    R = RelationIndependent(type,nprev,gpos);
                case {'start','end'}
                    attach_spot = randint(1,1,[1 nprev]);
                    R = RelationAttach(type,nprev,attach_spot);
                case 'mid'
                    attach_spot = randint(1,1,[1 nprev]);
                    nsub = previous_strokes{attach_spot}.nsub;
                    subid_spot = randint(1,1,[1 nsub]);                    
                    R = RelationAttachAlong(type,nprev,attach_spot,nsub,subid_spot,ncpt);
                    [~,lb,ub] = bspline_gen_s(ncpt,1);
                    R.eval_spot_type = lb + rand*(ub-lb);
                otherwise
                    error('invalid relation');
            end    
        end        
        
        function ll = score_relation_type(libclass,R)
            stroke_num = R.nprev+1;
            types = {'unihist','start','end','mid'};            
            logp = log(libclass.rel.mixprob);            
            if stroke_num > 1 % prior probability on a particular relation
                indx = strcmp(R.type,types);
                ll = logp(indx);
            else % there was only one option
                ll = 0; 
            end
            switch (R.type)
                case 'unihist'                    
                    ll = ll + libclass.Spatial.score(R.gpos,stroke_num);
                case {'start','end'}
                    ll = ll - log(R.nprev);
                case 'mid'
                    ll = ll - log(R.nprev);
                    ll = ll - log(R.nsub);
                    eval_missing = isempty(R.eval_spot_type);
                    if ~eval_missing
                        ncpt = libclass.ncpt;
                        [~,lb,ub] = bspline_gen_s(ncpt,1);
                        ll = ll - log(ub-lb);
                        if (R.eval_spot_type < lb || R.eval_spot_type > ub)
                            ll = -inf;
                        end
                    end
                otherwise
                    error('invalid relation');
            end 
        end
        
        function eval_spot_token = sample_relation_token(libclass,eval_spot_type)
        % sample an attachment, but within the bounds as defined by the model    
            eval_spot_token = eval_spot_type + libclass.tokenvar.sigma_attach.*randn;
            while isinf(CPD.score_relation_token(libclass,eval_spot_token,eval_spot_type))
                eval_spot_token = eval_spot_type + libclass.tokenvar.sigma_attach.*randn;
            end
        end
        
        function ll = score_relation_token_approx_marginalize(libclass,eval_spot_token)
        % score the token evaluation spot, for attached relations, assuming 
        % we are missing the type-level spot. 
            assert(isscalar(eval_spot_token));
            ncpt = libclass.ncpt;
            [~,lb,ub] = bspline_gen_s(ncpt,1);
            if (eval_spot_token < lb || eval_spot_token > ub)
                ll = -inf;
                return;
            end
            
            % cached structure that evaluates this complex integral 
            ll = libclass.score_eval_marg(eval_spot_token);            
            % ll = -log(ub-lb);
        end             
        
        function ll = score_relation_token(libclass,eval_spot_token,eval_spot_type)    
        % score the token evaluation spot, for attached relations, assuming 
        % the type-level spot is given
            assert(isscalar(eval_spot_token));
            assert(~isempty(eval_spot_type));
            ncpt = libclass.ncpt;
            [~,lb,ub] = bspline_gen_s(ncpt,1);
            if (eval_spot_token < lb || eval_spot_token > ub)
                ll = -inf;
                return
            end            
            ll = mvnormpdfln(eval_spot_token,eval_spot_type,libclass.tokenvar.sigma_attach);
            
            % correction for bounds
            p_within = normcdf(ub,eval_spot_type,libclass.tokenvar.sigma_attach) - normcdf(lb,eval_spot_type,libclass.tokenvar.sigma_attach); 
            ll = ll - log(p_within);
        end
        
        % ---
        % Local position model (L)
        % ---        
        function pos = sample_position(libclass,R,previous_strokes)
            % Given a relation R and the previous strokes,
            %  sample where the position of this stroke should be            
            base = getAttachPoint(R,previous_strokes);
            pos = [normrnd(base(1),libclass.rel.sigma_x) normrnd(base(2),libclass.rel.sigma_y)];
        end
        
        function ll = score_position(libclass,pos,R,previous_strokes)
            base = getAttachPoint(R,previous_strokes);
            ll = mvnormpdfln(pos(1),base(1),libclass.rel.sigma_x) + mvnormpdfln(pos(2),base(2),libclass.rel.sigma_y);
        end
        
        % ---
        % Image blur model (\sigma_b)
        % ---
        function blur_sigma = sample_image_blur(PM)
            % sample uniformly in [min_blur_sigma, max_blur_sigma]
            rg = PM.max_blur_sigma - PM.min_blur_sigma;
            blur_sigma = rg*rand + PM.min_blur_sigma;
        end
        
        function ll = score_image_blur(blur_sigma,PM)
            if (blur_sigma > PM.max_blur_sigma || blur_sigma < PM.min_blur_sigma)
                ll = -inf;
                return
            end                
            rg = PM.max_blur_sigma - PM.min_blur_sigma;
            ll = -log(rg);
        end        
        
        % ---
        % Pixel noise model (\epsilon)
        % ---
        function epsilon = sample_image_noise(PM)
            % sample uniformly in [min_epsilon, max_epsilon]
            rg = PM.max_epsilon - PM.min_epsilon;
            epsilon = rg*rand + PM.min_epsilon;
        end
        
        function ll = score_image_noise(epsilon,PM)
            if (epsilon > PM.max_epsilon || epsilon < PM.min_epsilon)
               ll = -inf;
               return
            end
            rg = PM.max_epsilon - PM.min_epsilon;
            ll = -log(rg);
        end
        
        % ---
        % Image model (I)
        % ---        
        function I = sample_image(pimg)
        %  sample a binary image
            I = binornd(1,pimg);
            I = logical(I);
        end
        
        function ll = score_image(I,pimg)
        %  score the image model
            assert(~isempty(I));
            if islogical(I)
                on = I;
                off = ~I;
            elseif isnumeric(I) % if we having missing data
                missing = isinf(I); 
                on = I > 0.5;
                off = I < 0.5;
                on(missing) = false;
                off(missing) = false;
            end  
            prob_on = pimg;
            ll_on = sum(log(prob_on(on)));
            ll_off = sum(log(1-prob_on(off)));
            ll = ll_on + ll_off;            
        end
        
%         function ll = score_image(I,pimg)
%         %  score the image model
%             assert(~isempty(I));
%             prob_on = pimg;
%             on = I;
%             off = ~I;
%             ll_on = sum(log(prob_on(on)));
%             ll_off = sum(log(1-prob_on(off)));
%             ll = ll_on + ll_off;            
%         end
        
        % ---
        % Affine model (A)
        % ---
        function sample_A = sample_affine(libclass,nsamp)
            % affine transformation [x-scale,y-scale,x-translate,y-translate]
            % the translation is relative to the center of mass
            % (x-scale and y-scale cannot be 0 or negative)
            %
            % sample_A: [nsamp x 4] each row is a sampled affine warp
            if ~exist('nsamp','var')
                nsamp = 1; 
            end
            
            sample_A = zeros(nsamp,4);
            
            % sample the image scale
            m_scale = libclass.affine.mu_scale;
            S_scale = libclass.affine.Sigma_scale;
            sample_A(:,1:2) = mvnrnd(m_scale,S_scale,nsamp);

            % sample the translation 
            m_x = libclass.affine.mu_xtranslate;
            m_y = libclass.affine.mu_ytranslate;
            s_x = libclass.affine.sigma_xtranslate;
            s_y = libclass.affine.sigma_ytranslate;
            sample_A(:,3) = normrnd(m_x,s_x,nsamp,1);
            sample_A(:,4) = normrnd(m_y,s_y,nsamp,1);

            % scales should not be zero or negative            
            if ~(all(vec(sample_A(:,1:2)) > 0))
               warning('sampled scale variable is less than zero'); 
            end
            
        end       
        
        % NOTE: doesn't bother to truncate Gaussian properly,
        %       since the variance is so tight it won't matter
        function ll = score_affine(libclass,A)
            if isempty(A)
               A = [1; 1; 0; 0]; 
            end
            % score the affine warp
            % (the number of elements is 4)            
            assert(isvector(A) && numel(A)==4);
            A = A(:);
            if (A(1) <= 0 || A(2) <= 0)
               ll = -inf;
               return
            end
            
            ll = 0;

            % sample the image scale
            m_scale = libclass.affine.mu_scale(:);
            S_scale = libclass.affine.Sigma_scale;
            ll = ll + mvnormpdfln(A(1:2),m_scale,[],S_scale);

            % sample the translation 
            m_x = libclass.affine.mu_xtranslate;
            m_y = libclass.affine.mu_ytranslate;
            s_x = libclass.affine.sigma_xtranslate;
            s_y = libclass.affine.sigma_ytranslate;
            ll = ll + mvnormpdfln(A(3),m_x,s_x);
            ll = ll + mvnormpdfln(A(4),m_y,s_y);
        end
        
    end    
end

function y = isunif(libclass)
    y = any(isnan(libclass.shape.mu));
end

% fast vectorized helper MVN likelihood
function ll = score_shape_type_helper(libclass,bspline_stack,subid,regcov)    
    mymu = libclass.shape.mu(subid,:);
    mysd = libclass.shape.vsd(subid,:);
    if exist('regcov','var')
       mysd = sqrt( (mysd.^2) + regcov );
    end    
    sz = size(mymu);    
    vbspline = reshape(bspline_stack,sz(2),sz(1));
    vbspline = vbspline';    
    llv = mvnormpdfln(vbspline(:)',mymu(:)',mysd(:)');       
    llv = reshape(llv,sz);
    ll = sum(llv,2);    
end