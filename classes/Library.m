classdef Library < BetterHandle
    % LIBRARY... hyper-parameters for the BPL model
    
    properties (SetAccess = private)
        shape = struct % .mu,Sigma,mixprob
        scale = struct % .theta,minfreq
        rel = struct   % .sigma_x,sigma_y,mixprob
        Spatial        % reference to class
        tokenvar = struct % .sigma_shape,sigma_invscale
        affine = struct   % .mu_scale,mu_xtranslate,mu_ytranslate,
                          %   Sigma_scale,sigma_xtranslate,sigma_ytranslate
        stat = struct     % .q_ink_emp,.q_ink_prior,.q_canvas_emp,.q_canvas_prior
        
        % misc. items
        logT % [N x N] log-transition probabilities
        logStart % [N x 1] log-prob of beginning in each state
        pkappa % [n x 1] probability of each number of strokes (starts at 1)
        pmat_nsub % [ns x nsub] each row is probability of each number of sub-strokes, for a charater with
                  % a certain number of strokes
        newscale
        smooth_bigrams
        diagSigma = true % Are the primivites diagonal Gaussians?
    end    
    
    properties (SetAccess = private, GetAccess = private)    
        % cached marignal likelihood for attach relation
        eval_int = 1e-2;
        int_eval_marg
        prob_eval_marg             
    end
    
    properties (Dependent = true)
        ncpt % number of control points
        N    % number of primitives
        endstate % the special endstate in the Markov process
    end
    
    methods
       
        % Constructor
        %        
        function this = Library(as)
            
            if ~exist('as','var') || isempty(as)
               return; 
            end
                                             
            % Learned model of stroke number
            load(as.number,'pkappa');
            this.pkappa = pkappa;
            
            % Learned model of number of sub-strokes
            load(as.nsub,'pmat_nsub');
            this.pmat_nsub = pmat_nsub;
                            
            % Learned primitives and transitions
            load(as.primitives,'vHMM');
            nv = numel(vHMM);
            bestHMM = [];
            for i=1:nv
                if (vHMM{i}.nclass == as.nprim)
                   bestHMM = vHMM{i}; 
                end
            end          
            if isempty(bestHMM)
               error('primitive model not found'); 
            end
            BH = bestHMM;
            this.shape.mu = BH.mu;
            this.shape.Sigma = BH.RawSigma;
            [N,dim] = size(this.shape.mu); % Compute diag. form of standard deviations
            this.shape.vsd = zeros(N,dim);
            for i=1:N
                this.shape.vsd(i,:) = diag(this.shape.Sigma(:,:,i))';
            end
            this.shape.vsd = sqrt(this.shape.vsd);
            freq = sum([BH.start_count(:)'; BH.bigrams]);
            this.shape.mixprob = freq ./ sum(freq);
            this.shape.freq = freq;
            this.scale.theta = BH.theta;
           
            % Learned transition probs.
            this.logT = bestHMM.logTmat;
            this.logStart = BH.logStart;
            this.smooth_bigrams = bestHMM.smooth_bigrams;
            
            % Learned relation model
            load(as.relations,'sigma_x','sigma_y','mixprob');
            this.rel.sigma_x = sigma_x;
            this.rel.sigma_y = sigma_y;
            this.rel.mixprob = mixprob;
            clear mixprob;
            
            % Learned position model
            load(as.position,'pModel');
            this.Spatial = pModel.copy();
            
            % Learned token-level variance
            load(as.tokenvar,'sigma_shape','sigma_invscale','sigma_seval');
            this.tokenvar.sigma_shape = sigma_shape;
            this.tokenvar.sigma_invscale = sigma_invscale;
            this.tokenvar.sigma_attach = sigma_seval;
            
            % Learned affine model
            load(as.affine,'affine');
            this.affine = affine;
            
            % Learn stats model
            if isfield(as,'stats')
                load(as.stats,'fit_ink_emp','fit_ink_M','fit_canvas_emp','fit_canvas_M');
                this.stat.q_ink_emp = fit_ink_emp;
                this.stat.q_ink_prior = fit_ink_M;
                this.stat.q_canvas_emp = fit_canvas_emp;
                this.stat.q_canvas_prior = fit_canvas_M;
            end
            
            % Learned ink model
            this.check_consistent();            
            ps = defaultps_clustering;
            this.newscale = ps.newscale_ss;
            
            % Caching structure
            this.create_eval_list;
        end
        
        function legacylib(this,oldlib)
            this.shape = oldlib.shape;
            this.scale = oldlib.scale;
            this.rel = oldlib.rel;
            this.Spatial = oldlib.Spatial.copy();
            this.tokenvar = oldlib.tokenvar;
            this.affine = oldlib.affine;
            this.logT = oldlib.logT;
            this.logStart = oldlib.logStart;
            this.pkappa = oldlib.pkappa;
            this.pmat_nsub = oldlib.pmat_nsub;
            this.smooth_bigrams = oldlib.smooth_bigrams;
            this.diagSigma = oldlib.diagSigma;            
            this.create_eval_list;
        end
        
        
        function restrict_library(this,keep)
        % remove primitives from library, except for those in "keep"
        % keep: [N x 1 logical] true if we want to keep primitive
            
            bool_keep = false(this.N,1);
            bool_keep(keep) = true;
            bool_rmv = ~bool_keep;
            frmv = find(bool_rmv);
            
            this.shape.mu(bool_rmv,:) = [];
            this.shape.vsd(bool_rmv,:) = [];
            this.shape.Sigma(:,:,bool_rmv) = [];
            this.shape.mixprob(bool_rmv) = [];
            this.shape.freq(bool_rmv) = [];
            this.shape.mixprob = this.shape.mixprob ./ sum(this.shape.mixprob);
            
            % renormalize transition probabilities
            T = exp(this.logT);
            T(frmv,:) = [];
            T(:,frmv) = [];
            nt = size(T,1);
            for i=1:nt
               row = T(i,:);
               ss = sum(row);
               if ss > 0 % case where we have no probability mass
                    T(i,:) = row ./ ss;
               else
                    T(i,:) = 0; 
               end
            end
            this.logT = log(T);
            
            % renomalize start prob.
            ST = exp(this.logStart);
            ST(frmv) = [];
            ST = ST ./ sum(ST);
            this.logStart = log(ST); 
            
            % select relevant scales
            this.scale.theta(bool_rmv,:) = [];
            
        end      
        
        function ncpt = get.ncpt(this)
        % get the number of control points    
           dim = size(this.shape.mu,2);
           ncpt = dim/2;
        end
        
        function N = get.N(this)
        % get the number of primitives    
           N = size(this.shape.mu,1);            
        end                
        
        function check_consistent(this)
        % check consistency of the number of primitives in the model
            N = this.N;
            ncpt = this.ncpt;
            assert( size(this.shape.Sigma,1) == ncpt*2);
            assert( size(this.logT,1) == N);
            assert( numel(this.logStart) == N);
            assert( numel(this.shape.mixprob) == N);
            assert( numel(this.shape.freq) == N);
            assert( size(this.shape.vsd,1) == N);
            assert( size(this.scale.theta,1) == N);              
            assert( aeq(sum(exp(this.logStart)),1) );
            for i=1:N
               assert( aeq(sum(this.pT(i)),1) ); 
            end
        end        
        
        function setfield(this,field,val)
        % since fields are private, we must do this to set them    
           this.(field) = val;  
        end
        
        function p = pT(this,prev_state)
        % get the probability of transitioning to a new state,
        %  given your current state is "prev_state"
            assert(isscalar(prev_state));
            logR = this.logT(prev_state,:);
            R = exp(logR);
            R = R(:)';
            p = R ./ sum(R);
        end        
        
        function ll = score_eval_marg(this,eval_spot_token)            
            % score the marginal likelihood of an attachment
            % position
            yi = interp1(this.int_eval_marg,this.prob_eval_marg,eval_spot_token);
            ll = log(yi);
        end
       
    end
    
    methods (Access = private)

        %% create caching structure for efficiently computing 
        % marginal likelihood of attachment
        function create_eval_list(this) 
            [~,lb,ub] = bspline_gen_s(this.ncpt,1); 
            x = lb:this.eval_int:ub;
            nint = numel(x);
            logy = zeros(nint,1);                        
            for i=1:nint
               logy(i) = this.score_relation_eval_marginalize_exact(x(i));
            end
            this.int_eval_marg = x;
            this.prob_eval_marg = exp(logy);            
        end               
        
        %% P(eval_spot_attach_token), marginalizing out the type level.
        % This is computed by approximate integration.
        function ll = score_relation_eval_marginalize_exact(this,eval_spot_token)
        
            assert(isscalar(eval_spot_token));
            ncpt = this.ncpt;
            [~,lb,ub] = bspline_gen_s(ncpt,1);                        
            fll = @(x) exp(CPD.score_relation_token(this,eval_spot_token,x) - log(ub-lb));
                       
            % approximate marginalize
            x = lb:(this.eval_int/100):ub;
            y = fll(x);
            Z = trapz(x,y);            
            if (eval_spot_token < lb || eval_spot_token > ub)
                ll = -inf;
                return;
            end
            ll = log(Z);
            
        end
        
        
    end
    
    
end

