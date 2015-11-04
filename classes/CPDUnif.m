classdef CPDUnif
   %
   % Replace CPDs that govern sub-stroke shape and scales
   % with just uniform distributions
   %  
   
   methods (Static)
              
        % ---
        % Shape model (X)
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
             assert(all(isnan(libclass.shape.mu)));
             assert(isvector(subid));
             k = length(subid);
             sz = [libclass.ncpt 2 k];
             lim = [-libclass.newscale libclass.newscale];
             rg = lim(2)-lim(1);
             bspline_stack = rg.*rand(sz) + lim(1);                   
        end
        
        function ll = score_shape_type(libclass,bspline_stack,subid)
        % 
        % Input
        %  subid [k x 1]: the id (index) of the sub-stroke
        %  bspline_stack: (ncpt x 2 x k) shapes of bsplines
        %
        % Output
        %  ll : [k x 1]: vector of scores
            assert(all(isnan(libclass.shape.mu)));
            assert(isvector(subid));            
            lim = [-libclass.newscale libclass.newscale];
            rg = lim(2)-lim(1);
            [ncpt,dim,k] = size(bspline_stack);
            nn = ncpt * dim;                       
            ll = -nn.*log(rg);
            ll = repmat(ll,[k 1]);
            for i=1:k
                below = bspline_stack(:,:,i) < lim(1);
                above = bspline_stack(:,:,i) > lim(2);
                out = below | above;
                if any(out(:))
                   ll(i) = -inf; 
                end
            end         
        end
       
        function ll = score_shape_marginalize(libclass,bspline_stack_token,subid)
        %
        % This score is only approximate...
        %
        % Input
        %  subid [k x 1]: the id (index) of the sub-stroke
        %  bspline_stack: (ncpt x 2 x k) shapes of bsplines
        %
        % Output
        %  ll : [k x 1]: vector of scores
            assert(all(isnan(libclass.shape.mu)));
            assert(isvector(subid));
            k = length(subid);
            assert(size(bspline_stack_token,3)==k);
            ll = CPDUnif.score_shape_type(libclass,bspline_stack_token,subid);
        end
        
        
        % --
        % Scale model (Y)
        % --        
        
        function invscales = sample_invscale_type(libclass,subid)
            % subid: [k x 1] vector of sub-stroke ids
            % 
            % sample from unif(0,1)
            %
            assert(all(isnan(libclass.shape.mu)));
            assert(isvector(subid));
            k = length(subid);            
            invscales = rand(k,1);                    
        end
        
        function lprob = score_invscale_type(libclass,invscales_type,subid)
            % subid: [k x 1] vector of sub-stroke ids
            assert(all(isnan(libclass.shape.mu)));
            assert(isvector(invscales_type));
            assert(isvector(subid));
            k = length(subid);
            assert(numel(invscales_type)==k);
            lprob = zeros(k,1); % -log(1), for uniform(0,1)    
            
            % check constraints
            out = invscales_type > 1 | invscales_type < 0;
            lprob(out) = -inf;            
            assert(numel(lprob)==k);
        end
        
       
       
   end
    
    
end