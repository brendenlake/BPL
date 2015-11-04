classdef SpatialHist
    % SpatialHist 
    %
    %   A 2D plane is divided into an evenly spaced grid,
    %   where a square is chosen randomly and then points
    %   are chosen uniformly from the square
    %    
    properties
      logpYX
      xlab
      ylab
      rg_bin % length, in pixels, of a side of a bin
      prior_count
    end
    
    methods
        
        %
        % Build a 2D histogram model of the data
        %
        % Input
        %   data: [n x 2 scalar] data to model
        %   xlim: [1 x 2] range of x-dimension
        %   ylim: [1 x 2] range of y-dimension
        %   nbin_per_side: number of bins per dimension 
        %   prior_count: prior counts in each cell (not added to edge cells)
        function  this = SpatialHist(data,xlim,ylim,nbin_per_side,prior_count)
          
            if ~exist('prior_count','var')
               prior_count = 0; 
            end
            [ndata,dim] = size(data);
            assert(numel(xlim)==2);
            assert(numel(ylim)==2);
            assert(dim==2);

            % Compute the "edges" of the histogram
            xrg = xlim(2)-xlim(1);
            yrg = ylim(2)-ylim(1);
            xint = xrg ./ nbin_per_side;
            yint = yrg ./ nbin_per_side;
            xtick = xlim(1):xint:xlim(2);
            ytick = ylim(1):yint:ylim(2);            
            assert(length(xtick)-1==nbin_per_side);
            assert(length(ytick)-1==nbin_per_side);
            edges = {xtick; ytick}; % for histogram function

            % store important information about the bins
            this.rg_bin(1) = xint;
            this.rg_bin(2) = yint;
            this.xlab = xtick;
            this.ylab = ytick;
            
            % Compute the histogram
            N = myhist3(data,edges);
            diff=ndata-sum(N(:));
            if diff>0
               fprintf(1,'Warning: %d position points are out of bounds\n',diff); 
            end

            % Add in the prior counts
            N = N';            
            N = N + prior_count;
            logN = log(N);

            % Convert to probability distribution
            logpN = logN - logsumexp(logN(:));
            assert(aeq(sum(exp(logpN(:))),1));

            this.logpYX = logpN;
            this.xlab = xtick;
            this.ylab = ytick;
            this.prior_count = prior_count;  
            
        end
        
        %
        % Sample from a 2D histogram model
        %
        % Input
        %  nsamp: number of samples
        %
        % Output
        %   samples: [n x 2 scalar] samples
        %   yi: [n x 1] y-bin index
        %   xi: [n x 1] x-bin index
        function [samples,yi,xi] = sample(this,nsamp)
            
            % Pick which bins the samples are from
            logpvec = this.logpYX(:);
            pvec = exp(logpvec);
            pvec = pvec ./ sum(pvec);
            lin = zeros(nsamp,1);
            for i=1:nsamp
               lin(i) = find(mnrnd(1,pvec));
            end

            % Retrieve the [y,x] indices of these bins
            [yi,xi] = ind2sub(size(this.logpYX),lin);

            % Retrieve the edges for each of these bins
            xmin = this.xlab(xi);
            ymin = this.ylab(yi);
            xmax = this.xlab(xi+1);
            ymax = this.ylab(yi+1);

            % Sample from a uniform distribution in each of the bins
            xsamp = vec(xmax-xmin).*rand(nsamp,1) + xmin(:);
            ysamp = vec(ymax-ymin).*rand(nsamp,1) + ymin(:);
            samples = [xsamp ysamp];
            
        end
        
        %
        % Compute the log-likelihood of data under 
        % a 2D histogram model
        %
        % Input
        %   data: [n x 2 scalar] data to model
        %   H: histogram model
        %
        % Output
        %   ll: [n x 1] log-likelihood scores
        %
        function ll = score(this,data)
            
            % Compute bin in histogram
            [n,dim] = size(data);
            edges = {this.xlab; this.ylab};
            
            mylogpYX = this.logpYX;
            
            % fast classification
            ll = fast_hclassif(data,mylogpYX,edges);
            
%             % incremental classification
%             ll2 = zeros(n,1);
%             xid = zeros(n,1);
%             yid = zeros(n,1);
%             parfor i=1:n
%                 [ll2(i),xid(i),yid(i)] = hclassif(data(i,:),mylogpYX,edges);
%             end            
%             assert(aeq(ll,sum(ll2)));
        
            % Adjust log-likelihoods to account
            % for the uniform component of the data
            ll = ll-n*log(this.rg_bin(1))-n*log(this.rg_bin(2));    
            assert(~any(isnan(ll)));
            
        end
        
        %
        % Compute the log-likelihood of data under 
        % a 2D histogram model
        %
        % Input
        %   data: [n x 2 scalar] data to model
        %   H: histogram model
        %
        % Output
        %   id: [n x 2] x and y id of each point in bins
        %   ll: [n x 1] log-likelihood of each point
        function [id,ll] = get_id(this,data)
            [n,dim] = size(data);
            edges = {this.xlab; this.ylab};
            ll = zeros(n,1);
            xid = zeros(n,1);
            yid = zeros(n,1);
            mylogpYX = this.logpYX;
            for i=1:n
                [ll(i),xid(i),yid(i)] = hclassif(data(i,:),mylogpYX,edges);
            end        
            id = [xid yid];
            ll = ll-log(this.rg_bin(1))-log(this.rg_bin(2)); 
            assert(~any(isnan(ll)));
        end
        
        % Visualize the learned position model
        function plot(this)
            pYX = exp(this.logpYX);
            img = pYX ./ max(pYX(:));            
            xlim = [this.xlab(1) this.xlab(end)];
            ylim = [this.ylab(1) this.ylab(end)];            
            image(xlim,[ylim(2) ylim(1)],repmat(img,[1 1 3]));
            set(gca,'XTick',[],'YTick',[]);
        end

        
    end
    
end

%
% Compute the log-likelihood of the point "pt"
%
% pt: [1 x 2]
% logpYX:  
% edges: 
%
function [logprob,xid,yid] = hclassif(pt,logpYX,edges)
    N = myhist3(pt,edges);
    N = N';    
    N = logical(N);
    [yid,xid] = find(N);
    if (sum(N(:))==0);
        logprob = -inf;
        xid = 1;
        yid = 1;
        return
    end
    logprob = logpYX(N);
    assert(~isnan(logprob));
end

%
% Vectorized version of hclassif
% 
function logprob = fast_hclassif(pt,logpYX,edges)
    npt = size(pt,1);
    N = myhist3(pt,edges);
    N = N';    
    sumN = sum(sum(N));
    MTPL = N.*logpYX;
    MTPL(N == 0) = 0; % incase a position is not valid
    logprob = sum(sum(MTPL));
    if (sumN < npt)
        logprob = -inf;
    end
    assert(~isnan(logprob));
end

% Modified histogram function, where
% datapoints on the edge are mapped to the last cell, 
% not their own cell
function N = myhist3(data,edges)
    
    % cluster with histogram function
    N = hist3(data,'Edges',edges);
       
    % move the last row/col to the second to last
    lastcol = N(:,end);
    lastrow = N(end,:);
    last = N(end,end);
    N(:,end-1) = N(:,end-1) + lastcol;
    N(end-1,:) = N(end-1,:) + lastrow;
    N(end-1,end-1) = N(end-1,end-1) + last;
    
    % Delete last row and column
    N(end,:) = [];
    N(:,end) = [];
end