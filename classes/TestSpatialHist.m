%
% testCase = TestSpatialHist;
% res = run(testCase);
%
classdef TestSpatialHist < matlab.unittest.TestCase
   % Tests for SpatialHist class
   %   Also looks at ability of the model to capture
   %    a mixture of Gaussians
   %
 
   properties
      H
      syndata
      xlim
      ylim
   end
    
   methods (TestMethodSetup)
       
        function createModels(this)            
            nbin_per_side = 50;
            prior_count = 0.1;
            xlm = [-10 10];
            ylm = [-10 10];

            mu1 = [-5 0];
            mu2 = [0 5];
            Sigma = eye(2);
            n = 1000;
            data1 = mvnrnd(mu1,Sigma,n);
            data2 = mvnrnd(mu2,Sigma,n);
            data = [data1; data2];

            % data = repmat([10 10],[1000 1]);

            figure(1);
            subplot(2,1,1);
            plot(data(:,1),data(:,2),'r.');
            xlim(xlm);
            ylim(ylm);
            title('original data');

            this.H = SpatialHist(data,xlm,ylm,nbin_per_side,prior_count);
            this.syndata = this.H.sample(1000);

            subplot(2,1,2);
            plot(this.syndata(:,1),this.syndata(:,2),'b.');
            xlim(xlm);
            ylim(ylm);
            title('reconstructed data');
            this.xlim = xlm;
            this.ylim = ylm;
        end
   end
    
   methods (Test)
       
       function dualMethodLL(this)
          % check two different ways of computing likelihood
          ll = this.H.score(this.syndata);
          [~,ll2] = this.H.get_id(this.syndata);
          ll2 = sum(ll2);
          this.verifyEqual(ll,ll2,'AbsTol',1e-2);
       end
       
       function validDensity(this)
           % numerically check the normalizing constant of the density
           
           nsamp = 1e4;
           
           xlim = this.xlim;
           ylim = this.ylim;
           xrg = xlim(2)-xlim(1);
           yrg = ylim(2)-ylim(1);
           area = xrg * yrg;

           x = xlim(1) + rand(nsamp,1)*xrg;
           y = ylim(1) + rand(nsamp,1)*yrg;
           D = [x y];

           [~,ll] = this.H.get_id(D);   
           ltot = logsumexp(ll(:));
           tsum = exp(ltot);
           tot = area * tsum ./ nsamp;

           fprintf(1,'  average score %s\n',num2str(tot,3));
           this.verifyEqual(tot,1,'AbsTol',1e-1);
       end
       
       
   end
   
   
end