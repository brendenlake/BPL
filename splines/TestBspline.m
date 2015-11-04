% To run unit tests:
%
% testCase = TestBspline;
% res = run(testCase);
%
classdef TestBspline < matlab.unittest.TestCase
    %
    % Test the spline scripts
    %
    
    properties
       cpts 
    end
    
    
    methods (TestClassSetup)
        
       function setup(this)
            cpts = [0 2;
                        0 4;
                        0 6;
                        0 8;
                        0 10];
           
            cpts = cpts + 10*rand*randn(size(cpts));
            this.cpts = cpts;
       end
       
    end
    
    methods (Test)
        
        function recover(this)
            nland = size(this.cpts,1);
            sval = bspline_gen_s(nland);
            X = bspline_eval(sval,this.cpts); 
            
            % Test fit of control points
            P = bspline_fit(sval,X,nland);
            if aeq(P,this.cpts)
               fprintf(1,'Landmarks recovered exactly! \n'); 
            else
               fprintf(1,'FAILED TO RECOVER landmarks \n'); 
            end
            
            % plot the spline
            figure(1);
            clf
            hold on
            x = X(:,1);
            y = X(:,2);
            plot(x,y,'b-');
            plot(this.cpts(:,1),this.cpts(:,2),'r.','MarkerSize',16);
            
        end
    end
    
    
end