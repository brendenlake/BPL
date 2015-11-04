% To run unit tests:
%
% testCase = TestMotorProgramFit;
% res = run(testCase);
%
classdef TestMotorProgramFit < matlab.unittest.TestCase
    
   properties
      M
      samples
   end
    
   methods (TestMethodSetup)
       
        function createModel(this)            
            
            ps = defaultps;
            load(ps.libname,'lib');
            motor_program = generate_all_rel_character(lib);
            M = motor_program();
            while M.ink_off_page
                motor_program = generate_all_rel_character(lib);
                M = motor_program();      
            end
            
            % make samples
            nsamp_mcmc = 5;
            samples = mcmc_all(M,lib,nsamp_mcmc,'type');
                        
            this.M = M;
            this.samples = samples;
        end
        
   end
    
   methods (Test)
       
       
       % returns an COPY of this class
       % for a particular type-level sample "indx", with all
        % of the variables instantiated
       function checkSampleTypeCopy(this)           
           Mfit = MotorProgramFit(this.samples);
           
           ns = Mfit.ns;
           nsamp = Mfit.nsamples;
           L = Mfit.TypeS;
           assert(size(L,1)==nsamp && size(L,2)==ns);
           
           val = true;
           for i=1:nsamp
               tic
               MS = Mfit.getSampleTypeCopy(i);
               t=toc;
               for sid=1:ns
                   this.samples{i}.S{sid}.motor();
                   MS.S{sid}.motor();
                   val = val && isequal(this.samples{i}.S{sid}, MS.S{sid});
               end
           end
           
           this.verifyTrue(val);
           fprintf(1,'time to copy %s\n',num2str(t,3));
       end
       
       % returns a REFERENCE of this class
       % for a particular type-level sample "indx", with all
       % of the variables instantiated
       function setSampleType(this)           
           Mfit = MotorProgramFit(this.samples);
           
           ns = Mfit.ns;
           nsamp = Mfit.nsamples;
           L = Mfit.TypeS;
           assert(size(L,1)==nsamp && size(L,2)==ns);
           
           val = true;
           for i=1:nsamp
               tic
               MS = Mfit.setSampleType(i);
               t=toc;
               for sid=1:ns
                   this.samples{i}.S{sid}.motor();
                   MS.S{sid}.motor();
                   val = val && isequal(this.samples{i}.S{sid}, MS.S{sid});
               end
               val = val && (Mfit.currentType == i);
               Mfit.clearSampleType;
               val = val && isnan(Mfit.currentType);
           end
           
           this.verifyTrue(val);
           fprintf(1,'time to switch reference %s\n',num2str(t,3));
       end
       
   end
   
   
end