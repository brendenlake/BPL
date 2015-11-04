classdef MotorProgramFit < MotorProgram
    %MOTORPROGRAMFIT Extension of the motor-program class
    % that handles multiple samples of "type-level" variables
    % and allows easy switching between them
    %
    properties
        TypeS
        currentType=nan;
    end
    
    properties (Dependent)
       nsamples 
    end
    
    properties (SetAccess = private)
       list_typeR % cell [ns x 1] type of each relation
    end
        
    methods
        
        function this = MotorProgramFit(TypeSamples)
            
           nsamp = length(TypeSamples);
           for i=1:nsamp
              assert(isa(TypeSamples{i},'MotorProgram'));
           end
           
           % superclass construtor
           anchor = TypeSamples{1}.copy();
           this = this@MotorProgram(anchor.ns);
           
           % copy over the token-level variables we want
           this.A = anchor.A;
           this.epsilon = anchor.epsilon;
           this.blur_sigma = anchor.blur_sigma;
           this.parameters = anchor.parameters;
           for sid=1:this.ns
              stk = anchor.S{sid}; 
              this.S{sid}.pos_token = stk.pos_token;
              this.S{sid}.shapes_token = stk.shapes_token;
              this.S{sid}.invscales_token = stk.invscales_token;
              if strcmp(stk.R.type,'mid')
                 this.S{sid}.R = stk.R;
              end
           end
            
           % copy over all of the type-level stroke parameters
           this.TypeS = cell(nsamp,this.ns);
           for i=1:nsamp
              for sid=1:this.ns
                  this.TypeS{i,sid} = TypeSamples{i}.S{sid}.myType.copy();
                  % already does not include myType.R.eval_spot_token
              end
           end   
           assert(checkSampleSameRelations(this));
           for sid=1:this.ns
              this.list_typeR{sid} = this.TypeS{1,sid}.R.type; 
           end
           
           % clear all of the extra type-level variables we might
           % have created
           clearSampleType(this);
            
        end
        
        function nsamp = get.nsamples(this)
            % returns the number of samples
            nsamp = size(this.TypeS,1);
        end
        
        function this = setSampleType(this,indx) 
            % returns a REFERENCE of this class
            % for a particular type-level sample "indx", with all
            % of the variables instantiated
            for sid=1:this.ns
                this.S{sid}.myType = this.TypeS{indx,sid};
            end
            this.currentType = indx;
        end
        
        function same = checkSampleSameRelations(this)
            % check to make sure all of the samples have the same
            % relation "type"
            typeR = cell(this.nsamples,this.ns);
            for indx=1:this.nsamples
                for sid=1:this.ns
                    typeR{indx,sid} = this.TypeS{indx,sid}.R.type;
                end
            end            
            same = true;
            for sid=1:this.ns
                same = same && all(strcmp(typeR{1,sid},typeR(:,sid)));
            end            
        end
        
        
        function clearSampleType(this)
            % clears the type-level variables in this model
            for sid=1:this.ns
                this.S{sid}.myType = StrokeType();
                R = struct;
                R.type = this.list_typeR{sid};
                if strcmp(R.type,'mid')
                    R.eval_spot_token = this.S{sid}.get_eval_spot_token();
                end
                this.S{sid}.R = R;
            end
            this.currentType = nan;
        end
        
        
        function out = isEmptyType(this)
            % if we haven't selected a specific type yet to use
            out = isnan(this.currentType);
        end
        
        function M = getSampleTypeCopy(this,indx) 
            % returns an COPY of this class
            % for a particular type-level sample "indx", with all
            % of the variables instantiated
            M = this.copy();
            for sid=1:M.ns
               M.S{sid}.myType = this.TypeS{indx,sid};
            end            
        end
                        
    end
    
end

