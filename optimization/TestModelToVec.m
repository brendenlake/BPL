% To run unit tests:
%
% testCase = TestModelToVec;
% res = run(testCase);
%
classdef TestModelToVec < matlab.unittest.TestCase
    
    properties
       M
    end
    
    methods (TestMethodSetup)
        
        function createModels(this)            
            lib = loadlib();
            
            motor_program = generate_all_rel_character(lib);
            M = motor_program();
            M.A = [1; 1; 0; 0];
            M.I = CPD.sample_image(M.pimg);
            for i=1:M.ns
                M.S{i}.R = [];
                M.S{i}.shapes_type = [];
            end
            this.M = M;
        end
        
    end
    
    methods (Test)
        
        function testType(this)
            % tests model_to_vec_fit_type function
            % and its inverse vec_to_model_fit_type
            M0 = this.M.copy();
            theta = model_to_vec_fit_type(this.M);
            
            this.M.blur_sigma = 0;
            this.M.epsilon = 0;
            for i=1:this.M.ns
                this.M.S{i}.pos_token = this.M.S{i}.pos_token*0;
                this.M.S{i}.R = [];
                this.M.S{i}.shapes_type = this.M.S{i}.shapes_type*0;
                this.M.S{i}.shapes_token = this.M.S{i}.shapes_token*0;
                this.M.S{i}.invscales_type = this.M.S{i}.invscales_type*0;
                this.M.S{i}.invscales_token = this.M.S{i}.invscales_token*0;
            end
           
            vec_to_model_fit_type(theta,this.M);           
            this.verifyTrue(isequal(M0.pimg,this.M.pimg));         
        end
        
        function testToken(this)
            % tests model_to_vec_fit_token for MotorProgram class
            % and its inverse vec_to_model_fit_token
            lib = loadlib();
            motor_program = generate_all_rel_character(lib);
            M = motor_program();
            this.M = M;
            
            M0 = this.M.copy();
            theta = model_to_vec_fit_token(this.M);
            
            this.M.blur_sigma = 0;
            this.M.epsilon = 0;
            this.M.A = zeros(4,1);
            for i=1:this.M.ns
                this.M.S{i}.pos_token = this.M.S{i}.pos_token*0;
                this.M.S{i}.shapes_token = this.M.S{i}.shapes_token*0;
                this.M.S{i}.invscales_token = this.M.S{i}.invscales_token*0;
            end
            
            vec_to_model_fit_token(theta,this.M);           
            this.verifyTrue(isequal(M0.pimg,this.M.pimg));         
        end
        
        function testAffine(this)
            % tests model_to_vec_fit_affine for MotorProgram class 
            lib = loadlib();
            this.M.A = CPD.sample_affine(lib);
            M0 = this.M.copy();
            theta = model_to_vec_fit_affine(this.M);
            
            this.M.blur_sigma = 0;
            this.M.epsilon = 0;
            this.M.A = 0*this.M.A;
           
            vec_to_model_fit_affine(theta,this.M);           
            this.verifyTrue(isequal(M0.pimg,this.M.pimg));         
        end
        
        function testTokenMfit(this)
            % tests model_to_vec_fit_token for MotorProgramFit class 
            lib = loadlib();
            motor_program = generate_all_rel_character(lib);
            M = motor_program();
            while M.ink_off_page
                motor_program = generate_all_rel_character(lib);
                M = motor_program();      
            end
            nsamp_mcmc = 5;
            samples = mcmc_all(M,lib,nsamp_mcmc,'type');
            clear M
            
            Mfit = MotorProgramFit(samples);
            
            M0 = Mfit.copy();
            theta = model_to_vec_fit_token(Mfit);
            
            Mfit.blur_sigma = 0;
            Mfit.epsilon = 0;
            Mfit.A = zeros(4,1);
            for i=1:Mfit.ns
                Mfit.S{i}.pos_token = Mfit.S{i}.pos_token*0;
                Mfit.S{i}.shapes_token = Mfit.S{i}.shapes_token*0;
                Mfit.S{i}.invscales_token = Mfit.S{i}.invscales_token*0;
            end
            
            vec_to_model_fit_token(theta,Mfit);           
            this.verifyTrue(isequal(M0.pimg,Mfit.pimg));         
        end
        
    end
    
end

