% To run unit tests:
%
% testCase = TestMotorProgram;
% res = run(testCase);
%
classdef TestMotorProgram < matlab.unittest.TestCase
    
   properties
      samples 
   end
    
   methods (TestMethodSetup)
       
        function createModels(this)            
            nsamp = 9;
            ps = defaultps;            
            lib = loadlib;
            motor_program = generate_all_rel_character(lib);
            this.samples = cell(nsamp,1);
            for i=1:nsamp
                this.samples{i} = motor_program();
            end
        end
   end
    
   methods (Test)
       
       function saveLoad(this)
       % saving/loading MotorProgram objects can be complex,
       % since there are listener properties
           M = this.samples{1}.copy();
           assert(M.onListener);
           save('temp','M');
           Q = M.copy();
           clear M
           load('temp','M');
           delete('temp.mat');
           this.verifyTrue(isequal(M,Q)); 
       end       
       
       function checkTied(this)
       % when we generate samples, make sure they are tied at the
       % type-level
           this.verifyTrue(MotorProgram.istied(this.samples)); 
       end
       
       % make sure relation functionality is working, 
       % in terms of a model indicating it has relations set or not.
       function checkRelations(this)
           out = true;
           M = this.samples{1};
           out = out && M.has_relations();
           M.S{1}.R = [];
           out = out && M.has_relations(2);
           M.clear_relations();
           out = out && ~M.has_relations();
           this.verifyTrue(out);           
       end
       
       function checkUntied(this)
           DeepCopy = this.samples{1}.copy();
           vsamples = [{DeepCopy}; this.samples(:)];
           this.verifyFalse(MotorProgram.istied(vsamples));           
       end
       
       % see if caching produces a speed-up when we modify a stroke
       % position
       function testCachingPosToken(this)
           
           MP = this.samples{randint(1,1,[1 length(this.samples)])};
           
           tic
           pimg = MP.pimg;
           t1=toc;
           fprintf(1,'\ntime for cached event: %s (s)\n',num2str(t1,3));

           MP.S{1}.pos_token = MP.S{1}.pos_token + [1 1];

           tic
           pimg = MP.pimg;
           t2=toc;
           fprintf(1,'time for recomputation: %s (s)\n',num2str(t2,3));

           tic
           pimg = MP.pimg;
           t3=toc;
           fprintf(1,'time for cached event: %s (s)\n\n',num2str(t3,3));
           
           this.verifyTrue(t2 > t1 && t2 > t3);
       end
       
       % 
       function testCachingMotorDisrupt(this)
           
           MP = this.samples{randint(1,1,[1 length(this.samples)])};
           
           MP.S{1}.pos_token = MP.S{1}.pos_token + [1 1];

           M1 = MP.copy();
           M2 = MP.copy();
           
           % update motor
           M1.motor;
           pimg1 = M1.pimg;
           
           % don't update motor
           pimg2 = M2.pimg;
           
           this.verifyEqual(pimg1,pimg2,'AbsTol',1e-4);
       end
       
       % see if caching produces a speed-up if we modify the stroke blur
       function testCachingBlur(this)
           
           MP = this.samples{randint(1,1,[1 length(this.samples)])};
           
           tic
           pimg = MP.pimg;
           t1=toc;
           fprintf(1,'\ntime for cached event: %s (s)\n',num2str(t1,3));

           MP.blur_sigma = MP.blur_sigma + 1;

           tic
           pimg = MP.pimg;
           t2=toc;
           fprintf(1,'time for recomputation: %s (s)\n',num2str(t2,3));

           tic
           pimg = MP.pimg;
           t3=toc;
           fprintf(1,'time for cached event: %s (s)\n\n',num2str(t3,3));
           
           this.verifyTrue(t2 > t1 && t2 > t3);
       end
       
       % Make sure that multiple examples of a given character type
       % are tied at the type-level, so changing it for one
       % changes it for all
       function testTiedSubid(this)           
           this.samples{1}.S{1}.ids = this.samples{1}.S{1}.ids + 1;           
           this.verifyTrue(all(this.samples{1}.S{1}.ids == this.samples{2}.S{1}.ids));
       end
       
       % make sure that for a "mid" relation, the eval_spot_type
       % is tied but the eval_spot_token is not
       function testNotTiedSubid(this)
           
           M = MotorProgram(2);
           M.S{2} = copy(this.samples{1}.S{1});
           M.S{2}.R = RelationAttachAlong('mid',1,1,M.S{2}.nsub,1,5);
           M.S{2}.R.eval_spot_type = 2;
           M.S{2}.R.eval_spot_token = 2;
           
           Q = MotorProgram(M);
           Q.S{2}.R.eval_spot_type = 3;
           Q.S{2}.R.eval_spot_token = 3;
           
           this.verifyTrue(M.S{2}.R.eval_spot_type==3 && M.S{2}.R.eval_spot_token == 2);        
       end
       
       % make sure copies are not tied at the type-level
       function testUntiedSubid(this)
           DeepCopy = this.samples{1}.copy();
           this.samples{1}.S{1}.ids = this.samples{1}.S{1}.ids + 1;           
           this.verifyTrue(all(this.samples{1}.S{1}.ids ~= DeepCopy.S{1}.ids));
       end
       
       
   end
   
   
end