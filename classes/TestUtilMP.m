% To run unit tests:
%
% testCase = TestUtilMP;
% res = run(testCase);
%
classdef TestUtilMP < matlab.unittest.TestCase
    % Test various UtilMP operations on the motor programs class
    
    properties
       M 
    end
    
    methods (TestMethodSetup)
       
        function createModels(this)            
            lib = loadlib;
            motor_program = generate_all_rel_character(lib);
            M = motor_program();        
            while M.ink_off_page
                motor_program = generate_all_rel_character(lib);
                M = motor_program();      
            end
            for i=1:M.ns
                M.S{i}.shapes_type = [];
            end
            this.M = M;
        end
        
    end
    
    methods (Test)
       
        function testFlipStroke(this)
            % flip the direction of a stroke
            
            this.M.S = this.M.S(1);
            O = this.M.copy();
            S = this.M.S{1};            
            S_before = S.copy();
            UtilMP.flip_stroke(S);
            S_after = S.copy();
            
            % visualize   
            figure(1)
            clf
            subplot(1,2,1);
            plot_motor_to_image(O.pimg,O.motor_warped);
            title('before flip');            
            subplot(1,2,2);
            plot_motor_to_image(this.M.pimg,this.M.motor_warped);            
            title('after flip');
            
            % Test flipping the direction of a stroke         
            start_before = S_before.motor{1}(1,:);
            start_after = S_after.motor{1}(1,:);
            end_before = S_before.motor{end}(end,:);
            end_after = S_after.motor{end}(end,:);
            
            valid = aeq(start_before,end_after) && aeq(start_after,end_before);
            this.verifyTrue(valid);
        end
        
        function testMergeStrokes(this)
            % Test merging two strokes together
            
            lib = loadlib;
            M = this.M;
            moves = UtilMP.all_merge_moves(M);
            if isempty(moves)
                fprintf(1,'no valid merge moves\n');
                return
            end           
            mv = moves{1};          
            
            S1 = M.S{mv.i1};
            S2 = M.S{mv.i2};
            
            M_before = M.copy();
            M_before.S{1} = S1.copy();
            M_before.S{2} = S2.copy();
            M_before.S(3:end) = [];
            argmax_relations(lib,M_before);
            
            Z = UtilMP.merge_strokes(S1,S2);
            M_after = M.copy();
            M_after.S{1} = Z;
            M_after.S(2:end) = [];
            argmax_relations(lib,M_after);
            
            figure(2)
            clf
            subplot(1,2,1);
            plot_motor_to_image(M_before.pimg,M_before.motor_warped);
            title('before merge');
            subplot(1,2,2);
            plot_motor_to_image(M_after.pimg,M_after.motor_warped);
            title('after merge');
            
            
        end
        
        function testSplitStrokes(this)
            % Test spliting a stroke apart at the sub-strokes
            
            lib = loadlib;
            M = this.M;          
            moves = UtilMP.all_split_moves(M);
            if isempty(moves)
                fprintf(1,'no valid split moves\n');
                return
            end           
            mv = moves{1};          
            
            Z = M.S{mv.sid}.copy();
            [S1,S2] = UtilMP.split_stroke(Z,mv.bid_start_of_second);            
            
            M_before = M.copy();
            M_before.S{1} = Z;
            M_before.S(2:end) = [];
            argmax_relations(lib,M_before);
            
            M_after = M.copy();
            M_after.S{1} = S1.copy();
            M_after.S{2} = S2.copy();
            M_after.S(3:4) = [];
            argmax_relations(lib,M_after);
            
            figure(3)
            clf
            subplot(1,2,1);
            plot_motor_to_image(M_before.pimg,M_before.motor_warped);
            title('before split');
            subplot(1,2,2);
            plot_motor_to_image(M_after.pimg,M_after.motor_warped);
            title('after split');
            
        end
        
        function testSplitMerge(this)
            % split a stroke, then merge it 
            % see if we can reconstruct the original character
            
            lib = loadlib;
            M = this.M;          
            moves = UtilMP.all_split_moves(M);
            if isempty(moves)
                fprintf(1,'no valid split moves\n');
                return
            end           
            mv = moves{1};          
            
            Z1 = M.S{mv.sid}.copy();
            [S1,S2] = UtilMP.split_stroke(Z1,mv.bid_start_of_second);            
            
            M_before = M.copy();
            M_before.S{1} = Z1;
            M_before.S(2:end) = [];
            argmax_relations(lib,M_before);
            
            M_after = M.copy();
            M_after.S{1} = S1.copy();
            M_after.S{2} = S2.copy();
            M_after.S(3:4) = [];
            argmax_relations(lib,M_after);
            
            
            moves = UtilMP.all_merge_moves(M_after);
            assert(numel(moves)==1);
            mv = moves{1};                      
            SS1 = M_after.S{mv.i1};
            SS2 = M_after.S{mv.i2};
            Z2 = UtilMP.merge_strokes(SS1,SS2);
            
            M_recon = M_before.copy();
            M_recon.S{1} = Z2.copy();
            argmax_relations(lib,M_recon);
            
            figure(4)
            clf
            subplot(1,2,1);
            plot_motor_to_image(M_before.pimg,M_before.motor_warped);
            title('original');
            subplot(1,2,2);
            plot_motor_to_image(M_recon.pimg,M_recon.motor_warped);
            title('reconstructed');

            this.verifyTrue(isequal(M_before,M_recon));            
            
        end
    end
end