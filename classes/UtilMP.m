classdef UtilMP
    % UtilMP Motor program utilities
    
    methods (Static)
        
        function clear_shape_type(M)
           % clear the type-level of the shape variable,
           % so that we marginalize over it
           for sid=1:M.ns
              M.S{sid}.shapes_type = [];               
           end            
        end
        
        function clear_invscales_type(M)
           % clear the type-level of the sub-stroke scale variables,
           % and just copyin the token level version
           for sid=1:M.ns
              M.S{sid}.invscales_type = M.S{sid}.invscales_token;               
           end            
        end
        
        
        function set_affine_to_image(M,I)
           % Revise the affine transformation
           % so that it fits the image as closely as possible
           target_com = UtilImage.img_com(I);
           target_range = UtilImage.img_range(I);
           missing = isinf(I);
           UtilMP.set_affine(M,target_com,target_range,missing); 
        end
        
        function set_affine(M,target_com,target_range,missing)
            % Revise the affine transformation
            % to a new center-of-mass and range
            %
            %  Base current com on image rather than motor
            %
            % Input
            %  target_com: [1 x 2] desired center-of-mass, motor coord
            %  target_range: [1 x 2] desired range in each dimension,
            %       motor coord
            %  missing: (optional) which pixels shoud we exclude?
                      
            bI = M.pimg > 0.5;
            if sum(sum(bI)) == 0 % special case when none of the pixels are inked
                rg = max(M.pimg(:)) - min(M.pimg(:));
                threshold = min(M.pimg(:)) + rg/2;
                bI = M.pimg > threshold;
            end
            if ~exist('missing','var')
                missing = false(size(bI));
            end
            curr_com = UtilImage.img_com(bI,missing);
            curr_range =  UtilImage.img_range(bI,missing);
                
            % compute the transformation
            scale = target_range ./ curr_range;
            shift = (target_com - curr_com);            
            if ~isempty(M.A)
                M.A(1:2) = vec(M.A(1:2)) .* scale(:);
                M.A(3:4) = vec(M.A(3:4)) + shift(:);
            else
                M.A = [scale(:); shift(:)];
            end
        end
        
        function set_affine_motor(M,target_com,target_range)
            % Revise the affine transformation
            % to a new center-of-mass and range.
            %
            % Current com/range is based on motor coordinates
            %
            % Input
            %  target_com: [1 x 2] desired center-of-mass, motor coord
            %  target_range: [1 x 2] desired range in each dimension,
            %       motor coord
            %  missing: (optional) which pixels shoud we exclude?            
            
            assert(numel(target_range)==2);
            
            % Base current com on the motor rather than image
            mymotor = M.motor_warped;
            mymotor = flatten_nested(mymotor);
            curr_com = com_char(mymotor);
            curr_range = range_char(mymotor);
                
            % compute the transformation
            scale = target_range ./ curr_range;
            scale(isnan(scale)) = 1; % if we divided by zero
            shift = (target_com - curr_com);            
            if ~isempty(M.A)
                M.A(1:2) = vec(M.A(1:2)) .* scale(:);
                M.A(3:4) = vec(M.A(3:4)) + shift(:);
            else
                M.A = [scale(:); shift(:)];
            end
        end            
            
        function nested = apply_each_substroke(nested,fnc)
        % apply a function handle "fnc" to each sub-stroke in the nested
        % cell array "nested"
            assert(isvector(nested)&&iscell(nested));
            ns = length(nested);
            for sid=1:ns
               nsub = length(nested{sid});
               for bid=1:nsub
                  nested{sid}{bid} = fnc(nested{sid}{bid}); 
               end
            end            
        end

        function  vcell = flatten_substrokes(nested)
        % flatten a nested cell array, to create
        % a cell array of [n x 1] elements
            assert(isvector(nested)&&iscell(nested));    

            % pre-allocate cell array
            ncell = 0;
            ns = length(nested);
            for sid=1:ns
               ncell = ncell + length(nested{sid});
            end
            vcell = cell(ncell,1);

            % fill in its elements
            count = 1;
            for sid=1:ns
               nsub = length(nested{sid});
               for bid=1:nsub
                  vcell{count} = nested{sid}{bid}; 
                  count = count + 1;
               end
            end
        end
        
        function flip_stroke(S)
            % given a Stroke object, flip it such that it
            % starts where it previously stopped            
            assert(isa(S,'Stroke'));            
            nsub = S.nsub;
            ncpt = size(S.shapes_token,1); % number of control points
            rev = nsub:-1:1;
            S.pos_token = S.motor{end}(end,:); % choose new position
            S.ids = S.ids(rev);
            if ~isempty(S.shapes_type)
                S.shapes_type = S.shapes_type(ncpt:-1:1,:,rev);
            end
            S.shapes_token = S.shapes_token(ncpt:-1:1,:,rev);
            S.invscales_type = S.invscales_type(rev);
            S.invscales_token = S.invscales_token(rev);        
        end
        
        function Z = merge_strokes(S1,S2)
            % create a new stroke "Z" by attaching S2 to the end of S1
            assert(isa(S1,'Stroke')); 
            assert(isa(S2,'Stroke')); 
            Z = Stroke();
            Z.pos_token = S1.pos_token;
            Z.ids = [S1.ids(:); S2.ids(:)];
            Z.invscales_type = [S1.invscales_type(:); S2.invscales_type(:)];
            Z.invscales_token = [S1.invscales_token(:); S2.invscales_token(:)];
            if ~isempty(S1.shapes_type)
                Z.shapes_type = cat(3,S1.shapes_type,S2.shapes_type);
            end
            Z.shapes_token = cat(3,S1.shapes_token,S2.shapes_token);            
        end        
        
        function [S1,S2] = split_stroke(Z,bid_start_of_second)
           % create two new strokes S1 and S2 by splitting stroke Z,
           % where S2 begins at the "bid_..." sub-stroke id 
           assert(isa(Z,'Stroke'));
           S1 = Stroke();
           S2 = Stroke();
           in1 = 1:bid_start_of_second-1;
           in2 = bid_start_of_second:Z.nsub;
           
           % fill in stroke 1
           S1.pos_token = Z.pos_token;
           S1.ids = Z.ids(in1);
           S1.invscales_type = Z.invscales_type(in1);
           S1.invscales_token = Z.invscales_token(in1);
           if ~isempty(Z.shapes_type);
               S1.shapes_type = Z.shapes_type(:,:,in1);
           end
           S1.shapes_token = Z.shapes_token(:,:,in1);
           
           % fill in stroke 2
           S2.pos_token = S1.motor{end}(end,:);
           S2.ids = Z.ids(in2);
           S2.invscales_type = Z.invscales_type(in2);
           S2.invscales_token = Z.invscales_token(in2);
           if ~isempty(Z.shapes_type);
               S2.shapes_type = Z.shapes_type(:,:,in2);
           end
           S2.shapes_token = Z.shapes_token(:,:,in2);
        end        
        
        function [moves,reverse_moves] = all_merge_moves(M)
            % M is a MotorProgram with instantiated relations.
            % Find all valid "merge" moves where a new stroke
            % begins at the end of a previous one.
            %
            % Output
            %  moves: [cell k x 1] list of possible merge structs
            %  reverse_moves: [cell k x 1] split moves that reverse them
            assert(isa(M,'MotorProgram'));  
            moves = [];
            reverse_moves = [];
            ns=M.ns;
            for sid=1:ns-1
                if valid_merge(M,sid)
                    
                   m = struct; % move
                   m.type = 'merge';
                   m.i1 = sid;
                   m.i2 = sid+1;                   
                   moves = [moves; {m}];
                   
                   r = struct; % reverse move
                   r.type = 'split';
                   r.sid = sid;
                   r.bid_start_of_second = M.S{sid}.nsub+1;
                   reverse_moves = [reverse_moves; {r}];
                   
                end                
            end
            moves = moves(:);
            reverse_moves = reverse_moves(:);
        end
        
        function [moves,reverse_moves] = all_split_moves(M)
            % M is a MotorProgram with instantiated relations
            % Find all valid "split" moves, where we have
            %   a sub-stroke break
            %
            % Output
            %  moves: [cell k x 1] list of possible split structs
            %  reverse_moves: [cell k x 1] merge moves that reverse them
            assert(isa(M,'MotorProgram'));
            moves = [];
            reverse_moves = [];
            ns = M.ns;
            for sid=1:ns
                nsub = M.S{sid}.nsub;
                for bid=2:nsub
                    
                    m = struct; % split move
                    m.type = 'split';
                    m.sid = sid;
                    m.bid_start_of_second = bid;
                    moves = [moves; {m}];
                    
                    r = struct; % reverse move
                    r.type = 'merge';
                    r.i1 = sid;
                    r.i2 = sid+1;
                    reverse_moves = [reverse_moves; {r}];
                end
            end
            moves = moves(:);
            reverse_moves = reverse_moves(:);
        end        
    end    
end

function val = valid_merge(M,sid)
    % Can we merge stroke "sid" and stroke "sid+1" in M?
    % return true/false
    S2 = M.S{sid+1};
    assert(~isempty(S2.R));
    val = strcmp(S2.R.type,'end');
    val = val && S2.R.attach_spot == sid; 
end