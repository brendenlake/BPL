classdef MotorProgram < BetterHandle
    %
    % MOTORPROGRAM Class defines a motor program object,
    %   which is the "concept" behind a character
    % 
    % MotorPrograms can either be independent objects,
    % or dependent ones that share type-level variables through
    % their array of Stroke classes.
    %
    % Constructor
    %  M = MotorProgram()  // creates independent object
    %  M = MotorProgram(Q) // creates new object M tied to Q at the
    %   type-level
    % 
    % Example, the function "generate_exemplar", creates multiple 
    %  objects tied at the type-level
    
    properties
       I % image (binary, true is inked)
       S % cell array of strokes       
       parameters % fixed parameters (not random variables)
       lh % listener handle
    end
    
    properties (SetObservable = true)
       epsilon % noise on pixels
       blur_sigma % image blur
       A = []; % affine transformation [x-scale,y-scale,x-translate,y-translate]
               %   the translation is relative to the center of mass
    end
    
    properties (Dependent = true)
       ns % number of strokes
       pimg % probability of inking each pixel
       ink_off_page % is any of the ink drawn off the page?
       motor % nested cell array of pen trajectories
       motor_warped % same as "motor", except after the affine warp
       cache_grand_current % is the cached "pimg" and "ink_off_page" still current?
    end
    
    properties (GetAccess = private, SetAccess = private)
       cache_noise_current = false; % the type-level parameters have not been update
       cache_pimg_motor % cache the state of the motor trajectory when pimg was updated
       cache_pimg
       cache_ink_off_page
       po = {'epsilon','blur_sigma','A'}; % tracked properties
    end
    
    methods (Access = protected)
        % make selective deep copy of the object
        function cpObj = copyElement(this)
            
            % make a shallow copy
            cpObj = copyElement@BetterHandle(this);
            
            % make a selective deep copy of the strokes
            for i=1:this.ns
                cpObj.S{i} = copy(this.S{i});
            end
            
            % add listener to the copies properties
            assert(~cpObj.onListener);
            cpObj.refresh_listener();
        end         
    end    
    
    methods        
     
        function this = MotorProgram(arg)
        % Constructor.    
        % Input
        %  arg: either the number of strokes (scalar),
        %     or a previous "MotorProgram" object that
        %     we want to link to this one at the type level
        % 
            this.refresh_listener();
            if isnumeric(arg)
               ns = arg;
               this.S = cell(ns,1);
               for i=1:ns
                   this.S{i} = Stroke();
               end                
            elseif isa(arg,'MotorProgram')
               Template = arg;
               this.S = cell(Template.ns,1);
               for i=1:Template.ns
                  this.S{i} = Stroke(Template.S{i});
               end
               this.parameters = Template.parameters;
            else
               error('invalid constructor'); 
            end
        end
        
        function load_legacy(this,oldMP)
        % Input
        %  oldMP: older version of MotorProgram class,where
        %  we may have updated it since
        %
        % This fills in all the variables of the new class
        % with the old ones
           this.I = oldMP.I;
           this.parameters = oldMP.parameters;
           this.epsilon = oldMP.epsilon;
           this.blur_sigma = oldMP.blur_sigma;
           this.A = oldMP.A;
           this.S = cell(size(oldMP.S));
           for sid=1:this.ns
              this.S{sid} = Stroke();
              this.S{sid}.load_legacy(oldMP.S{sid}); 
           end
        end        
        
        % get number of strokes
        function out = get.ns(this)
           out = length(this.S);            
        end
        
        % get pen trajectory (un-warped)
        function motor = get.motor(this) 
            ns = this.ns;
            motor = cell(ns,1);
            for i=1:ns
                motor{i} = this.S{i}.motor;
            end         
        end
        
        % get the pen trajectory (warped)
        function motor_warped = get.motor_warped(this)
            motor_warped = apply_warp(this);            
        end
        
        % get probability map of image
        function pimg = get.pimg(this)
            assert(this.onListener);
            if this.cache_grand_current            
                pimg = this.cache_pimg;
            else
                [pimg,ink_off_page] = apply_render(this);
                this.cache_pimg = pimg;
                this.cache_ink_off_page = ink_off_page;
                this.cache_noise_current = true;
                this.cache_pimg_motor = this.motor;
            end
        end
        
        % is there ink off of the page?
        function ink_off_page = get.ink_off_page(this)
            assert(this.onListener);
            if this.cache_grand_current            
                ink_off_page = this.cache_ink_off_page;
            else
                [pimg,ink_off_page] = apply_render(this);
                this.cache_pimg = pimg;
                this.cache_ink_off_page = ink_off_page;
                this.cache_noise_current = true;
                this.cache_pimg_motor = this.motor;
            end
        end
        
        % see if the cache is up-to-date
        function out = get.cache_grand_current(this)
            
            % make sure image-level parameters have not changed
            out = this.cache_noise_current;
            
            % make sure this.motor aligns with last time we
            % cached the pimg
            out = out && isequal(this.cache_pimg_motor,this.motor);
            
        end                  
        
        % return true if all the relations (in list_sid) are set and non-empty
        function out = has_relations(this,list_sid)
            if ~exist('list_sid','var')
               list_sid = 1:this.ns; 
            end
            present = false(size(list_sid));
            count = 1;
            for sid = list_sid
               present(count) = ~isempty(this.S{sid}.R);
               count = count + 1;
            end            
            if ~(all(present) || all(~present))
                error('all relations should be present or not');
            end            
            out = all(present);
        end
        
        % set all relations to the empty relation
        function clear_relations(this)
            for sid=1:this.ns
                this.S{sid}.R = [];
            end
        end
        
        % remove memory-heavy image matrices from the class
        function lightweight(this)
            this.cache_pimg = [];
            this.I = [];
            this.cache_noise_current = false;            
        end
        
        function clear_shapes_type(this)
        % remove the type-level shapes
            for sid=1:this.ns
                this.S{sid}.shapes_type = [];
            end
        end
        
        % Is the listener on and functional?
        function on = onListener(this)
            on = ~isempty(this.lh) && isvalid(this.lh) && this.lh.Enabled ...
                && eq(this,this.lh.Object{1});
        end   
        
        function Y = saveobj(this)
        % custom save method to delete listener    
           Y = this.copy();
           delete(Y.lh);
           for sid=1:Y.ns
              Y.S{sid} = Y.S{sid}.saveobj(); 
           end
        end              
        
    end
    
    methods (Access = private)
        
        
        function motor_warped = apply_warp(this)
        % apply affine transformation to the stroke trajectories,
        % and return as a nested cell array
            motor_unwarped = this.motor;
            if isempty(this.A)
                motor_warped = motor_unwarped;
                return
            end            
            cell_traj = UtilMP.flatten_substrokes(motor_unwarped);
            com = com_char(cell_traj);
            B = zeros(4,1);
            B(1:2) = this.A(1:2);
            B(3:4) = vec(this.A(3:4)) - (vec(this.A(1:2))-1) .* com';
            motor_warped = UtilMP.apply_each_substroke(motor_unwarped,@(stk)affine_warp(stk,B));            
        end
        
        function [pimg,ink_off_page] = apply_render(this)
        % apply affine warp and render the image    
            motor_warped = apply_warp(this);
            flat_warped = UtilMP.flatten_substrokes(motor_warped);
            [pimg,ink_off_page] = render_image(flat_warped,this.epsilon,this.blur_sigma,this.parameters);            
        end                
        
        function refresh_listener(this)
        % update the listener if it is invalid
            if ~this.onListener
                this.lh = addlistener(this,this.po,'PostSet',@MotorProgram.needs_update); 
            end  
        end
        
    end    
   
    methods (Static)      
        
        function Y = loadobj(X)
        % replaces the listener      
           Y = X.copy();
        end
        
        function needs_update(~,event)
        % notify that we have updated the noise parameters
            event.AffectedObject.cache_noise_current = false;
        end   
        
        function out = istied(varargin)
        % check whether list of MotorProgram objects reference then same
        % type object
        %
        % Input: varargin: cell array of MotorProgram objects
        % Output: out=true if they all co-reference the same type object       
            if numel(varargin)==1 && iscell(varargin{1})
                varargin = varargin{1};
            end
            n = length(varargin);
            base = varargin{1};
            ns = base.ns;
            out = true;
            for i=2:n
                item = varargin{i};
                if (item.ns == ns)
                   for j=1:ns
                       
                      % make sure the strokes have the same type
                      if (item.S{j}.myType ~= base.S{j}.myType)
                         out = false;
                         return
                      end
                      
                   end
                else
                   
                   % they don't have the same number of strokes
                   out = false;
                   return
                end
            end            
        end
        
    end
    
end
