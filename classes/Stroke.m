classdef Stroke < BetterHandle
    %
    % STROKETOKEN Random variables that define a continuous pen trajectory
    %   Utilized in the MotorProgram class
    %
    % Reference to StrokeType object "myType"
    %  This might be shared between multiple strokes (see MotorProgram).
    properties
       myType
       lh
    end
    
    % type-level parameters (which set/get accesses myType)
    properties (Dependent = true)
       R
       ids
       invscales_type                
       shapes_type
    end
    
    % token-level parameters
    properties (SetObservable = true)
        pos_token        
        invscales_token
        shapes_token
    end
    
    % variables that are fully-determined
    % given the others
    properties (Dependent = true)
        nsub
        motor
        motor_spline
    end
    
    properties (GetAccess = private, SetAccess = private)
       cache_current = false;
       cache_motor
       cache_motor_spline
       eval_spot_token  % token-level parameter for where
                        % this stroke attaches ona previous spline
       po = {'pos_token','invscales_token','shapes_token'}; % tracked properties
                        % we don't need to track "eval_spot_token" because
                        % it does not influence the image directly
    end   
    
    methods (Access = protected)
        
        % override the copy method to copy over the listeners
        function cpObj = copyElement(this)
            cpObj = copyElement@BetterHandle(this);
            
            % deep copy of the type-level variable
            cpObj.myType = copy(this.myType);         
            
            assert(~cpObj.onListener);
            cpObj.refresh_listener();
        end
        
    end    
    
    methods
        
        % constructor
        function this = Stroke(previousStroke)           
           this.refresh_listener();
           if exist('previousStroke','var')
              this.myType = previousStroke.myType;
           else
              this.myType = StrokeType(); 
           end           
        end
        
        function load_legacy(this,oldS)
        % Input
        %  oldMP: older version of Stroke class,where
        %  we may have updated it since
        %
        % This fills in all the variables of the new class
        % with the old ones
            this.myType = oldS.myType.copy();
            this.pos_token = oldS.pos_token;       
            this.invscales_token = oldS.invscales_token;
            this.shapes_token = oldS.shapes_token;           
        end
        
        % get methods for type-level variables
        function out = get.ids(this), out = this.myType.ids; end
        function out = get.invscales_type(this), out = this.myType.invscales_type; end
        function out = get.shapes_type(this), out = this.myType.shapes_type; end
        function out = get.R(this)
            % standard get, but for mid-stroke attachments, 
            % the token-level spline-attachment value is stored
            % in this object
            out = this.myType.R;
            if ~isempty(out) && strcmp(this.myType.R.type,'mid')
                out.eval_spot_token = this.eval_spot_token;
            end            
        end
        
        % set methods for type-level variables
        function set.ids(this,val), this.myType.ids = val; end
        function set.invscales_type(this,val), this.myType.invscales_type = val; end
        function set.shapes_type(this,val), this.myType.shapes_type = val; end
        function set.R(this,val)            
            this.myType.R = val;            
            if ~isempty(this.myType.R) && strcmp(this.myType.R.type,'mid')
                this.myType.R.eval_spot_token = [];
                this.eval_spot_token = val.eval_spot_token;
            end
        end        
        
        function out = get.nsub(this)
           % get the number of sub-strokes
           out = length(this.ids);
        end
        
        function motor = get.motor(this)
        % compute the [x,y,t] trajectory of this stroke, either from cached
        % item or from scratch
            assert(this.onListener);
            if this.cache_current
                 motor = this.cache_motor;
            else            
                [motor,motor_spline] = vanilla_to_motor(this.shapes_token,this.invscales_token,this.pos_token);
                this.cache_motor = motor;
                this.cache_motor_spline = motor_spline;
                this.cache_current = true;
            end            
        end
        
        function motor_spline = get.motor_spline(this)
        % compute the spline trajectory of this stroke, either from cached
        % item or from scratch
            assert(this.onListener);
            if this.cache_current
                motor_spline = this.cache_motor_spline;
            else            
                [motor,motor_spline] = vanilla_to_motor(this.shapes_token,this.invscales_token,this.pos_token);
                this.cache_motor = motor;
                this.cache_motor_spline = motor_spline;
                this.cache_current = true;
            end            
        end
        
        function out = get_eval_spot_token(this)
            % return value of private variable
            out = this.eval_spot_token;
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
        end
        
        
    end
    
    methods (Access = private)
       
        function refresh_listener(this)
        % update the listener if it is invalid    
           if ~this.onListener
               this.lh = addlistener(this,this.po,'PostSet',@Stroke.needs_update);
           end
        end
        
    end
    
    methods (Static)
        
        function Y = loadobj(X)
        % replaces the listener    
            Y = X.copy();
        end
        
        function needs_update(~,event)
            event.AffectedObject.cache_current = false;
        end
    end    
    
end

%
% Create the fine-motor trajectory of a stroke
% (denoted 'f()' in pseudocode)
%
%  Input (for k sub-strokes)
%   vanilla_shapes: [ncpt x 2 x k] spline points in normalized space
%   invscales: [k x 1] inverse scales for each sub-stroke
%   first_pos: [1 x 2] starting location of stroke
%
% Output
%   motor [k x 1 cell] fine motor sequence, in position and scale
%   motor_spline [k x 1 cell] same motor sequence, but in spline space
function [motor,motor_spline] = vanilla_to_motor(vanilla_shapes,invscales,first_pos)

    % make sure we have all of the variables we need
    if (isempty(vanilla_shapes) || isempty(invscales) || isempty(first_pos))
        motor = [];
        motor_spline = [];
        return
    end
    
    [ncpt,~,n] = size(vanilla_shapes);
    
    % re-scale the control points
    for i=1:n
       vanilla_shapes(:,:,i) = invscales(i) .* vanilla_shapes(:,:,i);
    end    

    % get trajectories from b-spline
    vanilla_traj = cell(n,1);
    for i=1:n
        vanilla_traj{i} = get_stk_from_bspline(vanilla_shapes(:,:,i));
    end
    
    % reposition
    motor = vanilla_traj;
    motor_spline = vanilla_shapes;
    for i=1:n
       if i==1
           offset = motor{i}(1,:) - first_pos;
       else
           offset = motor{i}(1,:) - motor{i-1}(end,:);
       end
       motor{i} = offset_stk(motor{i},offset);
       motor_spline(:,:,i) = motor_spline(:,:,i) - repmat(offset,[ncpt 1]);
    end

end