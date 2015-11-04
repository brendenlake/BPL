%
%  Map a sequence of pen points into a Motor Program
%
%  Input
%   S: [ns x 1 nested cell array] of trajectories in motor space [n x 2]
%   I: [binary image]
%   fclassify: [function handle] fclassify([n x 2] trajectory, [scalar]
%                                   scale)
%               this function classifies trajectories are belong to a
%               specific library component
%
%  Output
%   M: MotorProgram object
%
function M = sequence_to_MP(S,I,fclassify,lib)

    assert(iscell(S));
    assert(iscell(S{1}));

    % find the new scale of primitives
    assert(exist('lib','var')>0);
    newscale = lib.newscale;
    ncpt = lib.ncpt;
    
%     % find the new scale of primitives
%     ps_clustering = defaultps_clustering;
%     newscale = ps_clustering.newscale_ss;
%     ncpt = ps_clustering.ncpt;
%     if exist('lib','var')
%        newscale = lib.newscale;
%        ncpt = lib.ncpt;
%     end
    
    % get the normalized trajectories, the stroke centers, scales, etc.
    verbose = false;
    [S_norm,~,S_scales] = normalize_dataset(S,newscale,verbose);    
    
    % fit splines to the normalized trajectories, and assign them
    % indices in the library
    PM = defaultps;
    ps_bottomup = defaultps_bottomup;
    ns = length(S);
    S_splines = cell(ns,1);
    for sid=1:ns         
        nsub = length(S{sid});
        S_splines{sid} = cell(nsub,1);
        for b=1:nsub
            S_splines{sid}{b} = fclassify(S_norm{sid}{b},S_scales{sid}{b});
        end
    end

    % make motor program structure
    M = MotorProgram(ns);
    M.I = I; % set the image
    M.parameters = PM;
    for sid=1:ns % for each stroke
       nsub = length(S{sid}); 
       M.S{sid}.ids = []; %zeros(nsub,1);
       M.S{sid}.invscales_token = zeros(nsub,1);
       M.S{sid}.shapes_token = zeros(ncpt,2,nsub);
       
       for b=1:nsub % for each sub-stroke
          M.S{sid}.ids(b,:) = S_splines{sid}{b}.indx;
          M.S{sid}.invscales_token(b) = 1 ./ S_scales{sid}{b}(1);
          M.S{sid}.shapes_token(:,:,b) = S_splines{sid}{b}.bspline;
       end
       M.S{sid}.invscales_type = M.S{sid}.invscales_token;
       
       M.S{sid}.pos_token = S{sid}{1}(1,:); % start position in real trajectory     
    
    end    
    
    % initialize the noise parameters
    M.epsilon = ps_bottomup.init_epsilon;
    M.blur_sigma = ps_bottomup.init_blur_sigma;
    
end