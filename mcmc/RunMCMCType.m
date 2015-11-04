%
% Type-level resampling using MCMC
% 
% Sample a chain using ps.mcmc.nsamp_type_chain iterations
%
% Store a subset (ps.mcmc.nsamp_type_store) of these samples linearly
% spaced along the chain
%
% Return
%   samples_store: [nsamp_store x 1 cell]
%
function samples_store = RunMCMCType(M,lib)

    ps = defaultps;

    % fill-in type-level shape variables
    Q = M.copy();
    if isempty(Q.S{1}.shapes_type)
       for sid=1:M.ns
          Q.S{sid}.shapes_type = Q.S{sid}.shapes_token;
       end       
    end
    
    % fill-in type-level eval variable
    ns = Q.ns;
    for sid=1:ns
       R = Q.S{sid}.R;
       if strcmp(R.type,'mid') && isempty(Q.S{sid}.R.eval_spot_type)
           Q.S{sid}.R.eval_spot_type = Q.S{sid}.R.eval_spot_token;
       end       
    end
     
    % run mcmc
    samples_chain = mcmc_all(Q,lib,ps.mcmc.nsamp_type_chain,'type');
    
    % pick well-spaced set
    int = ps.mcmc.nsamp_type_chain ./ ps.mcmc.nsamp_type_store;
    indx = round(linspace(int,ps.mcmc.nsamp_type_chain,ps.mcmc.nsamp_type_store));
    samples_store = samples_chain(indx);
    
end