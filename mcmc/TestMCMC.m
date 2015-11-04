%
% Test to make sure the acceptance ratios are correct for the MH moves.
% This compares the ratios from the relevant conditional distributions
% to the overall generative score
%
% If this does not throw an assert, then it should be working.
%
ntry = 1;
verbose = true;
debug = true;
ps = defaultps;
load(ps.libname,'lib');
for t=1:ntry
    fprintf(1,'TRY %d\n',t);    
    motor_program = generate_all_rel_character(lib);
    M = motor_program();    
    mcmc_all(M,lib,1,'verbose','debug','type','token','relation');
end