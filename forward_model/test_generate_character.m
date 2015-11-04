%
% Generate new character types from the prior. 
%
% Note, this function does not generate characters as in
% the "Generating new concepts (unconstrained) task" which
% modeled additional global statistics.
%

ps = defaultps;
nsamp = 81;
lib = loadlib;

samples = cell(nsamp,1);
for i=1:nsamp
    [~,type] = generate_character(lib);    
    type = construct_full_type(type);
    my_center_mp(type);
    samples{i} = type;
end

figure(1)
clf;
nrow = ceil(sqrt(nsamp));
for i=1:nsamp
    sptight(nrow,nrow,i);
    simpleViz(samples{i});
end