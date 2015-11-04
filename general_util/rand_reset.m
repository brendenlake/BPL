%
% Reset the random number generator
%
function rand_reset
    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);
end

