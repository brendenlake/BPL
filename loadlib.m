%
% Load library from a file, 
% as specified in the function defaultps
%
function lib = loadlib()
    ps = defaultps;
    load(ps.libname,'lib');
    assert(isa(lib,'Library'));
end