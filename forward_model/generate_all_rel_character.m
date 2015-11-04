%
% Generate a character with 4 strokes, making sure
% to use all the different kinds of relations.
%
% Output
%   motor_program: function handle to produce exapmles
%   template: partial motor program at the type-level
function [motor_program,template] = generate_all_rel_character(libclass)
    
    libclass = libclass.copy();
    ns = 4;
    template = MotorProgram(ns);
    template.parameters = defaultps;
    
    % for each stroke, sample its template
    list = 2:ns;
    perm = randperm(numel(list));
    list = list(perm);
    list = [1; list(:)];
    
    % for each stroke, sample its template
    for i=1:ns
        
        if i>1
           rel = libclass.rel;
           rel.mixprob = zeros(4,1);
           rel.mixprob(list(i)) = 1;
           libclass.setfield('rel',rel);
        end
        
        template.S{i}.R = CPD.sample_relation_type(libclass,template.S(1:i-1));
        template.S{i}.ids = CPD.sample_sequence(libclass,ns);
        template.S{i}.shapes_type = CPD.sample_shape_type(libclass,template.S{i}.ids);
        template.S{i}.invscales_type = CPD.sample_invscale_type(libclass,template.S{i}.ids);
    end
    
    % return stochastic program that generates new images
    motor_program = @() generate_exemplar(template,libclass);
    
end