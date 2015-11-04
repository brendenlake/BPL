%
% Generate a character with "ns" strokes.
% This function samples from the prior.
%
% Output
%   motor_program: function handle to produce exapmles
%   template: partial motor program at the type-level
function [motor_program,template] = generate_character(libclass,ns)
    
    if ~exist('ns','var')
        ns = CPD.sample_number(libclass);
    end
    template = MotorProgram(ns);
    template.parameters = defaultps;
    
    % for each stroke, sample its template
    for i=1:ns
        template.S{i}.R = CPD.sample_relation_type(libclass,template.S(1:i-1));
        template.S{i}.ids = CPD.sample_sequence(libclass,ns);
        template.S{i}.shapes_type = CPD.sample_shape_type(libclass,template.S{i}.ids);
        template.S{i}.invscales_type = CPD.sample_invscale_type(libclass,template.S{i}.ids);
    end
    
    % return stochastic program that generates new images
    motor_program = @() generate_exemplar(template,libclass);
    
end