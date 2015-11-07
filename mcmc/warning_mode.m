function warning_mode(prompt)
    prompt = [prompt,' (PRESS "y" TO CONTINUE): '];
    u = '';
    while ~strcmpi(u,'y')
        u = input(prompt,'s');
    end
end