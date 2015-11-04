function [ pos ] = getAttachPoint(R,previous_strokes)
%GETATTACHPOINT Get the mean attachment point of 
%   where the start of the next stroke should be, given the previous
%   ones and its relations

    switch R.type
       case 'unihist'
            pos = R.gpos;
       case 'start'
            subtraj = previous_strokes{R.attach_spot}.motor{1};
            pos = subtraj(1,:);
       case 'end'
            subtraj = previous_strokes{R.attach_spot}.motor{end};
            pos = subtraj(end,:);
       case 'mid'
            bspline = previous_strokes{R.attach_spot}.motor_spline(:,:,R.subid_spot);
            pos = bspline_eval(R.eval_spot_token,bspline);
       otherwise
            error('invalid relation');
    end

end

