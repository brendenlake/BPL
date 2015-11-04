classdef UtilImage
   % Image utilities
    
   methods (Static)
       
    function out = check_black_is_true(I)
        % check whether the black pixels are labeled as 
        % "true" in the image format, since there should
        % be fewer black pixels
        missing = isinf(I);
        given = ~missing;
        out = sum(I(given)==true) < sum(I(given)==false);
    end
    
    function com = img_com(I,missing)
        %
        % Compute the center of mass from an image,
        % based on all of the "on" pixels.
        %
        % Return in MOTOR coordinates
        % 
        % Input
        %  missing: logical array the same size as I, but
        %    where "true" elements are "missing"
        %
        assert(UtilImage.check_black_is_true(I));
        if ~exist('missing','var')
            missing = isinf(I);
        end
        on = I==1;    
        [i,j] = find(on & ~missing);
        
        com = zeros(1,2);
        com(1) = mean(i(:));
        com(2) = mean(j(:));
        com = space_img_to_motor(com);
    end
    
    function rg = img_range(I,missing)
        %
        % Compute the range from an image,
        % based on all of the "on" pixels.
        %
        % Return in MOTOR coordinates
        %
        % Input
        %  missing: logical array the same size as I, but
        %    where "true" elements are "missing"
        %
        assert(UtilImage.check_black_is_true(I));    
        if ~exist('missing','var')
            missing = isinf(I);
        end
        on = I==1;
        [i,j] = find(on & ~missing);
        pt_img = [i j];

        pt_motor = space_img_to_motor(pt_img);    
        rg = zeros(1,2);
        rg(1) = range(pt_motor(:,1));
        rg(2) = range(pt_motor(:,2));        
    end
 
   end
        
end