classdef RelationAttachAlong < RelationAttach
    % Attach a stroke somewhere along a previous stroke
    
   properties
      %type
      %nprev
      %attach_spot % index from 1,...,nprev where this stroke attaches
      nsub
      subid_spot % index from 1,...,nsub, or the number of subids in the stroke(attach_spot)
      ncpt % number of control points
      eval_spot_type = [];
      eval_spot_token = [];
   end
   
   methods
       
       function this = RelationAttachAlong(type,nprev,attach_spot,nsub,subid_spot,ncpt)
          this = this@RelationAttach(type,nprev,attach_spot);
          this.nsub = nsub;
          this.subid_spot = subid_spot;
          this.ncpt = ncpt;
       end
       
   end
  
end