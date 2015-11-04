classdef RelationAttach < Relation
  % RELATIONATTACH Attach a stroke to either the "end" or "start" of a
  % previous stroke

   properties
      %type
      %nprev
      attach_spot % index from 1,...,nprev where this stroke attaches   
   end
   
   methods
       
       function this = RelationAttach(type,nprev,attach_spot)
          this.type = type;
          this.nprev = nprev;
          this.attach_spot = attach_spot;
          assert(this.validType);
       end 
       
   end
  
end