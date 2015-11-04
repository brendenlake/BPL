classdef RelationIndependent < Relation
    % RELATIONINDEPENDENT Stroke relation where there is no relation at all 
    
   properties
      %type
      %nprev
      gpos % global position      
   end
   
   methods 
       
       function this = RelationIndependent(type,nprev,gpos)
           this.type = type;
           this.nprev = nprev;
           if exist('gpos','var')
               this.gpos = gpos;
           end
           assert(this.validType);
       end
       
   end
  
end