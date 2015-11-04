classdef Relation
    
   properties
      type = ''; 
      nprev % number of previous strokes
      types_allowed = {'unihist','start','end','mid'};  
   end
  
   methods
       
       % is this relation of a valid type?
       function out = validType(this)
           out = any(strcmp(this.type,this.types_allowed));           
       end
       
   end
   
end