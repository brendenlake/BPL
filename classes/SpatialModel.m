classdef SpatialModel < BetterHandle
   % SPATIALMODEL
   %
   % Stores a set of SpatialHistograms, one for each stroke position,
   % and can evaluate the likelihood/sample new positions.
   % 
   properties
       list_SH        
   end
   
   properties (Dependent = true)
       last_model_id % all stroke ids of this value, or higher, use the last "group" model     
   end
   
   methods
       
       function this = SpatialModel(data_start,data_id,clump_id,xlim,ylim,nbin_per_side,prior_count)
       
            n = size(data_start,1);
            assert(numel(data_id)==n);
           
            % Learn specific spatial models
            this.list_SH = cell(clump_id,1);
            for i=1:clump_id-1
                data = data_start(data_id==i,:);
                this.list_SH{i} = SpatialHist(data,xlim,ylim,nbin_per_side,prior_count);
            end           
            
            % lump together datapoints from strokes after and includingclump_id
            sel = data_id >= clump_id;
            this.list_SH{end} = SpatialHist(data_start(sel,:),xlim,ylim,nbin_per_side,prior_count);
       
       end
       
       % stroke ids after this are given to the same model (inclusive)
       function out = get.last_model_id(this)
            out = numel(this.list_SH);
       end
       
       %
       % Compute log-likelihood of new points
       %
       % Input
       %  data_start: [n x 2] positions
       %  data_id: [n x 1] the stroke index of each position
       %
       % Output
       %  ll: [scalar] total log-likelihood
       function ll = score(this,data_start,data_id)
           
           new_id = this.map_indx(data_id);
           ndat = size(data_start,1);
          
           % for each stroke id
           ll = 0;
           for sid=1:this.last_model_id
             sel = new_id == sid;
             data = data_start(sel,:);             
             ll = ll + this.list_SH{sid}.score(data);
           end
              
           
       end
       
       %
       % Compute log-likelihood of new points,
       % and return breakdown for each one
       %
       % Input
       %  data_start: [n x 2] positions
       %  data_id: [n x 1] the stroke index of each position
       %
       % Output
       %  ll: [n x 1] the log-likelihood of each position
       function ll = score_vec(this,data_start,data_id)
           new_id = this.map_indx(data_id);
           ndat = size(data_start,1);
           ll = zeros(ndat,1);
           for sid=1:this.last_model_id
             sel = new_id == sid;
             data = data_start(sel,:);             
             [~,ll(sel)] = this.list_SH{sid}.get_id(data);
           end
       end
       
       %
       % Sample new stroke start positions
       %
       % Input
       %  data_id: [nsamp x 1] the stroke index of each position
       %
       % Output
       %  sampels: [nsamp x 2] positions drawn from the model
       function samples = sample(this,data_id)
           
           assert(isvector(data_id));
           nsamp = numel(data_id);
           new_id = this.map_indx(data_id);
           
           % for each stroke id
           samples = zeros(nsamp,2);
           for sid=1:this.last_model_id
             sel = new_id == sid;
             nsel = sum(sel);
             samples(sel,:) = this.list_SH{sid}.sample(nsel);
           end              
           
       end
       
       %
       % Plot the array of position models 
       %
       function plot(this)          
           n = this.last_model_id;
           nrow = ceil(sqrt(n));
           figure
           for sid=1:this.last_model_id
             subplot(nrow,nrow,sid);
             this.list_SH{sid}.plot();
             title(num2str(sid));
           end 
       end       
       
   end
   
   methods (Access = private)
       
       % map stroke ids to new ids
       function new_id = map_indx(this,old_id)
           new_id = old_id;
           new_id(new_id > this.last_model_id) = this.last_model_id; 
       end       
       
   end
    
end