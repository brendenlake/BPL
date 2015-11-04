% Same syntax as subplot(nrow,ncol,indx)
% except that everytying is plotted much tighter together.
%
function sptight(nrow,ncol,indx)
    if (indx > nrow*ncol)
       error('invalid index'); 
    end
    
    intcol = 1/ncol;
    introw = 1/nrow;

    % compute the row,col we are in
    myrow = ceil(indx/ncol);
    mycol = mod(indx,ncol);
    if mycol==0, mycol = ncol; end
    
    left = intcol*(mycol-1);    
    bottom = introw*(nrow-myrow);
    
    % make sub-plot
    subplot('Position',[left,bottom,intcol,introw]);
    set(gca,'XTick',[],'YTick',[]);    
end