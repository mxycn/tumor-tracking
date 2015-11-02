function pct_error=get_percentile_error(errorx,percentagex,n);
% errorx is an array.
nx=round(n*(1-percentagex)); % nx is the number of errors exceeding the threshold error
[initialx imin]=min(errorx(1:nx)); % choose the first nx errors, and find its minimum and corresponding index
        
for i=nx+1:n
    if errorx(i)>initialx %if a greater error is found, turn it into the collection.
        errorx(imin)=errorx(i);
    end
    [initialx imin]=min(errorx(1:nx)); % find the new minimum error in the collection.    
end
pct_error=initialx;
return   
    