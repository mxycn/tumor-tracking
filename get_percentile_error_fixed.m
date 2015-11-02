function pct_error=get_percentile_error_fixed(errorx,percentagex,j0,jx);
% errorx is an array.
% seperate the array with first (1-percentagex) part being the candidate large numbers.
newarray=zeros(jx-j0+1,1);
ik=1;
for ii=1:jx-j0+1
    if ~isnan(errorx(ii))
        newarray(ik)=errorx(ii);
        ik=ik+1;
    end
end        
pct_error=get_percentile_error(newarray,percentagex,1,ik-1);
return
    
    