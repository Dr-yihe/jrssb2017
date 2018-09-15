function restmat = genrest(datstruc,val,ifbeta)
% generate new radius estimates
% dastruc contains the generating function of estimated radius
% val: the p or the beta
% ifbeta: if the given value is beta

if nargin<3
    ifbeta=false;
end


restmat=zeros(size(datstruc.rest))*nan;
ix=1:size(restmat,2);
for i=ix(datstruc.ifEff)
    if ~ifbeta
        restmat(:,i)=datstruc.restgfp{i}(val);
    else
        restmat(:,i)=datstruc.restgfbeta{i}(val);
    end
end



end

