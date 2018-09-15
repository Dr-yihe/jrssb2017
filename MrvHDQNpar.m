function datstruc= MrvHDQNpar(X,u)
% Nonparametric estimation
%
% Beta known: only consider the case beta=1/n
%
dimX=size(u,2); S=size(X,2)/dimX;
if S~=fix(S)
    error('dimension of direction is not correct')
end
datstruc.rest=zeros(size(u,1),S)*nan;

display('Estimating using Non-parametric approach...')

for i=1:S
    
    if mod(i-1,S/10)==0
        fprintf([num2str(10*(i-1)/S),','])
    end
    % determine the extreme point
    k=convhulln(X(:,(i-1)*dimX+1:i*dimX));
    
    % revert the linear contstraints
    A=zeros(size(k,1),dimX)*nan;
    for j=1:size(k,1)
        A(j,:)=ones(1,dimX)/(X(k(j,:),(i-1)*dimX+1:i*dimX))';
    end
    datstruc.rest(:,i)=(1./max((A*u')))';
end
fprintf('10.\n')
datstruc.ifEff=(prod(datstruc.rest>0)>0);
datstruc.effN=sum(datstruc.ifEff);


end

