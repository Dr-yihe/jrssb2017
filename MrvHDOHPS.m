function datstruc = MrvHDOHPS(X,u,beta)
if nargin<3
    disp('Not enough input')
end
[M,dimX]=size(u);
S=size(X,2)/dimX;
if S~=fix(S)
    error('dimension of direction is not correct')
end
datstruc.rest=zeros(M,S)*nan;
display('Estimating using HPS approach...')

for i=1:S
    if mod(i-1,S/10)==0
        fprintf([num2str(10*(i-1)/S),','])
    end
    try
         [A,b] = hallinQ(X(:,(i-1)*dimX+1:i*dimX),beta,u);
         % A: M*dimX
         % b: M*1
    catch err
        continue
    end
    datstruc.rest(:,i)=(1./max((A*u')./(b*ones(1,M))))';
    
end
fprintf('10.\n')
datstruc.ifEff=(prod(datstruc.rest>0)>0);
datstruc.effN=sum(datstruc.ifEff);

end

