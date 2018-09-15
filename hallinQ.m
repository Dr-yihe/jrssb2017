function [A,b] = hallinQ(X,p,u)
%
% Regression in Hallin, Paindaveine and Siman (2010)
%
% X: n*dimX
% u: M*dimX
% p: scalar

if nargin<3
    error('Not enough inputs')
end
[M,dimX]=size(u);
S=size(X,2)/dimX;
if S~=fix(S)
    error('dimension of direction is not correct')
end

y=X*u'; % n*M;
q=zeros(M,dimX)*nan;
A=zeros(M,dimX)*nan;b=zeros(M,1)*nan;
for i=1:M
    gamu=null(u(i,:)); % dimX*(dimX-1)
    q(i,:)=rq([ones(size(X,1),1) X*gamu],y(:,i),p);
    A(i,:)=-(u(i,:)-q(i,2:end)*gamu');
    b(i)=-q(i,1);
end


end

