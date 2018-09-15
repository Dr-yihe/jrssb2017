function [X,datatype]= cloverrnd(n,ifDisplay)

% Simulate n iid sample from bivariate clover distribution

if nargin<2
    ifDisplay=false;
end
datatype='clover';
u=rand(n,1);
r=((1-u).^(-2)-1).^(1/6);
v=rand(n,1);
theta=zeros(n,1);
for i=1:n
    if ifDisplay
        display(['Randomizing the ',num2str(i),'-th Observation...'])
    end
    theta(i)=double(fzero(@(x)((3/(10*pi)*r(i)^5/(1+r(i)^6)^(3/2)*(5*x+sin(4*x)))...
        /(3*r(i).^5./(1+r(i).^6).^(3/2))-v(i)),0));
end
X=[r.*cos(theta) r.*sin(theta)];

end

