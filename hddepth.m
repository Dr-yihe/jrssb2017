function d = hddepth(x,datatype,A,b)

%
% Calculate the true halfspace depth value
%

n=size(x,1);
if nargin==4
    x=(x-ones(n,1)*b')/(A');
elseif mod(nargin,2)==1
    error('Not enough input')
end


switch datatype
    case 'bicauchy'
        d=1-tcdf(sqrt(sum(x.^2,2)),1);
    case 'bistudent'
        d=1-mvtcdf([sqrt(sum(x.^2,2)),inf*ones(n,1)],eye(2),3);
    case 'biellp'
        radius_x=sqrt(x(:,1).^2/4+x(:,2).^2);
        % the depth value is approximated value
        d=(1+radius_x.^6).^(-1/2)*2/(3*pi);
    case 'tricauchy'
        d=1-tcdf(sqrt(sum(x.^2,2)),1);
    otherwise
        error('Unknown datatype')
end
end

