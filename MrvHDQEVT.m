function datstruc = MrvHDQEVT(para,X,u,ifbeta,ifkdiff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimation of Extreme Depth-based Quantile Regions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ifbeta: ifbeta is known
% para: the simulation enviroment parameters
% datstrucX: data structure contains the observations
% datstruc.rtrue: the true radius (in columns) at all directions
% datstruc.betatrue: the true beta value
% datstruc.rest: the estimated radius for para.p(1) or para.beta(1)

dim=size(u,2);

%
% Default setting: beta is unknown
%

if nargin<4
    ifbeta=false;
end
if nargin<5
    ifkdiff=false;
end

%
% Initialization
%
datstruc.rest=zeros(size(u,1),para.S)*nan;
datstruc.restgfp=cell(para.S,1);
datstruc.restgfbeta=cell(para.S,1);
%
% Estimation
%

display('Estimating using EVT approach...')
for i=1:para.S
    if mod(i-1,para.S/10)==0
        fprintf([num2str(10*(i-1)/para.S),','])
    end
    if ~ifkdiff
    [datstruc.restgfp{i},datstruc.restgfbeta{i}]...
        =Q_hat(X(:,((i-1)*dim+1):(i*dim)),para.k,para.k,u);
    else
        [datstruc.restgfp{i},datstruc.restgfbeta{i}]...
        =Q_hat(X(:,((i-1)*dim+1):(i*dim)),para.k_gam,para.k_nu,u);
    end
    if ifbeta
        datstruc.rest(:,i)=datstruc.restgfbeta{i}(para.beta(1));
    else
        datstruc.rest(:,i)=datstruc.restgfp{i}(para.p(1));
    end
end
fprintf('10.\n')
datstruc.ifEff=(prod(datstruc.rest>0)>0);
datstruc.effN=sum(datstruc.ifEff);
end

