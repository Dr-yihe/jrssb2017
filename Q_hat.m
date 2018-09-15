function [restgfp,restgfbeta,gamma_est,nu_S_hat] =Q_hat(X,k_gam,k_nu,u)

% Extreme Estimator of Quantile Region
% X: N*d input data
% k_gam: choice of k for estimating gamma
% k_nu: choice of k for estimating nu
% u: M*d matrix contaning the directional vectors on unit sphere in rows


[n,dimX]=size(X); M=size(u,1);

if nargin<4
    error('not enough input');
end

norm_x=sqrt(sum(X.^2,2)); 
norm_x_sort=sort(norm_x);
tailData=log(norm_x_sort(end-k_gam:end,:));
gamma_est=mean(tailData(2:end,:)-ones(k_gam,1)*tailData(1,:));
alpha_est=1/gamma_est;
U_est=norm_x_sort(end-k_nu+1);
nu_hat_uniths=sum(u*X'>=U_est,2)/k_nu; % M*1 vector
hd_w_nuhatstar=min(((max((u*u'),0)).^(-alpha_est)).*(ones(M,1)*nu_hat_uniths'),[],2); % M*1 matrix

X_unit=X./(norm_x*ones(1,dimX));
[~,uix]=max(X_unit*u',[],2);

nu_S_hat=sum(norm_x/U_est>=hd_w_nuhatstar(uix).^gamma_est)/k_nu;

restgfbeta=@(b)(U_est*(k_nu/(n*b))^gamma_est*(hd_w_nuhatstar).^(1/alpha_est));
restgfp=@(p)(U_est*(k_nu*nu_S_hat/(n*p))^gamma_est*(hd_w_nuhatstar).^(1/alpha_est));
end