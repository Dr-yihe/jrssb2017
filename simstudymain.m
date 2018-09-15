disp('   Estimation of Extreme Depth-based Quantile Regions   ')
disp('       programmed by Yi He, Tilburg Univeristy       ')
disp('---------- Simulation Study -----------------')
%% Parameter Settings
disp('Initialzing envrioment parameters...')
para.n=5000; % sample size
%
% In the paper we also consider the case:
% para.n=1000;
% 
para.S=100; % number of scenarios 
para.Ac=[2 0.3;0.3 1]; % scaling and rotation
para.bc=[3;2]; % shift of data center
para.p=[1/5000,1/2000,1/10000]; 
para.beta=1/para.n;
para.M=1001; % nr of directions for depth computation
para.angle=(linspace(0,2*pi,para.M))'; % directions: 
para.u=[cos(para.angle) sin(para.angle)]; % points on the unit sphere
para.M3D=100+1; 
para.angle2=zeros(para.M3D,2)*nan;
para.angle2(:,1)=(linspace(0,pi,para.M3D))'; % inclination
para.angle2(:,2)=(linspace(0,2*pi,para.M3D))'; %azimuth
para.u3D=zeros(101*101,3);
for i=1:para.M3D
    para.u3D(((i-1)*para.M3D+1):(i*para.M3D),:)=[sin(para.angle2(i,1))*cos(para.angle2(:,2)),...
        sin(para.angle2(i,1))*sin(para.angle2(:,2)),cos(para.angle2(i,1))*ones(para.M3D,1)];
end
warning('off','all') % turn off system warning messages
rng(1) % set seed
clear i
%% True Radius Functiona and Depth Value Beta

disp('Computing the true radius function and the underlying beta...')
%
% First compute the true radius for all the distributions except bivariate
% clover
%

% Bivariate Cauchy Distribution
bicauchy.rtrue=sqrt(1./(para.p).^2-1);

% Bivariate Student-3 Distribution
bistudent.rtrue=sqrt(3*((para.p).^(-2/3)-1));

% Bivariate Elliptical
biellp.rtrue=zeros(length(para.angle),length(para.p))*nan;
for i=1:length(para.p)
    biellp.rtrue(:,i)=((para.p(i).^(-2)-1).^(1/6))./sqrt(cos(para.angle).^2/4+sin(para.angle).^2);
end

% Affine Bivariate Cauchy Distribution
bicauchyaff.rtrue=zeros(para.M,length(para.p))*nan;
for j=1:length(para.p)
    for i=1:para.M
        tempb=@(r)(para.Ac\[r*cos(para.angle(i))-para.bc(1);r*sin(para.angle(i))-para.bc(2)]);
        bicauchyaff.rtrue(i,j)=fzero(@(r)(tempb(r)'*tempb(r)-bicauchy.rtrue(j)^2),bicauchy.rtrue(j));
    end
end

% Trivariate Cauchy Distribution
tricauchy.rtrue=zeros(1,length(para.p))*nan;
for i=1:length(para.p)
    tricauchy.rtrue(i)=fzero(@(r)(2/pi*(atan(r)-r/(r^2+1))-(1-para.p(i))),1./para.p(i));
end

%
% Compute the underlying depth beta
%
bicauchy.betatrue=1-tcdf(bicauchy.rtrue,1);
bistudent.betatrue=1-mvtcdf([bistudent.rtrue',inf*ones(length(para.p),1)],eye(2),3);
biellp.betatrue=para.p*2/(3*pi);
bicauchyaff.betatrue=bicauchy.betatrue;
tricauchy.betatrue=1-mvtcdf([tricauchy.rtrue',inf*ones(length(para.p),2)],eye(3),1);

% Store the true probability p
bicauchy.ptrue=para.p;bistudent.ptrue=para.p;
biellp.ptrue=para.p;bicauchyaff.ptrue=para.p;
biclover.ptrue=para.p;tricauchy.ptrue=para.p;
%% Data Simulation
disp('Drawing random samples...')
% Initialization

bicauchy.dim=2;bistudent.dim=2;biellp.dim=2;bicauchyaff.dim=2;
tricauchy.dim=3;biclover.dim=2;

bicauchy.X=zeros(para.n,2*para.S);bistudent.X=zeros(para.n,2*para.S); 
biellp.X=zeros(para.n,2*para.S);bicauchyaff.X=zeros(size(bicauchy.X));
tricauchy.X=zeros(para.n,3*para.S);biclover.X=zeros(para.n,2*para.S);

% Data generating
for i=1:para.S
    if mod(i-1,para.S/10)==0
        fprintf([num2str(10*(i-1)/para.S),','])
    end
    bicauchy.X(:,[2*i-1 2*i])=mvtrnd(eye(2),1,para.n); % Bivariate Cauchy
    bistudent.X(:,[2*i-1 2*i])=mvtrnd(eye(2),3,para.n); % Bivariate student with nu=3
    biellp.X(:,[2*i-1 2*i])=elliprnd(para.n); % Bivariate Ellptical
    bicauchyaff.X(:,[2*i-1 2*i])=bicauchy.X(:,[2*i-1 2*i])*para.Ac'+ones(para.n,1)*para.bc'; % Aff Biv Cauchy
    tricauchy.X(:,[3*i-2 3*i-1 3*i])=mvtrnd(eye(3),1,para.n); %Trivariate Cauchy
    biclover.X(:,[2*i-1 2*i])=cloverrnd(para.n);  
end
fprintf('10.\n')

%% EVT Estimate and Relative Errors
para.k=400; % choice of k
%
% For n=5000 we also consider the cases k=200 and 800;
% For n=1000 we consider the case k=150.
% 

disp('-------- Bivariate Cauchy Distribution ----------')
bicauchy.evt = MrvHDQEVT(para,bicauchy.X,para.u);
bicauchy.evt = CmpErr(para,bicauchy.evt,bicauchy,'bicauchy');

disp('-------- Bivariate Student Distribution ----------')
bistudent.evt = MrvHDQEVT(para,bistudent.X,para.u);
bistudent.evt = CmpErr(para,bistudent.evt,bistudent,'bistudent');

disp('-------- Bivariate Ellptical Distribution ----------')
biellp.evt = MrvHDQEVT(para,biellp.X,para.u);
biellp.evt = CmpErr(para,biellp.evt,biellp,'biellp');

disp('-------- Bivariate Affine Cauchy Distribution ----------')
bicauchyaff.evt = MrvHDQEVT(para,bicauchyaff.X,para.u);
bicauchyaff.evt = CmpErr(para,bicauchyaff.evt,bicauchyaff,'bicauchy',true);

disp('-------- Bivariate Clover Distribution ----------')
biclover.evt = MrvHDQEVT(para,biclover.X,para.u);

disp('-------- Trivariate Cauchy Distribution ----------')
tricauchy.evt = MrvHDQEVT(para,tricauchy.X,para.u3D);
tricauchy.evt = CmpErr(para,tricauchy.evt,tricauchy,'tricauchy');

clear i j ans tempb
%% Approximation of the true radius function for Biclover distribution
biclover.rtrue=zeros(size(biclover.evt.rest,1),3);
biclover.rtrue(:,1)=exp(nanmean(log(biclover.evt.rest(:,biclover.evt.ifEff)),2));
biclover.rtrue(:,2)=exp(nanmean(log(genrest(biclover.evt,para.p(2))),2));
biclover.rtrue(:,3)=exp(nanmean(log(genrest(biclover.evt,para.p(3))),2));

tempmaxrlogerr=zeros(para.S,1);
for i=1:para.S
    tempmaxrlogerr(i)=max(abs(log(biclover.evt.rest(:,i))-log(biclover.rtrue(:,1))));
end
[~,biclover.evt.medianix]=min(abs(tempmaxrlogerr-nanmedian(tempmaxrlogerr)));
clear tempmaxrlogerr
%% Display Relative Errors: Table 1 (Manual change the choice of k and sample size n is needed to obtain the full table)
display('Relative Errors: Biv Cauchy, Biv Student-3, Elliptical, Affine Cauchy, Tri Cauchy');
display(['Prob RE: ', num2str(nanmedian([bicauchy.evt.pesterr,bistudent.evt.pesterr,biellp.evt.pesterr,...
    bicauchyaff.evt.pesterr,tricauchy.evt.pesterr]/para.p(1)))])
display(['Depth RE: ', num2str(nanmedian([max(abs(bicauchy.evt.betalogerr));max(abs(bistudent.evt.betalogerr));...
    max(abs(biellp.evt.betalogerr));max(abs(bicauchyaff.evt.betalogerr));...
    max(abs(tricauchy.evt.betalogerr))]'))])
%% Plot the EVT estimate examples: Figure 1

% Bivariate Cauchy
subplot(2,3,1)
drawMulCon(bicauchy,para,30,true,false)
xlim([-1.3 1.3]./min(para.p))
ylim([-1.3 1.3]./min(para.p))
title('Bivariate Cauchy')
axis_formatstring = '%3.0f';
set(gca, 'Units', 'normalized', 'Position', [0.05 0.6 0.33 0.33])
xtick = get(gca, 'xtick');
for i = 1:length(xtick)
     xticklabel{i} = sprintf(axis_formatstring, xtick(i));
end
set(gca, 'xticklabel', xticklabel)
ytick = get(gca, 'ytick');
for i = 1:length(ytick)
    yticklabel{i} = sprintf(axis_formatstring, ytick(i));
end
set(gca, 'yticklabel', yticklabel)

% Bivariate Student-3
subplot(2,3,2)
drawMulCon(bistudent,para,70,true,false)
 xlim([-5 5]*10)
 ylim([-5 5]*10)
title('Bivariate Student')
set(gca, 'Units', 'normalized', 'Position', [0.36 0.6 0.33 0.33])

% Bivariate Ellp
subplot(2,3,3)
drawMulCon(biellp,para,30,true,false)
xlim([-5 5]*10)
 ylim([-5 5]*10)
title('Bivariate Elliptical')
set(gca, 'Units', 'normalized', 'Position', [0.68 0.6 0.33 0.33])

% Affine Bivariate Cauchy
subplot(2,3,4)
drawMulCon(bicauchyaff,para,90,true,false)
axis equal
xlim([-2.8 2.8]*10^4)
ylim([-2.8 2.8]*10^4)
title('Aff. Biv. Cauchy')
set(gca, 'Units', 'normalized', 'Position', [0.18 0.15 0.33 0.33])

xtick = get(gca, 'xtick');
for i = 1:length(xtick)
     xticklabel{i} = sprintf(axis_formatstring, xtick(i));
end
set(gca, 'xticklabel', xticklabel)
ytick = get(gca, 'ytick');
for i = 1:length(ytick)
    yticklabel{i} = sprintf(axis_formatstring, ytick(i));
end
set(gca, 'yticklabel', yticklabel)

% Bivariate Clover
subplot(3,2,6)
drawMulCon(biclover,para,40,true,false)
xlim([-3 3]*10)
ylim([-3 3]*10)
title('Bivariate Clover')
set(gca, 'Units', 'normalized', 'Position', [0.54 0.15 0.33 0.33])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparsion with the estimation approach by Hallin et al. (2010) and
% Kong and Mizera (2012), when beta is known
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating the true radius for known beta=1/n, n=5000

disp('Computing the true radius function and the underlying beta...')

% Bivariate Cauchy Distribution
bicauchybeta.rtrue=tinv(1-para.beta,1);

% Bivariate Student-3 Distribution
bistudentbeta.rtrue=fzero(@(r)(mvtcdf([r,inf],eye(2),3)-(1-para.beta)),tinv(1-para.beta,3));

% Bivariate Elliptical
biellpbeta.rtrue=fzero(@(r)(hddepth([0,r],'biellp')-para.beta),para.beta^(-1/3))./sqrt(cos(para.angle).^2/4+sin(para.angle).^2);

% Affine Bivariate Cauchy Distribution
bicauchyaffbeta.rtrue=zeros(para.M,1)*nan;
for i=1:para.M
    tempb=@(r)(para.Ac\[r*cos(para.angle(i))-para.bc(1);r*sin(para.angle(i))-para.bc(2)]);
    bicauchyaffbeta.rtrue(i)=fzero(@(r)(tempb(r)'*tempb(r)-bicauchybeta.rtrue^2),bicauchybeta.rtrue);
end

%Trivariate Cauchy Distribution
tricauchybeta.rtrue=tinv(1-para.beta,1);

%
% Storing the depth value beta
%
bicauchybeta.betatrue=para.beta;bistudentbeta.betatrue=para.beta;
bicauchyaffbeta.betatrue=para.beta;binormalbeta.betatrue=para.beta;
biellpbeta.betatrue=para.beta;tricauchybeta.betatrue=para.beta;

%
% Calculate the true probability p
%
bicauchybeta.ptrue=(bicauchybeta.rtrue^2+1)^(-1/2);
bistudentbeta.ptrue=(bistudentbeta.rtrue^2/3+1)^(-3/2);
biellpbeta.ptrue=((biellpbeta.rtrue(1)/2)^6+1)^(-1/2);
bicauchyaffbeta.ptrue=bicauchybeta.ptrue;
tempf=@(x)((atan(x)-x/(x^2+1))*2/pi);
tricauchybeta.ptrue=1-tempf(tricauchybeta.rtrue);
clear tempf

%% EVT, HPS and KM estimate of Q
disp('Estimating using EVT, HPS and KM approaches...');

% Bivariate Cauchy
disp('-------- Bivariate Cauchy Distribution ----------')
bicauchybeta.evt.rest=genrest(bicauchy.evt,para.beta,true);
bicauchybeta.evt.ifEff=bicauchy.evt.ifEff;
bicauchybeta.dim=2;
bicauchybeta.evt = CmpErr(para,bicauchybeta.evt,bicauchybeta,'bicauchy');
bicauchybeta.Npar = MrvHDQNpar(bicauchy.X,para.u);
bicauchybeta.Npar = CmpErr(para,bicauchybeta.Npar,bicauchybeta,'bicauchy');
bicauchybeta.KM = MrvHDOKM(bicauchy.X,para.u,para.beta);
bicauchybeta.KM = CmpErr(para,bicauchybeta.KM,bicauchybeta,'bicauchy');
bicauchybeta.HPS=MrvHDOHPS(bicauchy.X,para.u,para.beta);
bicauchybeta.HPS = CmpErr(para,bicauchybeta.HPS,bicauchybeta,'bicauchy');

% Bivariate Student Distribution
disp('-------- Bivariate Student Distribution ----------')
bistudentbeta.evt.rest=genrest(bistudent.evt,para.beta,true);
bistudentbeta.evt.ifEff=bistudent.evt.ifEff;
bistudentbeta.dim=2;
bistudentbeta.evt = CmpErr(para,bistudentbeta.evt,bistudentbeta,'bistudent');
bistudentbeta.Npar = MrvHDQNpar(bistudent.X,para.u);
bistudentbeta.Npar = CmpErr(para,bistudentbeta.Npar,bistudentbeta,'bistudent');
bistudentbeta.KM = MrvHDOKM(bistudent.X,para.u,para.beta);
bistudentbeta.KM = CmpErr(para,bistudentbeta.KM,bistudentbeta,'bistudent');
bistudentbeta.HPS=MrvHDOHPS(bistudent.X,para.u,para.beta);
bistudentbeta.HPS = CmpErr(para,bistudentbeta.HPS,bistudentbeta,'bistudent');

% Bivariate Eilliptical Distribution
disp('-------- Bivariate Elliptical Distribution ----------')
biellpbeta.evt.rest=genrest(biellp.evt,para.beta,true);
biellpbeta.evt.ifEff=biellp.evt.ifEff;
biellpbeta.dim=2;
biellpbeta.evt = CmpErr(para,biellpbeta.evt,biellpbeta,'biellp');
biellpbeta.Npar = MrvHDQNpar(biellp.X,para.u);
biellpbeta.Npar = CmpErr(para,biellpbeta.Npar,biellpbeta,'biellp');
biellpbeta.KM = MrvHDOKM(biellp.X,para.u,para.beta);
biellpbeta.KM = CmpErr(para,biellpbeta.KM,biellpbeta,'biellp');
biellpbeta.HPS=MrvHDOHPS(biellp.X,para.u,para.beta);
biellpbeta.HPS = CmpErr(para,biellpbeta.HPS,biellpbeta,'biellp');

% Affine Bivariate Cauchy Distribution
disp('-------- Affined Bivariate Cauchy Distribution ----------')
bicauchyaffbeta.evt.rest=genrest(bicauchyaff.evt,para.beta,true);
bicauchyaffbeta.evt.ifEff=bicauchy.evt.ifEff;
bicauchyaffbeta.dim=2;
bicauchyaffbeta.evt = CmpErr(para,bicauchyaffbeta.evt,bicauchyaffbeta,'bicauchy',true);
bicauchyaffbeta.Npar = MrvHDQNpar(bicauchyaff.X,para.u);
bicauchyaffbeta.Npar = CmpErr(para,bicauchyaffbeta.Npar,bicauchyaffbeta,'bicauchy',true);
bicauchyaffbeta.KM = MrvHDOKM(bicauchyaff.X,para.u,para.beta);
bicauchyaffbeta.KM = CmpErr(para,bicauchyaffbeta.KM,bicauchyaffbeta,'bicauchy',true);
bicauchyaffbeta.HPS=MrvHDOHPS(bicauchyaff.X,para.u,para.beta);
bicauchyaffbeta.HPS = CmpErr(para,bicauchyaffbeta.HPS,bicauchyaffbeta,'bicauchy',true);

% Trivariate Cauchy Distribution
disp('-------- Trivariate Cauchy Distribution ----------')
tricauchybeta.evt.rest=genrest(tricauchy.evt,para.beta,true);
tricauchybeta.evt.ifEff=tricauchy.evt.ifEff;
tricauchybeta.dim=3;
tricauchybeta.evt = CmpErr(para,tricauchybeta.evt,tricauchybeta,'tricauchy');
tricauchybeta.Npar = MrvHDQNpar(tricauchy.X,para.u3D);
tricauchybeta.Npar = CmpErr(para,tricauchybeta.Npar,tricauchybeta,'tricauchy');
tricauchybeta.KM = MrvHDOKM(tricauchy.X,para.u3D,para.beta);
tricauchybeta.KM = CmpErr(para,tricauchybeta.KM,tricauchybeta,'tricauchy');
tricauchybeta.HPS=MrvHDOHPS(tricauchy.X,para.u3D,para.beta);
tricauchybeta.HPS = CmpErr(para,tricauchybeta.HPS,tricauchybeta,'tricauchy');

%% An illustrative example from Bivariate Cauchy Distribution: Figure 2
ix=13;
figure('Position',[200 200 800 800])
scatter(bicauchy.X(:,2*ix-1),bicauchy.X(:,2*ix),'k')
hold on
plot(bicauchybeta.rtrue.*cos(para.angle),bicauchybeta.rtrue.*sin(para.angle),'k')
plot(bicauchybeta.evt.rest(:,ix).*cos(para.angle),bicauchybeta.evt.rest(:,ix).*sin(para.angle),'--k')
plot(bicauchybeta.Npar.rest(:,ix).*cos(para.angle),bicauchybeta.Npar.rest(:,ix).*sin(para.angle),':k','LineWidth',1.5)
plot(bicauchybeta.KM.rest(:,ix).*cos(para.angle),bicauchybeta.KM.rest(:,ix).*sin(para.angle),'--k','LineWidth',1.8)
legend('Observation','True','EVT','NPar/HPS','KM',4)
axis equal
%axis tight
xlim([-2000 20000])
ylim([-2000 2500])
box on
hold off

%% Boxplots of Probability Relative Errors: Figure 3 part 1
boxname={'EVT','NPar','KM','HPS'};
% Bivariate Cauchy Distribution
subplot(2,5,1)
boxplot([bicauchybeta.evt.pesterr,bicauchybeta.Npar.pesterr,...
    bicauchybeta.KM.pesterr,bicauchybeta.HPS.pesterr]/bicauchybeta.ptrue,boxname,'color','k','symbol','k+');
ylim([0,4])
title('Biv. Cauchy')

% Bivariate Student-3 Distribution

subplot(2,5,2)
boxplot([bicauchybeta.evt.pesterr,bicauchybeta.Npar.pesterr,...
    bicauchybeta.KM.pesterr,bicauchybeta.HPS.pesterr]/bicauchybeta.ptrue,boxname,'color','k','symbol','k+');
ylim([0,4])
title('Biv. Student')

% Bivariate Elliptical Distribution
subplot(2,5,3)
boxplot([biellpbeta.evt.pesterr,biellpbeta.Npar.pesterr,...
    biellpbeta.KM.pesterr,biellpbeta.HPS.pesterr]/biellpbeta.ptrue,...
    boxname,'color','k','symbol','k+');
ylim([0,4])
title('Biv. Elliptical')

% Affine Bivariate Cauchy Distribution
subplot(2,5,4)
boxplot([bicauchyaffbeta.evt.pesterr,bicauchyaffbeta.Npar.pesterr,...
    bicauchyaffbeta.KM.pesterr,bicauchyaffbeta.HPS.pesterr]/bicauchyaffbeta.ptrue,...
    boxname,'color','k','symbol','k+');
ylim([0,4])
title('Aff. Cauchy')

% Trivariate Cauchy Distribution
subplot(2,5,5)
boxplot([tricauchybeta.evt.pesterr,tricauchybeta.Npar.pesterr,...
    tricauchybeta.KM.pesterr,tricauchybeta.HPS.pesterr]/tricauchybeta.ptrue,boxname,'color','k','symbol','k+');
ylim([0,4])
title('Triv. Cauchy')

%% Boxlplots of Depth Relative Errors: Figure 3 part 2
% Bivariate Cauchy Distribution
subplot(2,5,6)
boxplot([max(abs(bicauchybeta.evt.betalogerr))',max(abs(bicauchybeta.Npar.betalogerr))',...
    max(abs(bicauchybeta.KM.betalogerr))',max(abs(bicauchybeta.HPS.betalogerr))'],...
    boxname,'color','k','symbol','k+');
ylim([0,6])
title('Biv. Cauchy')

% Bivariate Student Distribution
subplot(2,5,7)
boxplot([max(abs(bistudentbeta.evt.betalogerr))',max(abs(bistudentbeta.Npar.betalogerr))',...
    max(abs(bistudentbeta.KM.betalogerr))',max(abs(bistudentbeta.HPS.betalogerr))'],...
    boxname,'color','k','symbol','k+');
ylim([0,6])
title('Biv. Student')

% Bivariate Elliptical Distribution
subplot(2,5,8)
boxplot([max(abs(biellpbeta.evt.betalogerr))',max(abs(biellpbeta.Npar.betalogerr))',...
    max(abs(biellpbeta.KM.betalogerr))',max(abs(biellpbeta.HPS.betalogerr))'],...
    boxname,'color','k','symbol','k+');
ylim([0,6])
title('Biv. Elliptical')

% Affine Bivariate Cauchy Distribution
subplot(2,5,9)
boxplot([max(abs(bicauchyaffbeta.evt.betalogerr))',max(abs(bicauchyaffbeta.Npar.betalogerr))',...
    max(abs(bicauchyaffbeta.KM.betalogerr))',max(abs(bicauchyaffbeta.HPS.betalogerr))'],...
    boxname,'color','k','symbol','k+');
ylim([0,6])
title('Aff. Cauchy')

% Trivariate Cauchy Distribution
subplot(2,5,10)
boxplot([max(abs(tricauchybeta.evt.betalogerr))',max(abs(tricauchybeta.Npar.betalogerr))',...
    max(abs(tricauchybeta.KM.betalogerr))',max(abs(tricauchybeta.HPS.betalogerr))'],...
    boxname,'color','k','symbol','k+');
ylim([0,6])
title('Triv. Cauchy')














