disp('   Estimation of Extreme Depth-based Quantile Regions   ')
disp('       programmed by Yi He, Tilburg Univeristy       ')
disp('---------- Application to Real Data -----------------')
%% Import daily price data
stockpricedata=importdata('stockpricedata.csv');
MarketIndex=stockpricedata.data;
datadate=stockpricedata.textdata(2:end,1);
ror=diff(log(MarketIndex(1:1565,:))); % calculate index return
clearvars stockpricedata
%% Filteration: EGARCH(1,1)
egarch11=egarch('Offset',NaN,'GARCHLags',1,'ARCHLags',1,...
    'LeverageLags',1,'Distribution','t');
USmdl=estimate(egarch11,ror(:,1));
USVar=infer(USmdl,ror(:,1));
USres=(ror(:,1)-USmdl.Offset)./sqrt(USVar);

UKmdl=estimate(egarch11,ror(:,2));
UKVar=infer(UKmdl,ror(:,2));
UKres=(ror(:,2)-UKmdl.Offset)./sqrt(UKVar);

JPNmdl=estimate(egarch11,ror(:,3));
JPNVar=infer(JPNmdl,ror(:,3));
JPNres=(ror(:,3)-JPNmdl.Offset)./sqrt(JPNVar);

aggres=[USres UKres JPNres];
%% Ljung-Box Test
lbqtesth=zeros(3,3); % row: rt,|rt|,rt^2   column: US,UK,JPN
for i=1:3
    lbqtesth(1,i)=lbqtest(aggres(:,i));
    lbqtesth(2,i)=lbqtest(abs(aggres(:,i)));
    lbqtesth(3,i)=lbqtest(aggres(:,i).^2);
end
%% Checking Tail Indices Equality
k=80;rng(1)
X=[sort(aggres) sort(-aggres)];
gamsort=sort(mean(log(X(end-k+1:end,:))-ones(k,1)*log(X(end-k,:))));
disp(['The tail indices in order: ',num2str(gamsort)]);
disp(['Maximal Difference: ',num2str(gamsort(end)-gamsort(1))]);
gamsimsort=sort(normrnd(mean(gamsort),mean(gamsort)/sqrt(k),1000,6),2);
disp(['Approximate p-value:',num2str(nanmean(gamsimsort(:,end)-gamsimsort(:,1)>gamsort(end)-gamsort(1)))]);
clear gamsimsort X

% The bivariate regular variation of the filtered returns are tested by
% Andrea Krajina, Georgia Augusta University Gottingen
%% Initializating Parameters
para.S=1;
para.p=[1/10000,1/5000,1/2000];
para.M3D=100+1;
para.angle2=zeros(para.M3D,2)*nan;
para.angle2(:,1)=(linspace(0,pi,para.M3D))'; % another directions: inclination
para.angle2(:,2)=(linspace(0,2*pi,para.M3D))'; %azimuth
para.u3D=zeros(101*101,3);
for i=1:para.M3D
    para.u3D(((i-1)*para.M3D+1):(i*para.M3D),:)=[sin(para.angle2(i,1))*cos(para.angle2(:,2)),...
        sin(para.angle2(i,1))*sin(para.angle2(:,2)),cos(para.angle2(i,1))*ones(para.M3D,1)];
end
para.M=1001; % nr of directions for depth computation
para.angle=(linspace(0,2*pi,para.M))'; % directions: 
para.u=[cos(para.angle) sin(para.angle)];
%% Predicted Bivariate Quantile Regions: Figure 4
predvar=zeros(3,1);predmu=zeros(3,1);
predvar(1)=forecast(USmdl,1);predmu(1)=USmdl.Offset;
predvar(2)=forecast(UKmdl,1);predmu(2)=UKmdl.Offset;
predvar(3)=forecast(JPNmdl,1);predmu(3)=JPNmdl.Offset;
para.k_gam=160;para.k_nu=160;
indexname={'S&P 500','FTSE 100','Nikkei 225'};
countryname={'U.S.','U.K.','Japan'};
bivix=[1,2;1,3;2,3];
axislim=[0.07;0.1;0.11];
mktbivQ=cell(3,1);
for i=1:3
    mktbivQ{i}=MrvHDQEVT(para,aggres(:,bivix(i,:)),para.u,false,true);
    subplot(1,3,i)
    hold on;
     scatter(aggres(:,bivix(i,1))*sqrt(predvar(bivix(i,1)))+predmu(bivix(i,1)),...
         aggres(:,bivix(i,2))*sqrt(predvar(bivix(i,2)))+predmu(bivix(i,2)),'k')
    for j=1:length(para.p)
        plot(mktbivQ{i}.restgfp{1}(para.p(j)).*para.u(:,1)*sqrt(predvar(bivix(i,1)))+predmu(bivix(i,1)),...
            mktbivQ{i}.restgfp{1}(para.p(j)).*para.u(:,2)*sqrt(predvar(bivix(i,2)))+predmu(bivix(i,2)),'k')
    end
    axis equal
    xlim([-axislim(i) axislim(i)])
    ylim([-axislim(i) axislim(i)])
    xlabel(indexname{bivix(i,1)})
    ylabel(indexname{bivix(i,2)})
    title([countryname{bivix(i,1)},' vs ',countryname{bivix(i,2)}])
    box on
    hold off
end
    
%% Trivariate Quantile Regions: Figure 5
para.k_gam=200;para.k_nu=300;
mktresQ = MrvHDQEVT(para,aggres,para.u3D,false,true);

colormap('gray')
figure('position',[300 200 800 600])
[meshX meshY]=meshgrid(para.angle2(:,1),para.angle2(:,2));
mesh(meshX,meshY,reshape(mktresQ.rest,para.M3D,para.M3D),'EdgeColor','k')
xlim([0 pi])
xlabel('Inclination','interpreter','latex')
ylim([0 2*pi])
ylabel('Azimuth','interpreter','latex')
zlabel('Radius','interpreter','latex')
hold on
[resaz,resele,resr]=cart2sph(aggres(:,1),aggres(:,2),aggres(:,3));
resaz(resaz<0)=resaz(resaz<0)+2*pi;
resele(resele<0)=resele(resele<0)+pi;
scatter3(resele,resaz,resr,'k','fill');

%% Finding Outliers
norm_aggres=sqrt(sum(aggres.^2,2));
aggres_unit=aggres./(norm_aggres*ones(1,size(aggres,2)));
[~,aggres_uix]=max(aggres_unit*para.u3D',[],2);
disp('----------------------------------')
disp('Outliers detected:')
datadate([false;norm_aggres>=mktresQ.rest(aggres_uix)])
clear norm_aggres aggres_unit aggres_uix





