function drawMulCon(datstruc,para,ix,iftrue,iflegend)
%
% Drawing the bivariate quantile contours
%

% datstruc: data and estimates
% para: enviromental parameters
% ix: the scenario index
% iftrue: if draw the true quantile contours
% iflegend: if display the legend

if nargin<4
    iftrue=true;iflegend=true;
elseif nargin<5
    iflegend=true;
end

scatter(datstruc.X(:,2*ix-1),datstruc.X(:,2*ix),10,'ok')
box on
hold on
for i=1:length(para.p)
    plot(max((datstruc.evt.restgfp{ix}(para.p(i))),0).*cos(para.angle),...
        max((datstruc.evt.restgfp{ix}(para.p(i))),0).*sin(para.angle),'--k','LineWidth',0.8)
    hold on
    if iftrue
        plot(max((datstruc.rtrue(:,i)),0).*cos(para.angle),...
            max((datstruc.rtrue(:,i)),0).*sin(para.angle),'k','LineWidth',0.8)
        hold on
    end
end
hold off
if iflegend
    if iftrue
        legend('Observation','EVT Estimate','True')
    else
        legend('Observation','EVT Estimate')
    end
end

axis equal

end

