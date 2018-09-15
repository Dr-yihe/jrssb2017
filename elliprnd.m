function [X,datatype] = elliprnd(n)
    % Simulate Elliptical Distribution in Cai, Einmahl,de Haan (2011)
    u=rand(n,1);
    r=((1-u).^(-2)-1).^(1/6);
    theta_r=rand(n,1)*2*pi;
    x1=2*r.*cos(theta_r);
    x2=r.*sin(theta_r);
    X=[x1 x2];
    datatype='Ellipse';
end

