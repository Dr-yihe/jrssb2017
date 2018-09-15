function datstruc = CmpErr(para,datstruc,datstrucX,datatype,ifAffine)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computation of Estimation Erros   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% para: the enviromental parameters
% datstruc: the data structure contatins estimates
% datstrucX: the data structure contains information of distribution
% ifAffine: if the data is a affine transformation of Cauchy distribution


% default setting
if nargin<4
    error('Not enough inputs')
elseif nargin<5
    ifAffine=false;
end
display('Computing estimation errors...')

%
% Initialization
%
datstrucX.rtrue(:,2:end)=[];
datstrucX.ptrue(2:end)=[];
datstrucX.betatrue(2:end)=[];

if ifAffine
    switch datatype
        case 'bicauchy'
            para.A=para.Ac;para.b=para.bc;
        case 'bistudent'
            para.A=para.At;para.b=para.bt;
    end
end

datstruc.rlogerr=zeros(size(datstruc.rest))*nan;
datstruc.betaest=zeros(size(datstruc.rest))*nan;
datstruc.betalogerr=zeros(size(datstruc.rest))*nan;
datstruc.pesterr=zeros(para.S,1)*nan;


%%
% the main procedures
%
sindex=1:para.S;

for i=sindex(datstruc.ifEff);
    
    % First compute the log relative error for estimated radius
    
    datstruc.rlogerr(:,i)=log(datstruc.rest(:,i))-log(datstrucX.rtrue);
    
   
    % Second, the log relative error at depth level
  if datstrucX.dim==2
      boundary_est=(datstruc.rest(:,i)*ones(1,2)).*...
            para.u;
  elseif datstrucX.dim==3
      boundary_est=(datstruc.rest(:,i)*ones(1,3)).*...
            para.u3D;
  end
    
    
    if ifAffine
        datstruc.betaest(:,i)=hddepth(boundary_est,datatype,para.A,para.b);
    else
        datstruc.betaest(:,i)=hddepth(boundary_est,datatype);
    end
    datstruc.betalogerr(:,i)=log(datstruc.betaest(:,i))-log(datstrucX.betatrue);
    
    % Third, the relative error at probability level
        switch datatype
            case 'bicauchy'
                f=@(r)((1-1./sqrt(r.^2+1))/(2*pi));
                if ~ifAffine
                    datstruc.pesterr(i)=mean(abs(f(datstruc.rest(2:end,i))-f(datstrucX.rtrue)))*(2*pi);
                else
                    % Transform back to standard bicauchy case
                    bestorigin=(boundary_est-ones(para.M,1)*para.b')/(para.A');
                    % calculate the corrsponding angle and radius
                    [angleorigin,restorigin]=cart2pol(bestorigin(:,1),bestorigin(:,2));
                    [~,sortix]=sort(angleorigin);
                    datstruc.pesterr(i)=sum(abs(f(restorigin(sortix(2:end)))-f(sqrt(1/(datstrucX.ptrue)^2-1)))...
                    .*diff(angleorigin(sortix)));

                end

            case 'bistudent'
                f=@(r)((1-1./((r.^2/3+1).^(3/2)))/(2*pi));
                if ~ifAffine
                    datstruc.pesterr(i)=mean(abs(f(datstruc.rest(2:end,i))-f(datstrucX.rtrue)))*(2*pi);
                    datstruc.pest(i)=mean(1-f(datstruc.rest(2:end,i))*(2*pi));
                else
                    bestorigin=(boundary_est-ones(para.M,1)*para.b')/(para.A');
                    [angleorigin,restorigin]=cart2pol(bestorigin(:,1),bestorigin(:,2));
                    [~,sortix]=sort(angleorigin);
                    datstruc.pesterr(i)=sum(abs(f(restorigin(sortix(2:end)))-f(sqrt(3*((datstrucX.ptrue)^(-2/3)-1))))...
                    .*diff(angleorigin(sortix)));
                end

            case 'biellp'
                bestorigin=boundary_est*[1/2 0;0 1];
                [angleorigin,restorigin]=cart2pol(bestorigin(:,1),bestorigin(:,2));
                f=@(r)((1-1./sqrt(r.^6+1))/(2*pi));
                [~,sortix]=sort(angleorigin);
                datstruc.pesterr(i)=sum(abs(f(restorigin(sortix(2:end)))-f((datstrucX.ptrue.^(-2)-1).^(1/6)))...
                    .*diff(angleorigin(sortix)));
            case 'tricauchy'
                f=@(r)((atan(r)-r./(r.^2+1))/(2*pi^2));
                firstint=zeros(para.M3D-1,1);
                for j=1:para.M3D-1
                    firstint(j)=mean(abs(f(datstruc.rest(((j-1)*para.M3D+2):(j*para.M3D),i))...
                        -f(datstrucX.rtrue)))*2*pi;
                end
                datstruc.pesterr(i)=mean(firstint.*sin(para.angle2(1:para.M3D-1,1)))*pi;
        end
    
end

end

