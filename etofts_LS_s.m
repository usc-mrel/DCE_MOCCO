function [cost,grad]=etofts_LS_s(x,Cp,Ct_m,tModel)
% Least-square equation with extended-tofts
% spatial version with spatial dimension stretch to one vector array, that
% is (kx*ky*kz,1) or (kx*ky*kz,nt)
% Yi Guo, 06/12/2014

%nt=length(Cp); % number of temporal points
Cp=Cp(:);
tModel=tModel(:);

[Ns,nt]=size(Ct_m);

% [kx,ky,kz,nt]=size(Ct_m);
% Ns=kx*ky*kz;
% Ct_m=reshape(Ct_m,[Ns,nt]);
%% Calculate cost function
Ktrans=x(1:Ns);
Kep=x(Ns+1:2*Ns);
vp=x(2*Ns+1:3*Ns);

dtime=diff(tModel);
dtime=dtime(20);

Ct = model_extended_tofts_s(Ktrans, Kep, vp, Cp, tModel);

cost_in=Ct-Ct_m;
cost=0.5*sum(abs(cost_in(:)).^2);

%% calculate gradient of cost function


dCdKt=zeros(Ns,nt);
dCdKep=zeros(Ns,nt);


for k = 1:nt 
    
    % The time for T
    Tc = tModel(1:k);
    Tc=repmat(Tc',[Ns,1]);
    
    Cp1= Cp(1:k);
    Cp1=repmat(Cp1',[Ns,1]);
    
    Kep1=repmat(Kep,[1,k]);
        
    Fkt = Cp1.*exp(-Kep1.*(Tc(end)-Tc));
    Fkep= -Cp1.*(Tc(end)-Tc).*exp(-Kep1.*(Tc(end)-Tc));
  %      dCdKt(k) = trapz(Tc,Fkt);
  %      dCdKep(k)=trapz(Tc,Fkep);
  dCdKt(:,k)=sum(Fkt,2)*dtime;
  dCdKep(:,k)=sum(Fkep,2).*dtime.*Ktrans;

end

grad1=sum((dCdKt).*cost_in,2);  % hermitian transpose needs this conjugate
grad2=sum((dCdKep).*cost_in,2);
grad3=sum(repmat(Cp',[Ns,1]).*cost_in,2);

grad=[grad1(:);grad2(:);grad3(:)];


end

