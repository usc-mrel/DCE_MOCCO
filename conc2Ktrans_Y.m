function [Kt,Vp,C_AIF]=conc2Ktrans_Y(CONC,time,C_AIF)
% forward model to calculate Ktrans from contrst concentration
% linear fitting using equation:
% Ct(t)=Kt*integral(0-t){Cp(tao)dtau}+vp*Cp(t)
% Yi Guo, 07/2014

% Simplify and verify to do only Patlak modelling, time unit in minutes
% Yi Guo, 10/2015

% Simplify to avoid loop, get rid of the mask
% Yi Guo, 04/2016

[np, nv, ns, nt] = size(CONC);
nL=np*nv*ns;

if nargin==2
    C_AIF=SAIF_p(time); % /3 moved to function
end
% this hematocrit correction should always be outside to be
% consistent with etofts estimation
% hct = 0.4;  
% C_AIF = C_AIF(:)./(1-hct);  

C_AIF=C_AIF(:);

if mean(diff(time))>0.5
time=time/60; % convert to minute unit
end

dtime=diff(time);
dtime=dtime(20); % just get one value for temporal resolution

[~,peak]=max(C_AIF);  %get rid of the part before peak
peak=peak-1;
C_aif=C_AIF(peak:nt);

nt1=nt-peak+1;

C_aifI = 1./C_aif;

A = [cumsum(C_aif).*dtime.*C_aifI ones(nt1,1)];

%   Compute TK parameters

CONC=CONC(:,:,:,peak:nt); % only use later points

CONC=reshape(CONC,[nL,nt1]); 
CONC=permute(CONC,[2 1]); % put nt on Column
CONC=CONC.*repmat(C_aifI,[1 nL]);  % divide by Cp(t), that is Ct(t)/Cp(t)

res=pinv(A)*CONC; % in this formulation, 2*nL matrix is returned, first row is ktrans, second row is Vp

Kt=reshape(res(1,:),[np nv ns]);
Vp=reshape(res(2,:),[np nv ns]);


%   Unit conversion [Kt] = 1/min 
% Time unit in minutes, will keep this consistent
% This way the Ktrans itself is intrincically small, make it possible to
% jointly estimate Ktrans and Vp together in direct recon
% Kt = 60*Kt; % No conversion if time is in minute

end

