function [CONC,A_out,AIF_out]=Ktrans2conc_Y(Kt,Vp,time,C_AIF)
% forward calculation of contrast consentration from Ktrans and Vp
% Yi Guo 04/30/2014

% Input C_AIF1 which is already converted from signal to contrast
% concentration for patient-specific AIF. Do not do any scalling before
% input

% get image size

if nargin==3
    C_AIF=SAIF_p(time);
end

% this hematocrit correction should always be outside to be
% consistent with etofts estimation
%  hct = 0.4; % move this hematocrit correction back here for consistency
% C_AIF = C_AIF(:)./(1-hct);   

if mean(diff(time))>0.5
time=time/60; % convert to minute unit
end
 
[nx,ny,nz]=size(Kt);
nt=length(time);

%   Precompute some variables
% ntind=[1:nt];
% %ntind = find(time>tlimits(1) & time<tlimits(2));
% ntind = ntind(ntind > find(C_AIF == max(C_AIF)));
% dtime = [0; diff(time(:))];
% dtime = dtime(ntind);
% C_aif = C_AIF(ntind);

% commented above is just to save time
dtime=diff(time);
dtime=dtime(20); % just get one value for temporal resolution
%dtime=[0;dtime'];

A=cumsum(C_AIF).*dtime;

%Kt=Kt/60; % do not do this conversion since time is already in minutes!

% for it=1:nt
% %CONC1=Kt.*A(count)+Vp;
% %CONC(:,:,:,it)=CONC1*C_aif(count);
% 
% %CONC(:,:,:,it)=Kt.*A(it)+Vp*C_AIF(it);
% 
% count=count+1;
% end

A_out=A;
AIF_out=C_AIF;
A=reshape(A,[1 1 1 nt]);
C_AIF=reshape(C_AIF,[1 1 1 nt]);

CONC=repmat(Kt,[1 1 1 nt]).*repmat(A,[nx,ny,nz,1])+repmat(Vp,[1 1 1 nt]).*repmat(C_AIF,[nx,ny,nz,1]);
% this avoids loop

end


