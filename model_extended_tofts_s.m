%% FXLStep1AAIF, vp
function Ct = model_extended_tofts_s(Ktrans, Kep, vp, Cp, tModel)

% spatial image version of extended-tofts modelling
% based off Sam Burnes' codes

% change input to Ktrans, Kep, Vp
% input is a vector array of length kx*ky*kz, that is, stretch the spatial
% coordinate
% output is Ct of dimension (kx*ky*kz, nt)
% Yi Guo, 10/15/2015

% make sure every input is a long vector array
Cp = Cp(:);
tModel = tModel(:);

Ktrans=Ktrans(:);
Kep=Kep(:);
vp=vp(:);

dtime=diff(tModel);
dtime=dtime(20);

Ns=length(Ktrans);
nt=length(Cp);
% Pre-alocate for speed
Ct = zeros(Ns,nt);


for k = 1:nt
    
    % The time for T
    T = tModel(1:k);
    T=repmat(T',[Ns,1]);
    
    Cp1= Cp(1:k);
    Cp1=repmat(Cp1',[Ns,1]);
    
    Kep1=repmat(Kep,[1 k]);

    F = Cp1.*exp((-Kep1).*(T(end)-T));
    
    M=sum(F,2)*dtime; % use sum instead of trapz
    
    Ct(:,k) = Ktrans.*M+vp.*Cp(k);
end


end
