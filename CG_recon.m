function [xout,count]=CG_recon(x0,v,sMaps,U1,kU,opt)
% this is a modified l2 constrained Conjugate Gradient SENSE to reconstruct image diff
% input 
%       x0: intial value
%        v: l2 norm term
%    sMaps: Sensitivity Map
%       U1: undersampling matrix(every coil)
%       kU: the undersampled data
%      opt: CG parameters
% output 
%     xout: reconstructed image diff
%    count: number of iterations

% 2013.04.22 Yi Guo

% 04/25/2016
% To test linear version of CG SENSE step in splitting model-based recon

[kx,ky,kz,nt,ncoil]=size(kU);

sMaps=repmat(sMaps,[1 1 1 nt 1]);

b=sum(conj(sMaps).*ifft2(kU-U1.*fft2(sMaps.*repmat(opt.Sb,[1 1 1 1 ncoil]))),5)+opt.beta1*v; 
% b=AHb=SHFHU(y-UFS*Sb)+beta*v;

count=1; % set initial parameters
xk=x0;
b=b(:);
rk=b-compQ(x0);
pk=rk;
%rk1=rk;

while count<opt.itermax && norm(rk)>opt.tol % CG iteration
    alpha=(rk'*rk)/(pk'*compQ(pk)); % modified CG equations
    xk=xk+alpha*pk;
    rk1=rk-alpha*compQ(pk);
    beta=(rk1'*rk1)/(rk'*rk);
    
    pk=rk1+beta*pk;
    count=count+1;
%    xD=reshape(xk1,kx,ky,kz,nt);       
   % imshow(real(xD(:,:,1,10)));drawnow;               
     rk=rk1; % only rk needs update
end

xout=xk;

    function Qx=compQ(x) % (A^H)Ax
        x=reshape(x,kx,ky,kz,nt);   
        Qx=sum(conj(sMaps).*ifft2(U1.*fft2(sMaps.*repmat(x,[1 1 1 1 ncoil]))),5);
        Qx=Qx+opt.beta1*x;
        Qx=Qx(:);      
    end
end
