%% load data
close all;
clc; clear all;
addpath(genpath('minFunc_2012'));
% direct reconstruction with model consistency constraint
% demo using digital referecen object (DRO), with the option to select 
% in-house gradient-based solver, or third-party Rocketship solver
% Both phantom and the solver are based on e-Tofts model

load ../Research_coding/Direct_split_phantom/phantom_data; % load phantom
% Phantom is generated based on R.J Bosca et al. Phys. Med. Biol, 2016 & Y Bliesener et al. ISMRM 2017, p1909
% please download demo phantom data from: 
% https://drive.google.com/file/d/0B4nLrDuviSiWT3ZKUmd0YjRwUEU/view?usp=sharing 

solver='G'; 
% 'R' to select Rocketship solver, 
% 'G' to select Gradient-based l-bfgs solver

ns=1; % take one slice of k-space for speed
k=k(:,:,ns,:,:);
opt.size=size(k);
[kx,ky,kz,nt,ncoil]=size(k);

sMaps=sMaps(:,:,ns,:,:);
sMaps=reshape(sMaps,[kx ky 1 1 ncoil]);

% coil combine to get the fully-sampled images
imgF=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(k),5);
opt.Sb=repmat(imgF(:,:,:,1),[1 1 1 nt]); % first frame is fully-sampled

%% set parameters
delay=8; % delay frames for contrast injection
tpres=5/60; % temporal resolution, unit in seconds!
opt.time=[zeros(1,delay),[1:(nt-delay)]*tpres];
opt.AIF=SAIF_p(opt.time); % get population-averaged AIF
hct=0.4;
opt.AIF = opt.AIF(:)./(1-hct); % correct for hematocrit

alpha=pi*15/180;
TR=0.006;
M0=5*ones(kx,ky,'single'); % use simulated M0, R1
R1=1*ones(kx,ky,'single');

options.MaxIter=10; % options for l-bfgs iteration for solving etofts maps
options.display = 'off';
options.Method ='lbfgs';
options.useMex=0;

opt.itermax=10; % options for CG iteration for solving image diff
opt.tol=0.001;
opt.beta1=0.1;

%% undersamping by RGR
R =40; %under-sampling rate
[~, U11] = genRGA(220, 220, kx,ky, round(kx*ky/R*nt), bin2dec('0001'), 0, nt);
U1=reshape(nshift(U11>0,[1 2]),[kx,ky,1,nt]);
U1=repmat(U1,[1 1 1 1 ncoil]);
U1(:,:,:,1,:)=1; % fully sample the first frame
kU=k.*U1;

%% splitting algorithm
e1_k=zeros(kx,ky,kz,nt);
v_k=e1_k;
w_k=e1_k;

Ns=sum(mask1(:)); % number of actual pixels in side the tumor ROI
res1=zeros(Ns*3,1);

iter=1;

while iter<101    
    % CG recon for image diff, warm start with previous result
    [w_k,count]=CG_recon(w_k(:),w_k,sMaps,U1,kU,opt); 
    w_k=reshape(w_k,[kx,ky,kz,nt]);
    
    % back-ward modelling from image difference
    w_k=sig2conc2D(w_k,R1,M0,alpha,TR);    
    
    % perform model fitting only on the known masked part
    CONC1=w_k(repmat(mask1,[1 1 1 nt]));
    CONC1=reshape(CONC1,[Ns,nt]);
    
    % get eTofts maps using selected solver
    if solver == 'G'   
        % warm start with previous result 
        [res1,f1,exitflag1,output1]=minFunc(@etofts_LS_s,res1,options,opt.AIF,CONC1,opt.time);
        KtS=res1(1:Ns); KepS=res1(Ns+1:2*Ns); VpS=res1(2*Ns+1:3*Ns);
    elseif solver == 'R' 
        % for Rocketship solver, please have the source codes ready 
        % from https://github.com/petmri/ROCKETSHIP, and add to path
        roi_data{1}.Cp = double(opt.AIF');
        roi_data{1}.timer = double(opt.time');
        roi_data{1}.Ct = double(permute(CONC1,[2 1]));
        [roi_results, roi_residuals] = FXLfit_generic(roi_data, Ns, 'ex_tofts', 0);
        KtS = roi_results(:,1); KepS = roi_results(:,2); VpS = roi_results(:,3);
    else
        disp('Unrecognized solver');
        break;
    end
    
    % get back the images
    Kt=zeros(kx,ky);Kep=Kt; Vp=Kt;
    Kt(mask1)=KtS;Kep(mask1)=KepS;Vp(mask1)=VpS;

    % forward modelling from TK parameters to image diff
    w_k = model_extended_tofts_s(Kt, Kep, Vp, opt.AIF, opt.time);
    w_k=reshape(w_k,[kx,ky,kz,nt]);
    w_k=conc2sigD(w_k,R1,M0,alpha,TR);
    
    % show intermediate result
    imagesc(real(cat(2,Kt,Kep/2,Vp)),[0 0.1]); title(['Iter: ',num2str(iter)]);
    colormap jet; axis off; axis image; drawnow;
    
    iter=iter+1;
end
 
%% display
% show original(true) and reconstructed etofts maps
multi_disp_e(cat(3,KtZ,Kt),cat(3,KepZ,Kep)/2,cat(3,VpZ,Vp),0.1)


