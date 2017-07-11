%% load data
close all;
clc; clear all;
addpath(genpath('minFunc_2012'));
% direct reconstruction with model consistency constraint
% demonstration of joint estimation of AIF and TK (Patlak model)

load ../Research_coding/Direct_recon/Direct_save_2/DCE50_0421.mat;
%please download demo in-vivo data from: 
%https://drive.google.com/file/d/0B4nLrDuviSiWXzJhLWFwN1c1ZG8/view?usp=sharing

ns=1; % take one slice of k-space for speed
k=k(:,:,ns,:,:);
opt.size=size(k);
[kx,ky,kz,nt,ncoil]=size(k);

sMaps=sMaps(:,:,ns,:,:);
sMaps=reshape(sMaps,[opt.size(1) opt.size(2) 1 1 ncoil]);

% coil combine to get the fully-sampled images
imgF=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(k),5);
opt.Sb=repmat(imgF(:,:,:,1),[1 1 1 nt]); % first frame is fully-sampled

%% set parameters
delay=8; % delay frames for contrast injection
tpres=5/60; % temporal resolution, unit in seconds!
opt.time=[zeros(1,delay),[1:(nt-delay)]*tpres];
opt.plot=1;

alpha=pi*15/180;
TR=0.006;
M0=5*ones(opt.size(1),opt.size(2),'single'); %use simulated M0, R1
R1=1*ones(opt.size(1),opt.size(2),'single');

opt.itermax=10; % options for CG iteration for solving image diff
opt.tol=0.01;
opt.beta1=0.1;

%% AIF ROI selection
CONCF=(imgF-opt.Sb);
CONCF=sig2conc2D(CONCF,R1,M0,alpha,TR); % fully-sampled concentration

AROIx=162:163; %0421 ROI selection
AROIy=119:121;

Cp0=squeeze(abs(CONCF(AROIx,AROIy,1,:))); % get AIF from the ROI
Cp0=reshape(Cp0,[numel(Cp0)/nt nt]);

Cp0=mean(Cp0,1);
hct=0.4;
Cp0 = Cp0(:)./(1-hct); % correct for hematocrit

[Ktf,Vpf]=conc2Ktrans_Y(CONCF,opt.time,Cp0); % get fully-sampled TK maps

%% undersamping by RGR
R = 20; %under-sampling rate
[~,U11] = genRGA(220.0, 220.0, kx, ky, round(kx*ky/R*nt), bin2dec('1111'), 0, nt, 0.3, 0.1, 10.0, 2500, 1, 1, 10.0);
U1=reshape(nshift(U11>0,[1 2]),[kx,ky,1,nt]);
U1=repmat(U1,[1 1 1 1 ncoil]);
U1(:,:,:,1,:)=1;  % first frame is fully-sampled
kU=k.*U1;

%% splitting algorithm
e1_k=zeros(opt.size(1),opt.size(2),opt.size(3),opt.size(4));
v_k=e1_k;
w_k=v_k;

iter=1;

while iter<51
    % CG recon for image diff, warm start with previous result
    [w_k,count]=CG_recon(w_k(:),w_k,sMaps,U1,kU,opt);    
    w_k=reshape(w_k,[kx,ky,kz,nt]);
    
    C_t=sig2conc2D(w_k,R1,M0,alpha,TR);
    
    % get AIF from ROI
    Cp_1=squeeze(abs(C_t(AROIx,AROIy,1,:)));
    Cp_1=reshape(Cp_1,[numel(Cp_1)/nt nt]);
    Cp_1=mean(Cp_1,1);
    Cp_1 = Cp_1(:)./(1-hct); % correct for hematocrit

    % get TK maps using patient-specific AIF
    [Kt,Vp]=conc2Ktrans_Y(C_t,opt.time,Cp_1);
    
    % forward modeling to get back concentration from TK maps
    C_t=Ktrans2conc_Y(Kt,Vp,opt.time,Cp_1);
    
    % get back image diff 
    w_k=conc2sigD(C_t,R1,M0,alpha,TR);
            
    if opt.plot % show intermediate result
        subplot(1,2,1);
        imagesc(real(cat(2,Kt,Vp)),[0 0.1]); axis image; axis off; title(num2str(iter));
        subplot(1,2,2);plot(1:nt,Cp0,1:nt,Cp_1);
        drawnow;
    end
    
    iter=iter+1;
end

%% display
figure;
subplot(1,2,1)
imagesc(real(cat(2,Ktf,Kt)),[0 0.1]); axis image; axis off; title('Fully/Recon Ktrans');
subplot(1,2,2);plot(1:nt,Cp0,1:nt,Cp_1); title('Population/Patient AIF');


