function C = sig2conc2D(imgD,R10,M0,alpha,TR)

% modified function that takes the difference of images as input
% 04/25/2016

% equation to calculate concentration from image intensity

% R1(t)=-1/TR*ln(1-((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*cos(alpha)))
%over 1-cos(alpha)*((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*sin(alpha)))

% where m=exp(-R10*TR)
% then C(t)=(R1(t)-R1(0))/r1

% Yi Guo, 06/12/2014
% some simplification, pay attention to R1, and alpha unit!!!

Rcs=4.39;
nt=size(imgD,4);

m=exp(-repmat(R10,[1 1 1 nt])*TR);
par2=(1-m)./(1-m*cos(alpha));

par1=(imgD)./repmat(M0+eps,[1 1 1 nt])/sin(alpha);

B=(par1+par2);

Rt=-1/TR*real(log((1-B)./(1-cos(alpha)*B+eps)));
% Rt=-1/TR*(log((1-B)./(1-cos(alpha)*B+eps)));

R1B=Rt(:,:,:,1); % baseline R1 is equal to R10 if img(1)==imgB
%R1B=R10;
C=Rt-repmat(R1B,[1 1 1 nt]);

C=C/Rcs;

C(C<0) = 0;
C(C>180)=0;

end

    
    



