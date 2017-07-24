function msekt=multi_disp_e(kt,kep,vp,intern)

if nargin==3
    intern=0.1;
end

nd=size(kt,3);
% kt=real(kt);
% kep=real(kep);
% vp=real(vp);

kt=abs(kt);
kep=abs(kep);
vp=abs(vp);

x1=1;y1=1;x2=256;y2=150;
% x1=46; y1=95; x2=80; y2=120;
% x1=155; y1=76; x2=202; y2=108;

 x1=100;x2=179;y1=52;y2=118; %phantom


ktd=kt(x1:x2,y1:y2,1);
kepd=kep(x1:x2,y1:y2,1);
vpd=vp(x1:x2,y1:y2,1);




for i=2:nd
    ktd=[ktd,kt(x1:x2,y1:y2,i)];
    kepd=[kepd,kep(x1:x2,y1:y2,i)];
    vpd=[vpd,vp(x1:x2,y1:y2,i)];
end

imagesc(cat(1,ktd,kepd,vpd),[0 intern]); axis image; axis off; colorbar; colormap jet;
% title(['MSE=',num2str(msekt)]);

end
