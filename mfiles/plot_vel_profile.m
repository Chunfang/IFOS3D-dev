% plots vertical velocity profile of real, inverted and starting model
% S. Butzer 2013

clear all;
%close all;

nx=320; ny=160; nz=320;

cuty1=1;
cuty2=320;
outx=1; outy=1; outz=1; 
dh=0.8;
nx=nx/outx;ny=ny/outy;nz=nz/outz;
fignum=44;

file_inp1='/data14/sdunkl/surface/results12/model/surface11_real.vs';
file_inp2='/data14/sdunkl/surface/results12/model/surface11_start.vs';
file_inp3='/data14/sdunkl/surface/results12/model/surface12.vs_it170';

px=160; %location of profile
pz=200;

%--------------------------------------------------------------------------

nxm=nx*dh+dh; nym=ny*dh+dh; nzm=nz*dh+dh;
xp1=dh; xp2=nx*dh; yp1=dh; yp2=ny*dh; zp1=dh; zp2=nz*dh;
x=xp1:dh*outx:xp2*outx;
y=yp1:dh*outy:yp2*outy;
z=zp1:dh*outz:zp2*outz;

fid_rot=fopen(file_inp1,'r','ieee-le');
rot1=zeros(nz/outz,nx/outx,ny/outy);
rot1=fread(fid_rot,(nx*ny*nz),'float');

MEAN=mean(rot1)
%rot1=rot1-div1;
%rot1=1./((rot1+5e-7));
    %min(min(rot))
    %max(max(rot1))
%rot1=rot1./max(max(rot1));
%
%norm222=norm(rot1)
%rot1=rot1.*div1;
%rot1=rot1./max(max(abs(rot1)));
    %rot1=log10(rot1);
MAX=max(max(rot1))
MIN=min(min(rot1))
rot=reshape(rot1,nz/outz,nx/outx,ny/outy);

depth=rot(pz,px,:);
depth1=depth(:);
A=size(depth1)
%--------------------------------------------------------------------------
fid_rot=fopen(file_inp2,'r','ieee-le');
rot1=zeros(nz/outz,nx/outx,ny/outy);
rot1=fread(fid_rot,(nx*ny*nz),'float');

MEAN=mean(rot1)
%rot1=rot1-div1;
%rot1=1./((rot1+5e-7));
    %min(min(rot))
    %max(max(rot1))
%rot1=rot1./max(max(rot1));
%
%norm222=norm(rot1)
%rot1=rot1.*div1;
%rot1=rot1./max(max(abs(rot1)));
    %rot1=log10(rot1);
MAX=max(max(rot1))
MIN=min(min(rot1))
rot=reshape(rot1,nz/outz,nx/outx,ny/outy);

depth=rot(pz,px,:);
depth2=depth(:);
A=size(depth2)
%--------------------------------------------------------------------------
fid_rot=fopen(file_inp3,'r','ieee-le');
rot1=zeros(nz/outz,nx/outx,ny/outy);
rot1=fread(fid_rot,(nx*ny*nz),'float');

MEAN=mean(rot1)
%rot1=rot1-div1;
%rot1=1./((rot1+5e-7));
    %min(min(rot))
    %max(max(rot1))
%rot1=rot1./max(max(rot1));
%
%norm222=norm(rot1)
%rot1=rot1.*div1;
%rot1=rot1./max(max(abs(rot1)));
    %rot1=log10(rot1);
MAX=max(max(rot1))
MIN=min(min(rot1))
rot=reshape(rot1,nz/outz,nx/outx,ny/outy);

depth=rot(pz,px,:);
depth3=depth(:);
A=size(depth3)


figure(fignum)
plot(depth1,y,'k-');hold on;
plot(depth2,y,'b-');hold on;
plot(depth3,y,'r-');

 xlabel('velocity in m/s');
    ylabel('depth in m');
          
    set(gca, 'xDir','normal')      
    set(gca, 'yDir','reverse')         
       
    set(get(gca,'Ylabel'),'FontSize',14);
   set(get(gca,'Ylabel'),'FontWeight','normal');
   set(get(gca,'Xlabel'),'FontSize',14);
   set(get(gca,'Xlabel'),'FontWeight','normal');
   set(get(gca,'Zlabel'),'FontSize',14);
   set(get(gca,'Zlabel'),'FontWeight','normal');
    set(gca,'FontSize',14);
    set(gca,'FontWeight','normal');
    set(gca,'Linewidth',1.0);
    
    xlim([350 1000]);
    ylim([1*dh 140*dh]);
    
