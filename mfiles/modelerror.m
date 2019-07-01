%calculates the error between the true and inverted model

clear all;

nx=300; ny=350; nz=200;

cuty1=40;
cuty2=290;
outx=1; outy=1; outz=1; 
dh=0.8;
nx=nx/outx;ny=ny/outy;nz=nz/outz;

file_inp1='/workag13/FWI/random_KTB/Model/random7.vs';
file_inp2='/workag13/FWI/random_KTB/2D_10sources_vert/model/random72D_vert.vs_it73';


fid_rot=fopen(file_inp1,'r','ieee-le');
rot1=zeros(1,(nx*ny*nz));
rot1=fread(fid_rot,(nx*ny*nz),'float');
%rot=reshape(rot1,(nx*ny*nz));
size(rot1)

fid_div=fopen(file_inp2,'r','ieee-le');
div=zeros(nz/outz,nx/outx,ny/outy);
div1=fread(fid_div,(nx*ny*nz),'float');
div=reshape(div1,nz/outz,nx/outx,ny/outy);


error=0.0;
for i=1:(nx*ny*nz)
               error=error+abs((rot1(i)-div1(i))/rot1(i));
 end
error=(error/(nx*ny*nz))

