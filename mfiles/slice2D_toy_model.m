%3D plot of model

clear all;
clc
close all

fname='toy.';
parameter='vp';
source_file='../par/sources/sources_toy.dat';
rec_file='../par/receiver/receiver_toy.dat';

iteration=60;

nx=160; ny=184; nz=160; %ny:vertical
outx=1; outy=1; outz=1; 
dh=0.8;
FW=10;
nx=nx/outx;ny=ny/outy;nz=nz/outz;

X=dh:dh*outx:nx*outx*dh;
Y=dh:dh*outy:ny*outy*dh;
Z=dh:dh*outz:nz*outz*dh;

% Slizes of the 3D-Model in m 
slx=56;
sly=24;
slz=52.8;

if strcmp(parameter,'vs')
caxis_value_1=3300;%vs
caxis_value_2=3900;%vs
end
if strcmp(parameter,'vp')
caxis_value_1=5700;%vp
caxis_value_2=6700;%vp
end


% No changes to be done after this line---------
%------------------------------------------------

%Input of receiver and source positions--

fp = fopen(source_file);
SRC=textscan(fp, '%f %f %f %f %f %f', 'delimiter','\t', 'commentstyle','%','MultipleDelimsAsOne',1);
fclose(fp);

fp = fopen(rec_file);
REC=textscan(fp, '%f %f %f ', 'delimiter','\t', 'commentstyle','#','MultipleDelimsAsOne',1);
fclose(fp);

% Lines to mark PML position
fontsize=14;
x_line=[FW FW nx-FW nx-FW FW]*dh;
y_line=[FW ny-FW ny-FW FW FW]*dh;
z_line=[FW nz-FW nz-FW FW FW]*dh;


%m --> model indices
slizex=int64(slx/dh/outx);
slizey=int64(sly/dh/outy);
slizez=int64(slz/dh/outz);


file_inp1=['../par/model/' fname parameter '.true'];
file_inp2=['../par/model/' fname parameter '_it' num2str(iteration)];




%Read 3D Models------------------------

fid=fopen(file_inp1,'r','ieee-le');
modelvec=zeros(ny/outy,nx/outx,nz/outz);
modelvec=fread(fid,(nx*ny*nz),'float');

model_true=reshape(modelvec,ny,nx,nz);

fid=fopen(file_inp2,'r','ieee-le');
modelvec=zeros(ny/outy,nx/outx,nz/outz);
modelvec=fread(fid,(nx*ny*nz),'float');

model=reshape(modelvec,ny,nx,nz);


%Creation of 2D Slices--------------

model_truex=zeros(ny,nz);
model_truey=zeros(nz,nx);
model_truez=zeros(ny,nx);

for y=1:1:ny
   for z=1:1:nz
    model_truex(y,z)=model_true(y,slizex,z);
   end
end

for x=1:1:nx
   for z=1:1:nz  
    model_truey(z,x)=model_true(slizey,x,z);
   end
end

for x=1:1:nx
   for y=1:1:ny  
    model_truez(y,x)=model_true(y,x,slizez);
   end
end


modelx=zeros(ny,nz);
modely=zeros(nz,nx);
modelz=zeros(ny,nx);

for y=1:1:ny
   for z=1:1:nz
    modelx(y,z)=model(y,slizex,z);
   end
end

for x=1:1:nx
   for z=1:1:nz  
    modely(z,x)=model(slizey,x,z);
   end
end

for x=1:1:nx
   for y=1:1:ny  
    modelz(y,x)=model(y,x,slizez);
   end
end

%Plotting-------------------------------------------

figure(1)
imagesc(Z,Y,model_truex);
line(x_line,y_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize);
colormap('jet');
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
%set(gca,'xdir','reverse');
xlabel('z in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['true ' parameter '-model. Slize at x=' num2str(slx) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{1}(ii)==slx;
        plot(SRC{3}(ii),SRC{2}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{1}(ii)==slx;
        plot(REC{3}(ii),REC{2}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off

figure(2)
imagesc(Z,Y,modelx);
line(x_line,y_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize); 
colormap('jet');
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
%set(gca,'xdir','reverse');
xlabel('z in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title([ parameter '-model. Iteration ' num2str(iteration) ' Slize at x=' num2str(slx) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{1}(ii)==slx;
        plot(SRC{3}(ii),SRC{2}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{1}(ii)==slx;
        plot(REC{3}(ii),REC{2}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off

figure(3)
imagesc(X,Z,model_truey);
line(x_line,z_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize); 
colormap('jet');
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
xlabel('x in m','FontSize',fontsize)
ylabel('z in m','FontSize',fontsize)
title(['true ' parameter '-model. Slize at x=' num2str(sly) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{2}(ii)==sly;
        plot(SRC{1}(ii),SRC{3}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{2}(ii)==sly;
        plot(REC{1}(ii),REC{3}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off


figure(4)
imagesc(X,Z,modely);
line(x_line,z_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize); 
colormap('jet');
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
xlabel('x in m','FontSize',fontsize)
ylabel('z in m','FontSize',fontsize)
title([parameter '-model. Iteration ' num2str(iteration) ' Slize at x=' num2str(sly) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{2}(ii)==sly;
        plot(SRC{1}(ii),SRC{3}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{2}(ii)==sly;
        plot(REC{1}(ii),REC{3}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off

figure(5)
imagesc(X,Y,model_truez);
line(x_line,y_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize); 
colormap('jet');
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
%set(gca,'xdir','reverse');
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title(['true ' parameter '-model. Slize at x=' num2str(slz) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{3}(ii)==slz;
        plot(SRC{1}(ii),SRC{2}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{3}(ii)==sly;
        plot(REC{1}(ii),REC{2}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off

figure(6)
imagesc(X,Y,modelz);
line(x_line,y_line,'LineStyle','--','Color','k')
colb=colorbar;
coll=get(colb,'xlabel');
colormap('jet');
set(coll,'String',[parameter ' in m/s'],'FontSize',fontsize); 
caxis([caxis_value_1 caxis_value_2])
set(gca,'ydir','normal');
%set(gca,'xdir','reverse');
xlabel('x in m','FontSize',fontsize)
ylabel('y in m','FontSize',fontsize)
title([parameter '-model. Iteration ' num2str(iteration) ' Slize at x=' num2str(slz) 'm'],'FontSize',fontsize)
axis tight
hold on
%Plot sources if in model slice
for ii=1:length(SRC{1})
    if SRC{3}(ii)==slz;
        plot(SRC{1}(ii),SRC{2}(ii),'*blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
%Plot receiver if in model slice
for ii=1:length(REC{1})
    if REC{3}(ii)==sly;
        plot(REC{1}(ii),REC{2}(ii),'.blue','MarkerSize',15,'LineWidth',2.0);
    end
    
end
hold off


 %exportfig(2, [fname parameter '.eps'],'bounds','tight', 'color','rgb', ...
 % 'preview','none', 'resolution',200, 'lockaxes',1);
