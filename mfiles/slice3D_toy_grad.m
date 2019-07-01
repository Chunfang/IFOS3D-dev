%3D plot of gradient

clear all;

nx=160; ny=186; nz=160; %ny:vertical
outx=1; outy=1; outz=1; 
dh=0.8;
nx=nx/outx;ny=ny/outy;nz=nz/outz;
fignum=71;


file_inp1='../par/grad/toy_grad.vs_160.00Hz_it1'; %preconditioned gradient: it1001
file_inp2='../par/grad/toy_grad.vp_160.00Hz_it1';    %"raw" gradient: it1



phi1=0; %rotation angles to x-z plane of first and second plane 
phi2=90;
rotaxes=[1,0,0]; %direction rotation axes [0,1,0] rotation of plane around vertical axis
                    % [1,0,0] rotation of plane around x-axis

                    %rotpoint=[120.8 140.8 80.8] ;       %turning point [x y z] in meter
                    rotpoint=[(nx/2)*dh (70)*dh (nz/2-13)*dh] ;
                    %rotpoint=[110 145 95] ;
                    rotpoint2=[110 130 70] ;
%viewpoint=[122,26];
viewpoint=[-20,12];
  
% caxis_value_1=-1; % preconditioned
% caxis_value_2=1; % preconditioned
% caxis_value_1 = -0.02; % raw, vp
% caxis_value_2 = 0.02; % raw, vp
caxis_value_1 = -0.03; % raw, vs
caxis_value_2 = 0.03; % raw, vs

xcut1=10; xcut2=149;
ycut1=10; ycut2=139;
zcut1=10; zcut2=149;


fid_rot=fopen(file_inp1,'r','ieee-le');
rot1=zeros(ny/outy,nx/outx,nz/outz);
rot1=fread(fid_rot,(nx*ny*nz),'float');

fid_div=fopen(file_inp2,'r','ieee-le');
div1=zeros(ny/outy,nx/outx,nz/outz);
div1=fread(fid_div,(nx*ny*nz),'float');


rot1=-rot1./max(max(rot1));
%rot1=rot1;
%rot1=1./((rot1+3e-10));
%rot1=rot1.*div1*50000;
%rot1=-rot1./max(max(abs(rot1)));
%rot1=log10(rot1);


rot=reshape(rot1,ny/outy,nx/outx,nz/outz);
rot=rot(ycut1:ycut2,xcut1:xcut2,zcut1:zcut2);

nx=xcut2-xcut1+1;
ny=ycut2-ycut1+1;
nz=zcut2-zcut1+1;

%size(rot)

xp1=xcut1*dh; xp2=xcut2*dh; yp1=ycut1*dh; yp2=ycut2*dh; zp1=zcut1*dh; zp2=zcut2*dh;
x=xp1:dh*outx:xp2*outx;
y=yp1:dh*outy:yp2*outy;
z=zp1:dh*outz:zp2*outz;



%x=xcut1:dh*outx:xcut2*outx;
%y=ycut1:dh*outy:ycut2*outy;
%z=zcut1:dh*outz:zcut2*outz;


% figure(12)
[Z,X,Y]=meshgrid(z,x,y);
xmin = min(X(:)); ymin = min(Y(:)); zmin = min(Z(:));
xmax = max(X(:)); ymax = max(Y(:)); zmax = max(Z(:));

%zplane=zeros(nx,nz);
%zplane(:,:)=rotpoint2(2);
%hslicez = surf(z,x,zplane);
%rotate(hslicez,rotaxes,phi1,rotpoint);
%xd1 = get(hslicez,'XData');
%yd1 = get(hslicez,'YData');
%zd1 = get(hslicez,'ZData');

%zplane=zeros(ny,nx);
%zplane(:,:)=rotpoint(3);
%hslicez = surf(x,y,zplane);
%rotate(hslicez,rotaxes,phi1,rotpoint);
%xd = get(hslicez,'XData');
%yd = get(hslicez,'YData');
%zd = get(hslicez,'ZData');

xd2 = rotpoint(1);
yd2 = rotpoint(2);
zd2 = rotpoint(3);

%zplane=zeros(nz,ny);
%zplane(:,:)=rotpoint(3);
%hslicez = surf(y,z,zplane);
%rotate(hslicez,rotaxes,phi2,rotpoint);
%xd2 = get(hslicez,'XData');
%yd2 = get(hslicez,'YData');
%zd2 = get(hslicez,'ZData');

%zplane=zeros(ny,nx);
%zplane(:,:)=rotpoint(3);
%hslicez = surf(x,y,zplane);
%rotate(hslicez,rotaxes,phi2,rotpoint);
%xd2 = get(hslicez,'XData');
%yd2 = get(hslicez,'YData');
%zd2 = get(hslicez,'ZData');

zplane=zeros(ny,nx);
zplane(:,:)=rotpoint(3);
hslicez = surf(x,y,zplane);
rotate(hslicez,[0 1 0],phi2,rotpoint);
xd3 = get(hslicez,'XData');
yd3 = get(hslicez,'YData');
zd3 = get(hslicez,'ZData');

%zplane=zeros(ny,nx);
%zplane(:,:)=rotpoint2(3);
%hslicez = surf(x,y,zplane);
%rotate(hslicez,[1 0 0],90,rotpoint2);
%xd4 = get(hslicez,'XData');
%yd4 = get(hslicez,'YData');
%zd4 = get(hslicez,'ZData');





rot=permute(rot,[2,3,1]);
%rott=size(rot)
%size(rot)
%size(Z)
figure(fignum)
%h3 = slice(Z,X,Y,rot,zd1,xd1,yd1);set(h3,'FaceColor','interp',...
%       'EdgeColor','none',...
%      'DiffuseStrength',.8) 
% alpha(.5)
%hold on
%h = slice(Z,X,Y,rot,zd,xd,yd);set(h,'FaceColor','interp',...
%        'EdgeColor','none',...
%        'DiffuseStrength',.8) 
%hold on
 h2 = slice(Z,X,Y,rot,zd2,xd2,yd2);set(h2,'FaceColor','interp',...
        'EdgeColor','none',...
        'DiffuseStrength',.8) 
   hold on
%h2 = slice(Z,X,Y,rot,zd3,xd3,yd3);set(h2,'FaceColor','interp',...
%        'EdgeColor','none',...
%       'DiffuseStrength',.8)  

    
hold on

  
             plot3([rotpoint(3) rotpoint(3)],[xcut1*dh xcut1*dh],[ycut1*dh ycut2*dh],'-black','LineWidth',0.5);
           plot3([rotpoint(3) rotpoint(3)],[xcut1*dh xcut2*dh],[ycut2*dh ycut2*dh],'-black','LineWidth',0.5);
           plot3([rotpoint(3) rotpoint(3)],[xcut2*dh xcut2*dh],[ycut1*dh ycut2*dh],'-black','LineWidth',0.5);
           plot3([rotpoint(3) rotpoint(3)],[xcut1*dh xcut2*dh],[ycut1*dh ycut1*dh],'-black','LineWidth',0.5);
     hold on             
            plot3([zcut1*dh zcut2*dh],[xcut1*dh xcut1*dh],[rotpoint(2) rotpoint(2)],'-black','LineWidth',0.5);
           plot3([zcut2*dh zcut2*dh],[xcut1*dh xcut2*dh],[rotpoint(2) rotpoint(2)],'-black','LineWidth',0.5);
            plot3([zcut1*dh zcut2*dh],[xcut2*dh xcut2*dh],[rotpoint(2) rotpoint(2)],'-black','LineWidth',0.5);
            plot3([zcut1*dh zcut1*dh],[xcut1*dh xcut2*dh],[rotpoint(2) rotpoint(2)],'-black','LineWidth',0.5);
       hold on      
            plot3([zcut1*dh zcut2*dh],[rotpoint(1) rotpoint(1)],[ycut1*dh ycut1*dh],'-black','LineWidth',0.5);
            plot3([zcut1*dh zcut2*dh],[rotpoint(1) rotpoint(1)],[ycut2*dh ycut2*dh],'-black','LineWidth',0.5);
            plot3([zcut1*dh zcut1*dh],[rotpoint(1) rotpoint(1)],[ycut1*dh ycut2*dh],'-black','LineWidth',0.5);
            plot3([zcut2*dh zcut2*dh],[rotpoint(1) rotpoint(1)],[ycut1*dh ycut2*dh],'-black','LineWidth',0.5);
        hold on    
            plot3([zcut1*dh zcut2*dh],[rotpoint(1) rotpoint(1)],[rotpoint(2) rotpoint(2)],'k--','LineWidth',0.5);
            plot3([rotpoint(3) rotpoint(3)],[xcut1*dh xcut2*dh],[rotpoint(2) rotpoint(2)],'k--','LineWidth',0.5);
            plot3([rotpoint(3) rotpoint(3)],[rotpoint(1) rotpoint(1)],[ycut1*dh ycut2*dh],'k--','LineWidth',0.5);
      
            
            % hold on
          plot3([32 32],[32 32],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
          plot3([52.8 52.8],[32 32],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
          plot3([74.4 74.4],[32 32],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
          plot3([96 96],[32 32],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);         
        plot3([32 32],[64 64],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([52.8 52.8],[64 64],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([74.4 74.4],[64 64],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([96 96],[64 64],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
      plot3([32 32],[96 96],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([52.8 52.8],[96 96],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([74.4 74.4],[96 96],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0);
     plot3([96 96],[96 96],[92 92],'*blue','MarkerSize',15,'LineWidth',2.0); 
         
         hold on
        for(i=1:13)
        for(kk=1:13)
        plot3([20*dh+(kk-1)*10*dh 20*dh+(kk-1)*10*dh],[20*dh+(i-1)*10*dh 20*dh+(i-1)*10*dh],[24 24],'+blue','MarkerSize',3,'LineWidth',2.0);
        end
        end
         hold on

       hold off   
         
    xlabel('z in m');
    ylabel('x in m');
    zlabel('y in m');        
    set(gca, 'xDir','reverse')      
    set(gca, 'yDir','normal')         
    set(gca, 'zDir','normal') 
    
   set(get(gca,'Ylabel'),'FontSize',14);
   set(get(gca,'Ylabel'),'FontWeight','normal');
   set(get(gca,'Xlabel'),'FontSize',14);
   set(get(gca,'Xlabel'),'FontWeight','normal');
   set(get(gca,'Zlabel'),'FontSize',14);
   set(get(gca,'Zlabel'),'FontWeight','normal');
    set(gca,'FontSize',14);
    set(gca,'FontWeight','normal');
    set(gca,'Linewidth',1.0);
    
    cb = colorbar('vert');
    xlabel(cb, 'norm. \newline gradient v_s');
    
        
      set(gca,'FontSize',14);      
   load('MyColormapsgrad','mycmap');
    set(figure(fignum),'Colormap',mycmap)
view(viewpoint);
    xlim([xcut1*dh xcut2*dh]);
    ylim([zcut1*dh zcut2*dh]);
    zlim([ycut1*dh ycut2*dh]);
    caxis([caxis_value_1 caxis_value_2])
    daspect([1,1,1]);
   % axis tight
   box on

   
   
   
 exportfig(fignum, 'rawgrad_it1_vs.eps','bounds','tight', 'color','rgb', ...
  'preview','none', 'resolution',200, 'lockaxes',1);
