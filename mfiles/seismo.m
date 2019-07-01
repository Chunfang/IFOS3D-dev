nt=3334; dt=0.000045;
nrec=459;
%nrec=130;
file_inp1='/workag13/FWI/random_KTB/3D/su_obs/obs_random7_vz_it1.bin.shot2filtLP320';
file_inp3='/workag13/FWI/random_KTB/3D/su/cal_random71_vz_it73.bin.shot2filtLP320';
file_inp2='/workag13/FWI/random_KTB/3D/su/cal_random7_vz_it1.bin.shot2filtLP320';
file_inp39='/workag13/FWI/random_KTB/2D_10sources/su/cal_random723_vz_it73.bin.shot2filtLP320';
file_inp29='/workag13/FWI/random_KTB/2D_10sources/su/cal_random723_vz_it1.bin.shot2filtLP320';
file_inp19='/workag13/FWI/random_KTB/2D_10sources/su_obs/obs_random72a_vz_it1.bin.shot2filtLP320';


fig=55;
%--------------------------------------------------------------------------

SEIS1=binread(file_inp1,nt,nrec);
SEIS2=binread(file_inp2,nt,nrec);

t=dt:dt:nt*dt;
trace=SEIS1(:,3);
size(trace)
figure(1)
plot(t,trace);
figure(fig);
imagesc(1:nrec, t, (SEIS1-SEIS2)/max(max(SEIS1)));
size(SEIS1)


fpal = f_colpal(77, 'linear', 4, 256);
colormap(fpal);

xlabel('receiver number');
ylabel('time in s');
set(get(gca,'Ylabel'),'FontSize',12);
set(get(gca,'Ylabel'),'FontWeight','normal');
set(get(gca,'Xlabel'),'FontSize',12);
set(get(gca,'Xlabel'),'FontWeight','normal');
set(gca,'FontSize',12);
set(gca,'FontWeight','normal');
set(gca,'Linewidth',1.0);

cb = colorbar('vertical');
zlab = get(cb,'ylabel');
set(zlab,'String','residuals normalized to observed data','fontsize',16)
set(gca,'ydir','reverse')


ylim([0.03 0.14]);

caxis([-0.1 0.1])