%plot of seismograms along two perpendicular receiver lines through shot
%position

clear all;
close all;

file_inp1='/data14/sdunkl/3DAWAIT/trunk_JURECA/results_toy/su/obs_toy_vy_it1.su.shot4';
file_inp2='/data14/sdunkl/3DAWAIT/trunk_JURECA/results_toy/su/obs_toy_vz_it1.su.shot4';
file_inp3='/data14/sdunkl/3DAWAIT/trunk_JURECA/results_toy/su/obs_toy_vy_it1.su.shot4';
%file_inp3='/data14/sdunkl/surface/results6/su/cal_surface5_vz_it70.su.shot10_int';

fignum=17;


shot1=su2matlab(file_inp1);
shot2=su2matlab(file_inp2);
shot3=su2matlab(file_inp3);

dt=shot1(1).dt*1e-6
nt=shot1(1).ns
t=dt:dt:nt*dt;
[shot1(1).sx]
[shot1(1).sy]
numgx=find([shot1.gx] == [shot1(1).sx])

for i=1:13
shot1cutx(i)=shot1(numgx(i)); 
shot2cutx(i)=shot2(numgx(i));
shot3cutx(i)=shot3(numgx(i));
offset(i)=([shot1cutx(i).gy]-[shot1cutx(i).sy])/100;
end
%max(max(shot1cutx(1).trace)) 
%max((shot2cutx(1).trace))
figure(fignum)
for i=1:13
plot(t,(shot1cutx(i).trace)/abs(max(shot2cutx(i).trace))+0.01*offset(i),'k-','LineWidth',1);hold on
%plot(t,(shot2cutx(i).trace)/abs(max(shot1cutx(i).trace))+0.05*offset(i),'b-','LineWidth',1);hold on
%plot(t,(shot3cutx(i).trace)/abs(max(shot1cutx(i).trace))+0.03*offset(i),'r-','LineWidth',1);
hold on
end
 xlabel('time in s');
 ylabel('offset in m');
 set(get(gca,'Ylabel'),'FontSize',14);
 set(get(gca,'Ylabel'),'FontWeight','normal');
 set(get(gca,'Xlabel'),'FontSize',14);
 set(get(gca,'Xlabel'),'FontWeight','normal');
 set(gca,'YTick',0:5:5)
 set(gca,'YTickLabel',{0,100})


numgy=find([shot1.gy] == [shot1(1).sy])

for i=1:13
shot1cuty(i)=shot1(numgy(i));
shot2cuty(i)=shot2(numgy(i));
shot3cuty(i)=shot3(numgy(i));
offset(i)=([shot1cuty(i).gx]-[shot1cuty(i).sx])/100;
end

figure(fignum+1)
for i=1:13
plot(t,(shot1cuty(i).trace)/abs(max(shot2cuty(i).trace))+0.01*offset(i),'k-','LineWidth',1);hold on
%plot(t,(shot2cuty(i).trace)/abs(max(shot1cuty(i).trace))+0.01*offset(i),'b-','LineWidth',1);hold on
%plot(t,(shot3cuty(i).trace)/abs(max(shot1cuty(i).trace))+0.03*offset(i),'r-','LineWidth',1);
hold on
end
 xlabel('time in s');
 ylabel('offset in m');
 set(get(gca,'Ylabel'),'FontSize',14);
 set(get(gca,'Ylabel'),'FontWeight','normal');
 set(get(gca,'Xlabel'),'FontSize',14);
 set(get(gca,'Xlabel'),'FontWeight','normal');
 set(gca,'YTick',0:3:3)
 set(gca,'YTickLabel',{0,100})






