%plots starting, observed and inverted seismograms for one source and
%receiver and x-, y- and z-component
%Input: binary format
clear all; clc;
close all

tracenum=50;

file_inp1='../par/su/cal_toy_vx_it1.su.shot4'; % initial data
file_inp4='../par/su_obs/obs_toy_vx.su.shot4.filt_200'; %filtered with lpb 200 Hz
file_inp10='../par/su_obs/obs_toy_vx.su.shot4.filt_320'; %filtered with lpb 320 Hz
file_inp7='../par/su/cal_toy_vx_it60.su.shot4'; % inverted data

file_inp2='../par/su/cal_toy_vy_it1.su.shot4';
file_inp5='../par/su_obs/obs_toy_vy.su.shot4.filt_200';
file_inp11='../par/su_obs/obs_toy_vy.su.shot4.filt_320'; %filtered with lpb 320 Hz
file_inp8='../par/su/cal_toy_vy_it60.su.shot4';

file_inp3='../par/su/cal_toy_vz_it1.su.shot4';
file_inp6='../par/su_obs/obs_toy_vz.su.shot4.filt_200';
file_inp12='../par/su_obs/obs_toy_vz.su.shot4.filt_320'; %filtered with lpb 320 Hz
file_inp9='../par/su/cal_toy_vz_it60.su.shot4';


tr1 = su2matlab(file_inp1);
tr2 = su2matlab(file_inp2);
tr3 = su2matlab(file_inp3);
tr4 = su2matlab(file_inp4);
tr5 = su2matlab(file_inp5);
tr6 = su2matlab(file_inp6);
tr7 = su2matlab(file_inp7);
tr8 = su2matlab(file_inp8);
tr9 = su2matlab(file_inp9);
tr10 = su2matlab(file_inp10);
tr11 = su2matlab(file_inp11);
tr12 = su2matlab(file_inp12);


nt = tr1.ns; % number of time samples
dt = tr1.dt; % sample interval in micro-seconds
ns = tr1.ns; % number of samples per trace


fig=55;
%--------------------------------------------------------------------------

recx=tr2(tracenum).gx
recy=tr2(tracenum).gy
recx=tr5(tracenum).gx
recy=tr5(tracenum).gy
soux=tr2(tracenum).sx
souy=tr2(tracenum).sy
soux=tr5(tracenum).sx
souy=tr5(tracenum).sy


trace1=tr1(tracenum).trace;
trace2=tr2(tracenum).trace;
trace3=tr3(tracenum).trace;
trace4=tr4(tracenum).trace;
trace5=tr5(tracenum).trace;
trace6=tr6(tracenum).trace;
trace7=tr7(tracenum).trace;
trace8=tr8(tracenum).trace;
trace9=tr9(tracenum).trace;
trace10=tr10(tracenum).trace;
trace11=tr11(tracenum).trace;
trace12=tr12(tracenum).trace;


t=(dt:dt:ns*dt)./10^6; % time in seconds

linewidth=1;

figure(fig)
subplot(1,2,1) 

plot(t,trace1/max(trace5)+4.8,'b-','LineWidth',linewidth);
hold on
plot(t,trace4/max(trace5)+4.8,'k-','LineWidth',linewidth);
hold on
plot(t,(trace2/max(trace5))+2.4,'b-','LineWidth',linewidth);
hold on
plot(t,(trace5/max(trace5))+2.4,'k-','LineWidth',linewidth);
hold on
plot(t,trace3/max(trace5),'b-','LineWidth',linewidth);
hold on
plot(t,trace6/max(trace5),'k-','LineWidth',linewidth);

xlabel('time in s');
xlim([0.015 0.05]);
ylabel('normalized amplitude');
ylim([-2 7]);
set(gca,'ytick',[])

set(get(gca,'Ylabel'),'FontSize',10);
set(get(gca,'Ylabel'),'FontWeight','normal');
set(get(gca,'Xlabel'),'FontSize',10);
set(get(gca,'Xlabel'),'FontWeight','normal');
set(gca,'FontSize',10);
set(gca,'FontWeight','normal');
set(gca,'Linewidth',1.0);

legend('starting','observed','Location','NorthWest');

txt1 = 'x-component';
text(0.038,5.3,txt1,'FontSize',8);
txt2 = 'y-component';
text(0.038,2.9,txt2,'FontSize',8);
txt3 = 'z-component';
text(0.038,0.5,txt3,'FontSize',8);


subplot(1,2,2) 
plot(t,trace10/max(trace8)+4.8,'k-','LineWidth',linewidth);
hold on
plot(t,trace7/max(trace8)+4.8,'r-','LineWidth',linewidth);
hold on
plot(t,(trace11/max(trace8))+2.4,'k-','LineWidth',linewidth);
hold on
plot(t,(trace8/max(trace8))+2.4,'r-','LineWidth',linewidth);
hold on
plot(t,trace12/max(trace8),'k-','LineWidth',linewidth);
hold on
plot(t,trace9/max(trace8),'r-','LineWidth',linewidth);

xlabel('time in s');
xlim([0.015 0.05]);
ylabel('normalized amplitude');
ylim([-2 7]);
set(get(gca,'Ylabel'),'FontSize',10);
set(get(gca,'Ylabel'),'FontWeight','normal');
set(get(gca,'Xlabel'),'FontSize',10);
set(get(gca,'Xlabel'),'FontWeight','normal');
set(gca,'FontSize',10);
set(gca,'FontWeight','normal');
set(gca,'Linewidth',1.0);
set(gca,'ytick',[]);
hleg1 = legend('observed','inverted','Location', 'NorthWest');

txt1 = 'x-component';
text(0.038,5.3,txt1,'FontSize',8);
txt2 = 'y-component';
text(0.038,2.9,txt2,'FontSize',8);
txt3 = 'z-component';
text(0.038,0.5,txt3,'FontSize',8);



% exportfig(fig, ['seismo_trace_toy.eps'],'bounds','tight', 'color','rgb', ...
%   'preview','none', 'resolution',200, 'lockaxes',1);
