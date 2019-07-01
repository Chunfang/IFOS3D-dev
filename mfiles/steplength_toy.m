%plot of steplength vs. iteration

clear all;
load /data14/sdunkl/3DAWAIT/trunk_JURECA/results_toy/steplength_toy.txt ;

misfit1=steplength_toy';

x=1:60;

figure(41)
plot(x,misfit1)
 xlabel('iteration');
 ylabel('steplength');
 set(get(gca,'Ylabel'),'FontSize',14);
 set(get(gca,'Ylabel'),'FontWeight','normal');
 set(get(gca,'Xlabel'),'FontSize',14);
 set(get(gca,'Xlabel'),'FontWeight','normal');

