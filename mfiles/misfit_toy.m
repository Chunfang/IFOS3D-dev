%plots misfit normalised to initial misfit value

clear all;

load /data14/sdunkl/3DAWAIT/trunk_JURECA/results_toy/misfit_toy.txt ;
misfit1=misfit_toy./misfit_toy(1)';



x=1:60;

figure(44)
plot(x,misfit1,'b')

 xlabel('iteration');
 ylabel('normalised misfit');
 set(get(gca,'Ylabel'),'FontSize',14);
 set(get(gca,'Ylabel'),'FontWeight','normal');
 set(get(gca,'Xlabel'),'FontSize',14);
 set(get(gca,'Xlabel'),'FontWeight','normal');

