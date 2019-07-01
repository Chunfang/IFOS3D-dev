%plot of taper funsctions for preconsitioning of gradients

r=-20:1:20;
w0=1./(1+1000*exp(-0.9*abs(r)))
w1=1./(1+1000*exp(-0.8*abs(r)))
w2=1./(1+1000*exp(-0.7*abs(r)))
w3=1./(1+1000*exp(-0.6*abs(r)))
w4=1./(1+1000*exp(-0.5*abs(r)))

%q=1./(1+4000*exp(-0.5*abs(r.*r)))

figure(1)
plot(r,w0)
hold on 
%plot(r,w1)
hold on 
plot(r,w2)
hold on 
%plot(r,w3)
hold on 
plot(r,w4)
hold off
xlabel('distance in m');
ylabel('D(r)');
set(get(gca,'Ylabel'),'FontSize',12);
set(get(gca,'Ylabel'),'FontWeight','normal');
set(get(gca,'Xlabel'),'FontSize',12);
set(get(gca,'Xlabel'),'FontWeight','normal');
set(gca,'FontSize',12);
set(gca,'FontWeight','normal');
set(gca,'Linewidth',1.0);

r1=-20:1:20;

q0=1./(1+1000*exp(-0.7*abs(r.*r)))
q1=1./(1+1000*exp(-0.2*abs(r.*r)))
q2=1./(1+1000*exp(-0.05*abs(r.*r)))

figure(2)
plot(r,q0)
hold on
plot(r,q1)
hold on
plot(r,q2)

xlabel('distance in m');
ylabel('D(r)');
set(get(gca,'Ylabel'),'FontSize',12);
set(get(gca,'Ylabel'),'FontWeight','normal');
set(get(gca,'Xlabel'),'FontSize',12);
set(get(gca,'Xlabel'),'FontWeight','normal');
set(gca,'FontSize',12);
set(gca,'FontWeight','normal');
set(gca,'Linewidth',1.0);