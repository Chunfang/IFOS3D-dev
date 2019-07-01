clear all
%close all

% Definition von Parametern
% Homogener Vollraum
npts=10000;
delta=0.0001;




t=(0:delta:(npts-1)*delta);

% Parameter: Homogener Vollraum
%T_d=0.004;
%T_d=0.0033;
T_d=1/300.0;
%T_d=0.002;
F0=1;

npts_start=1;
npts_end=T_d/delta;

wavelet=zeros(1,npts);


% Definition des Quellwavelets
wavelet(1:npts_end)=F0*(sin((pi*t(1:npts_end))/T_d)).^3;

%figure
%plot(t,wavelet)
%xlabel('Zeit in s')

WAVELET=delta*fft(wavelet);
maximum=max(max(abs(WAVELET)))
WAVELET=WAVELET./maximum;
frequenz=[0:fix(npts/2),-ceil(npts/2)+1:-1]'/(npts*delta);  %Definition des kompletten Frequenzvektors
f=[0:fix(npts/2)]'/(npts*delta);

%{
figure
subplot(211)
plot(f(1:fix(npts/2)+1),abs(WAVELET(1:fix(npts/2)+1)))
xlab=xlabel('Frequenz in Hz')
set(gca,'FontSize',14)
tit=title(['Laenge des Quellsignals: ' num2str(T_d) ' s']);
set(xlab,'FontSize',14)

subplot(212)
plot(f(1:fix(npts/2)+1),abs(WAVELET(1:fix(npts/2)+1)))
xlab=xlabel('Frequenz in Hz')
set(gca,'FontSize',14)
ylim([0 0.01])
%set(gca,'XTick',[0 30 50 70 90 110 150 200 250])
%set(gca,'XTickLabel',{'0';'30';'50';'70';'90';'110';'150';'200';'250'})
set(xlab,'FontSize',14)
%}


figure
plot(t,wavelet,'-k')
xlab=xlabel('Time in s');
ylab=ylabel('s(t) (dimensionless)');
set(gca,'FontSize',18)
set(xlab,'FontSize',20)
set(ylab,'FontSize',20)
xlim([0 1.5*T_d])


figure
plot(f(1:fix(npts/2)+1),abs(WAVELET(1:fix(npts/2)+1)),'-k')
xlab=xlabel('Frequency in Hz');
ylab=ylabel('|S(\omega)| in 1/Hz');
set(gca,'FontSize',18)
set(xlab,'FontSize',20)
set(ylab,'FontSize',20)
%tit=title(['Laenge des Quellsignals: ' num2str(T_d) ' s']);
%set(xlab,'FontSize',14)
xlim([0 600])
