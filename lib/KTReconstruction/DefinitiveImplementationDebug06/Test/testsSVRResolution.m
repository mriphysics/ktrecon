addpath(genpath('/home/lcg13/Work/DefinitiveImplementationDebug06'));

% k=0:0.01:0.99;
% 
% ord=[2 4 8 16 32 64];
% figure
% for n=1:length(ord)
%     v=sin(pi*k/2).^ord(n);
%     nv=norm(v)/sqrt(length(k));
%     plot(k,v/nv)
%     hold on
%     k0=2*acos(1/sqrt(ord(n)));
%     plot(k0/pi,(sin(k0/2).^ord(n))/nv,'*')
% end
% return


%n=cos(pi*k/2);
%n=n.^(-2);

%figure
%plot(k,n);

%return

%TIKH=[1 2 5 10 20];%[0 2 5 10 20];%[0 10 20];%[0 0.2 0.5 1 2 5];
%TIKH=[1 2];
TIKH=80;
RESO=1.1;%[1.2 1.3 1.4 1.5 1.6];

%WORST SNR
path='2018_04_09/MC_2731';
fil='mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4';

%GOOD SNR THIN CORTEX
%path='2018_04_06/SC_2630';
%fil='sc_06042018_1338179_8_2_sf7zt2mbzax0sensemi8s1V4';

%GOOD SNR
%path='2018_03_02/MI_43210';
%fil='mi_02032018_0943009_8_2_sf7zt2mbzax0sensemi8s1V4';

%RELATIVELY GOOD SNR MORE MOTION ARTIFACTED
%path='2018_02_23/MI_41210';
%fil='mi_23022018_1032032_7_2_sf7zt2mbzax0sensemi8s1V4';
fileRECON=sprintf('/home/lcg13/Data/rawDestin/ReconstructionsDebug06/%s/An-S2/%s',path,fil);

for o=1:length(TIKH)
    for p=1:length(RESO)
        fprintf('TIKH-%.3f\n',TIKH(o)/100);
        suff{1}=sprintf('20VR_TIKH-%.2f-RESO-%.2f',TIKH(o),RESO(p));
     
        [y,MS,MT]=readNII(fileRECON,suff,1);
        x=y{1};
        figure
        for di=1:3
            [f,k]=estimatePSD(x,di,1,'dB');
            hold on
            plot(k,f)
        end
        legend('LR','AP','IS')
        xlabel('mm')
        ylabel('dB')
     
        axis([-0.5 0.5 -30 50])
        v=[-1 -1.1 -1.2 -1.3 -1.4 -1.5 -2 -3 3 2 1.5 1.4 1.3 1.2 1.1 1];
        v=1./(2*v);
        set(gca,'xtick',v);
        set(gca,'xticklabel',{'-1','-1.1','-1.2','-1.3','-1.4','-1.5','-2','-3','3','2','1.5','1.4','1.3','1.2','1.1','1'})
        title(sprintf('Regularizer: %.2f / Resolution: %.2f',TIKH(o),RESO(p)))
        set(gcf, 'Position', get(0,'Screensize'))
    end
end

% 
% 
% fileRECON='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_AVR';
% suff={'Test1.00'};
% 
% [y,MS,MT]=readNII(fileRECON,suff,1);
% 
% x=y{1};
% res=0.5;
% 
% N=size(x);
% Nr=round(N/res);
% 
% x=real(resampling(x,Nr));
% MS{1}=MS{1}*res;
% for n=1:3
%     MT{1}(n,n)=MT{1}(n,n)*res;
% end
% 
% fileName='/home/lcg13/Data/rawDestin/ReconstructionsDebug06/2018_04_09/MC_2731/An-S2/mc_09042018_0942279_9_2_sf7zt2mbzax0sensemi8s1V4_AVR';
% 
% suff={'Test0.80Res'};
% 
% y{1}=x;
% 
% writeNII(fileName,suff,y,MS,MT); 
% 
% return
