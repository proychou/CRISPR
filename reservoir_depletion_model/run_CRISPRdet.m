%% plots
clear

%good colors
cz=[[31, 119, 180];[255, 127, 14];[44, 160, 44];[214, 39, 40];
    [148, 103, 189];[140, 86, 75];[227, 119, 194];[247, 182, 210]; 
[127, 127, 127];[199, 199, 199];[188, 189, 34];[219, 219, 141];[23, 190, 207]]/255.0;

num_doses=30;
tau=7;
tsim=0:(num_doses*tau); %30 weeks of therapy

L0=1e6; %reservoir size

thr=500; %Hill cure threshhold

%% figure 1 example with estimated parameters

lavg=100;
cov=0.5;
pot=0.5;
potx=0.5+0.24;
potn=0.5-0.24;

[Ls_0,thL,Lt,Lsum]=CRISPRdet(tsim,tau,L0,lavg,cov,pot);
%[Ls_0,thL,Lt,Lsumx]=CRISPRdet(tsim,tau,L0,lavg,cov,potx);
%[Ls_0,thL,Lt,Lsumn]=CRISPRdet(tsim,tau,L0,lavg,cov,potn);

figure(1); clf; hold on
semilogy(tsim/7,Lt,'LineWidth',1)
semilogy(tsim/7,Lsum,'Color','k','LineWidth',5)
%semilogy(tsim/7,Lsumn,'Color','b','LineWidth',2)
%semilogy(tsim/7,Lsumx,'Color','b','LineWidth',2)
semilogy(tsim/7,ones(length(tsim),1)*thr,'--k','LineWidth',2)
xlabel('time (weeks)')
ylabel('Reservoir size (cells)')
set(gca,'YScale','log')
set(gca,'YTick',10.^(0:7))
ylim([1,2e6])
xlim([0,num_doses])
title(['\rho =' num2str(cov) ', \epsilon =' num2str(pot) ', \mu_s = 100'])

%% figure 2 example with 100 coverage and varying efficacy
legpot={};

lavg=100;
cov=1;
figure(2); clf; hold on
for i=1:9
    eps=i/10;
    legpot{i}=['\epsilon = ' num2str(eps)];  
    [Ls_0,thL,Lt,Lsum]=CRISPRdet(tsim,tau,L0,lavg,cov,eps);
    semilogy(tsim/7,Lsum,'Color',cz(i,:),'LineWidth',2)
end
semilogy(tsim/7,ones(length(tsim),1)*thr,'--k','LineWidth',2)
xlabel('time (weeks)')
ylabel('Reservoir size (cells)')
set(gca,'YScale','log')
set(gca,'YTick',10.^(0:7))
ylim([1,2e6])
xlim([0,num_doses])
legpot{i+1}='cure threshold';
title('\rho = 1, mu_s = 100')
legend(legpot)

%% %figure S1, full sensitivity analysis
lavg=10.^(2:4); %avg clone size
cov=(1:4)/4; %coverage, rho
pot=(1:4)/5; %potency, epsilon

legpot={};
for i=1:length(lavg)
        legpot{i}=['\mu_s = ' num2str(lavg(i))];
end
legpot{i+1}='cure threshold';

figure(3); clf;
pind=1;
for i=1:length(cov)
    for j=1:length(pot)
        for k=1:length(lavg)
            [Ls_0,thL,Lt,Lsum]=CRISPRdet(tsim,tau,L0,lavg(k),cov(i),pot(j));
            subplot(4,4,pind)
            hold on
            %make the color
            semilogy(tsim/7,Lsum,'Color',cz(k,:))%,'k','LineWidth',5)
        end
        semilogy(tsim/7,ones(length(tsim),1)*thr,'--k','LineWidth',2)
        
        if pind>13
            xlabel('time (weeks)')
        end
        if sum(pind==[1,5,9,13])>0
            ylabel('Reservoir size (cells)')
        end
        set(gca,'YScale','log')
        set(gca,'YTick',10.^(0:7))
        ylim([1,1e7])
        xlim([0,num_doses])
        pind=pind+1;
        title(['\rho =' num2str(cov(i)) ', \epsilon =' num2str(pot(j))])
    end
end
legend(legpot)

%% %figure S2, histogram examples

figure(4); clf;
for k=1:length(lavg)
    subplot(130+k)
    [Ls_0,thL,Lt,Lsum]=CRISPRdet(tsim,tau,L0,lavg(k),cov(i),pot(j));
    histogram(log10(Ls_0),'BinMethod','integers','FaceColor',cz(k,:))
    xlim([-1,7])
    xlabel('log10 clone size')
    if k==1
        ylabel('counts')
    end
    set(gca,'YScale','log')
    ylim([1,2e4])
    set(gca,'XTick',0:6)
    set(gca,'YTick',10.^(0:6))
    title(['\mu_s = ' num2str(lavg(k))])
end
%%
