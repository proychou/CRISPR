%% deterministic multistrain simulation
%DBR May 2018
%function that simulates a multstrain reservoir decay
%inputs:
%   tsim, time series in days
%   L0, reservoir size
%   lavg, average clone size (from 1 to L0)
%   rho, coverage fraction of strains that CRISPR deletes
%   eps, potency of CRISPR per treatment
%outputs:
%   Ls_0, initial reservoir sizes
%   thL, all clearance rates
%   Lt, each strain over time, or 100 random if there are lots
%   Lt, each strain over time, or 100 random if there are lots
%   Lsum, total reservoir size over time


function [Ls_0,thL,Lt,Lsum]=CRISPRdet(tsim,tau,L0,lavg,rho,eps)
tic;

%make a reservoir with the approximate total size and distribution
Ls_0=round(lognrnd(log10(lavg),log10(lavg),[round(L0/lavg),1])); %initial clone sizes
Ls_0=Ls_0(Ls_0>0); %remove zeros

totl=sum(Ls_0);
R=length(Ls_0);
while totl<L0
    ls=round(lognrnd(log10(lavg),log10(lavg))); %initial clone sizes
    if ls>0
        R=R+1;
        Ls_0(R)=ls;
    end
    totl=sum(Ls_0);
end

num_clones=length(Ls_0); %number of clones

randintz=randperm(num_clones); %make a list of clones
covered=randintz(1:round(num_clones*rho)); %which clones are covered by CRISPR?

%clearance rates
%thL=zeros([num_clones,1]); %no clearance rate
thL=ones([num_clones,1])*5.2e-4; %single clearance rate for all

%hl=randn([num_clones,1])*1.5+3.6;%crooks halflife
%thL=log(2)./(hl*365);
%thL(thL>0)=0;

%%
int=1;
num_tracked=min(100,num_clones); %number to track, max is 100?
L=Ls_0; Lt=zeros(length(tsim),num_tracked); Lsum=zeros(length(tsim),1); %initialize
for ti=tsim
    %only keep max 100 lines over time
    Lt(int,:)=L(1:num_tracked);
    Lsum(int,:)=sum(L);
    L=L - thL.*L;
    if mod(ti,tau)==0
        L(covered)=L(covered)*(1-eps);
    end
    
    %other statistics if interested later
    %richness
    %Rt(int)=sum(L>1);
    %evenness
    %ra=L(L>1)/sum(L(L>1)); %relative abundances
    %St(int)=-sum(ra.*log(ra));
    %Et(int)=St(int)/log(Rt(int));
    
    int=int+1; %update index
end
    
toc;

