clear;close all;
%generate random vectors
%define size
size=10^6;
%define K-factors
Ks=[0.5 1 3 10];
%define signal,power and envelope matrices
s=zeros(length(Ks),size);
s_pow=zeros(length(Ks),size);
s_env=zeros(length(Ks),size);
%define number of bins
bins=75;

for j=1:length(Ks)
    K=Ks(j);
    %LOS components
    %a=sqrt(sqrt(K/2)/2);b=a;
    a=sqrt(K/2);b=a;
    %in-phase component
    i=a+randn([1 size]);
    %quadtature component
    q=b+randn([1 size]);
    %build fading component
    s=(i+1j*q);
    %calculate the fading envelope
    s_env=sqrt(s.*conj(s));
    
    %plots
    subplot(2,2,j);hold on;
    %plot histogram of samples
    histogram(s_env, bins, 'Normalization', 'pdf');
    %fit rician distribution to data
    env_dist=fitdist(s_env','Rician');
    %fit rayleigh for comparison
    rayleigh_dist=fitdist(s_env','Rayleigh');
    %generate x-axis values
    env_x=linspace(0,max(s_env),bins);
    %generate y-axis values
    env_y=pdf(env_dist,env_x);
    rayleigh_dist_y=pdf(rayleigh_dist,env_x);
    %plot rician distribution
    plot(env_x,env_y,'r--','LineWidth',2);
    %plot rayleigh distribution
    plot(env_x,rayleigh_dist_y,'g-.','LineWidth',2);
    %plot settings
    grid on;axis square;xlim([0 7]);ylim([0 1]);
    legend('Fading Envelope','Rician Dist.','Rayleigh Dist.');
    xlabel('h');ylabel('f_h(h)');axis square;
    title(strcat('K=',num2str(K)));
end
%set appropriate size
ss=get(0,'ScreenSize');
set(gcf,'Position',[0.2*ss(3),0.15*ss(4),0.66*ss(3),0.66*ss(4)]);
%mu represents m
%omega represents avg.power