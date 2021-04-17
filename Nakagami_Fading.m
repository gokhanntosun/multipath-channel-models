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
bins=150;

%create figures
env_fig=figure;
pow_fig=figure;

for j=1:length(Ks)
    K=Ks(j);
    %LOS components
    a=sqrt(K);b=a;
    %in-phase component
    i=a+randn([1 size]);
    %quadtature component
    q=b+randn([1 size]);
    %build fading component
    s=i+1j*q;
    %calculate the fading power
    s_pow=s.*conj(s);
    %calculate the fading envelope
    s_env=sqrt(s_pow);
    
    %plots
    %calculate nakagami-m distribution parameters
    m=power((K+1),2)/(2*K+1);
    omega=mean(s_pow);
    %approximate with nakagami-m distribution
    env_dist=makedist('Nakagami',m,omega);  %make distribution
    %generate x-axis values
    env_x=linspace(0,max(s_env),bins);
    pow_x=linspace(0,max(s_pow),bins);
    %generate y-axis values
    env_y=pdf(env_dist,env_x);
    pow_y=(power((m/omega),m).*power(pow_x,m-1)...
        .*exp(-pow_x*m/omega))./gamma(m);
    %plot the envelopes
    %set current figure to envelope 
    figure(env_fig);
    subplot(2,2,j);hold on;
    %plot histogram of samples
    histogram(s_env, bins, 'Normalization', 'pdf');
    plot(env_x,env_y,'r--','LineWidth',2);
    %plot settings
    grid on;axis square;xlim([0 7]);ylim([0 1]);
    legend('Fading Envelope','Nakagami-m Dist.');
    xlabel('h');ylabel('f_h(h)');axis square;
    title(strcat('m=',num2str(m)));
    
    %plot power
    figure(pow_fig);
    subplot(2,2,j);hold on;
    %plot power histogram
    histogram(s_pow,bins,'Normalization','pdf');
    plot(pow_x,pow_y,'r--','LineWidth',2);
    %plot settings
    grid on;axis square;xlim([0 20]);ylim([0 1]);
    legend('Fading Power','Nakagami-m Power Dist.');
    xlabel('v');ylabel('f_v(v)');axis square;
    title(strcat('\mu=',num2str(m),' | ','\Omega =',num2str(m/omega)));
end
%set appropriate size
figure(env_fig);
ss=get(0,'ScreenSize');
set(gcf,'Position',[0.2*ss(3),0.15*ss(4),0.66*ss(3),0.66*ss(4)]);
figure(pow_fig);
ss=get(0,'ScreenSize');
set(gcf,'Position',[0.2*ss(3),0.15*ss(4),0.66*ss(3),0.66*ss(4)]);
%mu represents m
%omega represents avg.power