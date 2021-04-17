clear;close all;
%generate random vectors
%in-phase component
i=randn([1 10^6]);
%quadtature component
q=randn([1 10^6]);
%build signal
s=(i+1j*q);

%calculate fading power
s_pow=s.*conj(s);
%calculate fading envelope
s_env=sqrt(s_pow);

%define number of bins
bins=150;

%plot envelope of fading with corresponding rayleigh distribution
%create a new figure
figure;subplot(1,2,1);hold on;
%create histogram
histogram(s_env,bins,'Normalization','pdf');
%fit rayleigh distribution to data
env_dist=fitdist(s_env','Rayleigh');
%generate x-axis values
env_x=linspace(0,max(s_env),bins);
%generate y-axis values
env_y=pdf(env_dist,env_x);
%plot the distribution
plot(env_x,env_y,'r--','LineWidth',1.5);
%show grid, add legend...
grid on;ylim([0 1]);xlim([0 5]);
legend('Fading Envelope','Rayleigh Distribution');
xlabel('h');ylabel('f_h(h)');axis square;
title('Envelope of Fading');

%plot power with corresponding exponential distribution
%create a new figure
subplot(1,2,2);hold on;
%create histogram
histogram(s_pow,bins,'Normalization','pdf');
%fit exponential dist to data
pow_dist=fitdist(s_pow','Exponential');
%generate x-axis values
pow_x=linspace(0,max(s_pow),bins);
%generate y_axis values
pow_y=pdf(pow_dist,pow_x);
%plot the distribution
plot(pow_x,pow_y,'r--','LineWidth',2);
%show grid, add legend...
grid on;ylim([0 1]);xlim([0 15]);
legend('Signal Power', 'Exponential Distribution');
xlabel('v');ylabel('f_v(v)');axis square;
title('Envelope-Squared, Power');