%clear;close all;
%defs
%define snr values in dBs
SNR_dB=1:30;
%define sample length
len=10^6;
%allocate output arrays
csi_rx=zeros(1,length(SNR_dB));
csi_trx=zeros(1,length(SNR_dB));
%convert to std deviations for scaling
n_pow=2./(10.^(SNR_dB./10));
%calculate awgn capacity, channel information irrelevant
awgn_c=log2(1+(2./n_pow));
%select m value
m=2;
%solve (K+1)^2/(2K+1)=m
f=@(K)(pow2(K+1)/(2*K+1))-m;
fzero(f,[1e-100 100]);

for j=1:length(SNR_dB)
    %select appropriate LOS values with K value
    a=sqrt(K/n_pow(j));b=a;
    %in-phase component
    i=a+randn([1 len])./sqrt(n_pow(j));
    %quadrature component
    q=b+randn([1 len])./sqrt(n_pow(j));
    %build fading component
    h=(i+1j*q);
    %calculate the fading power
    h_pow=h.*conj(h);
    %parameters for nakagami power distribution
    mu=m;omega=mean(h_pow);
    
    %csi at rx only
    %calculate the integral to obtain shannon capacity
    f=@(g)log2(1+g).*power(mu/omega,mu)...
        .*power(g,mu-1).*exp(-g.*mu/omega)/gamma(mu);
    c=integral(f,0,Inf);
    %store in array
    csi_rx(1,j)=c;
    
    %csi at both tx and rx
    %find the cut-off power for transmitter, set constraint
    f=@(gamma_0) integral(@(g)...
        (1./gamma_0-1./g).*power(mu/omega,mu)...
            .*power(g,mu-1).*exp(-g.*mu/omega)/gamma(mu),gamma_0,Inf)-1;
    %fins solution to constraint
    gamma_0=fzero(f,[1e-100,1000]);
    %calculate shannon capacity using the cut-off value
    c=integral(@(g)(log2(g/gamma_0)...
        .*power(mu/omega,mu).*power(g,mu-1)...
            .*exp(-g.*mu/omega)/gamma(mu)),gamma_0,Inf);
    %store in array
    csi_trx(1,j)=c;
end
%plot results
figure;hold on;
plot(SNR_dB,awgn_c,'r-.','LineWidth',2);
plot(SNR_dB,csi_rx,'b','LineWidth',2);
plot(SNR_dB,csi_trx,'g','LineWidth',2)
xlabel('SNR(dB)');ylabel('Shannon Capacity(Bits/Sec/Hz)');
title('Channel Capacity vs. SNR | Part II:Q2');
legend('AWGN Channel','Nakagami Fading Channel w/ CSI@RX',...
    'Nakagami Fading Channel w/ CSI@TXRX','Location','NorthWest');
grid on;axis square;