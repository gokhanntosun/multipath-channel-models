clear;close all;
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

for j=1:length(SNR_dB)
    %in-phase component
    i=randn([1 len])./sqrt(n_pow(j));
    %quadtature component
    q=randn([1 len])./sqrt(n_pow(j));
    %build fading component, unit variance
    h=(i+1j*q);
    %calculate the power of the fading
    h_pow=h.*conj(h);
    %extract mean value of the power
    gamma_bar=mean(h_pow);
    
    %csi at rx only
    %take integral to obtain shannon capacity
    c=integral(@(gamma)(log2(1+gamma))...
        .*(exp(-gamma/gamma_bar)/gamma_bar),0,Inf);
    %store in array
    csi_rx(1,j)=c;
        
    %csi at both tx and rx
    %find the cut-off power for transmitter, set constraint
    f=@(gamma_0) integral(@(gamma)...
        (1./gamma_0-1./gamma).*(exp(-gamma/gamma_bar)/gamma_bar),gamma_0,Inf)-1;
    %find the solution of the constraint
    gamma_0=fzero(f,[1e-100,100]);
    %calculate shannon capacity utilizing the cut-off value
    c=integral(@(gamma)(log2(gamma/gamma_0)...
        .*exp(-gamma/gamma_bar)/gamma_bar),gamma_0,Inf);
    %store in array
    csi_trx(1,j)=c;
end

%plot results
figure;hold on;
plot(SNR_dB,awgn_c,'r-.','LineWidth',2);
plot(SNR_dB,csi_rx,'b','LineWidth',2);
plot(SNR_dB,csi_trx,'g','LineWidth',2)
xlabel('SNR(dB)');ylabel('Shannon Capacity(Bits/Sec/Hz)');
title('Channel Capacity vs. SNR | Part II:Q1');
legend('AWGN Channel','Rayleigh Fading Channel w/ CSI@RX',...
    'Rayleigh Fading Channel w/ CSI@TXRX','Location','NorthWest');
grid on;axis square;

