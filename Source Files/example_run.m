% This examples shows how to estimate channels in MIMO-OFDM systems(LSE)

% for simplicity, this program only considers one OFDM symbol per
% transmission, the channel has two taps with power delay profile
% delay = [0,1](smaples), avarage power [0.8, 0.2]. Modultion: QPSK

% source paper: I. Barhumi, G. Leus and M. Moonen, 
% "Optimal training design for MIMO OFDM systems in mobile wireless channels," 
% in IEEE Transactions on Signal Processing, vol. 51, no. 6, pp. 1615-1624,
% June 2003.

nFFT = 64; % fft size
nCP= 16;   % cp length
Nt = 1;     % number of transmit antennas
Nr = 1;     % number of receive antenna
N0 = 0.2; %noise variance
pilotPos = 1:4:nFFT; % pilot positions

% Modulation QPSK
moduOrder = 2;
nOfBits = Nt* nFFT * moduOrder;
b =  randi([0,1],nOfBits,1);
symb = zeros(nOfBits/2,1);
for i = 0: nOfBits/2-1
    symb(i+1) = 1/sqrt(2)*((1-2*b(2*i+1))+1j*(1-2*b(2*i+2)));
end

% antenna mapping
symb = reshape(symb, nFFT, Nt);
% pilot
pilot = symb(pilotPos,:);

%ofdm + cp
txSig = zeros(nFFT+nCP,Nt);
for i = 1: Nt
    x = sqrt(nFFT)*ifft(symb(:,i));
    txSig(:,i) = [x(end-nCP+1:end);x];
end

% channel
nTaps = 2;
avgPow = [0.8;0.2];
delay = [0,1];
h = zeros(Nr,Nt,nTaps);
for i = 1: Nr
    for j = 1: Nt
         h(i,j,:) = avgPow.* (randn(nTaps,1) + 1i*randn(nTaps, 1));
    end
end

% recieved signal
rxSig = zeros(nFFT+nCP+1,Nr);
for i = 1: Nr
    
    for j =  1: Nt
        y = conv(txSig(:,j),squeeze(h(i,j,:)));
        noise = sqrt(N0/2) * (randn(size(y)) + 1j*randn(size(y)));
        rxSig(:,i) = rxSig(:,i) + y + noise;
    end
    
end

% remove delay spread
rxSig = rxSig(1:end-1,:);
% remove cp
rxSig = rxSig(nCP+1:end,:);
% FFT
RxSymbs = fft(rxSig)/sqrt(nFFT);

%---------------------------------------------------------------------------
%channel estimation LSE
[h_hat,H_hat] = mimoOfdmChannelEst(RxSymbs,pilot,pilotPos,Nt,Nr,nFFT,nTaps,N0,'lse');

% equalization with ZF
eqSymb = zeros(nFFT,Nt);
for iSC =  1: nFFT
    H = squeeze(H_hat(:,:,iSC));
    tmp = RxSymbs(iSC,:);
    eqSymb(iSC,:) = pinv(H)* tmp(:);        
end
mse_lse = norm(eqSymb(:)-symb(:))^2/length(symb(:));

close all
plot(symb(:),'r+');
hold on;
plot(eqSymb(:),'bo')

%---------------------------------------------------------------------------
%channel estimation MMSE
[h_hat,H_hat] = mimoOfdmChannelEst(RxSymbs,pilot,pilotPos,Nt,Nr,nFFT,nTaps,N0,'mmse');

% equalization with MMSE
eqSymb = zeros(nFFT,Nt);
for iSC =  1: nFFT
    H = squeeze(H_hat(:,:,iSC));
    tmp = RxSymbs(iSC,:);
    eqSymb(iSC,:) = inv(H' *H + N0*eye(Nt))*H'* tmp(:);        
end
mse_mmse = norm(eqSymb(:)-symb(:))^2/length(symb(:));

plot(eqSymb(:),'ks')
legend('ref symbols','rx symbs(LSE)','rx symbs(MMSE)');
xlim([-2,2]); ylim([-2,2]);
grid on;
xlabel("Inphase");
ylabel("Quadrature");
fprintf('Mean Square Error with noise variance %f: \n',N0);
fprintf('LSE: %f MMSE: %f \n',mse_lse, mse_mmse);



