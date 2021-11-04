function [h_hat,H_hat] = mimoOfdmChannelEst(rxSymbs,pilots,pilotPos,Nt,Nr,nFFT,nTaps,N0,estMethods)
% estimate channels in MIMO-OFDM systems(LSE/MMSE)
% rxSymbs - received symbols a  nFFT X Nr vector
% pilots - pilot symbols, a nP X 1 vector, nP isnumber of pilots
% pilotPos - positions of pilots a nP X 1 vector, 
% Nt -  number of transmit antennas
% Nr - number of receive antennas
% nFFT - fft size
% nTaps - the number of taps of the channel
% N0 -  noise variance
% estMethods: 'LSE' or 'MSE'
% h_hat -  estimated channel impulse response
% H_hat - estimated channel frequency response
% source paper: I. Barhumi, G. Leus and M. Moonen, 
% "Optimal training design for MIMO OFDM systems in mobile wireless channels," 
% in IEEE Transactions on Signal Processing, vol. 51, no. 6, pp. 1615-1624,
% June 2003.

    h_hat = zeros(Nr,Nt,nTaps);     % estimated channel impulse response
    H_hat = zeros(Nr,Nt,nFFT);      % estimated channel frequency response

    %Obtain F as submatrix of dft matrix
    F = dftmtx(nFFT);
    F = F(pilotPos,1:nTaps);

    % received pilots complex symbols
    P_hat = rxSymbs(pilotPos,:);
    
    % construct equivalent channel matrix
    A = zeros(length(pilotPos), nTaps*Nt);
    for j =1: Nt
        A(:,2*(j-1)+(1:2))= diag(pilots(:,j)) * F;
    end

    for i = 1: Nr

        % LSE or MMSE estimation
        if strcmpi(estMethods,'lse')
            hh = pinv(A)*P_hat(:,i);
        else
            hh = inv(A' * A + N0*eye(nTaps*Nt))*A'*P_hat(:,i);
        end
        
        % reshape to obtain channel impulse response and caclate channel 
        % frequence response
        for j = 1: Nt
            h_hat(i,j,:) = hh(nTaps*(j-1)+(1:nTaps));
            H_hat(i,j,:) = fft(squeeze(h_hat(i,j,:)),nFFT);
        end

    end

end

