% Multichannel Wind Noise Generation

% Author       : Daniele Mirabilii
% Date         : 21-05-2018
%
% Reference    : D. Mirabilii and E.A.P. Habets "Simulating multi-channel wind noise based on the Corcos model,” arXiv, 2018.

% Copyright (C) 2018 Friedrich-Alexander-Universität Erlangen-Nürnberg, Germany  
%  
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

clear
close all
clc

% Init parameters
Fs = 16000;                     % sampling frequency [Hz]
L =  30*Fs;                     % number of desired samples (seconds*sampling frequency)
type = 'const';                 % type of wind: constant or gusts
M = 2;                          % number of microphones

% Generate M uncorrelated wind noise signals with a shared long-term gain and individual short-term gains
uncorr_wind = generate_uncorr_wind(Fs, L, type, M); 

%% Wind noise coherence constraint

% Initialization
K = 2048;                       % FFT length
d = 0.004;                      % microphone spacing (assuming a Uniform Linear Array) [m]
theta = 0;                      % direction of arrival of the wind stream [rad]
d_m = d*cos(theta);             % DOA-dependent phase difference term
U = 1.8;                        % free-field wind speed [m/s]
alpha1 = -0.125;                % experimental longitudinal decay
alpha2 = -0.7;                  % experimental lateral decay
alpha = alpha1*cos(theta) + ...
    alpha2*sin(theta);          % DOA-dependent decay value
ww = 2*pi*Fs*(0:K/2)/K;         % angular frequency vector

% Uncorrelated wind noise signals as input (not yet filtered)
wind = uncorr_wind.';

%% Generate matrix with desired spatial coherence
DC = zeros(M,M,K/2+1); % init of the coherence matrix

% Build of the coherence matrix
for p = 1:M 
    for q = 1:M
        if p == q
            DC(p,q,:) = ones(1,1,K/2+1);
        else
            if p>q
                DC(p,q,:) = exp(alpha*ww*abs(p-q)*d/(0.8*U)).*exp(1i*ww*abs(p-q)*d_m/(0.8*U));
            else
                DC(p,q,:) = exp(alpha*ww*abs(p-q)*d/(0.8*U)).*exp(-1i*ww*abs(p-q)*d_m/(0.8*U));
            end
        end
    end
end

%% Generate wind noise signals with desired spatial coherence
x = mix_signals(wind,DC,'cholesky');

%% LPF AR filter (all-pole filter) for spectral characteristics matching with real wind noise
nf = zeros(size(x)); % init of the filtered wind noise samples

lpc_coeff = [2.4804   -2.0032    0.5610   -0.0794    0.0392]; % AR parameters from real wind recordings
lpc_order = 5; % AR order

for i=1:M
    % Filtering with the AR filter
    nf(:,i) = filter(1,[1 -lpc_coeff],x(:,i)); 
    
    % Scaling to avoid clippings
    nf(:,i) = nf(:,i)/max(abs(nf(:,i)))*0.95; 
end

%% Calculate theoretical and generated coherence
K_eval = 2048;                          % FFT length for the coherence analysis
ww = 2*pi*Fs*(0:K_eval/2)/K_eval;       % angular frequency vector for the analysis
sc_theory = zeros(M-1,K_eval/2+1);      % init of the theoretical coherence
sc_generated = zeros(M-1,K_eval/2+1);   % init of the generated coherence

for m = 1:M-1
    % Theoretical coherence
    sc_theory(m,:) = exp(alpha*ww*m*d/(0.8*U)).*exp(1i*ww*m*d_m/(0.8*U)); 
    
    % Generated coherence
    [sc_tmp, Freqs]= cohere_mod(nf(:,1),nf(:,m+1),K_eval,Fs,hanning(K_eval),0.75*K_eval); 
    sc_generated(m,:) = sc_tmp; 
end

%% Calculate MSE between theoretical and generated coherence for every microphone pair
nMSE = zeros(M-1,1);

for m = 1:M-1
    nMSE(m) = 10*log10(sum(abs(((sc_theory(m,:))-(sc_generated(m,:)))).^2)./sum(abs(sc_theory(m,:))).^2);
end

%% Plot complex coherence of each microphone pair and the related MSE
for m = 1:M-1
    figure(m);
    subplot(2,1,1);
    plot(Freqs/1000, real(sc_theory(m,:)), '-k', 'LineWidth', 1.5);
    hold on;
    plot(Freqs/1000, real(sc_generated(m,:)), '-.b', 'LineWidth', 1.5);
    hold off;
    xlabel('Frequency [kHz]');
    ylabel('Real Spatial Coherence'); axis([0 1 -1 1]);
    title(sprintf('Microphone distance %1.3f m ', m*d));
    legend('Theory', sprintf('Generated (nMSE = %2.1f dB)', nMSE(m)));
    grid on;
    subplot(2,1,2);
    plot(Freqs/1000, imag(sc_theory(m,:)), 'LineWidth', 1.5);
    hold on;
    plot(Freqs/1000, imag(sc_generated(m,:)), 'LineWidth',1.5), axis([0 1 -1 1]), xlabel('Frequency [kHz]');
    ylabel('Imaginary Spatial Coherence'), legend('Theory', sprintf('Generated (nMSE = %2.1f dB)', nMSE(m))); grid on;
end

%% Write .wav file with coherent wind noise
% audiowrite('Coherent wind noise crosswind 4mm.wav',nf,Fs);