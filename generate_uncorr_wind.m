function n = generate_uncorr_wind(Fs,L,type,M)
% generate_uncorr_wind: Function that generates M uncorrelated wind noise
% signals with specific temporal features. The AR filtering will be applied
% after the complex coherence constraint.

%   INPUT PARAMETERS:
%       Fs = sampling frequency
%       L = total number of generated samples (time*Fs)
%       type = typology of wind noise (constant or gusts)
%       M = number of microphones
%
%   OUTPUT PARAMETER:
%       n = matrix containing M uncorrelated wind noise signals

% This function is based on the single-channel generation algorithm
% developed in the related paper:
% C. Nelke, P. Vary: "Measurement, Analysis and Simulation of Wind Noise
% Signals for Mobile Communication Devices", International Workshop on
% Acoustic Signal Enhancement (IWAENC), September 2014

% The download is available at:
% https://www.iks.rwth-aachen.de/en/research/tools-downloads/databases/wind-noise-database/

%--------------------------------------------------------------------------
% Copyright (c) 2014, Christoph Nelke
% Institute of Communication Systems and Data Processing
% RWTH Aachen University, Germany
% Contact information: nelke@ind.rwth-aachen.de
%
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the RWTH Aachen University nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%--------------------------------------------------------------------------

% Parameters extracted from wind noise recordings

var_excitation_noise = 0.005874;                              % variance of the excitation signal
mean_state2 = 0.005;                                          % long term gain for low wind
mean_state3 = 0.25;                                           % long term gain for high wind
lpc_coeff = [2.4804   -2.0032    0.5610   -0.0794    0.0392]; % AR parameters
lpc_order = 5;                                                % AR order

% Switching factors of the excitation signals (control of the amount of white noise vs real excitation to use)

alpha1 = 0.15;                                                % low wind switching factor
alpha2 = 0.5;                                                 % high wind switching factor

% Load real excitation signals and generate white noise excitation signals (alpha1 or alpha2 will weight between the two)

data = open('exc_signals.mat');                               % matrix whose rows are the values of real excitation signals (last value=number of samples in the row)
exc_pulses = data.exc_pulses;                                 % the remaining elements of the row -> number contained in the last element are all zeros
nOfExcPulses = size(exc_pulses,1);                            % fixed value: 140 (used to randomly extract one excitation sequence from the matrix)

% Load transition probabilities for the Markov model (control in changing the states 'no wind' ,'low wind' and 'high wind')

switch type
    case 'gusts'
        load('transition_prob_gusts.mat');
    case 'const'
        load('transition_prob_const.mat');
end

% Compute cumulative transition probabilities (each row is associated with the state 1, 2 and 3)

transitionMatrix_cum = transitionMatrix;
for k=1:size(transitionMatrix_cum,1)
    
    for l=2:size(transitionMatrix_cum,1)
        transitionMatrix_cum(k,l) = transitionMatrix_cum(k,l) + transitionMatrix_cum(k,l-1);
    end
end

% Generate state sequence (create a vector of L samples, with flag of 1, 2 or 3 indicating in which state the wind noise is, for each sample)

state_act = 1;                                               % no wind
states_syn = zeros(1,L);                                     % vector that will contain the state flag

for k=lpc_order+1:L                                          % for every sample
    x = rand(1,1);                                           % generate a random number between 0 and 1 representing the probability of changing the state
    p = transitionMatrix_cum(state_act,:);                   % load the row of the cumulative prob. matrix associated to the state 1, 2 or 3
    if x<=p(1)      %no wind                                 % if the random number is less than the first probability: switch to state 1
        state_act = 1;
    elseif x>p(1) && x<=p(2)  %middle wind                   % if the random number is between the second and the third cum. probability: switch to low wind
        state_act = 2;
    else  %high wind                                         % if the random number is higher than the third cum. probability: switch to high wind
        state_act = 3;
    end
    states_syn(k) = state_act;                               % update the vector with the occured state in the current sample
end

% Generate gains for Long-Term behaviour from state sequence [SHARED BETWEEN THE SIGNALS](for no wind state 1, only the Short-Term gain is generated)
excite_sig_noise = zeros(M,L);                               % excitation signal for every microphone signal
g_apl = zeros(1,L);                                          % init of Long-term gain
g_apl(states_syn==2) = mean_state2;                          % low wind
g_apl(states_syn==3) = mean_state3;                          % high wind

% Generate gains for Short-Term behaviour from random processes (Gaussian White Process) [INDIVIDUAL FOR EACH SIGNAL] and apply a smoothing (hann filters)
g_apl_ST = zeros(M,length(g_apl));                           % init of Short-term gain

win1 = hanning(10e3);                                        % smoothing for the LT gain
win1 = win1/sum(win1);
g_apl = fftfilt(win1,g_apl);
g_apl_LT = abs(g_apl);

win2 = hanning(Fs*50e-3);                                    % smoothing for the ST gain
win2 = win2/sum(win2);
g_apl = zeros(M,L);                                          % init of LT*ST gain

for i=1:M
    excite_sig_noise(i,:) = randn(1,L)*sqrt(var_excitation_noise);   % Gaussian white noise for the i-th microphone signal
    g_apl_ST(i,:) = (randn(1,L));                               % Short-term gain as gaussian white process (can be replaced with a Weibull distribution process)
    g_apl_ST(i,:) = abs(fftfilt(win2,g_apl_ST(i,:)));           % Smoothing of the ST gain
    g_apl(i,:) = g_apl_LT.*g_apl_ST(i,:);                       % Combine LT and ST characteristic by modulation of gains (LT*ST)
end

% Init of vectors and parameters for the generation

hwb = waitbar(0,'Generating wind noise signals...');        % waitbar (for every signal contribution)

n = zeros(M,L);                                             % microphones wind signals init

exc_L = 0;                                                  % init of the number of samples in the chosen real excitation row
idx_exc = 1;                                                % init of the index of the current sample in the excitation row

states_syn(states_syn==0) = 1;

exc_pulse_cur = zeros(M,max(exc_pulses(:,end)));            % init of current excitation for the i-th signal

exc_sig = zeros(M,L);                                       % init of the excitation signal for the i-th signal

% Generate wind noise signals (the filtering will be performed after the complex coherence constraint step)
for i=1:M
    for k=lpc_order+1:L
        if ~mod(k,100)
            waitbar(k/L,hwb,sprintf('Generating wind noise signals: microphone %1d',i)); %update the bar
        end
        
        
        %Update excitation pulse position
        if states_syn(k)~=1
            if idx_exc<exc_L                                            % check if the current excitation signals has reached its end
                idx_exc = idx_exc+1;                                    % increment the index of the current sample in the current excitation
            else %end of current pulse -> load next pulse               % if the excitation reached its end load a new excitation
                
                r_pulse = ceil(rand(1,1)*nOfExcPulses);                 % pick a random row from the matrix of excitations for i-th signal
                exc_L = exc_pulses(r_pulse,end);                        % pick the last value of the row : the number of samples in that row
                exc_pulse_cur(i,1:exc_L) = exc_pulses(r_pulse,1:exc_L); % load the corresponding row (actual excitation for the i-th signal)
                idx_exc = 1;                                            % reset the index of the current sample of the excitation
            end
        end
        
        
        if states_syn(k)==1 % no wind                           % for state 1, use just the Gaussian white noise as excitation
            exc_sig(i,k) = excite_sig_noise(i,k)/2;             % GWN for the i-th signal
            
        elseif states_syn(k)==2 % low wind -> noise excitation  % for state 2, excitation signal as a weighted sum of GWN and real excitation with alpha1
            exc_sig(i,k) = alpha1*exc_pulse_cur(i,idx_exc) + (1-alpha1)*excite_sig_noise(i,k);
            
        else % middle/high wind -> turbulent pulse excitiation  % for state 3, excitation signal as a weighted sum of GWN and real excitation with alpha2
            exc_sig(i,k) = alpha2*exc_pulse_cur(i,idx_exc) + (1-alpha2)*excite_sig_noise(i,k);
        end
        
        % Generation of the M excitation signals modulated by long-term and short-term gains
        
        n(i,k) = (g_apl(i,k)* exc_sig(i,k) + g_apl_LT(k)*excite_sig_noise(i,k)*var_excitation_noise);
    end
end

close(hwb)

% Scale outputs
for i=1:M
    n(i,:) = n(i,:)/max(abs(n(i,:)))*0.95;
end

end