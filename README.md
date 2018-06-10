# Simulating multi-channel wind noise based on the Corcos model

This folder contains the proposed multi-channel artificial wind noise generator.

The main script is "generate_coherent_wind_noise.m".

The other files are:

- generate_uncorr_wind.m = Matlab function for the generations of M uncorrelated wind noise signals with a shared Long-Term gain and individual Short-Term gains

- cohere_mod.m = Matlab function for the computation of the complex coherence

- stft/istft/biorwin.m = Matlab functions for the computation of the STFT/iSTFT

- mix_signals.m = Matlab function for the decomposition of the coherence matrix to compute the mixing matrix
                 
- exc_signals.mat = Matrix whose rows are the excitation signals from real wind noise recordings
            
- transition_prob_const.mat = Matrix of the state transitions probabilities for no-wind,low-wind and high-wind for a constant flow

- transition_prob_gusts.mat = Matrix of the state transitions probabilities for no-wind,low-wind and high-wind for a gusty flow