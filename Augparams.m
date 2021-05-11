%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- Augparams()
% Version 1.0
%  
% ECGAug input parameter class for use with recombine_qrst(),
% transform_qrst(), and rhythmstrip()
%
% Copyright 2020 Jonathan W. Waks and Hans F. Stabenau
% Beth Israel Deaconess Medical Center, Boston MA
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%
% Generate a new Augparam class with default values by calling Augparams()
% Examples:
% ap = Augparams()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Augparams
        
properties
        
% *Augmentations*:
%   combine => swap QRS and T from 2 different ECGs.
%   shift => move the baseline up/down
%   stretch => stretch/shrink QRS and T in time dimension
%   scale => stretch/shrink QRS and T amplitude
%   invert => invert QRS or T
%   HF noise => add random Gaussian white noise to final augmented beat
%   LF noise => add low frequency sinusoidal noise to final augmented beat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% DATA INPUT/OUTPUT

% Specify folder and file extension for input
    input_folder = 'C:\ECGAug\templates';
    extension = '*.mat';

% Specify folder for output
    output_folder = 'C:\ECGAug\output';

% Max number of attempts to recombine before give up and move to next set of templates
% (Some combinations of QRS and TW may never be able to successfully combine)
	recombolimit = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSFORMATION PARAMETERS

% Shift: Linear translation in Y-axis
% Range = shift_min to shift_max
    shift_min = -0.1;
    shift_max = 0.1;

% Stretch: Deformation in X-axis
% Range = 1 + stretch_min to 1 + stretch_max
    stretch_min = -0.3;
    stretch_max = 0.3;

% Scale: Deformation in Y-axis
% Range =  = 1 + scale_min to 1 + scale_max
    scale_min = -0.3;
    scale_max = 0.3;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISE PARAMETERS

% HF Noise: Max amplitude of high frequency noise added to signal
% Units are same as signal units
% Range = noise_min to noise_max
% *Values of 0 disable HF noise*
    noise_min = 0.0
    noise_max = 0.03;

% LF Noise: Adds low frequency sinusoidal noise
% lf_noise = lf_amp * cos(2*pi*lf_freq/1000*(t-lf_offset))

% Amplitude range: lf_amp_min to lf_amp_max
    lf_amp_min = 0.0;
    lf_amp_max = 0.1;

% Frequency range: lf_freq_min to lf_freq_max
    lf_freq_min = 0;
    lf_freq_max = 3;
    
% Offset range: lf_offset_min to lf_offset_max
    lf_offset_min = 20;
    lf_offset_max = 100;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLIPPING QRS COMPLEX AND T WAVE

% Randomly flip QRS/T wave
% Values: 0 = disabled ; 1 = enabled
    random_qrsinvert = 1; 
    random_twinvert = 1;  

% If not randomly flipping QRS/T wave, do you want to flip manually?
% Values; 1 = no change; -1 = invert signal
    invert_qrs = 1;
    invert_tw = 1;

% Reverse signals appends QRS2 to TW1 instead of default appending QRS1 to TW2
    flip_signals = 0;
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING PARAMETERS

% Max difference that can be between end of QRS and start of TW to prevent
% rapid changes in signal.
    abs_tolerance = 0.04;       % Units are mV

% Max number of samples that can be interpolated before signal is rejected  
    k_tolerance = 10;

% Physiological limits on QRS and QT intervals (in ms)
    freq = 500;                 % Frequency of signal
    qrs_min = 50;               % Minimum acceptable QRS interval
    qrs_max = 300;              % Maximum acceptable QRS interval
    qt_min = 200;               % Minimum acceptable QT interval
    qt_max = 800;               % Maximum acceptable QT interval    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVING FIGURES AND SIGNALS

% Generate/Save figures? (0 = No, 1 = Yes)
    figures = 1;
      
% Save signals/parameters? (0 = No, 1 = Yes)
    save_data = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHYTHMSTRIP PARAMETERS

% hr: Estimated heart rate (HR) of rhythm strip in bpm.
% Due to variations in the RR intervals the final HR will not be exactly
% this value
    hr = 60;

% L: Length of rhythm strip in samples
    L = 5000; 

% minRR: Minimum RR interval with wobble on in milliseconds  
    minRR = 400;

% delta: max change in RR interval due to wobble in milliseconds    
    delta = 200;
    
% n_morph_max: maximum number of morphologies to include when
% rhythmstrip_script chooses number for you.
    n_morph_max = 6;
    
% n_morph: if want to manually specify number of different morphologies in 
% a rhythm strip.  Set to [] if want program to choose randomly between 1
% and n_morph_max
    n_morph = [];
    
% n_var: Total number of variations per QRST morphology to include
% Enter as vector eg [2 5 1 3] -- this will give 2 variations of 
% morphology 1, 5 variations of morphology 2, 1 variation of morphology 3 etc.    
% Set to [] if want program to choose value of 1-5 randomly for you 
    n_var = [];
    
% sig_ratio: Ratio of how prevalent each morphology is in the rhythmstrip  
% Enter as vector eg [0.7 0.2 0.1] that sums to 1 -- this will give you 70%
% morphology 1, 20% morphology 2, and 10% morphology 3.  Set to [] to want
% program to set this randomly, with most common beat being beat 1.
    sig_ratio = [];

end

methods
        
    % Constructor
    function [obj] = Augparams(varargin)
         
        % Create a 'default' Annoparam class with default values    
        if nargin == 0; return; end
            
    end  % End generation function         
end   % End methods    
end   % End class
