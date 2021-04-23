%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- Outparams()
% Version 1.0
%  
% ECGAug output parameter class for use with recombine_qrst(),
% transform_qrst(), and rhythmstrip()
%
% Copyright 2020 Jonathan W. Waks and Hans F. Stabeneau
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
% Generate a new Outparam class with default values by calling Outparams()
% Examples:
% op = Outparams()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef Outparams
    
    % Stores output parameters used in recombine_qrst(), transform_qrst(), and
    % rhythmstrip()
    
properties
    
    file1                       % Filename of sig 1
    file2                       % Filename of sig 2
    yshift                      % Y-axis shift of sig1 and sig2
    yscale_qrs                  % Y scaling of QRS template
    yscale_t                    % Y scaling of T template
    xscale_qrs                  % Temporal scaling of QRS template
    xscale_t                    % Temporal scaling of T template
    flip_qrs                    % T wave flipped (0 = No, 1 = Yes)
    flip_tw                     % QRS flipped (0 = No, 1 = Yes)
    high_freq_amp               % Max amplitude of high frequency noise
    low_freq_amp                % Amplitude of low frequncy noise
    low_freq_offset             % Offset of low frequency noise
    low_freq                    % Frequency of low frequency noise
    k1                          % Interpolation sample limit for QRST recombination
    k2                          % Interpolation sample limit for QRST transformation
    abs_tolerance               % Absolute voltage tolerance for recombining QRS and T templates
end
    
methods  
    
    % Constructor
    function [obj] = Outparams(varargin)
         
        % Create a 'default' Outparam class with default values    
        if nargin == 0; return; end
               
    end   % End generation function        
end   % End methods    
end   % End class
