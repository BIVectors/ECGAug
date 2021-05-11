function [rhythm_strip, Qloc, Sloc, Tloc, morph_selected, beat_selected] = ...
    rhythmstrip(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- rhythmstrip()
% Version 1.0
% [rhythm_strip, Qloc, Sloc, Tloc] = rhythmstrip(sigs, Q, S, Tend, n_var, sig_ratio, ap, label)
%
% Script to generate random, physiolgical, and annotated rhythm strips
%
% Uses ECGAug transform_qrst() function
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
% INPUT:
% sigs : 1 x n cell matrix of n beat template signals (n different
% morphologies).
% 
% Q : 1 x n matrix of the locations of QRS onset for each of the n beats in sigs.
% eg [34 56 33] -- sig{1} has QRS onset at sample 34, sig{2} has QRS onset
% at sample 56, and sig{3} has QRS onset at sample 33.
%
% S : 1 x n matrix of the locations of QRS offset for each of the n beat in sigs.
% Data structured the same way as in Q
%
% Tend : 1 x n matrix of the locations of T wave offset for each of the n beat in sigs.
% Data structured the same way as in Q
%
% n_var : Total number of variations per QRST morphology to include for  
% each of the different QRST morphologies.  
% Enter as vector eg [2 5 1] -- this will give 2 variations of 
% morphology 1, 5 variations of morphology 2, 1 variation of morphology 3 etc.
%
% sig_ratio : Ratio of how prevalent each morphology is in the rhythmstrip  
% Enter as vector eg [0.7 0.2 0.1] that sums to 1 -- this will give you 70%
% morphology 1, 20% morphology 2, and 10% morphology 3  
%
% ap : Augparams (passed in as Augparams class)
% See Augparams class comments for parameters that can be adjusted.  This also
% controls the location for loading/saving files, if figures are
% generated/saved, and if data is saved within execution of the function
% rather than saving the output manually after the function executes.
%
% label : how the output figures and data are named so that can keep them
% all together.  Will use the date/time because for rhythm strips with
% multiple beats/files used to create them, the file names get too long.
% The data output will save the files that were used to create the rhythm
% strip.  Format: rhythm_syn_'label'_.mat/.png where label is ideally the
% date and time the rhythmstrip function is called (see rhythmstrip_script)
% Enter as 'char' format.
%
%
% OUTPUT:
%
% rhythm_strip : rhythm strip signal
%
% QLoc : samples corresponding to Qon for all beats in rhythm strip
% SLoc : samples corresponding to Qoff for all beats in rhythm strip
% TLoc : sampels corresponding to Toff for all beats in rhythm strip:
%
% morph_selected : Which morphology was selected for each beat
%
% beat_selected : Which morphology variant was chosen for each beat
%
%
% Examples:
% * Generate your chosen signals and fiducial points outside of the
% rhythmstrip() function - see rhythmstrip_script() function for additional
% details and an examples of how to generate data in the appropriate format
% to pass into the  rhythmstrip() function
%
% (1)
% Generate a 10,000 sample rhythm strip with HR 70, 2 QRST templates stored 
% in cell matrix 'sigs', fiducial points for the 2 QRST templates stored in
% 'Q', 'S', and 'Tend', 3 variants of each QRST template, and a ratio of 
% 80:20 of template 1:template 2
% 
% aug_param = Augparams();
% n_var = [3 3];
% sig_ratio = [0.8 0.2];
% aug_param.hr = 70;
% aug_param.L = 10000;
% [rhythm_strip, Qloc, Sloc, Tloc, m, b] = ...
%       rhythmstrip(sigs, Q, S, Tend, n_var, sig_ratio, aug_param, 'label')
%
% (2)
% Generate a 5000 sample rhythm strip with HR 50, 3 QRST templates stored 
% in cell matrix 'sigs', fiducial points for the 3 QRST templates stored in
% 'Q', 'S', and 'Tend', 5 variants of the first 2 templates and 2 variants 
% of the 3rd template, and a ratio of 80:10:10 ratio of 
% template 1:template 2:template 3
% 
% aug_param = Augparams();
% n_var = [5 5 2];
% sig_ratio = [0.8 0.1 0.1];
% aug_param.hr = 50;
% aug_param.L = 5000;
% [rhythm_strip, Qloc, Sloc, Tloc, m, b] = ...
%       rhythmstrip(sigs, Q, S, Tend, n_var, sig_ratio, aug_param,'label')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

% Error checking input is appropriate
switch true
    case nargin > 8
        error("Wrong number of inputs")        
    case nargin < 8
        error("Wrong number of inputs")     
    case nargin == 8
        if ~isa(varargin{1},'cell') || ~isa(varargin{2},'double') ...
               || ~isa(varargin{3},'double') || ~isa(varargin{4},'double') ...
               || ~isa(varargin{5},'double') || ~isa(varargin{6},'double') ...
               || ~isa(varargin{7},'Augparams') || ~isa(varargin{8},'char')
            error("rhythmstrip input is (sig, Q, S, Tend, n_var, sig_ratio, Augparams, label")
        else
            sig = varargin{1}; 
            Q = varargin{2};
            S = varargin{3};
            T = varargin{4};
            n_var = varargin{5};
            sig_ratio = varargin{6};
            ap = varargin{7};    
            label = varargin{8};  
            n_morph = length(sig);
        end       
end

% Check that length of n_var = n_morph
    if length(n_var) ~= n_morph
        error("n_morp and n_vars should be equal")  
    end
    
% Check that sig_ratio sums to 1
    if round(sum(sig_ratio)) ~= 1
        error("Chcek that sig_ratio elements sum to 1")      
    end
    
 % Check that length of sig_ratio = length of n_var
    if length(sig_ratio) ~= length(n_var)
       error("sig_rato and n_var must have same number of elements")
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% GENERATE VARIANTS FOR SIGNALS PASSED INTO FUNCTION

% Re-adjust Rparams and Augparams if wanted
    ap.abs_tolerance = 0.2;     % Can be less agressive as not recombining

% Create variables
    sig_new = cell([n_morph, max(n_var)]);
    OP = cell([n_morph, max(n_var)]);           % Blank Outparams class
    Qnew = zeros(n_morph, max(n_var));
    Snew = zeros(n_morph, max(n_var));
    Tnew = zeros(n_morph, max(n_var));

% No figures or data saving when calling transform_qrst function
    old_figures = ap.figures;               % Store value for later
    ap.figures = 0;
    old_save_data = ap.save_data;        % Store value for later
    ap.save_data = 0;
    
% Add noise after create the rhythm strip to avoid "double noise"
% Disable noise during variant generation
    old_noise_min = ap.noise_min;
    old_noise_max = ap.noise_max;
    ap.noise_min = 0;                 
    ap.noise_max = 0;
    old_lf_amp_min = ap.lf_amp_min;
    old_lf_amp_max = ap.lf_amp_max;
    ap.lf_amp_min = 0;
    ap.lf_amp_max = 0;
    
for i = 1:n_morph  
    j = 1;
    
% Pass parameters into transform_qrst function to generate appropriate 
% number of variants of each QRST template morphology based on n_var

% In this case, we will not recombine different QRST templates, but 
% will vary the signals slightly to introduce physiological variation.
% If want to include new signals via recombination using recombine_qrst()
% can do this first and then pass in the recombined files to this script
% to generate rhythm strips using recombined QRST templates

bad_counter = 0;        % For keeping track of bad combinations

while j <= n_var(i)
    [sig_new{i,j}, Qnew(i,j), Snew(i,j), Tnew(i,j), bad_combo, OP{i,j}] = transform_qrst(sig{i}, Q(i), S(i), T(i), '', ap,1);
        if bad_combo == 0 
            j = j + 1;
        else
            bad_counter = bad_counter + 1;
                
            if bad_counter > 10
                    % If after 100 iterations it cant find a way to
                    % successfully transform the QRST template, start over
                    % This is usually due to the transformation parameters 
                    % in Augparams being too agressive
                    error("Can't find a good combination with these parameters - Try running again with different parameters")
            end    
        end
    end
end
         
% Restore prior noise values and figure flags
    ap.noise_min = old_noise_min;
    ap.noise_max = old_noise_max;
    ap.lf_amp_min = old_lf_amp_min;
    ap.lf_amp_max = old_lf_amp_max;
    ap.figures = old_figures;
    ap.save_data = old_save_data;

% We now have multiple versions of each QRST morphology stored in sig_new
% and fiducial points for each QRST signal are stored in correponding row, col
% of Qnew, Snew, Tnew


% View all signals and fiducial points
if ap.figures   
    figure('visible','off')

for f = 1:size(sig_new,1)
    row = floor(sqrt(size(sig_new,1)));
    col = ceil(size(sig_new,1)/row);
    subplot(row, col, f)
    %subplot(1,2,f)
    hold on
    rand_color = rand(1,3);
    sgtitle("Included Beats")
    for j = 1:size(sig_new,2)
        if Qnew(f,j) > 0
            plot(sig_new{f,j},'color','k');
            qq = scatter(Qnew(f,j),sig_new{f,j}(Qnew(f,j)),30,'r', 'filled');
            ss = scatter(Tnew(f,j),sig_new{f,j}(Tnew(f,j)),30,'b', 'filled');
            tt =scatter(Snew(f,j),sig_new{f,j}(Snew(f,j)),30,'MarkerEdgeColor',[0 .6 .6], 'MarkerFaceColor',[0 .6 .6]);
            title(sprintf("Morphology %i - %i Variations",f,j));
            xlabel("Samples")
            ylabel("mV")
            legend([qq ss tt],{'Qon', 'Qoff','Toff'})
            
        end
    end
        set(gcf,'Position',[200 200 300*col 300*row]);
end
        saveas(gcf,fullfile(ap.output_folder,strcat('rhythm_syn_',label,'_beats_.png')));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE LOCATIONS OF START OF EACH BEAT ALONG THE RHYTHM STRIP

% RR interval (distance between subsequent QRST complexes) based on heart
% rate (hr)
    rr = round((60000*(ap.freq/1000))/ap.hr);
 
% Calculate QRS and QT intervals
% Generate variables    
    qrs_dur =  zeros(n_morph, max(n_var)); 
    qt_int=  zeros(n_morph, max(n_var));  

% Iterate through to get QRS and QT intervals for each beat and variant
    for i = 1:size(sig_new,1)
        for j = 1:size(sig_new,2)
            qrs_dur(i,j) = Snew(i,j)-Qnew(i,j);
            qt_int(i,j) = Tnew(i,j)-Qnew(i,j);
        end 
    end
    
% Multiply Lx3 so can deal with edge effects.
    L = 3*ap.L;
    
% Number of beats in L samples
    num_beats = L/rr;

% start_seed is the sample number where first QRS will be placed.
% Nominally this defaults to a random number between 100 and 500.
    start_seed = round(100 + (500-100).*rand(1,1));

% Space out QRS complexes based on RR intervals starting at start_seed
    QRS = start_seed:rr:L;

% Placeholder for rhythm strip signal    
    rhythm = nan(1,L);

% Max negative and positive range for random number generator to 
% wobble QRS complex locations (in samples)   
    delta_neg = -ap.delta/(1000/ap.freq); 
    delta_pos = ap.delta/(1000/ap.freq);  

% Change Minimum RR interval from milliseconds to samples
     minRR = ap.minRR/(1000/ap.freq);          

% Add random variation into RR inervals
    Rshift = round((delta_neg + (delta_pos-delta_neg).*rand(1,length(QRS))));
    QRS = QRS + Rshift;

% Check if any RR are below the minimum RR distance and adjust if needed
for i = 1:length(QRS)-1
    if QRS(i+1)-QRS(i) < minRR
       QRS(i+1) = QRS(i)+minRR; 
    end
end

% Check that all QRS are within the signal length and adjust if % needed
for i = 1:length(QRS)
    if QRS(i)>L
        QRS(i) = L-10;
    elseif QRS(i)<1
        QRS(i) = 10;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECT DIFFERENT BEATS/VARIANTS TO ADD TO EACH BEAT LOCATION

% QRS is a list of locations where individual beats will be inserted into
% the rhythm strip signal
    beat_selected = [];
    morph_selected = [];
    c = [];

% Now choose a beat/morphology at random to insert into each QRS location

% Choose which morphology to select from:
% c = matrix of morphologies with 1s representing morphology 1 and 2s
% representing morphology 2 etc. with relative number of each morphology
% based on sig_ratio -- 
% eg [1 1 1 1 1 1 1 2 2 3] is 70% morphology 1, 20% morphology 2, and 10%
% morphology 3.  This is specified based on sig_ratio [0.7 0.2 0.1]
    sig_ratio = sig_ratio * 10;
    
    for i = 1:n_morph
        c = [c ones(1,round(sig_ratio(i)))*i];
    end
    
    % Fix c if sig_ratio isnt ideal based on rounding
    if length(c) < 10
        c = [c ones(1,(10 - length(c)))];
    end

    if length(c) > 10
        c = c(1:10);
    end

        
% Select a beat morphology at random from c, and then select a signal 
% variant of that beat morphology at random    

for i=1:length(QRS)
        
    % Choose a random number between 1 and 10:
    morph = c(round(1 + (10-1).*rand(1,1)));
    morph_selected = [morph_selected morph];
    b = round(1 + (n_var(morph)-1).*rand(1,1));
    beat_selected = [beat_selected b];

    % Check if new inserted beat will be too close to the prior
    % beat, and if so push forward the QRS locations for all subsequent
    % beats

    if i > 1 
        start_pos = QRS(i);
        end_pos = QRS(i-1) + length(sig_new{morph_selected(i-1),beat_selected(i-1)}) - 1;

        if end_pos - start_pos < 5
            QRS([i:end]) =  QRS([i:end]) + 25;  
        end
    end
        
    rhythm(QRS(i):QRS(i)+length(sig_new{morph,b})-1) = sig_new{morph,b}; 
    
end
 
% Now assign fiducial points based on which QRST morphology and which 
% variant was selected for each location within the rhythm strip

% Start by looking at matrix morph_selected to see which QRST morphology
% was used:

% QRS is the locations of the start of each template

for i = 1:length(QRS)
    Qloc(i) = QRS(i) + Qnew(morph_selected(i),beat_selected(i));    
    Sloc(i) = Qloc(i) + qrs_dur(morph_selected(i),beat_selected(i));
    Tloc(i) = Qloc(i) + qt_int(morph_selected(i),beat_selected(i));
    start_sig(i) = QRS(i);
    end_sig(i) = QRS(i) + length(sig_new{morph_selected(i),beat_selected(i)});
end

% Take middle third of rhythm strip
    startpt = L/3+1;
    endpt = 2*L/3;

% Delete anything after start/end of rhythm strip
    start_ind = find(Qloc<=startpt);
    end_ind = find(Tloc>=endpt);

    Qloc(end_ind) = NaN;
    Sloc(end_ind) = NaN;
    Tloc(end_ind) = NaN;
    start_sig(end_ind) = NaN;
    end_sig(end_ind) = NaN;
    morph_selected(end_ind) = NaN;
    beat_selected(end_ind) = NaN;

    Qloc(start_ind) = NaN;
    Sloc(start_ind) = NaN;
    Tloc(start_ind) = NaN;
    start_sig(start_ind) = NaN;
    end_sig(start_ind) = NaN;
    morph_selected(start_ind) = NaN;
    beat_selected(start_ind) = NaN;

    Qloc = Qloc(~isnan(Qloc));
    Sloc = Sloc(~isnan(Sloc));
    Tloc = Tloc(~isnan(Tloc));
    start_sig = start_sig(~isnan(start_sig));
    end_sig = end_sig(~isnan(end_sig));
    morph_selected = morph_selected(~isnan(morph_selected));
    beat_selected = beat_selected(~isnan(beat_selected));

% 'rhythm' now contains all of the beats, but the TP segments are missing and
% will need to be interpolated.  'rhythm2' is the interpolated signal
    [rhythm2,~] = fillmissing(rhythm,'pchip');

% Add random noise to signal
% If noise_min or noise_max = 0, do not add any noise
if ap.noise_min == 0 || ap.noise_max == 0
    % Do nothing
    noise = 0;
else
    noise = abs((ap.noise_max-ap.noise_min).*rand(1,1) + ap.noise_min);
    rhythm2 = rhythm2 + noise * rand(1, length(rhythm2));
end
    
% Add sinusoidal LF noise to signal
% Generate low frequency cosine and alter freq, amplitude, and phase
% X axis samples
    lf_noise_t = 1:1:length(rhythm2);

% Random values
    lf_amp = ap.lf_amp_min + (ap.lf_amp_max-ap.lf_amp_min).*rand(1,1);
    lf_offset = floor(ap.lf_offset_min + (ap.lf_offset_max-ap.lf_offset_min).*rand(1,1));
    lf_freq = ap.lf_freq_min + (ap.lf_freq_max-ap.lf_freq_min).*rand(1,1);

% LF noise signal    
    lf_noise_sig = lf_amp * (cos(0.001*lf_freq*(lf_noise_t-lf_offset)));
    rhythm_strip = rhythm2 + lf_noise_sig;
    
% Cut out middle L samples and adjust locations of fiducial points 
    rhythm_strip = rhythm_strip(startpt:endpt);
    Qloc = Qloc - L/3;
    Sloc = Sloc - L/3;
    Tloc = Tloc - L/3;
    start_sig = start_sig - L/3;
    end_sig = end_sig - L/3;
    end_sig(end_sig > L/3) = [];
 
% Display final generated rhythm strip and annotations   
if ap.figures
    % Plot results
    figure('visible','off')
    plot(rhythm_strip,'k')
    hold on
    q=scatter(Qloc,rhythm_strip(Qloc),15,'r', 'filled');
    s=scatter(Tloc,rhythm_strip(Tloc),15,'b', 'filled');
    t=scatter(Sloc,rhythm_strip(Sloc),15,'MarkerEdgeColor',[0 .6 .6], 'MarkerFaceColor',[0 .6 .6]);
    %st = scatter(start_sig,rhythm2(start_sig),15);
    %ed = scatter(end_sig,rhythm2(end_sig),15);
    %legend([q s t st ed],{'Qon', 'Qoff','Toff','St', 'End'})
    legend([q s t],{'Qon', 'Qoff','Toff'})
    xlim([0 L/3])
    set(gcf,'Position',[100 100 1700 150]);
    xlabel("Samples")
    ylabel("mV")
    
    saveas(gcf,fullfile(ap.output_folder,strcat('rhythm_syn_',label,'_.png')));
    
end
    
    

        
