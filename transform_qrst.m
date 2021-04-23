function [sig_new, Qnew, Snew, Tnew, bad_combo, op] = ...
    transform_qrst(sig, Q, S, T, fname, ap, syn_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- transform_qrst()
% Version 1.0
%  
% Takes a QRST template and annotations (Qon, Qoff, and Toff) and performs
% random transformation of the QRS template and T template indepdently and
% then recombines into a new signal with adjusted annotations
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
% INPUT:
%
% sig : QRST template
%
% Q : Qon of sig
% S : Qoff of sig
% T : Toff of sig
%
% fname : filename of sig (for renaming files)
%
% ap : Augparams (passed in as Augparams class)
%
% syn_num : number to differentiate multiple transformations of same QRST
% template
%
%
% OUTPUT:
%
% sig_new : recombined QRST template
%
% Qnew : Qon of sig_new
% Snew : Qoff of sig_new
% Tnew : Toff od sig_new
%
% bad_combo : flag if cannot make a good transformation
%
% op : Outparams (passed out as Outparams class)
%
%
% Examples:
% (1)
% Transform QRST template 'sig' with filename 'sig_name' and its associated 
% fiducial points 'Qon', 'Qoff', and 'Toff' by using the default 
% transformation parameters from Augparams
% (see documentation for Augparams for additional information on variables)
%
% [sig_new, Qnew, Snew, Tnew, bad_combo, op] = ...
%       transform_qrst(sig, Qon, Qoff, Toff, fname, Augparams(),1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start blank Outparams class
    op = Outparams();

% Bad combo flag initializes = 0
    bad_combo = 0;

% Reshape to row vectors if entered as column vectors
    if iscolumn(sig)
        sig = sig';
    end 
    
% Define "QRS" template as everything before Qoff (1:S)
% Define "TW" template as everything after Qoff+1 (S+1:end)
    sig_orig = sig;
    Q_orig = Q;
    S_orig = S;
    T_orig = T;

    QRS = sig(1:S);
    TW = sig(S+1:end);   
    TW_alone_fid = T-S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Further augment QRST complex with distotion

% *Augmentations*:
%   shift => move the baseline up/down
%   stretch => stretch/shrink QRS and T in time dimension
%   scale => stretch/shrink QRS and T amplitude
%   noise => add hf and LF noise to final augmented beat

% Main loop generates 'num_aug_beats' of transformed variants for each 
% QRST template that is input

% display("Start Transform")         % Debug
% Start with shifting entire QRST signal equally in Y axis
    new_shift = ap.shift_min + (ap.shift_max-ap.shift_min).*rand(1,1);
    QRS_shifted = QRS + new_shift;
    TW_shifted = TW + new_shift;
    
% xscales work better when restrict to 0.1 unit increments
% Generate QRS and T stretch/scale factors
% _1 = QRS
% _2 = TW
    xscale_1 = round(1+((ap.stretch_max-ap.stretch_min).*rand(1,1) + ap.stretch_min),1); 
    xscale_2 = round(1+((ap.stretch_max-ap.stretch_min).*rand(1,1) + ap.stretch_min),1);

    yscale_1 = round(1+((ap.scale_max-ap.scale_min).*rand(1,1) + ap.scale_min),2); 
    yscale_2 = round(1+((ap.scale_max-ap.scale_min).*rand(1,1) + ap.scale_min),2);     
    
% Multiply by yscale to create new QRS and TW
    QRS_scaled = QRS_shifted * yscale_1;
    TW_scaled = TW_shifted * yscale_2;

% Stretch QRS - interpolate based on new samples
    delta_1 = 1/xscale_1;
    interp_t_1 = 1:delta_1:length(QRS_scaled);
    QRS_stretched = interp1(QRS_scaled,interp_t_1,'spline');
 %display("QRS INTERP")
% Adjust fiducial points for stretch - floor to adjust for discrete samples
    Qnew = floor(Q_orig * xscale_1);
    Snew = floor(S_orig * xscale_1);

% Stretch TW
    delta_2 = 1/xscale_2;
    interp_t_2 = 1:delta_2:length(TW_scaled);
    TW_stretched = interp1(TW_scaled,interp_t_2,'spline');

% Adjust fiducial points for stretch
    Tnew = floor(TW_alone_fid*xscale_2);
    
% At this point, QRS is QRS_stretched, TW is TW_stretched, and 
% fiducial points are Qnew, Snew, and Tnew    

% Stretching in the Y axis changes the voltages of QRS(end) and TW(1).
% Need to now reinterpolate to re-join QRS and TW:

% Find where QRS_stretched is > TW_stretched(1) and take the final signal that satisfies this
% See recombine_qrst.m for additonal details of why calculate for signal
% and absolute value of the signal
    
    loc_ind1 = find(abs(QRS_stretched) > abs(TW_stretched(1)));
    loc_ind2 = find(QRS_stretched > TW_stretched(1));
        
    if isempty(loc_ind1)
        p_1 = NaN;    % Deal with empty possibility - set to NaN so will be ignored
    else
        p_1 = loc_ind1(end);
    end
        
    if isempty(loc_ind2)
        p_2 = NaN;    % Deal with empty possibility - set to NaN so will be ignored
    else
        p_2 = loc_ind2(end);
    end
  
    p = min([p_1 p_2]);
    
    if isnan(p)
       p = 1    % If p_1 and p_2 are NaN then set k2 = length(QRS_stretched) so it will fail filtering 
    end

% Open k2 sample window at end of QRS1 to fill with NaN for later interpolation
    k2 = length(QRS_stretched)-p+1;

% Set conditions for rejecting the combination:

if k2 > ap.k_tolerance || abs(QRS_stretched(end-k2) - TW_stretched(1)) > ap.abs_tolerance 
    
    % If bad combination set bad_combo = 1 and set output for signal and
    % fiducial points to 0
        bad_combo = 1;
        sig_new = 0;
        Qnew = 0;
        Snew = 0;
        Tnew = 0;

    % % Can display information if wanted on why combination failed
    % display("Bad Combination During Transformation!")
    % display(k2);
    % display(abs(QRS_stretched(end-k2) - TW_stretched(1)));
    
else
    
% If rejection criteria are not satisfied, will continue....    
% Replace signal to be interpolated with NaN
    QRS_stretched(end-k2:end) = NaN;

% Set QRS1(end) = TW2(start) as this sets the knot for interpolation
    QRS_stretched(end) = TW_stretched(1);

% Interpolate NaN values
    [QRS_stretched2,~] = fillmissing(QRS_stretched,'pchip');

% Do you want to flip QRS or T wave randomly (or for all cases)?
    if ap.random_qrsinvert == 1
            flip_qrs = 0;
                while flip_qrs == 0
                    flip_qrs = round(-1 + 2.*rand(1,1)); % Random -1 or +1 chosen
                end
        else
            flip_qrs = ap.invert_qrs;   % can manually set flip if you want (1 = no change)
    end
    
    if ap.random_twinvert == 1
            flip_tw = 0;
                while flip_tw == 0
                    flip_tw = round(-1 + 2.*rand(1,1)); % Random -1 or +1 chosen
                end
        else
            flip_tw = ap.invert_tw;   % can manually set flip if you want (1 = no change)
    end
    
% If flip TW or QRS have to make sure that the first sample is at 0,
% (rotate around y = 0) otherwise get really jagged ECGs.

% Will therfore shift the entire signal so that TW(1) = 0 or QRS(end) = 0 
% and then shift back to baseline after invert T wave or QRS complex

    shift_val_tw = TW_stretched(1);
    TW_new_stretched = TW_stretched - shift_val_tw;

    % Flip around y=0
    TW_new_stretched = flip_tw*TW_new_stretched;

    % reshift
    TW_stretched = TW_new_stretched + shift_val_tw;

% Can similarly flip QRS if you want:

    shift_val_qrs = QRS_stretched2(end);
    QRS_new_stretched = QRS_stretched2 - shift_val_qrs;

    % Flip around y=QRS(end)
    QRS_new_stretched = flip_qrs*QRS_new_stretched;

    % reshift
    QRS_stretched2 = QRS_new_stretched + shift_val_qrs;
           
% Join new QRS and other Twave
    sig_new = [QRS_stretched2 TW_stretched];
    Tnew = Tnew + Snew;
        
    
% Add random noise to signal
% If noise_min or noise_max = 0, do not add any noise

if ap.noise_min == 0 || ap.noise_max == 0
    % Do nothing
    noise = 0;
else
    noise = abs((ap.noise_max-ap.noise_min).*rand(1,1) + ap.noise_min);
    sig_new = sig_new + noise * rand(1, length(sig_new));
end
    

% Add LF noise to signal
% Generate low frequency cosine and alter freq, amplitude, and phase

% X axis smples
    lf_noise_t = 1:1:length(sig_new);

% Random values
    lf_amp = ap.lf_amp_min + (ap.lf_amp_max-ap.lf_amp_min).*rand(1,1);
    lf_offset = floor(ap.lf_offset_min + (ap.lf_offset_max-ap.lf_offset_min).*rand(1,1));
    lf_freq = ap.lf_freq_min + (ap.lf_freq_max-ap.lf_freq_min).*rand(1,1);

% LF noise signal    
    lf_noise_sig = lf_amp * (cos(2*pi*lf_freq/1000*(lf_noise_t-lf_offset)));
    sig_new = sig_new + lf_noise_sig;
    

end        % End loop that occurs only if have good beat to combine


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second checkpoint to reject augmented beat if QRS or QT outside of limits
% Set limits in samples as Augparams has user input QRS and QT limits in
% millisconds

samples = 1000/ap.freq;

if (Snew - Qnew > ap.qrs_max/samples | Snew - Qnew < ap.qrs_min/samples | Tnew - Qnew > ap.qt_max/samples | Tnew - Qnew < ap.qt_min/samples) | sig_new == 0     
    bad_combo = 1;
    sig_new = 0;
    Qnew = 0;
    Snew = 0;
    Tnew = 0;
    return
    % End here if checkpoint picks up a bad beat
else

% If figures checked off continue adding new signal to figures
if ap.figures
    
    figure('visible','off')
    
    p1 = plot(sig_orig, 'linewidth',2,'color', '[0 0 0]','linestyle','-');
    hold on
    xlabel("Samples")
    ylabel("mV")
    scatter(Q_orig,sig_orig(Q_orig), 60, 'k')
    scatter(S_orig,sig_orig(S_orig), 60, 'k', 'filled')
    scatter(T_orig,sig_orig(T_orig), 60, 'MarkerEdgeColor','[0 0 0]','MarkerFaceColor','[0.7 0.7 0.7]')
    
    if ~bad_combo
        title("Augmented Beat Comparison",'fontsize',12)
    else
        title("Augmented Beat Comparison - Bad Combination",'fontsize',12)
    end
    
    xlabel("Samples")
    ylabel("mV")
    p2 = plot(sig_new, 'linewidth',2,'color', '[1 0 0]','linestyle','-');
    scatter(Qnew,sig_new(Qnew), 60, 'k')
    scatter(Snew,sig_new(Snew), 60, 'k', 'filled')
    scatter(Tnew,sig_new(Tnew), 60, 'MarkerEdgeColor','[0 0 0]','MarkerFaceColor','[0.7 0.7 0.7]')

    %set(gcf,'Position',[100 100 800 650]);
    legend([p1 p2],{'Combined' 'Augmented'},'fontsize',12);

    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    %set(gcf,'name',strcat(fname,'_syn_',num2str(syn_num)),'numbertitle','off');

    saveas(gcf,fullfile(ap.output_folder,strcat(fname,'_syn_',num2str(syn_num),'.png')));

end        % End if figures


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put Outparams into Outparams class
    op.file1 = fname;
    op.yshift = new_shift;
    op.yscale_qrs = yscale_1;
    op.yscale_t = yscale_2;
    op.xscale_qrs = xscale_1;
    op.xscale_t = xscale_2;
    op.flip_qrs = flip_qrs;
    op.flip_tw = flip_tw;
    op.high_freq_amp = noise;
    op.low_freq_amp = lf_amp;
    op.low_freq_offset = lf_offset;
    op.low_freq = lf_freq;
    op.k2 = k2;
    op.abs_tolerance = ap.abs_tolerance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% display("End Transform")          % Debug

end        % End loop for second checkpoing using QRS and QT limits

	close all  % Close figures
    





