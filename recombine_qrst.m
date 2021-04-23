function [sig_new, Qnew, Snew, Tnew, bad_combo, op] = ...
    recombine_qrst(sig1, sig2, Q1, S1, T1, Q2, S2, T2, fname1, fname2, ap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- recombine_qrst()
% Version 1.0
%  
% Takes 2 QRST templates and annotations (Qon, Qoff, and Toff) and 
% performs recombination of QRS1 and TW2 to create a novel QRST template
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
% sig1 : first QRST template  (time series)
% sig2 : second QRST template (time series)
%
% Q1 : Qon of sig1
% S1 : Qoff of sig1
% T1 : Toff of sig1
%
% Q2 : Qon of sig2
% S2 : Qoff of sig2
% T2 : Toff of sig2
%
% fname1 : filename of sig1
% fname2 : filename of sig2
%
% ap : Augparams (passed in as Annoparams class)
%
%
% OUTPUT:
%
% sig_new : recombined signal
%
% Qnew : Qon of sig_new
% Snew : Qoff of sig_new
% Tnew : Toff od sig_new
%
% bad_combo : flag if cannot make a good recombination
%
% op : Outparams (passed out as Outparams class)
%
% Examples:
% (1)
% Combine QRS from 'sig1' which has filename 'sig_name1', and its associated 
% fiducial points 'Qon1', 'Qoff1', and 'Toff1' and the T wave of
% 'sig2' which has filename 'sig_name2' and its associated fiducial points 
% 'Qon2', 'Qoff2', and 'Toff2' by using the default parameters from Augparams
% (see documentation for Augparams for additional information on variables)
%
% [sig_new, Qnew, Snew, Tnew, bad_combo, op] = ...
%   recombine_qrst(sig1, sig2, Qon1, Qoff1, Toff1, Qon2, Qoff2, Tof2, sig_name1, sig_name2, Augparams())
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start blank Outparams class
    op = Outparams();

% Bad combo flag initializes = 0
    bad_combo = 0;

% Reshape to row vectors if entered as column vectors
    if iscolumn(sig1)
        sig1 = sig1';
    end
    if iscolumn(sig2)
        sig2 = sig2';
    end
   
% Default is to join QRS1 to TW2.  If flip_signals = 1 will instead join QRS2 to TW1
% If flip_signals = 1, swap sig1 and sig2
    if ap.flip_signals == 1
        sig1_old = sig1;
        sig2_old = sig2;
        Q1_old = Q1;
        Q2_old = Q2;
        S1_old = S1;
        S2_old = S2;
        T1_old = T1;
        T2_old = T2;

        sig1 = sig2_old;
        sig2 = sig1_old;
        Q1 = Q2_old;
        Q2 = Q1_old;
        S1 = S2_old;
        S2 = S1_old;
        T1 = T2_old;
        T2 = T1_old;
        
        clear sig2_old sig2_old
    end
    
% Define "QRS" as everything before Qend (1:S)
% Define "TW' as everything after Qend+1 (S+1:end)
    QRS1 = sig1(1:S1);
    QRS1_orig = QRS1;
    TW1 = sig1(S1+1:end);

    QRS2 = sig2(1:S2);
    QRS2_orig = QRS2;
    TW2 = sig2(S2+1:end);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Will join QRS1 to TW2

% To accomplish this without large discontinuities will need to interpolate 
% at the end of the QRS and link it to the start of the TW

% Want to avoid "wiggles" or "jumps" at the end of the QRS which occur if 
% the QRS1(end) is much lower/higher in voltage than TW2(1)

% Find where the end of QRS1 is > TW2(1) and take the final signal that 
% satisfies this condition.  This will prevent large discontinuities in the
% signal once combined if QRS(end) is < TW2(1)

% Because the end of QRS1 can be positive or negative, and greater than or
% less than TW2 sample 1, we look for where QRS1 is > TW2(1) and where 
% abs(QRS1) is > abs(TW2(1)) and take the smaller value

    loc_ind1 = find(abs(QRS1) > abs(TW2(1)));
    loc_ind2 = find(QRS1 > TW2(1));
    
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
       p = 1    % If p_1 and p_2 are NaN then set k1 = length(QRS1) so it will fail filtering 
    end

% Open k1 sample window at end of QRS1 and fill with NaN
% In this case, k1 = number of samples that QRS1 is < TW2(1)

    k1 = length(QRS1)-p+1;    

% if k1 is very large, it may not be a good merge because the start of TW2
% is at a significantly different voltage than the end of QRS1.  Likewise
% if there is a large voltage difference between QRS(end) and TW(1) there
% will be a large/rapid deflection to connect the QRS and TW and this is
% non-physiologic.

% Set condutions for rejecting the combination:
% Will use a large interpolation window (suggests that QRS1 is
% significantly below TW2, and if the difference between the end of the QRS
% and TW is > the user set absolute tolerance (in mV)

if k1 > ap.k_tolerance || abs(QRS1(end-k1) - TW2(1)) > ap.abs_tolerance 
    
    bad_combo = 1;      % Set this flag to know if combination failed

    % Zero output because no good combination
    sig_new = 0;
    Qnew = 0;
    Snew = 0;
    Tnew = 0;
    
    % % Can display information if wanted on why combination failed
    % display("Bad Combination of QRS1 and TW2!")
    % display(k1);
    % display(abs(QRS1(end-k1) - TW2(1)));
    
    return;    % Break out of function
       
else
    bad_combo = 0;      % Do not set flag equal to 1
end

% Replace signal to be interpolated with NaN
    QRS1(end-k1:end) = NaN;

% Set QRS1(end) = TW2(start) as this sets the knot for interpolation
    QRS1(end) = TW2(1);

% Interpolate NaN values using piecewise cubic hermite interpolating polynomial
% QRS with interpolation is now QRS1_new
    [QRS1_new,~] = fillmissing(QRS1,'pchip');

% Join new QRS and other Twave
    sig_new = [QRS1_new TW2];

% Update locations of fiducial points:
% Q1 and S1 do not change because we do not change anything related to Q
% and we know the location of S based on where the T wave starts (unchanged)
% Tnew is shifted forward by the difference between S2 and S1
    Qnew = Q1;
    Snew = S1;
    Tnew = S1+(T2-S2);

% Generate and save figures if wanted
if ap.figures

    % Specify limits of X and Y axes
    maxy = 1.1*max([max(QRS1) max(QRS2) max(TW1) max(TW2)]);
    miny = min([min(QRS1) min(QRS2) min(TW1) min(TW2)])-0.1;
    maxx = 1.1*max([length(sig1) length(sig2)]);

    figure('visible','off')

        subplot(1,3,1)
        plot(sig1, 'linewidth',1.5,'color', '[0.7 0.7 0.7]','linestyle','-')
        hold on
        thick = plot(QRS1_orig, 'linewidth',3,'color', 'k');
        ylim([miny maxy]);
        xlim([0 maxx])
        xlabel("Samples")
        ylabel("mV")
        ylim2 = ylim;
        scatter(Q1,sig1(Q1), 60, 'k')
        scatter(S1,sig1(S1), 60, 'k', 'filled')
        scatter(T1,sig1(T1), 60, 'MarkerEdgeColor','[0 0 0]','MarkerFaceColor','[0.7 0.7 0.7]')
        title("Beat 1 - QRS",'fontsize',12)

        subplot(1,3,2)
        plot(sig2, 'linewidth',1.5,'color', '[0.7 0.7 0.7]','linestyle','-');
        hold on
        plot([nan(1,length(QRS2)) TW2], 'linewidth',3,'color', 'k')
        ylim([miny maxy]);
        xlabel("Samples")
        ylabel("mV")
        ylim2 = ylim;
        xlim([0 maxx])
        scatter(Q2,sig2(Q2), 60, 'k')
        scatter(S2,sig2(S2), 60, 'k', 'filled')
        scatter(T2,sig2(T2), 60, 'MarkerEdgeColor','[0 0 0]','MarkerFaceColor','[0.7 0.7 0.7]')
        title("Beat 2 - T Wave",'fontsize',12)

        subplot(1,3,3)
        plot(sig_new, 'linewidth',3,'color', 'k','linestyle','-')
        hold on
        ylim([miny maxy]);
        xlabel("Samples")
        ylabel("mV")
        ylim2 = ylim;
        xlim([0 maxx])
        q =scatter(Q1,sig_new(Q1), 60, 'k');
        s = scatter(S1,sig_new(S1), 60, 'k', 'filled');
        t =scatter(Tnew,sig_new(Tnew), 60, 'MarkerEdgeColor','[0 0 0]','MarkerFaceColor','[0.7 0.7 0.7]');
        legend([q s t],{'Qon','Qoff' 'Toff'},'fontsize',12)
        title("New Combined Beat",'fontsize',12)

    sgtitle("Creating a New QRST Complex", 'fontsize',14,'fontweight', 'bold')

    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[100 100 1000 400]);
    
    saveas(gcf,fullfile(ap.output_folder,strcat(fname1,'__',fname2,'_recombined.png')));
    
    % Close figures automatically
    close(gcf)

end     % End if statement for generating figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outparams
    op.file1 = fname1;
    op.file2 = fname2;
    op.k1 = k1;
    op.abs_tolerance = ap.abs_tolerance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\








