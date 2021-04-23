%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- ecgaug_script
% Version 1.0
%  
% Script for loading files and running recombine_qrst() and transform_qrst()
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
% NOTE; you may need to modify the way signals and fiducial
% points are extracted from files, as the format in this function may not
% match the format for your files.  The section where files are loaded 
% is between %----% markers.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Load parameters from Augparams class
    ap = Augparams();

% NUMBER OF RECOMBINED & AUGMENTED BEATS TO CREATE
% Total number of combined ECGs to generate
    tot_aug = 10;

% Number of augmented beats per combined ECG
    num_aug_beats = 3;

% Total number of ECGs = tot_aug * num_aug_beats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% Load list of files in load directory
    file_list_struct = dir(fullfile(ap.input_folder, ap.extension));
    num_files = length(file_list_struct);
    file_list = cell(num_files, 1);
    for i = 1:num_files
        file_list{i} = fullfile(file_list_struct(i).folder, file_list_struct(i).name);
    end
    orig_directory = file_list_struct(i).folder;
    file_list = file_list';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables to keep track of which ECGs have been used
    ind_list = [ -1 -1 ];
   
% Reset indexing to deal with bad combinations
    i = 1;

while i <= tot_aug              % Generate tot_aug number of combined ECGs

    % display(i)    % Uncomment this to follow progress if you want
    
    % Choose 2 random files
    % This first selection of random files wont be used, but is needed to 
    % initialize how account for avoiding repeats
    ind1 = round(1 + (num_files-1).*rand(1,1));
    ind2 = round(1 + (num_files-1).*rand(1,1)); 
    
        % Set ind1 and ind2 here if you want to use specific files
		% If want to go through all files in sequence, can add code that
		% includes something like ind1 = i etc.
       
        %%%%%%%%%%%%%
        % ind1 = 1;
        % ind2 = 2;
        %%%%%%%%%%%%%
    
    % Loop to check if the combination of ind1 *AND* ind2 has already been used
    % and if so choose another ECG to avoid over representation of a single ECG in
    % the final augmented dataset
    
    w = 0;
    while sum(ismember(ind_list,[ind1 ind2],'rows'))
    
    if w > 50
        display("No Combinations Found");
        return;
    end
        
        
    % display("Combination Already Used!")

        % Generate 2 new random numbers between 1 and number of signal files
        ind1 = round(1 + (num_files-1).*rand(1,1));
        ind2 = round(1 + (num_files-1).*rand(1,1));    
    
        % Set ind1 and ind2 here if you want to use specific files
		% If want to go through all files in sequence, can add code that
		% includes something like ind1 = i etc.
       
        %%%%%%%%%%%%%
        % ind1 = 1;
        % ind2 = 2;
        %%%%%%%%%%%%%
        
     % This loop continues until the combination of ind1 and ind2 has
     % not already been used.  If a novel combination of ind1 and ind2
     % is found, it breaks out of thile while loop
     
     w = w + 1; 
     
     end
   
% Keep tally of which indices have been chosen to avoid duplication
    ind_list = [ind_list ; ind1 ind2];
    
% Pull filenames out of file list
    filename1 = file_list{ind1};
    filename2 = file_list{ind2};

    [~,fname1,~] = fileparts(filename1);
    [~,fname2,~] = fileparts(filename2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code section should be modified as needed to extract the
% signals and fiducial points from the relevant data files.  Will need to
% modify to adjust for different variable naming conventions

%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this code we have stored the signals and fiducial points in a
% structure called vm_sig where the signal is vm_sig.vm, and the fiducial
% points are vm_sig.Q, vm_sig.S, and vm_sig.Tend -- Modify as needed!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load both ECGs/fiducial points in and save original signals/fiducial points
load(filename1)
    sig1 = vm_sig.vm;
    sig1_orig = sig1;
    Q_orig1 = vm_sig.Q;
    S_orig1 = vm_sig.S;
    Tend_orig1 = vm_sig.Tend;

load(filename2)
    sig2 = vm_sig.vm;
    sig2_orig = sig2;
    Q_orig2 = vm_sig.Q;
    S_orig2 = vm_sig.S;
    Tend_orig2 = vm_sig.Tend;
%-------------------------------------------------------------------------%
    
% Rename fiducial points for use and save old fiducial points
    Q1 = Q_orig1;
    S1 = S_orig1;
    T1 = Tend_orig1;

    Q2 = Q_orig2;
    S2 = S_orig2;
    T2 = Tend_orig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust what is saved during pass through recombine_qrst()
% ap.figures = 1;
% ap.save_data = 0;

% Pass filenames and augmentation parameters into generative function

[sig_new, Qnew, Snew, Tnew, bad_combo, op] = recombine_qrst(sig1, sig2, Q1, S1, T1, Q2, S2, T2, fname1, fname2, ap);

if bad_combo == 0
    i = i+1; 
 
    fname = strcat(fname1,'__',fname2);
    
    j = 1;
	c = 0;
    
    while j <= num_aug_beats  
	
	if c > ap.recombolimit
	    % break out and go to next recombination if cant find any combination that works 
		j = num_aug_beats + 1	
	end
    
        [sig_new2, Qnew2, Snew2, Tnew2, bad_combo2, op2] = transform_qrst(sig_new, Qnew, Snew, Tnew, fname, ap, j);

    % Retry if bad combination
    
        if bad_combo2 == 1
			c = c + 1;
            % display("Bad Transformation")         % Debug
        else
        % Export signal and other information
            if ap.save_data
                export_filename1 = fullfile(ap.output_folder,strcat(fname1,'__',fname2,'_recombined.mat')); 
                recombo = struct('sig_new', sig_new, 'Q', Qnew, 'S', Snew, 'Tend', Tnew, 'ap', ap, 'op', op, 'fname1', fname1, 'fname2', fname2);    
                save(export_filename1,'recombo');
                
                export_filename2 = fullfile(ap.output_folder,strcat(fname,'_syn_',num2str(j),'.mat'));
                sig = struct('file1', fname, 'sig',sig_new2, 'Q', Qnew2, 'S',Snew2, 'Tend', Tnew2, 'ap', ap, 'op', op2);
                save(export_filename2,'sig');
            end

            j = j+1;    % increment successful transformation counter
        end
    end
end
    
end
    
disp(" ")
fprintf('Complete! Generated %i combined ECGs each with %i augmentations = %i Total ECGs\n', i-1, num_aug_beats, (i-1)*num_aug_beats)
fprintf('Output directory = %s\n', ap.output_folder)

close all
    
