%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- rhythmstrip_script
% Version 1.0
%  
% Loads QRST signals and fidicual points and puts data into correct format
% for then passing into the rhythmstrip() function which generates random,
% annotated rhythm strips.  Sections between %----% will likely 
% need to be modified based on your specific files and data structure.
%
% The script loops ('num_strips' times) to generate multiple rhythm strips 
% as needed.  Adjust Augparams between loop iterations to generate 
% additional diversity.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% Number of different rhythm strips to create
    num_strips = 4;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Create Augparams classes with default values.  Edit as needed
    ap = Augparams();   
    
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

for k = 1:num_strips         % Loop num_strips times

% Loop to catch errors and continue on
    loopflag = 1;        
    while loopflag == 1
        
% Specify a different Augparams for different iterations if wanted
% Could create an array of Augparams and load in different ones for each
% loop here
    % ap = Augparams();
        
% Reset sigs and fiducial points for each loop
    clearvars -except loopflag num_strips ap file_list num_files k i
   
    try       
     
% n_morph : Total number of QRST morphologies to include in the rhythm strip.
% Choose if want to specify manually or choose random number (if >9
% will only see morphology 1 due to how morpholgoies are chosen based on
% sig_ratio)
    
    if isempty(ap.n_morph)
       n_morph = round(1 + (ap.n_morph_max-1).*rand(1,1)); % Random number 1-n_morph_max
    else
       n_morph = ap.n_morph;
    end
    
% n_var : Total number of variations per QRST morphology to include
% Enter as vector eg [2 5 1 3] -- this will give 2 variations of 
% morphology 1, 5 variations of morphology 2, 1 variation of morphology 3 etc.

% If n_var is empty the script will randomly generate the appropriate 
% number of values (based on n_morph) with range 1-5 (see below)
    n_var = ap.n_var;

% sig_ratio : Ratio of how prevalent each morphology is in the rhythmstrip  
% Enter as vector eg [0.7 0.2 0.1] that sums to 1 -- this will give you 70%
% morphology 1, 20% morphology 2, and 10% morphology 3    

% If sig_ratio is empty the script will randomly generate appropriate number 
% of values (based on n_morph) such that the first signal is most frequent, 
% and subsequent signals are less frequent (see below)
    sig_ratio = ap.sig_ratio;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate n new random numbers between 1 and n_morph to use
% for the different beats to be included in rhythm strip

% Use random QRST templates
    ind = zeros(1,n_morph);
    for i = 1:n_morph 
        ind(i) = round(1 + (num_files-1).*rand(1,1));     
    end
    
% Or can manually specify files based on their index (ind) after loading
    % ind = [1 2 4];          
    % ind                   % Uncomment if want to see indices chosen
    
% Extract filenames for loading signals
    for i = 1:length(ind) 
        filename{i} = file_list{ind(i)};
        [~,fname{i},~] = fileparts(filename{i}); 
    end
 
% Note: unlike with ecgaug_script, this script as provieded does NOT prevent 
% reusing the same combinations of QRST morphologies.  Modify if needed.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load ECGs/fiducial points in and save original signals/fiducial points

%-------------------------------------------------------------------------%
% In this section we have stored the signals and fiducial points in a
% structure called vm_sig where the signal is vm_sig.vm, and the fiducial
% points are vm_sig.Q, vm_sig.S, and vm_sig.Tend -- Modify as needed!

for i = 1:n_morph 
    load(filename{i})
    sig{i} = vm_sig.vm;
    Q_orig(i) = vm_sig.Q;
    S_orig(i) = vm_sig.S;
    Tend_orig(i) = vm_sig.Tend; 
end
%-------------------------------------------------------------------------%

% Rename fiducial points for use and save old fiducial points for use later
for i = 1:n_morph 
    Q(i) = Q_orig(i);
    S(i) = S_orig(i);
    T(i) = Tend_orig(i);
end

% Adjust minRR based on signals that are loaded
    for i = 1:length(sig)
        sig_len(i) = length(sig{i});
    end
    ap.minRR = max(sig_len)*(1000/ap.freq); 

% fprintf("%i QRST morphologies loaded",n_morph)    % Uncomment to view

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
% If n_var is empty or does not exist, will choose random integers 
% between 1 and 5 for each signal morphology's number of variations 
if ~exist('n_var','var') || isempty(n_var) 

    n_var = zeros(1,length(sig));

    for i = 1:length(sig)
       n_var(i) = round(1 + (5-1).*rand(1,1)); 
    end 
end

if length(n_var) ~= length(sig)
   error("n_var must be same length as number of signals (n_morph)");
end

% If sig_ratio is empty or does not exist, will choose random numbers 
% between 0 and 1 for each signal morphology, with the first signal never
% being < 0.5 (presumed to be dominant beat morphology)
% Note: sum of sig_ratio must = 1

if ~exist('sig_ratio','var') || isempty(sig_ratio) 

    if length(sig) == 1
        sig_ratio(1) = 1;

    else
        % Choose number between 0.5 and 0.9 for sig_ratio(1)
        sig_ratio(1) = round(0.5 + (0.9-0.5).*rand(1,1),1);

        r = rand(1,length(sig)-1);
        r = ((r * (1-sig_ratio(1))) / sum(r));

        sig_ratio = [sig_ratio r];
    end 
end


% view filenames, n_var, and sig_ratio if wanted.  Comment out if not wanted.
    %fname
    %n_var
    %sig_ratio

% Date/time for labeling files
    dt = datestr(now, 'yyyymmdd_HHMMSS');
    
% Call rhythmstrip function.   
    [rhythm_strip, Qloc, Sloc, Tloc, m, b] = rhythmstrip(sig, Q, S, T, n_var, sig_ratio, ap, dt);
    loopflag = 0;
    
    catch
        loopflag = 1;   % Restart loop to generate rhythmstrip if error
        % display("Rhythmstrip Creation Failed")
    end       % end try/catch loop
    
end    % end main for loop

% Save data
% Saves with date/time as filename because would be too unweidly to have
% filename made up of every possible beat in the rhythm strip
    if ap.save_data
        export_filename = fullfile(ap.output_folder,strcat('rhythm_syn_',dt,'_.mat'));
        rhythm = struct('rhythm', rhythm_strip, 'n_morph', length(n_var), 'n_variants',n_var, 'Qloc', Qloc, 'Sloc', Sloc, 'Tendloc', Tloc, ...
            'morph', m, 'beat', b, 'ap', ap, 'files',string(fname'), 'date_time', dt);
        save(export_filename,'rhythm');
    end

end       % End main while loop

disp(" ")
fprintf('Complete! Generated %i rhythm strips\n', num_strips)
fprintf('Output directory = %s\n', ap.output_folder)


