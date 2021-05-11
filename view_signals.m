%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECGAUG -- view_signals
% Version 1.0
%  
% Views template signals provided with ECGAug software
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Folder where signals are present
input_folder = 'C:\ECGAug\templates';

% Load list of files in load directory
file_list_struct = dir(fullfile(input_folder, '*.mat'));
num_files = length(file_list_struct);
file_list = cell(num_files, 1);
for i = 1:num_files
    file_list{i} = fullfile(file_list_struct(i).folder, file_list_struct(i).name);
end
orig_directory = file_list_struct(i).folder;
file_list = sort(file_list');

for i=1:length(file_list)
    load(file_list{i}) 
    sig{i} = vm_sig.vm;
    Q(i) = vm_sig.Q;
    S(i) = vm_sig.S;
    Tend(i) = vm_sig.Tend;
end

% Number of figures with 5x5 subplots
pgs = ceil(length(file_list)/25)

% Plot
for j=1:pgs
    if j<pgs
    figure(j)
        for i=1+(25*(j-1)):25+(25*(j-1))
            subplot(5,5,i-(25*(j-1)))
            plot(sig{i})
            hold on
            scatter(Q(i),sig{i}(Q(i)))
            scatter(S(i),sig{i}(S(i)))
            scatter(Tend(i),sig{i}(Tend(i)))
            line([0 length(sig{i})],[0 0])
            [~,fname,~] = fileparts(file_list{i});
            title(fname)
        end
    else
        
    remainder =  length(file_list)-(25*(pgs-1))
    figure(j)   
        for i=1+(25*(j-1)):remainder+(25*(j-1))
            subplot(5,5,i-(25*(j-1)))
            plot(sig{i})
            hold on
            scatter(Q(i),sig{i}(Q(i)))
            scatter(S(i),sig{i}(S(i)))
            scatter(Tend(i),sig{i}(Tend(i)))
            line([0 length(sig{i})],[0 0])
            [~,fname,~] = fileparts(file_list{i});
            title(fname)
        end
    end
end




