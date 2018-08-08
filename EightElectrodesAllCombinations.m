%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stim = EightElectrodesAllCombinations()
% all possible stimulus and measurement combinations in an eight electrode
% configuration

% initialize stim2
stim2.stimulation(1).stimulation = 'Amps';
stim2.stimulation(1).stim_pattern = [];
stim2.stimulation(1).meas_pattern = [];

% loop through all possible combinations
idx = 1;
for i = 1:8
    for k = 1:8
        if i == k
            % do nothing
        else
            stim2.stimulation(idx).stimulation = 'Amps';
            % stimulate current from i to k
            stim2.stimulation(idx).stim_pattern = stimulates(i,k,8);
            % measure all pairs from remaining (not i and k) electrodes
            for j = 1:8
                if j ~= i && j ~= k
                    % all eight electrodes
                    pairs = 1:8;
                    % remove i,k,j to get the list of agains which is
                    % measured
                    pairs = pairs(pairs ~= i & pairs ~= j & pairs ~= k);
                    stim2.stimulation(idx).meas_pattern = ...
                        [ stim2.stimulation(idx).meas_pattern;
                          measures(j,pairs,8)];
                end
            end
            idx = idx + 1;
        end
    end
end
stim = stim2.stimulation;

end

% create a stimulation vector from pos to neg with model that has n
% electrodes. This is a wrapper for patternmatrix()
function s = stimulates(pos_electrode,neg_electrodes,n_electrodes)
    if numel(pos_electrode) ~= 1
        error('Number of positive electrodes not equal to 1.');
    end
    if numel(neg_electrodes) ~= 1
        error('Number of negative electrodes not equal to 1.');
    end
    s = patternmatrix(pos_electrode,neg_electrodes,n_electrodes)';
end

% create a measurement vector from pos to neg with model that has n
% electrodes. This is a wrapper for patternmatrix()
function m = measures(pos_electrode,neg_electrodes,n_electrodes)
    if numel(pos_electrode) ~= 1
        error('Number of positive electrodes not equal to 1.');
    end
    m = patternmatrix(pos_electrode,neg_electrodes,n_electrodes);
end


% Sub-function that creates an EIDORS pattern matrix for either stimulation
% of measurement pattern.
function pm = patternmatrix(pos_electrode,neg_electrodes,n_electrodes)
    if pos_electrode == 0
        pos_electrode = n_electrodes;
        %warning('pos_electrode = 0 parameter given, setting to n_electrodes = %i.\n',n_electrodes);
    end
    if sum(neg_electrodes == 0) > 0
        neg_electrodes(neg_electrodes == 0) = n_electrodes;
        %warning('neg_electrode = 0 parameter given, setting to n_electrodes = %i.\n',n_electrodes);
    end
    pm = zeros(numel(neg_electrodes),n_electrodes);
    pm(:,pos_electrode) = 1;
    for i = 1:numel(neg_electrodes)
        pm(i,neg_electrodes(i)) = -1;
    end
    pm = sparse(pm);
end