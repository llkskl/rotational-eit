%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Finite Element Mapping for
% Efficient Image Reconstruction in Rotational Electrical Impedance 
% Tomography".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ data ] = ClickingRotationalMeasurement( eidors_img, n_rotational_positions )
%ROTATIONALEITMEASUREMENT Returns rotational measurement data for given
%EIDORS image object.
%
% Returned data is matrix of NxM where N is the amount measurements in
% total for all stimulations and M = numel(thetas).
%
% Stimulation patterns are to be set in eidors_img

% get the initial data
data1 = fwd_solve(eidors_img);
data1 = data1.meas;

% allocate
data = zeros(numel(data1),n_rotational_positions);
data(:,1) = data1;

% rotational procedure
for i = 2:n_rotational_positions
    % rotate
    rotated_img = ClickRotationalMesh(eidors_img,i);
    % measure
    data2 = fwd_solve(rotated_img);
    % save
    data(:,i) = data2.meas;

end

end

