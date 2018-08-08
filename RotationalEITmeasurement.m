%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ data ] = RotationalEITmeasurement( eidors_img, thetas )
%ROTATIONALEITMEASUREMENT Returns rotational measurement data for given
%EIDORS image object with given angles (thetas) and stimulation patterns
%(EIDORS structure).
%
% Returned data is matrix of NxM where N is the amount measurements in
% total for all stimulations and M = numel(thetas).
%
% Stimulation patterns are to be set in eidors_img

% get the initial data
data1 = fwd_solve(eidors_img);
data1 = data1.meas;

% allocate
data = zeros(numel(data1),numel(thetas));
data(:,1) = data1;

% rotational procedure
for i = 2:numel(thetas)
    % rotate
    rotated_img = rotationalmapping(eidors_img,thetas(i));
    % measure
    data2 = fwd_solve(rotated_img);
    % save
    data(:,i) = data2.meas;

end

end


function [ rotatedimg ] = rotationalmapping( meshimg, theta )
%ROTATIONALMAPPING Rotates the given EIDORS mesh data by given angle theta.
%The mesh elements stay as they are, only data in the elements is rotated.

rotatedimg = meshimg;
rotatedimg.elem_data(:) = 0;

% rotate mesh
for i = 1:size(meshimg.fwd_model.elems,1)
    % get nodes
    iN1 = meshimg.fwd_model.elems(i,1);
    iN2 = meshimg.fwd_model.elems(i,2);
    iN3 = meshimg.fwd_model.elems(i,3);
    Nxy = [meshimg.fwd_model.nodes(iN1,:);
        meshimg.fwd_model.nodes(iN2,:);
        meshimg.fwd_model.nodes(iN3,:)];
    
    % rotate nodes
    Nrotated = ElemNodesPolarRotation(Nxy,theta);
    
    % get elem weights under rotation
    %   first define the new elem as polygon
    Npoly = [Nrotated; Nrotated(1,:)];
    %   define the membership selection function
    %   here points lying on edges of polygon are defined out, need to use
    %   addition on output of inpolygon if desired otherwise (see MATLAB
    %   doc)
    select_fcn = @(x,y,z) inpolygon(x,y,Npoly(:,1),Npoly(:,2));
    %   now get the elem membership element weights
    memb_weights = elem_select( meshimg.fwd_model, select_fcn );
    
    % sum elem data to new elem with weights
    rotatedimg.elem_data = plus(rotatedimg.elem_data, ...
        memb_weights.*meshimg.elem_data(i));
    
end


end


