%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ M ] = CreateRotationMatrixCW( eidors_img, thetas )
%CREATEROTATIONMATRIXCW creates rotation matrix for clockwise rotation.

M = CreateRotationMatrix( eidors_img , -thetas );

end

