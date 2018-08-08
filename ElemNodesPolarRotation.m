%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ nodesrotated ] = ElemNodesPolarRotation( Nxy, theta )
%ELEMNODESPOLARROTATION 
% Rotates given element nodes Nxy with given angle theta using polar
% coordinates transform

    %   first to polar coordinates
    %       radius
    R = sqrt(Nxy(:,1).^2+Nxy(:,2).^2);
    %       angle (note that we need to check 90 deg and above 180 deg
    %       separately, we get these respective to x-coordinate)
    iposx = Nxy(:,1) > 0;
    inegx = Nxy(:,1) < 0;
    izerox = Nxy(:,1) == 0;
    alpha = zeros(size(Nxy,1),1);
    alpha(iposx) = atand(Nxy(iposx,2)./Nxy(iposx,1));
    alpha(inegx) = 180+atand(Nxy(inegx,2)./Nxy(inegx,1));
    alpha(izerox) = sign(Nxy(izerox,2))*90;
    % now rotate
    nodesrotated = [R.*cosd(alpha+theta) R.*sind(alpha+theta)];

end

