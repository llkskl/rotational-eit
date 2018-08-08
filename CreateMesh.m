%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ img ] = CreateMesh( diam, elect_size, electrode_positions, ...
    max_mesh_elem_size, max_mesh_elem_size_rim, max_mesh_elem_size_elec, ...
    outerlayer_width_percent )
%MESHFUNC Returns EIDORS forward mesh with phantom for simulations or inverse mesh
%
% INPUT
%   diam:       diameter of the sample
%   elect_size: size of electrodes
%   electrode_positions:    positions of electrodes around the cylinder, in
%                           degrees
%   max_mesh_elem_size
%   max_mesh_elem_size_rim
%   max_mesh_elem_size_elec
%               :   maximum element sizes of sample, rim, and near electrode
%                   respectively
%   outerlayer_width_percent : width of the outer layer, in per cents of
%                   the radius of the sample
%   phantom:    name of the phantom if wanted
% OUTPUT
%   mesh (including phantom if requested)
%
% EXAMPLE INPUT
% diam = 3*3.14;
% elect_size = 0.5;
% electrode_positions = [0 20 40 60 180 200 220 240];
% max_mesh_elem_size = 1;
% max_mesh_elem_size_rim = 1;
% max_mesh_elem_size_elec = 0.2;
% outerlayer_width_percent = 0.05;
% phantom = 'one_brick';

% conversions
radius = diam/2;
sample_r = radius*(1-outerlayer_width_percent);

% 0 for 2d circular, radius, max size of mesh elems
mdl_shape = [0,radius,max_mesh_elem_size_rim];
% degrees of position, z level
elect_pos = horzcat(electrode_positions', zeros(numel(electrode_positions),1));
% shape of circular elecs
% radius, 0 (for circular), max mesh size near electrode
elect_shape = [elect_size, 0, max_mesh_elem_size_elec];

extra={'ball',sprintf('solid ball = sphere(0,0,0;%f) and orthobrick(-%f,-%f,0;%f,%f,0) -maxh=%f;',sample_r,sample_r,sample_r,sample_r,sample_r,max_mesh_elem_size)}; %-maxh=0.5
fmdl= ng_mk_cyl_models(mdl_shape,elect_pos,elect_shape,extra);
img = mk_image( fmdl, 0);


end

