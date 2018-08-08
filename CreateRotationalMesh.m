%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Finite Element Mapping for
% Efficient Image Reconstruction in Rotational Electrical Impedance 
% Tomography".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rotationalMesh ] = CreateRotationalMesh( n_rotpos, electrode_size, ...
    electrode_centres, surface_r, object_r, max_vertex_length)
% This function produces rotational mesh that has a rotationally invariant
% boundary, such as presented in the paper (number of rotationally
% invariant positions is given with n_rotpos and the radius of the
% rotational boundary with object_r). Total domain radius is given with
% surface_r.

% Constraction in the further code is performed using polar coordinates
% (radius, alpha).

% in-function modified inputs
electrode_centres_polarcoord = sort(mod(90 - electrode_centres,360));
EleNodesInPolarCoord_Ralpha = [];

% Insert electrodes and in-between boundary nodes.

for e = 1:numel(electrode_centres)
    % Centre node
    EleNodesInPolarCoord_Ralpha = [EleNodesInPolarCoord_Ralpha; ...
                                surface_r electrode_centres_polarcoord(e)];
    % Angular difference between edge and the center, number of nodes
    % in between and the nodes
    HalfElecAlpha = ( electrode_size / (2*2*pi*surface_r) ) * 360;
    NHalfElecNodes = ceil(HalfElecAlpha/MaxAngleAtR(surface_r,max_vertex_length));
    HalfElecAlphaCoords = linspace(0,HalfElecAlpha,NHalfElecNodes+1);
    HalfElecAlphaCoords(1) = []; % centre coordinate
    size_before_addition = size(EleNodesInPolarCoord_Ralpha,1);
    EleNodesInPolarCoord_Ralpha = [EleNodesInPolarCoord_Ralpha; ...
                                repmat(surface_r,NHalfElecNodes,1) electrode_centres_polarcoord(e)-HalfElecAlphaCoords'; ...
                                repmat(surface_r,NHalfElecNodes,1) electrode_centres_polarcoord(e)+HalfElecAlphaCoords' ];
	electrodeNodes(e).nodes = size_before_addition:size(EleNodesInPolarCoord_Ralpha,1);
    electrodeNodes(e).z_contact = 1;
	% nodes between electrodes
    between_begin = electrode_centres_polarcoord(e) + HalfElecAlpha;
    if e == numel(electrode_centres)
        between_end = electrode_centres_polarcoord(1) - HalfElecAlpha;
    else
        between_end = electrode_centres_polarcoord(e+1) - HalfElecAlpha;
    end
    n_between_nodes = ceil(mod(between_end-between_begin,360)/MaxAngleAtR(surface_r,max_vertex_length));
    % create nodes and check if they go over 0 alpha
    if between_begin > between_end; c = 360; else; c = 0; end
    between_nodes = linspace(between_begin,between_end+c,n_between_nodes);
    % remove first and last, they are electrode nodes
    between_nodes(end) = []; between_nodes(1) = [];
    % add nodes
    EleNodesInPolarCoord_Ralpha = [EleNodesInPolarCoord_Ralpha; ...
                                repmat(surface_r,numel(between_nodes),1) between_nodes']; 
	
end

% Insert nodes between rotational boundary and domain boundary.

points_between = [];
% between boundaries nodes, minimum one layer of nodes
r_steps_between = linspace(object_r, surface_r, max(3,ceil((surface_r - object_r)/max_vertex_length)+1));
% remove first and last r_steps, they are already included elsewhere
r_steps_between(1) = []; r_steps_between(end) = [];
for i = 1:numel(r_steps_between)
    n_nodes = floor(2*pi*r_steps_between(i)/max_vertex_length)+1;
    alphas = linspace(0,360,n_nodes+1);
    alphas(end) = [];
    if numel(alphas) > 2 % mod(i,2) == 1 && 
        alphas = alphas + abs(alphas(2)-alphas(1))/2;
    end
    points_between = [points_between;
              ones(n_nodes,1)*r_steps_between(i) alphas'];
end

% insert nodes on rotational boundary and inside it

points = [];
rotational_boundary = [];
% korjaa side length t‰sm‰‰m‰‰n electrodijaon sattumaa sek‰ object_r eik‰
% surface_r
%   lis‰ksi etsitt‰v‰ kirjallisuudesta, voiko mesh olla liian tasainen? vai
%   onks se just hyv‰ tasaisena
r_steps = linspace(0,object_r,ceil(object_r/max_vertex_length));
% TODO: r steps lis‰tt‰v‰ min edellisest‰ ja thetas jaosta
for i = 1:numel(r_steps)
    if i < numel(r_steps)
        n_nodes = floor(2*pi*r_steps(i)/max_vertex_length)+1;
        alphas = linspace(0,360,n_nodes+1);
        alphas(end) = [];
        if numel(alphas) > 2 % mod(i,2) == 1 && 
            alphas = alphas + abs(alphas(2)-alphas(1))/2;
        end
        points = [points;
                  ones(n_nodes,1)*r_steps(i) alphas'];
    else % i == numel(r_steps)
        boundary_refinement = floor((2*pi*r_steps(i)/n_rotpos)/max_vertex_length);
        n_nodes_boundary = n_rotpos*(boundary_refinement+1);
        alphas = linspace(0,360,n_nodes_boundary+1);
        alphas(end) = [];
        rotational_boundary = [ones(n_nodes_boundary,1)*r_steps(i) alphas'];
    end
end


DT = delaunayTriangulation(Polar2XY([EleNodesInPolarCoord_Ralpha; points_between; rotational_boundary; points]));

% rotational parameters
rotparams.n_rotational_positions = n_rotpos;
rotparams.outernodes = ismember(DT.Points,Polar2XY([EleNodesInPolarCoord_Ralpha; points_between]),'rows');
rotparams.rotational_boundary = ismember(DT.Points,Polar2XY(rotational_boundary),'rows');
rotparams.innernodes = ismember(DT.Points,Polar2XY(points),'rows');

% Start the EIDORS img
rotationalMesh = [];
rotationalMesh.type = 'image';
rotationalMesh.name = 'Rotational mesh';
rotationalMesh.fwd_model.type = 'fwd_model';
rotationalMesh.fwd_model.name = 'Rotational mesh';
rotationalMesh.fwd_model.nodes = DT.Points;
rotationalMesh.fwd_model.elems = DT.ConnectivityList;
rotationalMesh.fwd_model.electrode = electrodeNodes;
rotationalMesh.rotational_parameters = rotparams;
rotationalMesh.elem_data = zeros(size(rotationalMesh.fwd_model.elems,1),1);

rotationalMesh.fwd_model.solve = 'eidors_default';
rotationalMesh.fwd_model.jacobian = 'eidors_default';
rotationalMesh.fwd_model.system_mat = 'eidors_default';
rotationalMesh.fwd_model.normalize_measurement = 0;
rotationalMesh.fwd_model.gnd_node = 1;

end


%** SUBFUNCTIONS


function alpha = MaxAngleAtR(R, max_vertex_length)
    % Function returns the max element side length in polar coordinates.
    alpha = asind( (max_vertex_length/2) / R ) * 2;
end

function xy = Polar2XY(Ralpha)
    % Polar coordinates to cartesian coordinates
    xy = Ralpha(:,1).*[cosd(Ralpha(:,2)) sind(Ralpha(:,2))];
end


