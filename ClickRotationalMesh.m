%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Finite Element Mapping for
% Efficient Image Reconstruction in Rotational Electrical Impedance 
% Tomography".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EIDORS_img = ClickRotationalMesh(rEIDORS_img, click_steps)
% Rotates the mesh for a number of steps defined in the click_steps
% argument. The total number of possible steps is defined in the mesh
% construction. For click_steps values higher than the total, result will
% be modulo of the total possibilities.
%
% click_steps == 0 -> no rotation

EIDORS_img = rEIDORS_img;

thetas = linspace(0,360,rEIDORS_img.rotational_parameters.n_rotational_positions+1);
thetas(end) = [];
theta = thetas(mod(click_steps,numel(thetas))+1);

% [x' y']' = [cos(t) -sin(t); sin(t) cos(t)]*[x y]'

R = [cosd(theta) -sind(theta);
     sind(theta) cosd(theta)];
 
idx_rotnodes = rEIDORS_img.rotational_parameters.innernodes | rEIDORS_img.rotational_parameters.rotational_boundary;
rotated_nodes = R*rEIDORS_img.fwd_model.nodes(idx_rotnodes,:)';
EIDORS_img.fwd_model.nodes(idx_rotnodes,:) =  rotated_nodes';

idx_boundary_nodes = rEIDORS_img.rotational_parameters.rotational_boundary;
% K = (click_steps+1)*sum(idx_boundary_nodes)/rEIDORS_img.rotational_parameters.n_rotational_positions;
boundary_nodes = 1:size(rEIDORS_img.fwd_model.nodes,1);
boundary_nodes = find(boundary_nodes'.*idx_boundary_nodes);
% boundary_mapping = circshift(boundary_nodes,K);
boundary_mapping = CreateBoundaryMapping(rEIDORS_img,EIDORS_img,boundary_nodes);

idx_outer_nodes = rEIDORS_img.rotational_parameters.outernodes;
outer_nodes = 1:size(rEIDORS_img.fwd_model.nodes,1);
outer_nodes = find(outer_nodes'.*idx_outer_nodes);


% for each segment
% line = 1;
for line = 1:size(EIDORS_img.fwd_model.elems,1)
    % check if one end outer node and the other boundary node
    % (1,2)
    node1 = EIDORS_img.fwd_model.elems(line,1);
    node2 = EIDORS_img.fwd_model.elems(line,2);
    node3 = EIDORS_img.fwd_model.elems(line,3);
    [b12, b_idx12,b_map_idx12] = isOuterAndBoundary(node1,node2,outer_nodes,boundary_nodes);
    [b23, b_idx23,b_map_idx23] = isOuterAndBoundary(node2,node3,outer_nodes,boundary_nodes);
    [b13, b_idx13,b_map_idx13] = isOuterAndBoundary(node1,node3,outer_nodes,boundary_nodes);
    if b12
%         EIDORS_img.fwd_model.elems(line,:) = [];
%         EIDORS_img.elem_data(line) = [];
        EIDORS_img.fwd_model.elems(line,b_idx12) = boundary_mapping(b_map_idx12);
    end
    if b23
%         EIDORS_img.fwd_model.elems(line,:) = [];
%         EIDORS_img.elem_data(line) = [];
        EIDORS_img.fwd_model.elems(line,b_idx23+1) = boundary_mapping(b_map_idx23);
    end
    if b13
%         EIDORS_img.fwd_model.elems(line,:) = [];
%         EIDORS_img.elem_data(line) = [];
        if b_idx13 == 2; b_idx13 = 3; end
        EIDORS_img.fwd_model.elems(line,b_idx13) = boundary_mapping(b_map_idx13);
    end
%     line = line + 1;
    
end


end

function boundary_mapping = CreateBoundaryMapping(orig_img, rot_img, boundary_nodes)
    % for each node on boundary
    for i = 1:numel(boundary_nodes)
        % check original coordinates
        orig_xy = orig_img.fwd_model.nodes(boundary_nodes(i),:);
        % find row of new node with same coordinates
        d = sqrt( (rot_img.fwd_model.nodes(:,1)-orig_xy(1)).^2 + (rot_img.fwd_model.nodes(:,2)-orig_xy(2)).^2 );
        [~,row_idx] = min(d);;
        % save mapping
        boundary_mapping(i) = row_idx;
    end
end

function [b_isOuterAndBoundary,boundary_node_arg_idx,boundary_node_idx] = ...
    isOuterAndBoundary(node1,node2,outer_nodes,boundary_nodes)
    
    if ismember(node1,outer_nodes) && ismember(node2,boundary_nodes)
        b_isOuterAndBoundary = true;
        boundary_node_arg_idx = 2;
        [~,boundary_node_idx] = ismember(node2,boundary_nodes);
    elseif ismember(node2,outer_nodes) && ismember(node1,boundary_nodes)
        b_isOuterAndBoundary = true;
        boundary_node_arg_idx = 1;
        [~,boundary_node_idx] = ismember(node1,boundary_nodes);
    else
        b_isOuterAndBoundary = false;
        boundary_node_arg_idx = 0;
        boundary_node_idx = 0;
    end
end


