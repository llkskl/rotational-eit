%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Finite Element Mapping for
% Efficient Image Reconstruction in Rotational Electrical Impedance 
% Tomography".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ J ] = ClickingRotationalJacobian( mdl, img )
% Computes jacobian based on rotational mapping.

    NROTPOS = mdl.rotational_parameters.n_rotational_positions;
    % first compute the initial jacobian
    J1 = calc_jacobian( img );
    
    % allocate
    J = zeros(size(J1,1)*NROTPOS, size(J1,2));
       
	% save initial jacobian
    J(1:size(J1,1),:) = J1;
          
	% compute rest of the jacobians
	for i = 2:NROTPOS
        rotimg = ClickRotationalMesh(mdl,i);
        rotimg = data_mapper(calc_jacobian_bkgnd( rotimg ));
        Ji = calc_jacobian(rotimg);
        J( (i-1)*size(J1,1) + (1:size(Ji,1)) , :) = Ji;
    end

end