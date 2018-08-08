%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our article "Rotational electrical 
% impedance tomography using electrodes with limited boundary coverage
% provides window for multimodal sensing".
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ M ] = CreateRotationMatrix( eidors_img, thetas )
%CREATEROTATIONMATRIX Returns a set of sparse mapping matrices such that 
%M = N x N x D where N is number of elements in EIDORS mesh eidors_img and
%D is the number of rotation angles i.e. numel(thetas)
%
%Function computes the mapping matrix by rotating each element separately
%for all thetas. Results such that if Mi is the i'th matrix in M, then new
%element positions e1 can be computed from old elements e by
% e1 = M*e
%Here e and e1 are both in vectoriced form.
%
%Rotation direction of the mapping matrix is counterclockwise.

% DEVELOPMENT
%eidors_img = img; thetas = 180; %linspace(0,360,3);

%%
% allocate M
N = num_elems(eidors_img.fwd_model);
M = zeros(N,N,numel(thetas));

for t = 1:numel(thetas)
for i = 1:size(eidors_img.fwd_model.elems,1)
    
    % Pick the coordinates of the nodes of i'th element
    iN1 = eidors_img.fwd_model.elems(i,1);
    iN2 = eidors_img.fwd_model.elems(i,2);
    iN3 = eidors_img.fwd_model.elems(i,3);
    Nxy = [eidors_img.fwd_model.nodes(iN1,:);
        eidors_img.fwd_model.nodes(iN2,:);
        eidors_img.fwd_model.nodes(iN3,:)];
    
    % Rotate the coordinates
    Nrotated = ElemNodesPolarRotation(Nxy,thetas(t));
    
    % Get elem weights under rotation
    %   first define the new elem as polygon
    Npoly = [Nrotated; Nrotated(1,:)];
    %   define the membership selection function
    %   here points lying on edges of polygon are defined out, need to use
    %   addition on output of inpolygon if desired otherwise (see MATLAB
    %   doc)
    select_fcn = @(x,y,z) inpolygon(x,y,Npoly(:,1),Npoly(:,2));
    %   now get the elem membership element weights
    memb_weights = elem_select( eidors_img.fwd_model, select_fcn );
    
    % Normalize the weights to 1 (Only setting data to one element with
    % this M. Conversely: think about situation where element in e1 is
    % covering two elements in e.)
    memb_weights = memb_weights / sum(memb_weights);
    
    % Insert the element weights to corresponding row in M
    M(i,:,t) = memb_weights';
    
end
end

end

