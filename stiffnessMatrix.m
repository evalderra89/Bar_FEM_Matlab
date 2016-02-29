%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function created the global stiffness matrix.
% Written By: Esteban Valderrama, 11/25/2015
%             University of Wisconsin at Platteville
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, Ke] = stiffnessMatrix(A, E, Xy, elems, nelems, nodes)

% Create Element Stiffness Matrix
Ke = zeros(2,2,nelems);
%  Create Stiffness Matrix
K = zeros(nodes,nodes);
for i = 1:nelems
    stiffness = (A(i) * E(i))/ (Xy(i,2)-Xy(i,1));
    Ke(:,:,i) = [stiffness -stiffness;-stiffness stiffness];
    dof = elems(i,:);            % Select the element
    K(dof,dof) = K(dof,dof) + Ke(:,:,1); % Fill K
end