%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This code finds the displacements, forces and nodal forces for a
%          bar, using the finite element method.
% Written By: Esteban Valderrama, 11/25/2015
%             University of Wisconsin at Platteville
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%%  Read Input file
file = 'problem1.txt'; % problem2.txt
[nodes, Xy, Ndof, dof, fixedDof, Nglob, nelems, elems, A, E, F, nF] = readFiles(file);

%% Create Element and Global Stiffness Matrices
[K, Ke] = stiffnessMatrix(A, E, Xy, elems, nelems, nodes);

%%  Solve for Displacements
u = zeros(nodes,1);
u(fixedDof) = 0; 
free = setdiff([1:Ndof]', [fixedDof]);
u(free,1) = K(free,free)\F(free);

%%  Solve for Force
F = K * u;

%%  Solve for node forces
% Solve for f when the problem with different K
for i = 1:nelems
    dof = elems(i,:);      % Select the element
    f(i,:) = Ke(:,:,1) * u(dof,1); % Solve for f
    show_f(i,:) = [i f(i,:)];     % Accumulate results
end
%%  Display results
disp('Local Stiffness Matrix:')
for i = 1:nelems
    fprintf('Element: %d\n', i)
    disp(Ke(:,:,i))
    
end
disp('Global Stiffness Matrix:')
disp(K)
disp('Displacements:')
n = 1:nodes; format
disp([n' u])
disp('Reactions:')
disp([fixedDof F(fixedDof)])
disp('Element Forces:')
disp(show_f)