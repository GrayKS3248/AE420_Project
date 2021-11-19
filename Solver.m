%% Setup
close all
clear


%% Parameter definition
E = 200.0;
v = 0.28;


%% Load mesh
load('elements.mat')
load('type.mat')
nodes = ELEMENT_TRIANGULATION.Points;
connectivity = ELEMENT_TRIANGULATION.ConnectivityList;


%% Calculate E matrix
E_mat = E/((1+v)*(1-2*v)) * [1-v v v 0 0 0; v 1-v v 0 0 0; v v 1-v 0 0 0; 0 0 0 0.5*(1-2*v) 0 0; 0 0 0 0 0.5*(1-2*v) 0; 0 0 0 0 0 0.5*(1-2*v)];


%% Elementwise stiffness
n_nodes = length(nodes(:,1));
n_elements = length(connectivity(:,1));
K_elements = zeros(n_elements, 12, 12);
for ele = 1:n_elements
    % Load element nodes
    P1 = nodes(connectivity(ele,1),:)';
    P2 = nodes(connectivity(ele,2),:)';
    P3 = nodes(connectivity(ele,3),:)';
    P4 = nodes(connectivity(ele,4),:)';
    
    % Define element normals
    N1 = cross( (P4-P2), (P4-P3) )';
    N2 = cross( (P4-P1), (P4-P3) )';
    N3 = cross( (P4-P1), (P4-P2) )';
    N4 = cross( (P3-P1), (P3-P2) )';
    N = [N1; N2; N3; N4];
    
    % Solve for shape function weights
    c1 = 1.0 / dot( N1, (P1-P4) );
    c2 = 1.0 / dot( N2, (P2-P4) );
    c3 = 1.0 / dot( N3, (P3-P4) );
    c4 = 1.0 / dot( N4, (P4-P3) );
    c = [c1, c2, c3, c4];
    
    % Solve for dN
    dN = zeros(6,12);
    for node=1:4
        dN(:,3*(node-1)+1) = [c(node)*N(node,1); 0.0; 0.0; c(node)*N(node,2); 0.0; c(node)*N(node,3)];
        dN(:,3*(node-1)+2) = [0.0; c(node)*N(node,2); 0.0; c(node)*N(node,1); c(node)*N(node,3); 0.0];
        dN(:,3*(node-1)+3) = [0.0; 0.0; c(node)*N(node,3); 0.0; c(node)*N(node,2); c(node)*N(node,1)];
    end
    dN(dN<1e-6) = 0.0;

    % Solve for the stiffness
    K_elements(ele,:,:) = dN'*E*dN;

end

%% Stiffness assembly
K_global = zeros(n_nodes*3, n_nodes*3);
for ele = 1:n_elements

    K_ele = squeeze(K_elements(ele,:,:));
    ele_nodes = connectivity(ele,:);
    node_start_ind = (ele_nodes-1)*3+1;
    node_end_ind = node_start_ind + 2;

    for node_x = 1:4
        for node_y =1:4
            node_x_start = (node_x-1)*3+1;
            node_x_end = node_x_start + 2;
            node_y_start = (node_y-1)*3+1;
            node_y_end = node_y_start + 2;

            K_global(node_start_ind(node_x):node_end_ind(node_x),node_start_ind(node_y):node_end_ind(node_y)) = ...
            K_global(node_start_ind(node_x):node_end_ind(node_x),node_start_ind(node_y):node_end_ind(node_y)) + ...
            K_ele(node_x_start:node_x_end, node_y_start:node_y_end);

        end
    end

end


%% Elementwise loading



%% Load assembly



%% Solution generation



%% Post processing



%% Visualization

