%% Setup
close all
clear


%% Problem definition
% Material properties
E = 200.0e9;
v = 0.28;

% Loading
tension_load_on_top_of_head = 350.0e6; % Pascals
tension_load_on_tip = 0.0;             % Newtons

% Boundary conditions
% 1 - Body
% 2 - Approach face of thread
% 3 - Crest of thread
% 4 - Receeding face of thread
% 5 - Shank
% 6 - Shank to head transition
% 7 - Head
node_types_that_are_fixed_in_x = [];
node_types_that_are_fixed_in_y = [];
node_types_that_are_fixed_in_z = [];
node_indices_that_are_fixed = linspace(1,87,87);


%% Load mesh
load('elements.mat')
load('type.mat')
nodes = 1.0e-3 * ELEMENT_TRIANGULATION.Points;
connectivity = ELEMENT_TRIANGULATION.ConnectivityList;


%% Calculate E matrix
E_mat = E/((1+v)*(1-2*v)) * [1-v, v, v, 0, 0, 0; 
                             v, 1-v, v, 0, 0, 0; 
                             v, v, 1-v, 0, 0, 0; 
                             0, 0, 0, 0.5*(1-2*v), 0, 0; 
                             0, 0, 0, 0, 0.5*(1-2*v), 0; 
                             0, 0, 0, 0, 0, 0.5*(1-2*v)];


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
    dN(abs(dN)<1e-6) = 0.0;

    % Solve for the stiffness
    K_elements(ele,:,:) = dN'*E_mat*dN;
end

%% Stiffness assembly
K_global = zeros(n_nodes*3, n_nodes*3);
for ele = 1:n_elements

    K_ele = squeeze(K_elements(ele,:,:));
    ele_nodes = connectivity(ele,:);
    nodes_start = (ele_nodes-1)*3+1;
    nodes_end = nodes_start + 2;

    for i = 1:4
        i_start = (i-1)*3+1;
        i_end = i_start + 2;

        for j =1:4
            j_start = (j-1)*3+1;
            j_end = j_start + 2;

            K_global(nodes_start(i):nodes_end(i),nodes_start(j):nodes_end(j)) = ...
            K_global(nodes_start(i):nodes_end(i),nodes_start(j):nodes_end(j)) + ...
            K_ele(i_start:i_end, j_start:j_end);
        end
    end
end
K_global = sparse(K_global);


%% Elementwise loading

% Top of head traction
top_of_head_height = max(nodes(:,3));
R_T_elements = zeros(12, n_elements);
for ele=1:n_elements

    % If we are on an element who has a normal surface on the top of the head of the screw
    if sum(nodes(connectivity(ele,:),3) == top_of_head_height) == 3

        % Calculate the equally distributed load per node on top of head
        % elements
        top_of_head_nodes = connectivity(ele,(nodes(connectivity(ele,:),3))==top_of_head_height);
        coords = nodes(top_of_head_nodes,1:2);
        area = 0.5*abs(coords(1,1)*(coords(2,2)-coords(3,2)) + coords(2,1)*(coords(3,2)-coords(1,2)) + coords(3,1)*(coords(1,2)-coords(2,2)));
        total_load = area*tension_load_on_top_of_head;
        load_per_node = total_load / 3.0;

        % Apply calculated load to element-wise load vector
        for node=1:4
            if sum(connectivity(ele,node)==top_of_head_nodes)==1
                R_T_elements((node-1)*3+3, ele) = load_per_node;
            end
        end
    end
end

% TODO: IMPLEMENT BODY FORCES
% TODO: IMPLEMENT POINT LOADS
% TODO: IMPLEMENT TORSION
% TODO: IMPLEMENT BENDING


%% Load assembly
R_T_global = zeros(3*n_nodes,1);
for ele = 1:n_elements

    R_T_ele = R_T_elements(:, ele);
    ele_nodes = connectivity(ele,:);
    nodes_start = (ele_nodes-1)*3+1;
    nodes_end = nodes_start + 2;

    for i = 1:4
        i_start = (i-1)*3+1;
        i_end = i_start + 2;
        R_T_global(nodes_start(i):nodes_end(i),1) = R_T_global(nodes_start(i):nodes_end(i),1) + R_T_ele(i_start:i_end, 1);
    end
end

% TODO: IMPLEMENT BODY FORCES
% TODO: IMPLEMENT POINT LOADS
% TODO: IMPLEMENT TORSION
% TODO: IMPLEMENT BENDING

R_global = R_T_global;


%% Application of fixed BC
% Determine node fixing type
x_fixed = [];
y_fixed = [];
z_fixed = [];
for node=1:n_nodes
    if sum(TYPE(node,1) == node_types_that_are_fixed_in_x)==1 || sum(node == node_indices_that_are_fixed)==1
        x_fixed(end+1)=node;
    end
    if sum(TYPE(node,1) == node_types_that_are_fixed_in_y)==1 || sum(node == node_indices_that_are_fixed)==1
        y_fixed(end+1)=node;
    end
    if sum(TYPE(node,1) == node_types_that_are_fixed_in_z)==1 || sum(node == node_indices_that_are_fixed)==1
        z_fixed(end+1)=node;
    end
end
indices_2_remove = sort([(x_fixed-1)*3+1, (y_fixed-1)*3+2, (z_fixed-1)*3+3])';

% Trim global stiffness matrix and global load vector
K_global_trimmed = K_global;
K_global_trimmed((indices_2_remove),:) = [];
K_global_trimmed(:,(indices_2_remove)) = [];
R_global_trimmed = R_global;
R_global_trimmed((indices_2_remove),:) = [];

%% Solution generation
displacement_trimmed = K_global_trimmed\R_global_trimmed;
displacement = zeros(3*n_nodes, 1);
num_added = 0;
for i=1:3*n_nodes
    if sum(i==indices_2_remove)==1
        displacement(i,1) = 0.0;
        num_added = num_added + 1;
    else
        displacement(i,1) = displacement_trimmed(i-num_added,1);
    end
end


%% Post processing
displaced_nodes = 1000.0*(nodes + reshape(displacement,3,[])');


%% Visualization
trisurf(ELEMENT_TRIANGULATION.freeBoundary, displaced_nodes(:,1), displaced_nodes(:,2), displaced_nodes(:,3), 'EdgeAlpha', '0.1');
title("Displaced Mesh")
xlabel('mm')
ylabel('mm')
zlabel('mm')
axis equal;
view(30,30);
saveas(gcf,'displaced.png',"png")