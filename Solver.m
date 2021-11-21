%% Setup
close all
clear


%% Problem definition
% Material properties
E = 200e9;    % Pascals
v = 0.30;     % Unitless

% Loading
top_traction = 100e4;  % Pascals

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
node_types_that_are_fixed_in_z = [1, 2, 3, 4];
node_indices_that_are_fixed_in_x = [];
node_indices_that_are_fixed_in_y = [];
node_indices_that_are_fixed_in_z = [];
fixed_node_z_range = [5.0, 15.9];  %mm


%% Simulation and visualization options
make_sparse = true;
show_fixed = true;
show_load = true;
visualize = true;
cross_section = true;
save_data = true;


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


%% Elementwise stiffness and loading
n_nodes = length(nodes(:,1));
n_elements = length(connectivity(:,1));
K_elements = zeros(n_elements, 12, 12);
top_z = max(nodes(:,3));
R_T_elements = zeros(12, n_elements);
dN_elements = zeros(n_elements, 6,12);

for ele = 1:n_elements

    % Get element nodes
    P1 = nodes(connectivity(ele,1),:)';
    P2 = nodes(connectivity(ele,2),:)';
    P3 = nodes(connectivity(ele,3),:)';
    P4 = nodes(connectivity(ele,4),:)';
    
    % Calculate element face normals
    N1 = cross( (P4-P2), (P4-P3) )' / norm(cross( (P4-P2), (P4-P3) ));
    N2 = cross( (P4-P1), (P4-P3) )' / norm(cross( (P4-P1), (P4-P3) ));
    N3 = cross( (P4-P1), (P4-P2) )' / norm(cross( (P4-P1), (P4-P2) ));
    N4 = cross( (P3-P1), (P3-P2) )' / norm(cross( (P3-P1), (P3-P2) ));
    N = [N1; N2; N3; N4];
    
    % Solve for shape function weights
    k1 = 1.0 / dot( N1, (P1-P4) );
    k2 = 1.0 / dot( N2, (P2-P4) );
    k3 = 1.0 / dot( N3, (P3-P4) );
    k4 = 1.0 / dot( N4, (P4-P3) );
    k = [k1, k2, k3, k4];

    % Volume calculation
    volume = abs(det([1 P1'; 1 P2'; 1 P3'; 1 P4']) / 6.0);

    % Solve for dN
    for node=1:4
        dN_elements(ele,:,3*(node-1)+1) = [k(node)*N(node,1); 0.0; 0.0; k(node)*N(node,2); 0.0; k(node)*N(node,3)];
        dN_elements(ele,:,3*(node-1)+2) = [0.0; k(node)*N(node,2); 0.0; k(node)*N(node,1); k(node)*N(node,3); 0.0];
        dN_elements(ele,:,3*(node-1)+3) = [0.0; 0.0; k(node)*N(node,3); 0.0; k(node)*N(node,2); k(node)*N(node,1)];
    end

    % Solve for the element stiffness
    K_elements(ele,:,:) = (squeeze(dN_elements(ele,:,:))'*E_mat*squeeze(dN_elements(ele,:,:)))*volume;

    % Calculate top of screw head traction loading
    if sum(nodes(connectivity(ele,:),3) == top_z) == 3

        % Determine which element nodes are on the top of the head
        top_of_head_nodes = connectivity(ele,(nodes(connectivity(ele,:),3))==top_z);
        coords = nodes(top_of_head_nodes,1:2);

        % Calculate the element top surface area of resultant load
        area = 0.5*abs(coords(1,1)*(coords(2,2)-coords(3,2)) + coords(2,1)*(coords(3,2)-coords(1,2)) + coords(3,1)*(coords(1,2)-coords(2,2)));
        total_load = area*top_traction;
        load_per_node = total_load / 3.0;

        % Apply calculated load to element-wise load vector
        for node=1:4
            if sum(connectivity(ele,node)==top_of_head_nodes)==1
                R_T_elements((node-1)*3+3, ele) = load_per_node;
            end
        end
    end
end

clear P1 P2 P3 P4 N1 N2 N3 N4 N k1 k2 k3 k4 k ele node top_of_head_nodes coords area total_load load_per_node top_z volume


%% Stiffness and load assembly
if make_sparse
    sparse_i_indices = zeros(12*12*n_elements,1);
    sparse_j_indices = zeros(12*12*n_elements,1);
    sparse_data = zeros(12*12*n_elements,1);
    sparse_index = 0;
else
    K_global = zeros(n_nodes*3, n_nodes*3);
end
R_T_global = zeros(3*n_nodes,1);

for ele = 1:n_elements

    K_ele = squeeze(K_elements(ele,:,:));
    R_T_ele = R_T_elements(:, ele);
    ele_nodes = connectivity(ele,:);
    global_ind_start = (ele_nodes-1)*3+1;
    global_ind_end = global_ind_start + 2;

    for i = 1:4
        i_ind_start = (i-1)*3+1;
        i_ind_end = i_ind_start + 2;

        % Assemble load vector
        R_T_global(global_ind_start(i):global_ind_end(i),1) = ...
        R_T_global(global_ind_start(i):global_ind_end(i),1) + ...
        R_T_ele(i_ind_start:i_ind_end, 1);

        for j =1:4
            j_ind_start = (j-1)*3+1;
            j_ind_end = j_ind_start + 2;

            if make_sparse
                % Assemble sparse stiffness matrix
                for sparse_i = i_ind_start:i_ind_end
                    for sparse_j = j_ind_start:j_ind_end
                        sparse_index = sparse_index + 1;
                        sparse_i_indices(sparse_index,:) = global_ind_start(i)+sparse_i-i_ind_start;
                        sparse_j_indices(sparse_index,:) = global_ind_start(j)+sparse_j-j_ind_start;
                        sparse_data(sparse_index,:) = K_ele(sparse_i, sparse_j);
                    end
                end


            else
                % Assemble stiffness matrix
                K_global(global_ind_start(i):global_ind_end(i),global_ind_start(j):global_ind_end(j)) = ...
                K_global(global_ind_start(i):global_ind_end(i),global_ind_start(j):global_ind_end(j)) + ...
                K_ele(i_ind_start:i_ind_end, j_ind_start:j_ind_end);
            end

        end
    end
end
if make_sparse
    K_global = sparse(sparse_i_indices, sparse_j_indices, sparse_data);
else
    K_global = sparse(K_global);
end
R_global = R_T_global;

clear K_elements R_T_elements R_T_global K_ele R_T_ele ele_nodes global_ind_start global_ind_end i_ind_start i_ind_end j_ind_start j_ind_end ele i j;
clear sparse_i sparse_j sparse_index sparse_i_indices sparse_j_indices sparse_data


%% Application of fixed BC
% Determine node fixing type
x_fixed = [];
y_fixed = [];
z_fixed = [];
for node=1:n_nodes
    if nodes(node,3) <= (1e-3*fixed_node_z_range(2) + min(nodes(:,3))) && nodes(node,3) >= (1e-3*fixed_node_z_range(1) + min(nodes(:,3)))
        if sum(TYPE(node,1) == node_types_that_are_fixed_in_x)==1 || sum(node == node_indices_that_are_fixed_in_x)==1
            x_fixed(end+1)=node;
        end
        if sum(TYPE(node,1) == node_types_that_are_fixed_in_y)==1 || sum(node == node_indices_that_are_fixed_in_y)==1
            y_fixed(end+1)=node;
        end
        if sum(TYPE(node,1) == node_types_that_are_fixed_in_z)==1 || sum(node == node_indices_that_are_fixed_in_z)==1
            z_fixed(end+1)=node;
        end
    end
end
indices_2_remove = sort([(x_fixed-1)*3+1, (y_fixed-1)*3+2, (z_fixed-1)*3+3])';

% Trim global stiffness matrix and global load vector
K_global_trimmed = K_global;
K_global_trimmed((indices_2_remove),:) = [];
K_global_trimmed(:,(indices_2_remove)) = [];
R_global_trimmed = R_global;
R_global_trimmed((indices_2_remove),:) = [];

clear node


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

clear K_global_trimmed R_global_trimmed displacement_trimmed num_added indices_2_remove i


%% Post processing
% Displacement
dimensional_displacement = reshape(displacement,3,[])';
displacement_mag = vecnorm(dimensional_displacement,2,2);

% Node positions
displaced_nodes = 1000.0*(nodes + dimensional_displacement);
scale_factor = 0.25*(max(nodes(:,3)) - min(nodes(:,3))) / max(abs(displacement));
scaled_displaced_nodes = 1000.0*(nodes + scale_factor*dimensional_displacement);

% Reaction forces
reaction_load = K_global*displacement - R_global;
reaction_load(reaction_load<1.0e-5) = 0.0;
reaction_load = reshape(reaction_load,3,[])';
reaction_load_mag = vecnorm(reaction_load,2,2);

% Mesh wide element stress and strain
elements_strain = zeros(n_elements, 6);
elements_stress = zeros(n_elements, 6);
elements_von_mises = zeros(n_elements, 1);
elements_principal = zeros(n_elements, 3);
for ele=1:n_elements
    % Calculate element stress and strain
    global_ind_start = 3*(connectivity(ele,:)-1) + 1;
    global_ind_end = 3*(connectivity(ele,:)-1) + 3;
    displacment_sub_vec = [displacement(global_ind_start(1):global_ind_end(1));
                           displacement(global_ind_start(2):global_ind_end(2));
                           displacement(global_ind_start(3):global_ind_end(3)); 
                           displacement(global_ind_start(4):global_ind_end(4))];
    elements_strain(ele,:) = squeeze(dN_elements(1,:,:)) * displacment_sub_vec;
    elements_stress(ele,:) = E_mat * elements_strain(ele,:)';
    xx = elements_stress(ele,1);
    yy = elements_stress(ele,2);
    zz = elements_stress(ele,3);
    xy = elements_stress(ele,4);
    yz = elements_stress(ele,5);
    xz = elements_stress(ele,6);
    elements_von_mises(ele,:) = sqrt(0.5*((xx-yy)^2.0+(yy-zz)^2+(zz-xx)^2+6.0*((xy^2)+(yz^2)+(xz^2))));
    elements_principal(ele,:) = flip(eig([xx xy xz; xy yy yz; xz yz zz])');
end

if cross_section
    % Take a cross section of the mesh
    cross_section_connectivity = zeros(0,4);
    cross_section_orig_ele_ind = zeros(0,1);
    for ele = 1:n_elements
        ele_nodes = nodes(connectivity(ele,:),2);
    
        % Tet has at least one node on the XZ plane
        if min(abs(ele_nodes)) < 1.0e-10
            nodes_on_xz = abs(ele_nodes) < 1.0e-10;
    
            % Tet has a face on the XZ plane
            if sum(nodes_on_xz) == 3
                cross_section_connectivity(end+1,:) = connectivity(ele,:);
                cross_section_orig_ele_ind(end+1,:) = ele;
            
            % Tet has a segment on the XZ plane
            elseif sum(nodes_on_xz) == 2
                nodes_off_xz = ele_nodes(~nodes_on_xz);
    
                % Tet stradles XZ
                if min(nodes_off_xz) < 1.0e-10 && max(nodes_off_xz) > 1.0e-10
                    cross_section_connectivity(end+1,:) = connectivity(ele,:);
                    cross_section_orig_ele_ind(end+1,:) = ele;
                end
    
            % Tet has a point on the XZ plane
            elseif sum(nodes_on_xz) == 1
                nodes_off_xz = ele_nodes(~nodes_on_xz);
    
                % Tet stradles XZ
                if min(nodes_off_xz) < 1.0e-10 && max(nodes_off_xz) > 1.0e-10
                    cross_section_connectivity(end+1,:) = connectivity(ele,:);
                    cross_section_orig_ele_ind(end+1,:) = ele;
                end
    
            end
    
        % Tet intersects the XZ plane
        elseif min(nodes(connectivity(ele,:),2)) < 1.0e-10 && max(nodes(connectivity(ele,:),2)) > -1.0e-10
            cross_section_connectivity(end+1,:) = connectivity(ele,:);
            cross_section_orig_ele_ind(end+1,:) = ele;
        end
    
    end
    cross_section_orig_node_ind = reshape(cross_section_connectivity,[],1);
    cross_section_orig_node_ind = unique(cross_section_orig_node_ind);
    
    % Gather node-wise stress data in mesh cross section
    cross_section_node_stress = zeros(length(cross_section_orig_node_ind(:,1)), 6);
    cross_section_node_von_mises = zeros(length(cross_section_orig_node_ind(:,1)), 1);
    cross_section_node_principal = zeros(length(cross_section_orig_node_ind(:,1)), 3);
    for node_ind=1:length(cross_section_orig_node_ind(:,1))
        for ele = 1:length(cross_section_connectivity(:,1))
        
            if sum(cross_section_connectivity(ele,:)==cross_section_orig_node_ind(node_ind))==1
                
                % XX
                if abs(elements_stress(cross_section_orig_ele_ind(ele),1)) > abs(cross_section_node_stress(node_ind, 1))
                    cross_section_node_stress(node_ind, 1) = elements_stress(cross_section_orig_ele_ind(ele),1);
                end
                
                % YY
                if abs(elements_stress(cross_section_orig_ele_ind(ele),2)) > abs(cross_section_node_stress(node_ind, 2))
                    cross_section_node_stress(node_ind, 2) = elements_stress(cross_section_orig_ele_ind(ele),2);
                end
        
                % ZZ
                if abs(elements_stress(cross_section_orig_ele_ind(ele),3)) > abs(cross_section_node_stress(node_ind, 3))
                    cross_section_node_stress(node_ind, 3) = elements_stress(cross_section_orig_ele_ind(ele),3);
                end
        
                % XY
                if abs(elements_stress(cross_section_orig_ele_ind(ele),4)) > abs(cross_section_node_stress(node_ind, 4))
                    cross_section_node_stress(node_ind, 4) = elements_stress(cross_section_orig_ele_ind(ele),4);
                end
        
                % YZ
                if abs(elements_stress(cross_section_orig_ele_ind(ele),5)) > abs(cross_section_node_stress(node_ind, 5))
                    cross_section_node_stress(node_ind, 5) = elements_stress(cross_section_orig_ele_ind(ele),5);
                end
        
                % XZ
                if abs(elements_stress(cross_section_orig_ele_ind(ele),6)) > abs(cross_section_node_stress(node_ind, 6))
                    cross_section_node_stress(node_ind, 6) = elements_stress(cross_section_orig_ele_ind(ele),6);
                end
        
                % Von mises
                if abs(elements_von_mises(cross_section_orig_ele_ind(ele))) > abs(cross_section_node_von_mises(node_ind))
                    cross_section_node_von_mises(node_ind) = elements_von_mises(cross_section_orig_ele_ind(ele));
                end
    
                % 11
                if abs(elements_principal(cross_section_orig_ele_ind(ele),1)) > abs(cross_section_node_principal(node_ind, 1))
                    cross_section_node_principal(node_ind, 1) = elements_principal(cross_section_orig_ele_ind(ele),1);
                end
                
                % 22
                if abs(elements_principal(cross_section_orig_ele_ind(ele),2)) > abs(cross_section_node_principal(node_ind, 2))
                    cross_section_node_principal(node_ind, 2) = elements_principal(cross_section_orig_ele_ind(ele),2);
                end
        
                % 33
                if abs(elements_principal(cross_section_orig_ele_ind(ele),3)) > abs(cross_section_node_principal(node_ind, 3))
                    cross_section_node_principal(node_ind, 3) = elements_principal(cross_section_orig_ele_ind(ele),3);
                end
        
            end
        end
    end
    
    % Create vector of cross section nodes and use it to update connectivity indices
    cross_section_nodes = zeros(length(cross_section_orig_node_ind(:,1)), 3);
    for node = 1:length(cross_section_orig_node_ind(:,1))
        cross_section_nodes(node,:) = nodes(cross_section_orig_node_ind(node),:);
    end
    for ele = 1:length(cross_section_connectivity(:,1))
        for node = 1:4
            [~,arg_min] = min(abs(cross_section_orig_node_ind - cross_section_connectivity(ele,node)));
            cross_section_connectivity(ele,node) = arg_min;
        end
    end
    
    % Used for y projection plotting the cross section
    cross_section_mean_y_val = zeros(length(cross_section_connectivity(:,1)),1);
    cross_section_min_y_node = zeros(length(cross_section_connectivity(:,1)),1);
    for ele=1:length(cross_section_connectivity(:,1))
        cross_section_mean_y_val(ele) = mean(cross_section_nodes(cross_section_connectivity(ele,:),2));
        [~,arg_min] = min(cross_section_nodes(cross_section_connectivity(ele,:),2));
        cross_section_min_y_node(ele) = arg_min;
    end
    [~,cross_section_draw_order] = sort(cross_section_mean_y_val);
end

clear arg_min cross_section_mean_y_val cross_section_orig_ele_ind cross_section_orig_node_ind displacment_sub_vec ele ele_nodes
clear global_ind_end global_ind_start node node_ind nodes_off_xz nodes_on_xz scale_factor xx yy zz xy yz xz


%% Visualization
if visualize
    % Problem definition
    figure(1)
    colormap jet
    if show_fixed
        x_fixed_scatter = scatter3(1000.0*nodes(x_fixed',1), 1000.0*nodes(x_fixed',2), 1000.0*nodes(x_fixed',3), 30, 'o', 'm');
        hold on
        y_fixed_scatter = scatter3(1000.0*nodes(y_fixed',1), 1000.0*nodes(y_fixed',2), 1000.0*nodes(y_fixed',3), 30, '+', 'm');
        z_fixed_scatter = scatter3(1000.0*nodes(z_fixed',1), 1000.0*nodes(z_fixed',2), 1000.0*nodes(z_fixed',3), 60, '.', 'm');
    end
    if show_load
        norm_node_loading = reshape(R_global,3,[])';
        norm_node_loading = norm_node_loading / max(vecnorm(norm_node_loading,2,2));
        min_z = 1000.0*min(nodes(:,3));
        max_z = 1000.0*max(nodes(:,3));
        for node=1:n_nodes
            if norm(norm_node_loading(node,:))~=0.0
                p0 = displaced_nodes(node,:);
                p1 = p0 + norm_node_loading(node,:)*0.05*(max_z - min_z);
                vectarrow(p0,p1,0.005*(max_z - min_z))
                hold on
            end
        end
    end
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, displaced_nodes(:,1), displaced_nodes(:,2), displaced_nodes(:,3), 1000.0*displacement_mag, 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    hold off
    title("Problem Definition",'FontSize',16)
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    c = colorbar('eastoutside');
    c.Label.String = 'Displacement Magnitude (mm)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    if show_fixed
        legend([x_fixed_scatter, y_fixed_scatter, z_fixed_scatter], {'X Fixed', 'Y Fixed', 'Z Fixed'}, 'Location', 'westoutside','FontSize',12)
    end
    axis equal;
    view(30,30);
    set(gcf,'Position',[200 55 750 1000])
    saveas(gcf,'problem_def.png',"png")
    
    % Scaled displacement magnitude
    figure(2)
    colormap jet
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, scaled_displaced_nodes(:,1), scaled_displaced_nodes(:,2), scaled_displaced_nodes(:,3), 1000.0*displacement_mag, 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    title("Scaled Displaced Mesh",'FontSize',16)
    c = colorbar('eastoutside');
    c.Label.String = 'Displacement Magnitude (mm)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    h.ZAxis.TickLength = [0 0];
    axis equal;
    view(30,30);
    set(gcf,'Position',[1098 55 750 1000])
    saveas(gcf,'scaled_displaced.png',"png")
    
    % X displacement
    figure(3)
    colormap jet
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1000.0*dimensional_displacement(:,1), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    title("X Displacement",'FontSize',16)
    c = colorbar('eastoutside');
    c.Label.String = 'X Displacement (mm)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    h.ZAxis.TickLength = [0 0];
    axis equal;
    view(30,30);
    set(gcf,'Position',[200 55 750 1000])
    saveas(gcf,'x_displacement.png',"png")
    
    % Y displacement
    figure(4)
    colormap jet
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1000.0*dimensional_displacement(:,2), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    title("Y Displacement",'FontSize',16)
    c = colorbar('eastoutside');
    c.Label.String = 'Y Displacement (mm)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    h.ZAxis.TickLength = [0 0];
    axis equal;
    view(30,30);
    set(gcf,'Position',[1098 55 750 1000])
    saveas(gcf,'y_displacement.png',"png")
    
    % Z displacement
    figure(5)
    colormap jet
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1000.0*dimensional_displacement(:,3), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    title("Z Displacement",'FontSize',16)
    c = colorbar('eastoutside');
    c.Label.String = 'Z Displacement (mm)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    h.ZAxis.TickLength = [0 0];
    axis equal;
    view(30,30);
    set(gcf,'Position',[200 55 750 1000])
    saveas(gcf,'z_displacement.png',"png")
    
    % Reaction load magnitude
    figure(6)
    colormap jet
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), reaction_load_mag, 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
    title("Reaction Load Magnitude",'FontSize',16)
    c = colorbar('eastoutside');
    c.Label.String = 'Load (N)';
    c.Label.FontSize = 12;
    c.FontSize = 12;
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    h.ZAxis.TickLength = [0 0];
    axis equal;
    view(30,30);
    set(gcf,'Position',[1098 55 750 1000])
    saveas(gcf,'reaction_mag.png',"png")
    
    if cross_section
        % Von Mises
        figure(7)
        colormap jet
        for tri =1:length(cross_section_draw_order(:,1))
            ele = cross_section_draw_order(tri);
            dropped_node = [1,2,3,4] ~= cross_section_min_y_node(ele,:);
            nodes = cross_section_connectivity(ele,dropped_node);
            trisurf(nodes, cross_section_nodes(:,1), zeros(length(cross_section_nodes(:,1)), 1), cross_section_nodes(:,3), 1.0e-6*cross_section_node_von_mises);
            hold on
        end
        title("Von Mises Stress Y=0",'FontSize',16)
        c = colorbar('eastoutside');
        c.Label.String = 'σv (MPA)';
        c.Label.FontSize = 12;
        c.FontSize = 12;
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        h.ZAxis.TickLength = [0 0];
        axis equal;
        view(0,0);
        shading interp
        set(gcf,'Position',[200 55 750 1000])
        saveas(gcf,'von_mises.png',"png")
        
        % First principal
        figure(8)
        colormap jet
        for tri =1:length(cross_section_draw_order(:,1))
            ele = cross_section_draw_order(tri);
            dropped_node = [1,2,3,4] ~= cross_section_min_y_node(ele,:);
            nodes = cross_section_connectivity(ele,dropped_node);
            trisurf(nodes, cross_section_nodes(:,1), zeros(length(cross_section_nodes(:,1)), 1), cross_section_nodes(:,3), 1.0e-6*cross_section_node_principal(:,1));
            hold on
        end
        title("First Principal Stress Y=0",'FontSize',16)
        c = colorbar('eastoutside');
        c.Label.String = 'σ1 (MPA)';
        c.Label.FontSize = 12;
        c.FontSize = 12;
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        h.ZAxis.TickLength = [0 0];
        axis equal;
        view(0,0);
        shading interp
        set(gcf,'Position',[1098 55 750 1000])
        saveas(gcf,'first_principal.png',"png")
        
        % Third principal
        figure(9)
        colormap jet
        for tri =1:length(cross_section_draw_order(:,1))
            ele = cross_section_draw_order(tri);
            dropped_node = [1,2,3,4] ~= cross_section_min_y_node(ele,:);
            nodes = cross_section_connectivity(ele,dropped_node);
            trisurf(nodes, cross_section_nodes(:,1), zeros(length(cross_section_nodes(:,1)), 1), cross_section_nodes(:,3), 1.0e-6*cross_section_node_principal(:,3));
            hold on
        end
        title("Third Principal Stress Y=0",'FontSize',16)
        c = colorbar('eastoutside');
        c.Label.String = 'σ3 (MPA)';
        c.Label.FontSize = 12;
        c.FontSize = 12;
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        h.ZAxis.TickLength = [0 0];
        
        axis equal;
        view(0,0);
        shading interp
        set(gcf,'Position',[200 55 750 1000])
        saveas(gcf,'third_principal.png',"png")
    end
end

% Save
if save_data
    save('strain.mat','elements_strain');
    save('stress.mat','elements_stress');
    save('von_mises.mat','elements_von_mises');
    save('principal.mat','elements_principal');
    save('displacement.mat','dimensional_displacement');
    save('reaction.mat','reaction_load');
end

% Readout
n_elements = n_elements
max_von_mises = 1.0e-6*max(elements_von_mises)
max_first = 1.0e-6*max(elements_principal(:,1))
min_third = 1.0e-6*min(elements_principal(:,3))
max_displacement = max(displacement_mag)*1e6

clear dropped_node nodes tri ele min_z max_z node p0 p1 c h x_fixed_scatter y_fixed_scatter z_fixed_scatter norm_node_loading s show_load show_fixed make_sparse cross_section save


%% Helper functions
function vectarrow(p0,p1,head)
  if max(size(p0))==3
      if max(size(p1))==3
          x0 = p0(1);
          y0 = p0(2);
          z0 = p0(3);
          x1 = p1(1);
          y1 = p1(2);
          z1 = p1(3);
          plot3([x0;x1],[y0;y1],[z0;z1],'k','LineWidth',2.0);   % Draw a line between p0 and p1
          [X,Y,Z]=cylinder([head 0],10);
          M=makehgtform('translate',p1);
          surf(X,Y,Z,zeros(1,length(X(1,:))),'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',1.0,'FaceColor','k');

      else
          error('p0 and p1 must have the same dimension')
      end
  else
      error('this function only accepts 2D or 3D vector')
  end
end