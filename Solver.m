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
node_types_that_are_fixed_in_z = [2, 3, 4];
node_indices_that_are_fixed_in_x = [];
node_indices_that_are_fixed_in_y = [];
node_indices_that_are_fixed_in_z = [];
fixed_node_z_range = [5.0, 10.5];  %mm


%% Visualization options
show_fixed = true;
show_load = true;


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
K_global = zeros(n_nodes*3, n_nodes*3);
R_T_global = zeros(3*n_nodes,1);

for ele = 1:n_elements

    K_ele = squeeze(K_elements(ele,:,:));
    R_T_ele = R_T_elements(:, ele);
    ele_x_coords = connectivity(ele,:);
    global_ind_start = (ele_x_coords-1)*3+1;
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

            % Assemble stiffness matrix
            K_global(global_ind_start(i):global_ind_end(i),global_ind_start(j):global_ind_end(j)) = ...
            K_global(global_ind_start(i):global_ind_end(i),global_ind_start(j):global_ind_end(j)) + ...
            K_ele(i_ind_start:i_ind_end, j_ind_start:j_ind_end);

        end
    end
end
K_global = sparse(K_global);
R_global = R_T_global;

clear K_elements R_T_elements R_T_global K_ele R_T_ele ele_x_coords global_ind_start global_ind_end i_ind_start i_ind_end j_ind_start j_ind_end ele i j;


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
scaled_displaced_nodes = 1000.0*(nodes + (0.25 * ((max(nodes(:,3)) - min(nodes(:,3))) / max(abs(displacement))))*dimensional_displacement);

% Reaction forces
reaction_load = K_global*displacement - R_global;
reaction_load(reaction_load<1.0e-5) = 0.0;
reaction_load = reshape(reaction_load,3,[])';
reaction_load_mag = vecnorm(reaction_load,2,2);

% Stress and strain
elements_strain = zeros(n_elements, 6);
elements_stress = zeros(n_elements, 6);
elements_von_mises = zeros(n_elements, 1);
nodes_strain = zeros(n_nodes, 7);
nodes_stress = zeros(n_nodes, 7);
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

    % Assign average stress and strain of all connected elements to each node
    for node=1:4
        nodes_strain(connectivity(ele,node),1:6) = elements_strain(ele,:);
        nodes_strain(connectivity(ele,node),7) = nodes_strain(connectivity(ele,1),7) + 1;
        nodes_stress(connectivity(ele,node),1:6) = elements_stress(ele,:);
        nodes_stress(connectivity(ele,node),7) = nodes_stress(connectivity(ele,1),7) + 1;
    end
end
nodes_strain = nodes_strain(:,1:6)./nodes_strain(:,7);
nodes_stress = nodes_stress(:,1:6)./nodes_stress(:,7);
nodes_von_mises = zeros(n_nodes, 1);
for node=1:n_nodes
    xx = nodes_stress(node,1);
    yy = nodes_stress(node,2);
    zz = nodes_stress(node,3);
    xy = nodes_stress(node,4);
    yz = nodes_stress(node,5);
    xz = nodes_stress(node,6);
    nodes_von_mises(node,:) = sqrt(0.5*((xx-yy)^2.0+(yy-zz)^2+(zz-xx)^2+6.0*((xy^2)+(yz^2)+(xz^2))));
end

% Stress and strain cross section data
precision = 750;
delta_x = max(nodes(:,1))-min(nodes(:,1));
delta_z = max(nodes(:,3))-min(nodes(:,3));
norm_delta_x = delta_x / (delta_x+delta_z);
norm_delta_z = delta_z / (delta_x+delta_z);
norm_precision = precision / max([norm_delta_x, norm_delta_z]);
x_precision = round(norm_precision * norm_delta_x);
z_precision = round(norm_precision * norm_delta_z);
x_space = linspace( min(nodes(:,1)), max(nodes(:,1)), x_precision);
z_space = linspace( max(nodes(:,3)), min(nodes(:,3)), z_precision);
delta_x = abs(x_space(2) - x_space(1));
delta_z = abs(z_space(2) - z_space(1));
[X,Z] = meshgrid(x_space, z_space);
cross_section_stress = zeros(z_precision, x_precision, 6);
cross_section_strain = zeros(z_precision, x_precision, 6);
cross_section_von_mises = zeros(z_precision, x_precision, 1);
x_ind = 0;

% Find those elements who are around the y axis
mid_plane_elements = [];
for ele=1:n_elements
    if min(nodes(connectivity(ele,:),2)) <= 0.0 && max(nodes(connectivity(ele,:),2)) >= 0.0
        mid_plane_elements(end+1) = ele;
    end
end
n_mid_plane_elements = length(mid_plane_elements);

% Populate the cross sectional space
for x=x_space
    z_ind = 0;
    x_ind = x_ind + 1;

    % Search for those elements whose x range surround the current x point
    mid_plane_elements_with_good_x = [];
    for ele_ind=1:n_mid_plane_elements
        ele_x_coords = nodes(connectivity(mid_plane_elements(ele_ind),:),1);
        if min(ele_x_coords) <= x && max(ele_x_coords) >= x
            mid_plane_elements_with_good_x(end+1) = mid_plane_elements(ele_ind);
        end
    end
    n_mid_plane_elements_with_good_x = length(mid_plane_elements_with_good_x);

    for z=z_space
        z_ind = z_ind + 1;

        found = false;
        ele_ind = 1;
        coord = [x 0.0 z];

        while ~found && ele_ind <= n_mid_plane_elements_with_good_x
            ele = mid_plane_elements_with_good_x(ele_ind);

            if min(nodes(connectivity(ele,:),3)) <= z && max(nodes(connectivity(ele,:),3)) >= z
                if pt_in_tet(nodes(connectivity(ele,:),:), coord) == 1
                    found = true;
                end

            end

            ele_ind = ele_ind + 1;
        end
        
        if found
            cross_section_stress(z_ind, x_ind, :) = elements_stress(ele,:);
            cross_section_strain(z_ind, x_ind, :) = elements_strain(ele,:);
            cross_section_von_mises(z_ind, x_ind, :) = elements_von_mises(ele,:);
        else
            cross_section_stress(z_ind, x_ind, :) = [nan, nan, nan, nan, nan, nan];
            cross_section_strain(z_ind, x_ind, :) = [nan, nan, nan, nan, nan, nan];
            cross_section_von_mises(z_ind, x_ind, :) = nan;
        end
    end
end

clear ele global_ind_start global_ind_end displacment_sub_vec node xx yy zz xy yz xz precision delta_x delta_z n_mid_plane_elements ele_x_coords
clear norm_delta_x norm_delta_z norm_precision x_precision z_precision x_space z_space x z x_ind z_ind coord found ele ele_ind mid_plane_elements


%% Visualization
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

% Stress XX
figure(7)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,1), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("X Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxx (MPa)';
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
saveas(gcf,'stress_xx.png',"png")

% Stress YY
figure(8)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,2), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Y Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σyy (MPa)';
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
saveas(gcf,'stress_yy.png',"png")

% Stress ZZ
figure(9)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,3), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Z Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σzz (MPa)';
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
saveas(gcf,'stress_zz.png',"png")

% Stress XY
figure(10)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,4), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("XY Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxy (MPa)';
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
saveas(gcf,'stress_xy.png',"png")

% Stress YZ
figure(11)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,5), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("YZ Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σyz (MPa)';
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
saveas(gcf,'stress_yz.png',"png")

% Stress XZ
figure(12)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,6), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("XZ Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxz (MPa)';
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
saveas(gcf,'stress_xz.png',"png")

% Von Mises Stress
figure(13)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_von_mises(:,1), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Von Mises Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σv (MPa)';
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
saveas(gcf,'von_mises.png',"png")

% Von Mises Stress cross-section
figure(14)
colormap jet
s = surf(1000.0*X,1000.0*Z,1.0e-6*cross_section_von_mises);
s.EdgeAlpha = 0.0;
s.FaceColor = 'interp';
title("Von Mises Stress Y=0",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σv (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
xlabel('x [mm]')
ylabel('z [mm]')
set(gcf,'Position',[1098 55 750 1000])
ylim(1000.0*[min(nodes(:,3)) - 0.05*(max(nodes(:,3)) - min(nodes(:,3))), max(nodes(:,3)) + 0.05*(max(nodes(:,3)) - min(nodes(:,3)))])
view(0,90);
axis equal
saveas(gcf,'von_mises_cross_section.png',"png")

clear min_z max_z node p0 p1 c h x_fixed_scatter y_fixed_scatter z_fixed_scatter norm_node_loading s


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

function result = pt_in_tet(tet,pt)
        
    % Check if in interior
    result = 0;
    node_list = [1 2 3 4; 2 3 4 1; 1 2 4 3; 1 3 4 2];
    ind = 1;
    while result==0 && ind <= 4
        a = tet(node_list(ind,1),:);
        b = tet(node_list(ind,2),:);
        c = tet(node_list(ind,3),:);
        d = tet(node_list(ind,4),:);
        if same_side_3D(pt,a,b,c,d) && same_side_3D(pt,b,c,d,a) && same_side_3D(pt,c,d,a,b) && same_side_3D(pt,d,a,b,c)
            result = 1;
        end
        ind = ind + 1;
    end

end

function val = same_side_3D(p1,p2,a,b,c)
    normal = cross(b-a, c-a);
    d1 = dot(normal, (p1-a));
    d2 = dot(normal, (p2-a));
    if sign(d1) == sign(d2)
        val = 1;
    else
        val = 0;
    end
end