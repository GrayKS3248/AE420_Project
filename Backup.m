%%Get mean circumference as a function of height and convert to a CDF
mean_circ_h = zeros(floor((h_end-h_start)/0.01)+1,2);
height_index = 1;
for h = h_start:0.01:h_end
    mean_r_at_h = 0.0;
    count = 0.0;
    for theta = 0.0:0.01:2.0*pi
        mean_r_at_h = mean_r_at_h + get_r(h, theta, arg);
        count = count + 1.0;
    end
    mean_r_at_h = mean_r_at_h / count;
    mean_circ_h(height_index,1) = h;
    mean_circ_h(height_index,2) = 2.0*pi*mean_r_at_h;
    height_index = height_index + 1;
end
mean_circ_h(:,1) = mean_circ_h(:,1);
mean_circ_h(:,2) = mean_circ_h(:,2) / trapz(mean_circ_h(:,1),mean_circ_h(:,2));
cdf = cumtrapz(mean_circ_h(:,1),mean_circ_h(:,2));

%% Triangulate surface
nodes_per_height = length(theta_track);

% Tip
triangle = 1;
for tri = 1:nodes_per_height
    if tri == nodes_per_height
        surface_triangulation(tri,:) = [1,2+tri-1,2];
        surface_element_type(tri,:) = [type_ext(1), type_ext(2+tri-1), type_ext(2)];
    else
        surface_triangulation(tri,:) = [1,2+tri-1,3+tri-1];
        surface_element_type(tri,:) = [type_ext(1), type_ext(2+tri-1), type_ext(3+tri-1)];
    end
    triangle = triangle + 1;
end

% Body
for level = 2:(length(0.0:precision:1.0)-2)
    % Upper facing triangles
    for tri = 1:nodes_per_height
        ind_of_right_angle = tri+1+(level-2)*nodes_per_height;
        if tri == nodes_per_height
            surface_triangulation(triangle,:) = [ind_of_right_angle,ind_of_right_angle-nodes_per_height+1,ind_of_right_angle+nodes_per_height];
            surface_element_type(triangle,:) = [type_ext(ind_of_right_angle), type_ext(ind_of_right_angle-nodes_per_height+1), type_ext(ind_of_right_angle+nodes_per_height)];
        else
            surface_triangulation(triangle,:) = [ind_of_right_angle,ind_of_right_angle+1,ind_of_right_angle+nodes_per_height];
            surface_element_type(triangle,:) = [type_ext(ind_of_right_angle), type_ext(ind_of_right_angle+1), type_ext(ind_of_right_angle+nodes_per_height)];
        end
        triangle = triangle + 1;
    end
    
    % Lower facing triangles
    for tri = 1:nodes_per_height
        if tri == nodes_per_height
            ind_of_right_angle = tri+2+(level-1)*nodes_per_height - nodes_per_height;
            surface_triangulation(triangle,:) = [ind_of_right_angle, ind_of_right_angle-nodes_per_height, ind_of_right_angle-1+nodes_per_height];
            surface_element_type(triangle,:) = [type_ext(ind_of_right_angle), type_ext(ind_of_right_angle-nodes_per_height), type_ext(ind_of_right_angle-1+nodes_per_height)];
        else
            ind_of_right_angle = tri+2+(level-1)*nodes_per_height;
            surface_triangulation(triangle,:) = [ind_of_right_angle, ind_of_right_angle-nodes_per_height, ind_of_right_angle-1];
            surface_element_type(triangle,:) = [type_ext(ind_of_right_angle), type_ext(ind_of_right_angle-nodes_per_height), type_ext(ind_of_right_angle-1)];
        end
        triangle = triangle + 1;
    end
end

% Head taper
for tri = 1:nodes_per_height
    if tri == nodes_per_height
        surface_triangulation(triangle,:) = [length(nodes_ext),length(nodes_ext)-tri,length(nodes_ext)-1];
        surface_element_type(triangle,:) = [type_ext(length(nodes_ext)), type_ext(length(nodes_ext)-tri), type_ext(length(nodes_ext)-1)];
    else
        surface_triangulation(triangle,:) = [length(nodes_ext),length(nodes_ext)-tri,length(nodes_ext)-tri-1];
        surface_element_type(triangle,:) = [type_ext(length(nodes_ext)), type_ext(length(nodes_ext)-tri), type_ext(length(nodes_ext)-tri-1)];
    end
    triangle = triangle + 1;
end

%% Determine coloration of exterior nodes
C_ext = zeros(length(nodes_ext),1);
for i = 1:length(nodes_ext)
     [theta, r, h] = cart2pol(nodes_ext(i,1), nodes_ext(i,2), nodes_ext(i,3));
     if theta < 0.0
         theta = theta + 2.0 * pi;
     end
     nodes_ext_pol(i,:) = [r, rad2deg(theta), h];
     C_ext(i) = type_ext(i)/8.0;
end