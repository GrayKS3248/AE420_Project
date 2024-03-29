%% Setup
close all
clear


%% Parameters
% Radii 
r_head = 6.5;
r_shank = 3.0;
r_body = 2.5;

% Lengths
l_head = 2.75;
l_shank_transition = 3.0;
l_shank = 1.5;
l_body = 42.0;
l_tip = 5.0;

% Thread
thread_number = 0.5;
thread_depth = 0.75;
thread_width = 1.25;
thread_crest = 0.5;

% Options
fillet_shank = false;

% Mesh
num_height = 138;
num_ring = 6;
num_theta = 43;

% Render options
render_mesh = false;
rotations = 1.0;
color_cycles = 10.0;
    
% Visualization options
vis = true;


%% Compile parameters
arg = struct('rh',r_head,...
    'rs',r_shank,...
    'rb',r_body,...
    'lh',l_head,...
    'lst',l_shank_transition,...
    'ls',l_shank,...
    'lb',l_body,...
    'lt',l_tip,...
    'tn',thread_number,...
    'td',thread_depth,...
    'tw',thread_width,...
    'tc',thread_crest,...
    'fillet_shank',fillet_shank);


%% Generate nodes
% Calculate general parameters
total_length = arg.lt + arg.lb + arg.ls + arg.lst + arg.lh;
num_nodes = num_height*num_ring*num_theta;

% Calculate height, radius, and theta space
height_list = linspace(0.0,total_length,num_height);
radius_list = linspace(1.0,0.0,num_ring+2);
theta_list = linspace(0.0,2.0*pi,num_theta+1);
theta_list = theta_list(1:end-1);

X = zeros(num_nodes, 1);
Y = zeros(num_nodes, 1);
Z = zeros(num_nodes, 1);
R = zeros(num_nodes, 1);
T = zeros(num_nodes, 1);
H = zeros(num_nodes, 1);
TYPE = zeros(num_nodes, 1);
height_level = 0;
node = 1;
for height = height_list
    height_level = height_level + 1;
    for radius = radius_list
        for theta = theta_list

            % Get radius, convert to cartesian
            [r, node_type] = get_r(theta,height,arg);
            r = r * radius;
            if radius ~= 1.0 && radius ~= 0.0
                [~,inner_ring_number] = min(abs(radius_list-radius));
                inner_ring_number = inner_ring_number - 1;
                if mod(inner_ring_number,2)~=0
                    t = theta+pi/num_theta;
                else
                    t = theta;
                end
            else
                t = theta;
            end
            if t > 2.0*pi
                t = t - 2.0*pi;
            elseif t < 0.0
                t = t + 2.0*pi;
            end
            [x,y,z] = pol2cart(t, r, height);
            X(node) = x;
            Y(node) = y;
            Z(node) = z;
            R(node) = r;
            T(node) = t;
            H(node) = height;
            
            % Assign type to surface nodes only
            if radius == 1.0
                TYPE(node) = node_type;
            else
                TYPE(node) = -1;
            end
            node = node + 1;
            
        end
    end
end

% Remove duplicate nodes
NODES = [X Y Z];
[NODES, ia] = unique(NODES,'rows','stable');
X = NODES(:,1);
Y = NODES(:,2);
Z = NODES(:,3);
R = R(ia);
T = T(ia);
H = H(ia);
CYL_NODES = [R T H];
TYPE = TYPE(ia);
num_nodes = length(NODES);


%% Element generation
element = 1;

% Tip
for i = 1:num_ring
    for j = 1:num_theta
        
        % Inward pointing elements
        node_1 = 1 + j + (i-1)*num_theta;
        if j ~= num_theta
            node_2 = node_1 + 1;
        else
            node_2 = node_1 + 1 - num_theta;
        end
        if mod(i,2) == 0
            node_3 = node_2 + num_theta;
        else
            node_3 = node_1 + num_theta;
        end
        ELEMENTS(element,:) = [1,node_1,node_2,node_3];
        element = element + 1;
        
        % Outward pointing elements
        node_4 = node_1;
        if j==1
            node_5 = node_3 - 1 + num_theta - mod(i+1,2)*num_theta;
        elseif j==num_theta
            node_5 = node_3 - 1 + num_theta - mod(i,2)*num_theta;
        else
            node_5 = node_3 - 1;
        end
        node_6 = node_3;
        ELEMENTS(element,:) = [1,node_4,node_5,node_6];
        element = element + 1;
        
        % Core elements
        if i==num_ring
            node_7 = node_5;
            node_8 = node_6;
            ELEMENTS(element,:) = [1,node_7,node_8,(num_ring+1)*num_theta+2];
            element = element + 1;
        end
        
    end
end

% Body and head taper
for k = 1:num_height-2
    for i = 1:num_ring
        
        % Calculate relevant element indices
        ele_h_1 = (num_ring+1)*num_theta+3 + (k-2)*((num_ring+1)*num_theta + 1) + (i-1)*(num_theta);
        ele_h_2 = (num_ring+1)*num_theta+3 + (k-1)*((num_ring+1)*num_theta + 1) + (i-1)*(num_theta);
        
        for j = 1:num_theta

            % Upward pointing elements of inward pointing origin
            node_1 = ele_h_1 + j - 1;
            if j ~= num_theta
                node_2 = node_1 + 1;
            else
                node_2 = node_1 + 1 - num_theta;
            end
            node_3 = node_1 + num_theta + mod(i+1,2);
            if mod(i,2)==0 && j == num_theta
                node_3 = node_3 - num_theta;
            end
            node_4 = ele_h_2 + j - 1;
            ELEMENTS(element,:) = [node_1,node_2,node_3,node_4];
            element = element + 1;

            % Downward pointing elements of inward pointing origin
            node_5 = node_2;
            if j~= num_theta
                node_6 = ele_h_2+j;
            else
                node_6 = ele_h_2;
            end
            node_7 = node_4;
            node_8 = ele_h_2+num_theta+j-1+mod(i+1,2);
            if mod(i,2)==0 && j==num_theta
                node_8 = node_8-num_theta;
            end
            ELEMENTS(element,:) = [node_5,node_6,node_7,node_8];
            element = element + 1;
            
            % Inner elements of inward pointing origin
            node_9 = node_2;
            node_10 = node_3;
            node_11 = node_8;
            node_12 = node_4;
            ELEMENTS(element,:) = [node_9,node_10,node_11,node_12];
            element = element + 1;

            % Upward pointing elements of outward pointing origin
            node_13 = node_1;
            node_14 = node_3;
            node_15 = node_4;
            if j==1
                node_16 = node_8 - 1 + num_theta - num_theta*mod(i+1,2);
            elseif j==num_theta && mod(i,2)==0
                node_16 = node_8 - 1 + num_theta;
            else
                node_16 = node_8 - 1;
            end
            ELEMENTS(element,:) = [node_13,node_14,node_15,node_16];
            element = element + 1;
            
            % Downward pointing elements of outward pointing origin
            node_17 = node_13;
            node_18 = node_14;
            node_19 = node_16;
            if j==1
                node_20 = node_14 - 1 + num_theta - num_theta*mod(i+1,2);
            elseif j==num_theta && mod(i,2)==0
                node_20 = node_14 - 1 + num_theta;
            else
                node_20 = node_14 - 1;
            end
            ELEMENTS(element,:) = [node_17,node_18,node_19,node_20];
            element = element + 1;
            
            % Inner elements of outward pointing origin
            node_21 = node_14;
            node_22 = node_15;
            node_23 = node_16;
            if j==num_theta && mod(i,2)==0
                node_24 = node_15 + mod(i+1,2);
            else
                node_24 = node_15 + num_theta + mod(i+1,2);
            end
            ELEMENTS(element,:) = [node_21,node_22,node_23,node_24];
            element = element + 1;
            
            if i==num_ring
                %Core elements 1
                node_25 = node_18;
                node_26 = node_19;
                node_27 = node_20;
                node_28 = ele_h_1 + 2*num_theta;
                ELEMENTS(element,:) = [node_25,node_26,node_27,node_28];
                element = element + 1;

                %Core elements 2
                node_29 = node_25;
                node_30 = node_26;
                node_31 = node_24;
                node_32 = ele_h_1 + 2*num_theta;
                ELEMENTS(element,:) = [node_29,node_30,node_31,node_32];
                element = element + 1;
                
                %Core elements 3
                node_33 = ele_h_2 + 2*num_theta;
                node_34 = node_30;
                node_35 = node_31;
                node_36 = node_32;
                ELEMENTS(element,:) = [node_33,node_34,node_35,node_36];
                element = element + 1;
            end
        end
    end
end


%% Post-Process elements 
num_elements = length(ELEMENTS);
ELEMENT_TRIANGULATION = triangulation(ELEMENTS, X, Y, Z);
[ELEMENT_CIRCUMCENTER, ELEMENT_CIRCUMRADIUS] = circumcenter(ELEMENT_TRIANGULATION);
ELEMENT_VOLUME = zeros(num_elements, 1);
ELEMENT_INRADIUS = zeros(num_elements, 1);
ELEMENT_INCENTER = zeros(num_elements, 3);
ELEMENT_AR = zeros(num_elements, 1);
for i = 1:num_elements
    % Volume calculation
    ele_node = [X(ELEMENTS(i,:)), Y(ELEMENTS(i,:)), Z(ELEMENTS(i,:))];
    ele_vol  = [ele_node(1,:) - ele_node(4,:); ele_node(2,:) - ele_node(4,:); ele_node(3,:) - ele_node(4,:)];
    ELEMENT_VOLUME(i) = abs(det(ele_vol)) / 6.0;
    
    % Inradius calculation
    ele_norm_dirns = [ele_node(4,:)-ele_node(1,:);
                      ele_node(3,:)-ele_node(1,:);
                      ele_node(1,:)-ele_node(2,:);
                      ele_node(2,:)-ele_node(1,:)];
    ele_norms = [cross(ele_node(1,:) - ele_node(2,:), ele_node(2,:) - ele_node(3,:));
                 cross(ele_node(1,:) - ele_node(2,:), ele_node(2,:) - ele_node(4,:));
                 cross(ele_node(2,:) - ele_node(3,:), ele_node(3,:) - ele_node(4,:));
                 cross(ele_node(1,:) - ele_node(4,:), ele_node(4,:) - ele_node(3,:))];
    ele_norms = (ele_norms ./ vecnorm(ele_norms,2,2)) .* (-1*sign(dot(ele_norm_dirns,ele_norms,2)));
    A = [1.0, ele_norms(1,:);
         1.0, ele_norms(2,:);
         1.0, ele_norms(3,:);
         1.0, ele_norms(4,:);];
    B = [ele_norms(1,:)*ele_node(1,:)';
         ele_norms(2,:)*ele_node(2,:)';
         ele_norms(3,:)*ele_node(3,:)';
         ele_norms(4,:)*ele_node(4,:)';];
    soln = A\B;
    ELEMENT_INRADIUS(i) = soln(1);
    ELEMENT_INCENTER(i,:) = soln(2:end)';
    
    % Aspect ratio calculation
    ELEMENT_AR(i) = 3.0*ELEMENT_INRADIUS(i) / ELEMENT_CIRCUMRADIUS(i);
end
MESH_AR = [mean(ELEMENT_AR), std(ELEMENT_AR), median(ELEMENT_AR), min(ELEMENT_AR), max(ELEMENT_AR)]
save('elements.mat','ELEMENT_TRIANGULATION');
save('type.mat','TYPE')
save('quality.mat','MESH_AR');

%% Visualization
if render_mesh
    colormap hsv;
    cmap = colormap;
    color = zeros(num_elements,3);
    for i=1:num_elements
        color(i,:) = cmap(floor(mod((256*color_cycles*((i-1)/(num_elements-1))),256))+1,:);
    end
    set(gcf,'Position',[10 50 1100 1010])
    for i=1:num_elements
        tetramesh(ELEMENTS(i,:),NODES,'FaceColor',color(i,:),'FaceAlpha','1.0','EdgeAlpha', '1.0','EdgeColor','k');
        hold on
        xlim([-1.05*max([r_head,r_shank,r_body]), 1.05*max([r_head,r_shank,r_body])])
        ylim([-1.05*max([r_head,r_shank,r_body]), 1.05*max([r_head,r_shank,r_body])])
        zlim([-0.05*total_length, 1.05*total_length])
        view(rad2deg(2.0*pi - rotations*2.0*pi*((i-1)/(num_elements-1)))+5.0,50+30*((i-1)/(num_elements-1)));
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'ZTickLabel',[]);
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        h.ZAxis.TickLength = [0 0];
        grid on
        axis vis3d
        saveas(gcf,strcat("mesh_",num2str(i),'.png'),"png")
    end
    for i=num_elements+1:num_elements+601
        xlim([-1.05*max([r_head,r_shank,r_body]), 1.05*max([r_head,r_shank,r_body])])
        ylim([-1.05*max([r_head,r_shank,r_body]), 1.05*max([r_head,r_shank,r_body])])
        zlim([-0.05*total_length, 1.05*total_length])
        view(rad2deg(2.0*pi - rotations*2.0*pi*((i-1)/(num_elements-1)))+5.0,80-50*((i-num_elements-1)/600));
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'ZTickLabel',[]);
        h=gca;
        h.XAxis.TickLength = [0 0];
        h.YAxis.TickLength = [0 0];
        h.ZAxis.TickLength = [0 0];
        grid on
        axis vis3d
        saveas(gcf,strcat("mesh_",num2str(i),'.png'),"png")
    end
end
close all

if vis
    figure(1)
    trisurf(ELEMENT_TRIANGULATION.freeBoundary, X, Y, Z, 'EdgeAlpha', '0.25', 'FaceColor', [0.75, 0.75, 0.75]);
    title("Surface Mesh")
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    axis equal;
    view(30,30);
    set(gcf,'Position',[50 55 750 1000])
    saveas(gcf,'surf.png',"png")

    figure(2)
    histogram(ELEMENT_AR)
    title("Mesh AR")
    xlabel("Aspect Ratio")
    ylabel("Number of Elements")
    xlim([0.0,1.0])
    set(gcf,'Position',[825 55 1200 1000])
    saveas(gcf,'condition.png',"png")
end


%% Define body geometry
function [r, node_type] = get_r(theta,h,arg)
    
    % Calculate relevant heights
    tip_2_body = arg.lt;
    body_2_shank = tip_2_body + arg.lb;
    shank_2_shank_transition = body_2_shank + arg.ls;
    shank_transition_2_head = shank_2_shank_transition + arg.lst;
    total_length = shank_transition_2_head + arg.lh;
        
    % Head
    if h>=shank_transition_2_head && h<=total_length
        node_type = 7;
        r = arg.rh;
        
    % Shank transition
    elseif h>=shank_2_shank_transition && h<=shank_transition_2_head
        node_type = 6;
        loc_in_section = (h - shank_2_shank_transition) / (shank_transition_2_head - shank_2_shank_transition);
        if arg.fillet_shank
            loc_in_section = 3.0*loc_in_section*loc_in_section - 2.0*loc_in_section*loc_in_section*loc_in_section;
        end
        r = arg.rs*(1.0-loc_in_section) + arg.rh*loc_in_section;
        
    % Shank
    elseif h>=body_2_shank && h<=shank_2_shank_transition
        node_type = 5;
        r = arg.rs;
        
    % Body
    elseif h>=tip_2_body && h<=body_2_shank
        [thread_height, node_type] = get_thread_height(theta,h,arg);
        r = arg.rb + thread_height;
        
    % Tip
    elseif h>=0.0 && h<=tip_2_body
        [thread_height, node_type] = get_thread_height(theta,h,arg);
        loc_in_section = h/tip_2_body;
        r = (arg.rb+thread_height)*loc_in_section;
        
    % Outside of body
    else
        node_type = -1;
        r = 0.0;
        
    end
    
end

%% Define thread geometry
function [thread_height, node_type] = get_thread_height(theta,h,arg)

    % Calculate parameters
    nearest_thread = round( (2.0*pi*h*arg.tn - theta) / (2.0*pi) );
    mid_height = theta/(2.0*pi*arg.tn) + nearest_thread/arg.tn;
    min_height = mid_height - 0.5*arg.tw;
    max_height = mid_height + 0.5*arg.tw;
    crest_min_height = mid_height - 0.5*arg.tc;
    crest_max_height = mid_height + 0.5*arg.tc;
    
    % Approach
    if h>=min_height && h<crest_min_height
        node_type = 2;
        loc_in_section = (h - min_height) / (crest_min_height - min_height);
        thread_height = arg.td*loc_in_section;
    
    % Crest
    elseif h>=crest_min_height && h<=crest_max_height
        node_type = 3;
        thread_height = arg.td;
        
    
    % Receeding
    elseif h>crest_max_height && h<=max_height
        node_type = 4;
        loc_in_section = (h - crest_max_height) / (max_height - crest_max_height);
        thread_height = arg.td*(1.0-loc_in_section);
    
    % Body
    else
        node_type = 1;
        thread_height = 0.0;
    end
end