%% Setup
close all;
clear all;


%% Parameters
% Radii
r_head = 6.0;
r_shank = 2.5;
r_body = 2.0;

% Lengths
l_head_taper = 1.0;
l_head = 4.0;
l_shank_transition = 1.0;
l_shank = 4.0;
l_body = 40.0;
l_tip = 5.0;

% Thread
thread_number = 0.0;
thread_depth = 0.0;
thread_width = 0.0;
thread_crest = 0.0;

% Mesh
num_height = 20;
num_ring = 2;
num_theta = 5;


%% Compile parameters
arg = struct('rh',r_head,...
    'rs',r_shank,...
    'rb',r_body,...
    'lht',l_head_taper,...
    'lh',l_head,...
    'lst',l_shank_transition,...
    'ls',l_shank,...
    'lb',l_body,...
    'lt',l_tip,...
    'tn',thread_number,...
    'td',thread_depth,...
    'tw',thread_width,...
    'tc',thread_crest);


%% Generate nodes
% Calculate general parameters
total_length = arg.lt + arg.lb + arg.ls + arg.lst + arg.lh + arg.lht;
num_nodes = num_height*num_ring*num_theta;

% Calculate height, ring, and theta space
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
node = 1;
for height = height_list
    for radius = radius_list
        for theta = theta_list
            
            % Get radius, convert to cartesian
            [r, node_type] = get_r(theta,height,arg);
            r = r * radius;
            if radius ~= 1.0 && radius ~= 0.0
                [~,inner_ring_number] = min(abs(radius_list-radius));
                inner_ring_number = inner_ring_number - 1;
                t = inner_ring_number*(theta+pi/num_theta);
            else
                t = theta;
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
TYPE = TYPE(ia);
num_nodes = length(NODES);


%% Element generation
ELEMENTS = [1,2,3,4];


%% Visualization
%scatter3(X,Y,Z,'.','k');
tetramesh(ELEMENTS, NODES);
axis equal;
view(30,15);

%% Define body geometry
function [r, node_type] = get_r(theta,h,arg)
    
    % Calculate relevant heights
    tip_2_body = arg.lt;
    body_2_shank = tip_2_body + arg.lb;
    shank_2_shank_transition = body_2_shank + arg.ls;
    shank_transition_2_head = shank_2_shank_transition + arg.lst;
    head_2_head_taper = shank_transition_2_head + arg.lh;
    total_length = head_2_head_taper + arg.lht;

    % Head taper
    if h>=head_2_head_taper && h<=total_length
        node_type = 8;
        loc_in_section = (h - head_2_head_taper) / (total_length - head_2_head_taper);
        r = arg.rh*(1.0-loc_in_section);
        
    % Head
    elseif h>=shank_transition_2_head && h<=head_2_head_taper
        node_type = 7;
        r = arg.rh;
        
    % Shank transition
    elseif h>=shank_2_shank_transition && h<=shank_transition_2_head
        node_type = 6;
        loc_in_section = (h - shank_2_shank_transition) / (shank_transition_2_head - shank_2_shank_transition);
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
    thread_height = 0.0;
    node_type = 0;
end
