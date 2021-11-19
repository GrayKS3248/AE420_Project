%% Load mesh
load('elements.mat')
nodes = ELEMENT_TRIANGULATION.Points;
connectivity = ELEMENT_TRIANGULATION.ConnectivityList;

%% Elementwise stiffness
