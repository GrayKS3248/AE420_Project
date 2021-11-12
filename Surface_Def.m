%% Setup
close all 
clear all

%% Parameters

r_head = 0.0;
r_shank = 0.0;
r_body = 0.0;

l_head_taper = 0.0;
l_head = 0.0;
l_shank_transition = 0.0;
l_shank = 0.0;
l_body = 0.0;
l_tip = 0.0;

thread_number = 0.0;
thread_depth = 0.0;
thread_width = 0.0;
thread_crest = 0.0;

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