% Get relevant min coordinates for integration
[min_x,argmin_x] = min(coords(:,1));
min_y = coords(argmin_x,2);

% Get relevant max coordinates for integration
[max_x,argmax_x] = max(coords(:,1));
max_y = coords(argmax_x,2);

% Get relevant mid coordinates for integration
argmid_x = ((argmin_x ~= [1,2,3]) + (argmax_x ~= [1,2,3])-1)*[1;2;3];
mid_x = coords(argmid_x,1);
mid_y = coords(argmid_x,2);

% Get relevant slopes for integration
m1 = (mid_y - min_y) / (mid_x - min_x);
m2 = (max_y - min_y) / (max_x - min_x);
m3 = (max_y - mid_y) / (max_x - mid_x);

% Integrate shape functions in two parts (1/2)
min_2_mid_N1_int = (-c1*(m1-m2)*(mid_x-min_x)^2*(3*N1(3)*top_z-N1(1)*(3*P4(1)-2*mid_x-min_x)+N1(2)*(m1*(mid_x-min_x)+m2*(mid_x-min_x)-3*(P4(2)-min_y))-3*N1(3)*P4(3)))/6.0;
min_2_mid_N2_int = (-c2*(m1-m2)*(mid_x-min_x)^2*(3*N2(3)*top_z-N2(1)*(3*P4(1)-2*mid_x-min_x)+N2(2)*(m1*(mid_x-min_x)+m2*(mid_x-min_x)-3*(P4(2)-min_y))-3*N2(3)*P4(3)))/6.0;
min_2_mid_N3_int = (-c3*(m1-m2)*(mid_x-min_x)^2*(3*N3(3)*top_z-N3(1)*(3*P4(1)-2*mid_x-min_x)+N3(2)*(m1*(mid_x-min_x)+m2*(mid_x-min_x)-3*(P4(2)-min_y))-3*N3(3)*P4(3)))/6.0;
min_2_mid_N4_int = (-c4*(m1-m2)*(mid_x-min_x)^2*(3*N4(3)*top_z-N4(1)*(3*P3(1)-2*mid_x-min_x)+N4(2)*(m1*(mid_x-min_x)+m2*(mid_x-min_x)-3*(P3(2)-min_y))-3*N4(3)*P3(3)))/6.0;

% Integrate shape functions in two parts (2/2)
mid_2_max_N1_int = (-c1*(m2-m3)*(max_x-mid_x)^2*(3*N1(3)*top_z-N1(1)*(3*P4(1)-max_x-2*mid_x)-N1(2)*(m2*(max_x-mid_x)+m3*(max_x-mid_x)+3*(P4(2)-max_y))-3*N1(3)*P4(3)))/6.0;
mid_2_max_N2_int = (-c2*(m2-m3)*(max_x-mid_x)^2*(3*N2(3)*top_z-N2(1)*(3*P4(1)-max_x-2*mid_x)-N2(2)*(m2*(max_x-mid_x)+m3*(max_x-mid_x)+3*(P4(2)-max_y))-3*N2(3)*P4(3)))/6.0;
mid_2_max_N3_int = (-c3*(m2-m3)*(max_x-mid_x)^2*(3*N3(3)*top_z-N3(1)*(3*P4(1)-max_x-2*mid_x)-N3(2)*(m2*(max_x-mid_x)+m3*(max_x-mid_x)+3*(P4(2)-max_y))-3*N3(3)*P4(3)))/6.0;
mid_2_max_N4_int = (-c4*(m2-m3)*(max_x-mid_x)^2*(3*N4(3)*top_z-N4(1)*(3*P3(1)-max_x-2*mid_x)-N4(2)*(m2*(max_x-mid_x)+m3*(max_x-mid_x)+3*(P3(2)-max_y))-3*N4(3)*P3(3)))/6.0;

% Sum integration parts and assemble
N1_int = abs(min_2_mid_N1_int) + abs(mid_2_max_N1_int);
N2_int = abs(min_2_mid_N2_int) + abs(mid_2_max_N2_int);
N3_int = abs(min_2_mid_N3_int) + abs(mid_2_max_N3_int);
N4_int = abs(min_2_mid_N4_int) + abs(mid_2_max_N4_int);
N_int = [N1_int 0 0 N2_int 0 0 N3_int 0 0 N4_int 0 0;
         0 N1_int 0 0 N2_int 0 0 N3_int 0 0 N4_int 0;
         0 0 N1_int 0 0 N2_int 0 0 N3_int 0 0 N4_int]';

% Calculate the element load vector given constant traction
load = N_int * [ 0.0; 0.0; top_traction ];