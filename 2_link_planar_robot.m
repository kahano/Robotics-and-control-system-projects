%% 2-link planar Robot %%
syms L1 L2
syms T1 T2
syms pi

% Dynamic symbols %
syms g              % Gravity
syms q1_dot q2_dot  % Generalized coordinates
syms m1 m2          % Mass for each links 
%% Forward-kinematics %%
A_1 = A_matrix(T1,0,L1,0);
A_2 = A_matrix(T2,0,L2,0);
fk = simplify(A_1*A_2);
disp(fk);
%% Jacobian betwen coordinate 0 and coordinate 1 %%
z_0 = [0 0 1];
z_1 = A_1(1:3,3).';

O_0 = [0 0 0];
O_1 = A_1(1:3,4).';

Jv_1 = cross(z_0, (O_1-O_0));

Jw_1 = z_0;

Jv = [Jv_1];
Jw = [Jw_1];

J_01 = simplify([Jv.'; Jw.';]);

disp(J_01);

%% Jacobian betwen coordinate 0 and coordinate 2 %%
z_0 = [0 0 1];
z_1 = A_1(1:3,3).';

O_0 = [0 0 0];
O_1 = A_1(1:3,4).';
O_2 = fk(1:3,4).';

Jv_1 = cross(z_0, (O_2-O_0));
Jv_2 = cross(z_1, (O_2-O_1));

Jw_1 = z_0;
Jw_2 = z_1;

Jv = [Jv_1; Jv_2;];
Jw = [Jw_1; Jw_2;];

J_02 = simplify([Jv.'; Jw.';]);

disp(J_02);

%% Inertia Tensors %%
I_1 = [0 0 0;
       0 0 0;
       0 0 0;];
   
I_2 = [0 0 0;
       0 0 0;
       0 0 0;];

%% Kinetic Energy %%
q_dot = [q1_dot q2_dot];
v_1 = J_01(1:3)*q1_dot;
v_2 = J_02(1:3,:)*[q1_dot; q2_dot;];

w_1 = J_01(4:end)*q1_dot;
w_2 = J_02(4:end,:)*[q1_dot; q2_dot;]

K_1 = simplify(((1/2)*m1)*(v_1).'*(v_1)+(w_1).'*(I_1)*w_1);
K_2 = simplify(((1/2)*m2)*(v_2).'*(v_2)+(w_2).'*(I_2)*w_2);

K = simplify(K_1+K_2);
disp("K_1: ")
disp(K_1);
disp("K_2: ")
disp(K_2);
disp("Kinetic Energy-(K_1+K_2): ")
disp(K);

%% Potential Energy %%
% Assume uniform-distributed density mass %
P_1 = simplify(m1*g*((L1*sin(T1))/2));
P_2 = simplify(m2*g*(L1*sin(T1) + (L2*sin(T1+T2))/2));

P = simplify(P_1+P_2);
disp("P_1: ");
disp(P_1);
disp("P_2: ");
disp(P_2);
disp("Potential Energy-(P_1+P_2): ");
disp(P);
%% Lagrangian term %%
L = K-P;
disp("L: ");
disp(L);

%% END %%