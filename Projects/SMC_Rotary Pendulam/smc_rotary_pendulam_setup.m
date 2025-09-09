%% Motor
% Resistance
Rm = 8.4;
% Current-torque (N-m/A)
kt = 0.042;
% Back-emf constant (V-s/rad)
km = 0.042;
%
%% Rotary Arm
% Mass (kg)
mr = 0.095;
% Total length (m)
r = 0.085;
% Moment of inertia about pivot (kg-m^2)
Jr = mr*r^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
br = 1e-3; % damping tuned heuristically to match QUBE-Sero 2 response
%
%% Pendulum Link
% Mass (kg)
mp = 0.024;
% Total length (m)
Lp = 0.129;
% Pendulum center of mass (m)
l = Lp/2;
% Moment of inertia about pivot (kg-m^2)
Jp = mp*Lp^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
bp = 5e-5;
% Gravity Constant
g = 9.81;

beta = 1;

delta = 1;

syms x_1 x_2 x_3 x_4 vm;

a_1 = Jr + (Jp*sin((x_2*x_2)));
a_2 = mp * l* r*cos(x_2);

b_1 = ((km*vm)/(Rm))  + (mp*l*r*(sin(x_2))*x_4*x_4) - (2*Jp*(sin(x_2))*(cos(x_2))*x_3*x_4) - (br*x_3) - (   (km*km*x_3)/Rm);

a_4 = Jp;
a_3 = mp*l*r*cos(x_2);
b_2 = (Jp*(sin(x_2))*(cos(x_2))*x_3*x_3) - (mp*g*l*sin(x_2)) - (bp*x_4 );

X = [a_1,a_2;a_3,a_4];
Y = [b_1;b_2];

D = linsolve(X,Y);

f_1 = x_3;
f_2 = x_4;
f_3 = D(1);
f_4 = D(2);

J = jacobian([f_1,f_2,f_3,f_4],[x_1,x_2,x_3,x_4,vm]);

x_1_val = 0;
x_2_val = pi;
x_3_val = 0;
x_4_val = 0;


J_subs = vpa(subs(J,[x_1, x_2, x_3, x_4,vm], [x_1_val, x_2_val, x_3_val, x_4_val,0]),4);

A = double(J_subs(:,1:4));

B = double(J_subs(:,5));

Q_1 = [2.5,0,0,0,;0,2,0,0;0,0,0.5,0;0,0,0,0.5;]

R_1= 2.5 ;


K_lin=lqr(A,B,Q_1,R_1);%
C = [1,0,0,0;0,1,0,0];

D = [0,0,0,0];


Closed_loop_sys = (A-(B*K_lin));

[eig_vec,eig_val] = eig(transpose(Closed_loop_sys));

% On the sliding surface the dynamics becomes

sliding_surface = transpose(eig_vec(:,2))/(transpose(eig_vec(:,2))*B);

dynamics = A - (B*sliding_surface*A);
eig_va = eig(dynamics);

K_nom = -(sliding_surface*A);





