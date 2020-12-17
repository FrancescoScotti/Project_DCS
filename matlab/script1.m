% DCS progetto: "Controllo della glicemia nei pazienti diabetici
%                di tipo 1"
%
% A.A. 2020/2021
%
% Sara Rossi
% Francesco Scotti 

clear

% syms theta1 theta2 theta3 k s
% 
% assume(theta1,'real')
% assume(theta1,'positive')
% 
% assume(theta2,'real')
% assume(theta2,'positive')
% 
% assume(theta3,'real')
% assume(theta3,'positive')

s = tf('s');

theta1 = 0.21;
theta2 = 17.5;
theta3 = 70;

T = 5; % sampling time
k = 0.005;
A = [0 -theta2 0; 0 -1/theta3 1/theta3; 0 0 -1/theta3];
B = [0; 0; 1/theta3];
C = [1 0 0];

% [V,J] = jordan(A);

F_k = k*[1/theta2 -theta3 -theta3];



