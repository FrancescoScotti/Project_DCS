% DCS progetto: "Controllo della glicemia nei pazienti diabetici
%                di tipo 1"
%
% A.A. 2020/2021
%
% Sara Rossi
% Francesco Scotti 

clear


%%%%%%% Continous Time %%%%%%

s = tf('s');

theta1 = 0.21;
theta2 = 17.5;
theta3 = 70;

Gr = 120;

A = [0 -theta2 0; 0 -1/theta3 1/theta3; 0 0 -1/theta3];
B = [0; 0; 1/theta3];
C = [1 0 0];

% [V,J] = jordan(A);

SYSC = ss(A,B,C,0);

P_s = tf(SYSC);

% SYSC = chgTimeUnit(SYSC,'minutes');

%%%%%% State Feedback %%%%%%

k = 0.001;
F_k = k*[1/theta2 -theta3 -theta3];



%%%%%% Discrete Time %%%%%%

Ts = 300; % sampling time (problema con scratch_luenberger)

SYSD_zoh = c2d(SYSC,Ts,'zoh'); % zoh

SYSD_foh = c2d(SYSC,Ts,'foh'); % foh

SYSD_tustin = c2d(SYSC,Ts,'tustin'); % tustin

SYSD_matched = c2d(SYSC,Ts,'matched'); %forward;



rank(obsv(A,C)) %check observability



%%%%%% Luenberger Continous %%%%%%

disp("Eigenvalues of original A")
eig(A)
disp(' ')

disp("Eigenvalues of original A+B*F_k")
eig(A+B*F_k)
disp(' ')

disp("Eigenvalues I want to impose to my observer")
lambda = [-10*k -1/theta3 -1/theta3] 
disp(' ')

disp("Observer gain")
L = acker((A+B*F_k)',C',[lambda(1), lambda(2), lambda(3)])'
disp(' ')

disp("Check")
eig((A+B*F_k)-L*C)



%%%%%%%%%%%% ZOH %%%%%%%%%%

A_zoh = SYSD_zoh.A;
B_zoh = SYSD_zoh.B;
C_zoh = SYSD_zoh.C;

disp("Eigenvalues of original Ad")
eig(A_zoh)
disp(' ')

disp("Eigenvalues of original Ad+Bd*F_k")
eig(A_zoh+B_zoh*F_k)
disp(' ')

disp("Eigenvalues I want to impose to my observer")
lambda_zoh = exp(lambda*Ts)
disp(' ')

disp("Observer gain")
L_zoh = acker((A_zoh+B_zoh*F_k)',C_zoh',[lambda_zoh(1), lambda_zoh(2), lambda_zoh(3)])'
disp(' ')

disp("Check")
eig((A_zoh+B_zoh*F_k)-L_zoh*C_zoh)

%%%%%%%%%%%% TUSTIN %%%%%%%%%%

A_tustin = SYSD_tustin.A;
B_tustin = SYSD_tustin.B;
C_tustin = SYSD_tustin.C;

disp("Eigenvalues of original Ad_tustin")
eig(A_tustin)
disp(' ')

disp("Eigenvalues of original Ad+Bd*F_k")
eig(A_tustin+B_tustin*F_k)
disp(' ')

disp("Eigenvalues I want to impose to my observer")
lambda_tustin = exp(lambda*Ts)
disp(' ')

disp("Observer gain")
L_tustin = acker((A_tustin+B_tustin*F_k)',C_tustin',[lambda_tustin(1), lambda_tustin(2), lambda_tustin(3)])'
disp(' ')

disp("Check")
eig((A_tustin+B_tustin*F_k)-L_tustin*C_tustin)

%%%%%%%%%%%% FOH %%%%%%%%%%

A_foh = SYSD_foh.A;
B_foh = SYSD_foh.B;
C_foh = SYSD_foh.C;

disp("Eigenvalues of original Ad")
eig(A_foh)
disp(' ')

disp("Eigenvalues of original Ad+Bd*F_k")
eig(A_foh+B_foh*F_k)
disp(' ')

disp("Eigenvalues I want to impose to my observer")
lambda_foh = exp(lambda*Ts)
disp(' ')

disp("Observer gain")
L_foh = acker((A_foh+B_zoh*F_k)',C_zoh',[lambda_foh(1), lambda_foh(2), lambda_foh(3)])'
disp(' ')

disp("Check")
eig((A_foh+B_foh*F_k)-L_foh*C_foh)


