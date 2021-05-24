%% ASSIGNMENT 1 : NCDS

% A= [ 5 9; 0 -3] B= [ 0; 1]
%
%% Question 1 

% In the first question we will assume a constant h and no delays tau
% 1. Place poles to -2 and -3

A= [ 5 8.5; 0 -3];
B= [ 0; 1];
K= place(A,B,[-2 -3]);

%2. Construct discrete-time model

% The restulting expressions for the discrete model:
% xk+1= e^(Ah)*xk + limintegral[h,0]e^(As)*B*ds*uk
%     = F(h)xk + G(h)uk
% Closing the loop:
% xk+1= (F(h)-G(h)*K) xk

syms h

A_h=[5*h 8.5*h; 0 -3*h];
F_h= expm(A_h);
G_h=int(expm(A_h),[0 h])*B;

A_bar= F_h-G_h*K;
A_bar= matlabFunction(A_bar);


%3. Analyze stability (check maximum eigenvalue norm) < 1
max_eig=[];

for i= 1:150
    h=0.005*(i-1)+0.0001;    
    max_eig(i)=max(abs(eig(A_bar(h))));
end

xaxis= linspace(0,0.005*150,150);
plot(xaxis,max_eig,'b', 'linewidth', 1.2)
yline(1,'r', 'linewidth', 1.2)
title('Spectral Radius of $\bar{A}$','Interpreter','latex')
legend('max($|\lambda_i|$)', '$\rho(\bar{A})$=1','Interpreter','latex')
ylabel('max($|\lambda_i|$)','Interpreter','latex')
xlabel('h','Interpreter','latex')

%% Question 2

% 1. Expression for the NCS with delays

% The restulting expressions for the discrete model:
% xk+1= e^(Ah)*xk + limintegral[h,h-tau]e^(As)*B*ds*uk-1+limintegral[h-tau,0]e^(As)*B*ds*uk
%     = Fx(h)xk +Fu(h,tau)uk-1 + G1(h,tau)uk
% Closing the loop and extending xek=[xk' uk-1']
% xek+1= (F(h,tau)-G(h,tau)*K) xke

syms tau
syms h 

Fx= expm(A_h);
Fu= int(expm(A_h),[h-tau h])*B;
G1= int(expm(A_h),[0 h-tau])*B;

F_tau= [Fx Fu; zeros(1,3)];
G_tau= [G1;eye(1)];
A_tau=F_tau-G_tau*[K 0];
F_tau_f= matlabFunction(F_tau);
G_tau_f= matlabFunction(G_tau);
A_tau= matlabFunction(A_tau);

% 2. Study the stability with the varying parameters tau and h and plot the
% results

max_eig_tau=[];

for i= 1:500
    for j=1:500
        h=0.001*(i-1);
        tau=0.001*(j-1);
        max_eig_tau(i,j)=max(abs(eig(A_tau(h,tau))));
    end
end

%% Plotting Q2.2
max_eig_tau=max_eig_tau';
xaxis= linspace(0,0.001*500,500);
yaxis= linspace(0,0.001*500,500);
copy_stable=zeros(500,500);
copy_unstable=zeros(200,200);

for i=1:500
    for j=1:500
        if xaxis(i) < yaxis(j)
            if max_eig_tau(i,j)<1
                copy_stable(j,i)=max_eig_tau(i,j);
            else
                copy_unstable(j,i)=max_eig_tau(i,j);
            end
        end
    end
end

%subplot(2,1,1)
copy_stable(copy_stable==0) = NaN;
surf(xaxis,yaxis,copy_stable')
colormap autumn
colorbar
ylabel('$\tau$','Interpreter','latex')
xlabel('h','Interpreter','latex')

title('Stable Combinations of $\{h,\tau\}$','Interpreter','latex')
legend('max($|\lambda_i|$)','Interpreter','latex')
%subplot(2,1,2)
%copy_unstable(copy_unstable==0) = NaN;
%surf(xaxis,yaxis,copy_unstable')
%title('Unstable Combinations of $\{h,\tau\}$','Interpreter','latex')
%ylabel('$\tau$','Interpreter','latex')
%xlabel('h','Interpreter','latex')
%legend('max($|\lambda_i|$)','Interpreter','latex')
%colorbar

%%
% 3. Tune the controller K in order to gain robustness, so the range of
% taus with a fixed h that makes the system stable becomes larger. Our
% fixed h will be h=0.2

syms tau
max_eig_h02=[];

h=0.2;
stability_h02=max(abs(eig(A_tau(h,0))));
True=stability_h02<1;

for i= 1:150
    tau=0.0025*(i-1);    
    max_eig_h02(i)=max(abs(eig(A_tau(h,tau))));
end
%%
xaxis_r= linspace(0,0.0025*150,150);
plot(xaxis_r,max_eig_h02,'b', 'linewidth', 1.2)
yline(1,'r', 'linewidth', 1.2)
title('Spectral Radius of $\bar{A}$','Interpreter','latex')
legend('max($|\lambda_i|$)', '$\rho(\bar{A})$=1','Interpreter','latex')
ylabel('max($|\lambda_i|$)','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')
%%
% With h=0.2 the maxixmum possible tau is 0.0612, now we will explore the
% tuning of the controller K so we can extend the admissible tau to higher
% values and make the controller more robust K = [6.222222 7 0]
syms tau
syms h
Q=eye(3);
R=eye(1);
K_lqr = dlqr(F_tau_f(0.2,0),G_tau_f(0.2,0),Q,R);
h=0.2;
K_lqr(3)=0.35;
A_tau_rob=F_tau-G_tau*K_lqr;
A_tau_rob= matlabFunction(A_tau_rob);

for i= 1:150
    tau=0.0025*(i-1)+0.0001;    
    max_eig_h02_rob(i)=max(abs(eig(A_tau_rob(h,tau))));
end

hold on
xaxis_r= linspace(0,0.0025*150,150);
plot(xaxis_r,max_eig_h02,'g', 'linewidth', 1.2)
plot(xaxis_r,max_eig_h02_rob,'b', 'linewidth', 1.2)
yline(1,'r', 'linewidth', 1.2)
title('Spectral Radius of $\bar{A}$','Interpreter','latex')
legend('max($|\lambda_i|$) Static Controller','max($|\lambda_i|$) Dynamic Controller','$\rho(\bar{A})$=1','Interpreter','latex')
ylabel('max($|\lambda_i|$)','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')

%% QUESTION 3

% The delay now is higher than h but bounded for tau=[h,2h)
% The system now will have to store an extra state and we have to rebuild
% the matrices

syms h
syms tau

Fx= expm(A_h);
Fk1= int(expm(A_h),[0 2*h-tau])*B;
Fk2= int(expm(A_h),[2*h-tau h])*B;

F_2h= [Fx Fk1 Fk2;0 0 0 0 ; 0 0 eye(1) 0];
G_2h= [zeros(2,1);eye(1);0];
A_2h=F_2h-G_2h*[K 0 0];
F_2h_f= matlabFunction(F_2h);
A_2h= matlabFunction(A_2h);

for i= 1:500
    for j=1:500
        h=0.001*(i-1)+0.0001;
        tau=0.001*(j-1)+0.0001;
        max_eig_tau_2h(i,j)=max(abs(eig(A_2h(h,tau))));
    end
end

%% Plotting Q3.3
max_eig_tau_2h=max_eig_tau_2h';
xaxis= linspace(0,0.001*500,500);
yaxis= linspace(0,0.001*500,500);
copy_stable_2h=zeros(500,500);
copy_unstable_2h=zeros(500,500);

for i=1:500
    for j=1:500
        if 2*j>i && i>= j
            if max_eig_tau_2h(i,j)<1
                copy_stable_2h(j,i)=max_eig_tau_2h(i,j);
            else
                copy_unstable_2h(j,i)=max_eig_tau_2h(i,j);
            end
        end
    end
end


copy_stable_2h(copy_stable_2h==0) = NaN;
surf(xaxis,yaxis,copy_stable_2h')
colormap spring
colorbar
ylabel('$\tau$','Interpreter','latex')
xlabel('h','Interpreter','latex')
title('Stable Combinations of $\{h,\tau\}$','Interpreter','latex')
legend('max($|\lambda_i|$)','Interpreter','latex')
hold on
surf(xaxis,yaxis,copy_stable')
colormap spring
colorbar
ylabel('$\tau$','Interpreter','latex')
xlabel('h','Interpreter','latex')
title('Stable Combinations of $\{h,\tau\}$','Interpreter','latex')
legend('max($|\lambda_i|$)','Interpreter','latex')

%% Question 3.3.

% Tuning K to increase robustness for unstable h=tau=0.1

h=0.1;
stability_005=max(abs(eig(A_2h(h,tau))));
True=stability_005<1;

%% Initial stable range

% The system is AS for h=0.1 and tau between 0.17 and 0.2 (MIRAR ESTO,
% ESTA MAL)
max_eig_005=[];

for i= 1:400
    tau=0.001*(i-1)+0.1;    
    max_eig_005(i)=max(abs(eig(A_2h(h,tau))));
end

xaxis= linspace(0.1,0.001*400+0.1,400);
plot(xaxis,max_eig_005,'b', 'linewidth', 1.2)
yline(1,'r', 'linewidth', 1.2)
xlim([0.1 0.2])
title('Spectral Radius of $\bar{A}$','Interpreter','latex')
legend('max($|\lambda_i|$)', '$\rho(\bar{A})$=1','Interpreter','latex')
ylabel('max($|\lambda_i|$)','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')

%%
% With h=0.1 the range of possible tau is {0.17,0.2}, now we will explore the
% tuning of the controller K so we can extend the admissible tau to higher
% values and make the controller more robust
syms tau
syms h
Q=eye(4);
R=eye(1);
h=0.1;
K_lqr = dlqr(F_2h_f(h,0.1),G_2h,Q,R);
K_lqr(3)=0.6;
K_lqr(4)=0.45;
A_2h_rob=F_2h-G_2h*[K 0.19 0.19];
A_2h_rob= matlabFunction(A_2h_rob);

for i= 1:400
    tau=0.001*(i-1)+0.1;    
    max_eig_01_rob(i)=max(abs(eig(A_2h_rob(h,tau))));
end

xaxis= linspace(0.1,0.001*400+0.1,400);
plot(xaxis,max_eig_01_rob,'b', 'linewidth', 1.2)
hold on
plot(xaxis,max_eig_005,'g', 'linewidth', 1.2)
yline(1,'r', 'linewidth', 1.2)
xlim([0.1 0.2])
title('Spectral Radius of $\bar{A}$','Interpreter','latex')
legend('max($|\lambda_i|$) Static Controller','max($|\lambda_i|$) Dynamic Controller','$\rho(\bar{A})$=1','Interpreter','latex')
ylabel('max($|\lambda_i|$)','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')

%% QUESTION 4

% Definition of four system matrices depending on h: 
syms tau
syms h 
F_2h= [Fx Fk1 Fk2;0 0 0 0 ; 0 0 eye(1) 0];
G_2h= [zeros(2,1);eye(1);0];

% tau = 0.2h
Fx1= expm(A_h);
Fk11= int(expm(A_h),[h-0.2*h h])*B;
Gk21= int(expm(A_h),[0 h-0.2*h])*B;

F1= [Fx1 Fk11 zeros(2,1) ; 0 0 0 0 ; 0 0 eye(1) 0];
G1= [Gk21;eye(1);0];
F1= matlabFunction(F1);
G1= matlabFunction(G1);

% tau = 0.5h
Fx2= expm(A_h);
Fk12= int(expm(A_h),[h-0.5*h h])*B;
Gk22= int(expm(A_h),[0 h-0.5*h])*B;

F2= [Fx2 Fk12 zeros(2,1) ; 0 0 0 0 ; 0 0 eye(1) 0];
G2= [Gk22;eye(1);0];
F2= matlabFunction(F2);
G2= matlabFunction(G2);

%tau=h

Fx3= expm(A_h);
Fk13= int(expm(A_h),[h-h h])*B;
Gk23= int(expm(A_h),[0 h-h])*B;

F3= [Fx3 Fk13 zeros(2,1);0 0 0 0 ; 0 0 eye(1) 0];
G3= [zeros(2,1);eye(1);0];
F3= matlabFunction(F3);

%tau=1.5h

Fx4= expm(A_h);
Fk14= int(expm(A_h),[2*h-1.5*h h])*B;
Fk24= int(expm(A_h),[0 2*h-1.5*h])*B;

F4= [Fx4 Fk14 Fk24;0 0 0 0 ; 0 0 eye(1) 0];
G4= [zeros(2,1);eye(1);0];
F4= matlabFunction(F4);

%% Defining LMIs 4.1 and Solving LMI for K 4.2
valid_h=0;

for i= 1:200    
    h=0.01+0.0125*(i-1);
    F1_h=F1(h);
    G1_h=G1(h);
    F2_h=F2(h);
    G2_h=G2(h);
    F3_h=F3(h);
    G3_h=G3;
    F4_h=F4(h);
    G4_h=G4;

    cvx_begin sdp
    cvx_expert true
        variable M(4,4) symmetric semidefinite 
        variable Y(1,4)       
        [-M, M*F1_h'-Y'*G1_h'; F1_h*M-G1_h*Y, -M]<=-eye(8);
        [-M, M*F2_h'-Y'*G2_h'; F2_h*M-G2_h*Y, -M]<=-eye(8);
        [-M, M*F3_h'-Y'*G3_h'; F3_h*M-G3_h*Y, -M]<=-eye(8);
        [-M, M*F4_h'-Y'*G4_h'; F4_h*M-G4_h*Y, -M]<=-eye(8);        
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        valid_h=h;
        K_cvx=Y*inv(M);
        P=inv(M);
    end
    if strcmp(cvx_status,'Infeasible')
        break
    end
end

%% Define automata for tau=(0.2h,h,0.5h)^w 

% We will find the closed loop matrices A_i and multiply them, check
% Lyapunov decrease for the resulting matrix

valid_h_w=0;

for i= 1:200    
    h=0.01+0.0125*(i-1);
    F1_h=F1(h);
    G1_h=G1(h);
    F2_h=F2(h);
    G2_h=G2(h);
    F3_h=F3(h);
    G3_h=G3;    
    
    cvx_begin sdp
    cvx_expert true
        variable P_i(4,4) symmetric semidefinite 
                
        Ai1=F1_h-G1_h*K_cvx;
        Ai2=F2_h-G2_h*K_cvx;
        Ai3=F3_h-G3_h*K_cvx;
        
        Ai=Ai1*Ai2*Ai3;
        Ai'*P_i*Ai-P_i<=-eye(4);                
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        valid_h_w=h;        
        P_i_b=P_i;
    end  
    if strcmp(cvx_status,'Failed')
        break
    end
    if strcmp(cvx_status,'Infeasible')
        break
    end
end

