%% ASSIGNMENT 2 : NCDS

% A= [ 5 8.5; 0 -3] B= [ 0; 1]
%
%% Question 3

% Define the static controller derived in ASS1 and define the system
% matrices
A= [ 5 8.5; 0 -3];
B= [ 0; 1];
K= place(A,B,[-2 -3]);

syms h

A_h=[5*h 8.5*h; 0 -3*h];
F_h= expm(A_h);
G_h=int(expm(A_h),[0 h])*B;

A_cl= F_h-G_h*K;
A_cl= matlabFunction(A_cl);
F_h= matlabFunction(F_h);

%% Question 3.1 

%using the to zero mechanism

h_zero=linspace(0.001,0.4,40);
deltas=ones(1,40)*-1;

for i= 1:40    
    A1=A_cl(h_zero(i));
    A0=F_h(h_zero(i));
    cvx_begin sdp    
        variable P(2,2) symmetric semidefinite 
        subject to 
        P>=eye(2)*0.0000001;
        A1'*P*A1-P<=-eye(2);        
    cvx_end    
    if strcmp(cvx_status,'Solved')
        deltas(i)=0;
        cvx_begin sdp    
            variable P(2,2) symmetric semidefinite 
            subject to 
            P>=eye(2)*0.0000001;
            A1'*P*A1-P<=-0.001*eye(2);
            (A1*A0)'*P*(A1*A0)-P<=-0.001*eye(2)
        cvx_end
        if strcmp(cvx_status,'Solved')
            deltas(i)=1;
            cvx_begin sdp    
                variable P(2,2) symmetric semidefinite 
                subject to 
                P>=eye(2)*0.0000001;
                A1'*P*A1-P<=-0.001*eye(2);
                (A1*A0)'*P*(A1*A0)-P<=0.001*eye(2)
                (A1*A0*A0)'*P*(A1*A0*A0)-P<=-0.001*eye(2)
            cvx_end
            if strcmp(cvx_status,'Solved')
                deltas(i)=2;
                cvx_begin sdp    
                    variable P(2,2) symmetric semidefinite 
                    subject to 
                    P>=eye(2)*0.0000001;
                    A1'*P*A1-P<=-0.001*eye(2);
                    (A1*A0)'*P*(A1*A0)-P<=-0.001*eye(2)
                    (A1*A0*A0)'*P*(A1*A0*A0)-P<=-0.001*eye(2)
                    (A1*A0*A0*A0)'*P*(A1*A0*A0*A0)-P<=-0.001*eye(2)
                cvx_end
                if strcmp(cvx_status,'Solved')
                    deltas(i)=3;
                    cvx_begin sdp    
                        variable P(2,2) symmetric semidefinite 
                        subject to 
                        P>=eye(2)*0.0000001;
                        A1'*P*A1-P<=-0.001*eye(2);
                        (A1*A0)'*P*(A1*A0)-P<=-0.001*eye(2)
                        (A1*A0*A0)'*P*(A1*A0*A0)-P<=-0.001*eye(2)
                        (A1*A0*A0*A0)'*P*(A1*A0*A0*A0)-P<=-0.001*eye(2)
                        (A1*A0*A0*A0*A0)'*P*(A1*A0*A0*A0*A0)-P<=-0.001*eye(2)
                    cvx_end
                end
            end
        end
    end  
end

deltas(deltas==-1) = NaN;
%%
scatter(h_zero,deltas,70,'x','b')
grid on 
ylabel('$\bar{\delta}$','Interpreter','latex')
xlabel('h','Interpreter','latex')
title('Maximum number of admisible $\bar{\delta}$ for stable system depending on $h$','Interpreter','latex')
legend('$\bar{\delta}_{max}$','Interpreter','latex')
ylim([0 2.2])

%% Question 3.2

%MSS guaranteed for Bernoulli if: P-(1-p)*F_cl^t*P*Fcl0-p*Fcl1^t*P*Fcl1 > 0
%ojo porque esto es derivado de una simplificacion despues de ver que las
%possibilidades son iguales Lecture 3 slide 36

h=0.2;

F_cl0=A_cl(h);
F_cl1=F_h(h);
p_mss=linspace(0.01,0.5,50);

for i= 1:50    
    p=p_mss(i);
    cvx_begin sdp    
        variable P(2,2) symmetric semidefinite 
        subject to 
        P>=eye(2)*0.0000001;
        P-(1-p)*F_cl0'*P*F_cl0-p*F_cl1'*P*F_cl1>=eye(2);
        cvx_quiet true
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        p_bmss=p;
        P_bmss=P;        
    end  
end

%% Question 3.3.

time=linspace(0.0001,3,15);
state=0;
p_bmss=0.12;

h=0.2;

F_cl0=A_cl(h);
F_cl1=F_h(h);

x0= randn(2,1);

for i=1:15
    a=rand;
    if state==1
        if a<p_bmss
            xk1=F_cl1*x0;
            states(i,1)=xk1(1);
            states(i,2)=xk1(2);
            state=0;
        elseif a>=p_bmss
            xk1=F_cl0*x0;
            states(i,1)=xk1(1);
            states(i,2)=xk1(2);
            state=1;
        end
        
    elseif state==0
        if a<p_bmss
            xk1=F_cl1*x0;
            states(i,1)=xk1(1);
            states(i,2)=xk1(2);
            state=0;
        elseif a>=p_bmss
            xk1=F_cl0*x0;
            states(i,1)=xk1(1);
            states(i,2)=xk1(2);
            state=1;
        end
    end 
    x0=xk1;    
end

plot(time,states,'k-','LineWidth',1.2)
%hold on
%plot(time,states(:,2),'k-','LineWidth',1.2)
xlabel('Time (s)','Interpreter','latex');
ylabel('$\xi$(s)','Interpreter','latex');
legend('$\xi(s_0)$ ','$\xi(s_1)$','Interpreter','latex')
title('Evolution of the states with $p<p^*$','Interpreter','latex')
hold off 


%% Question 3.4 

% This is going to be funny xd we have 3 equations to satisfy, however the
% last is parametric lambda^(1-p)*L^p < 1. So we will find the largest
% possible L in advance for each p and lambda and then solve for the two
% other LMI V(Fcl0x)<= lambda* V(x) and V(Fcl1x) <= L* V(x)

h=0.2;

F_cl0=A_cl(h);
F_cl1=F_h(h);

p_ass=linspace(0.01,0.5,20);
lambda=linspace(0.001,1-0.001,10);

for i= 1:length(p_ass)
    for j=1:length(lambda)
        
        p=p_ass(i);
        l=lambda(j);

        L=((1-0.0001)/(l^(1-p)))^(1/p);

        cvx_begin sdp        
            variable P(2,2) symmetric semidefinite 
            subject to 
            P>=eye(2)*0.0001;
            F_cl0'*P*F_cl0<=l*P;
            F_cl1'*P*F_cl1<=L*P;
            cvx_quiet true
        cvx_end

        if strcmp(cvx_status,'Solved')
            p_bass=p;
            break
        end 
    end
end

%% Question 3.5 Gilbert-Elliot Model

% For assessing MSS with a gilbert-elliot model we have to solve two LMIs
% P0-p*F_cl0^t*P0*Fcl0-(1-p)*Fcl1^t*P1*Fcl1 > 0 and
% P1-q*F_cl1^t*P1*Fcl1-(1-q)*Fcl0^t*P0*Fcl0 > 0

p_gmss=linspace(0.01,1,50);
q_gmss=linspace(0.01,1,50);
valid_pq=zeros(length(p_gmss),length(q_gmss));

for i= 1:length(p_gmss)
    for j=1:length(q_gmss)
        
        cvx_begin sdp        
            variable P0(2,2) symmetric semidefinite
            variable P1(2,2) symmetric semidefinite            
            subject to
                P0>=eye(2)*0.0000000001;
                P1>=eye(2)*0.0000000001;
                P0-p_gmss(i)*F_cl0'*P0*F_cl0-(1-p_gmss(i))*F_cl1'*P1*F_cl1 >= 0.0001*eye(2);
                P1-q_gmss(j)*F_cl1'*P1*F_cl1-(1-q_gmss(j))*F_cl0'*P0*F_cl0 >= 0.0001*eye(2);
                cvx_quiet true
            cvx_end

        if strcmp(cvx_status,'Solved')
            valid_pq(i,j)=1;            
        end
        
        if strcmp(cvx_status,'Infeasible')
            break;            
        end 
    end
end
%% Plot this tomorrow
valid_pq(valid_pq==0) = NaN;
surf(p_gmss,q_gmss,valid_pq')
ylabel('$q$','Interpreter','latex')
xlabel('$p$','Interpreter','latex')
title('Stable Combinations of $p$ and $q$ for fixed $h=0.2$','Interpreter','latex')
ylim([0 1])
xlim([0 1])
%% Question 4

% Jordan approach, decompose to find F0, F1 , F2 and G0 , G1 and G2 also
% alpha1 and alpha2, now we have tau = [0,h) so we have to change the
% matrices

syms tau
syms h 

Fx= expm(A_h);
Fu= int(expm(A_h),[h-tau h])*B;
G1= int(expm(A_h),[0 h-tau])*B;

F_tau= [Fx Fu; zeros(1,3)];
G_tau= [G1;eye(1)];

F_tau_f= matlabFunction(F_tau);
G_tau_f= matlabFunction(G_tau);

%%

% Based on the observed F_tau and G_tau:

alpha1=exp(-3*h)*exp(3*tau);
alpha2=exp(5*h)*exp(-5*tau);
alpha1= matlabFunction(alpha1);
alpha2= matlabFunction(alpha2);

F0=[exp(5*h) (17*exp(5*h))/16-(17*exp(-3*h))/16 (17*exp(5*h))/80+17*exp(-3*h)/48;0 exp(-3*h) -exp(-3*h)/3; 0 0 0];
F1=[0 0 -17/48;0 0 1/3;0 0 0];
F2=[0 0 -17/80;0 0 0;0 0 0];
F0= matlabFunction(F0);

G0=[-17/30;1/3;1];
G1=[17/48; -1/3;0];
G2=[17/80;0;0];

% We can rewrite the system as xek+1= F0+a1*F1+a2*F2+G0+a1*G1+a2*G2
%%

% here we find the range of values for tau=h so maximum range
tau_min=0.001*h;
tau_max=(1-0.001)*h;
tau_max= matlabFunction(tau_max);
tau_min= matlabFunction(tau_min);
h_test=linspace(0.001,0.5,200);
valid_hmax=zeros(1,length(h_test));

%% Vamo a jugalllll LEGO 4 vertices

taus=linspace(0.0001,0.15,60);
h_test=linspace(0.001,0.5,60);
combis=zeros(60,60);

for j=1:length(taus)
    for i=1:length(h_test)
        if taus(i)<h_test(j)
            tau_mx=taus(i);
            tau_mn=taus(1);
            alpha1_min=alpha1(h_test(j),tau_mn);
            alpha1_max=alpha1(h_test(j),tau_mx);
            alpha2_min=alpha2(h_test(j),tau_mn);
            alpha2_max=alpha2(h_test(j),tau_mx);

            Ha1min2max=F0(h_test(j))+F1*alpha1_min+F2*alpha2_max;
            Ha1max2min=F0(h_test(j))+F1*alpha1_max+F2*alpha2_min;
            Ha1max2max=F0(h_test(j))+F1*alpha1_max+F2*alpha2_max;
            Ha1min2min=F0(h_test(j))+F1*alpha1_min+F2*alpha2_min;

            Ga1min2max=G0+G1*alpha1_min+G2*alpha2_max;
            Ga1max2min=G0+G1*alpha1_max+G2*alpha2_min;
            Ga1max2max=G0+G1*alpha1_max+G2*alpha2_max;
            Ga1min2min=G0+G1*alpha1_min+G2*alpha2_min;
            cvx_begin sdp        
                variable P(3,3) symmetric semidefinite            
                P>=eye(3)*0.00001;
                (Ha1min2max-Ga1min2max*[K 0])'*P*(Ha1min2max-Ga1min2max*[K 0])-P <= -eye(3);
                (Ha1max2min-Ga1max2min*[K 0])'*P*(Ha1max2min-Ga1max2min*[K 0])-P <= -eye(3);
                (Ha1max2max-Ga1max2max*[K 0])'*P*(Ha1max2max-Ga1max2max*[K 0])-P <= -eye(3);
                (Ha1min2min-Ga1min2min*[K 0])'*P*(Ha1min2min-Ga1min2min*[K 0])-P <= -eye(3);
                cvx_quiet true
            cvx_end

            if strcmp(cvx_status,'Solved')
                combis(i,j)=1;            
            end
        end
    end
end

%%
combis(combis==0) = NaN;
surf(h_test,taus,combis)
ylabel('$\tau$','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
title('Stable Combinations of $h$ and $\tau$','Interpreter','latex')
ylim([0 0.15])
xlim([0 0.5])

%%
for j=(2/3)*length(taus):length(taus)
    for i=3:7
        for p=1:3
            if taus(i)<h_test(j)
                tau_mx=taus(i);
                tau_mn=taus(3-(p-1));
                alpha1_min=alpha1(h_test(j),tau_mn);
                alpha1_max=alpha1(h_test(j),tau_mx);
                alpha2_min=alpha2(h_test(j),tau_mn);
                alpha2_max=alpha2(h_test(j),tau_mx);

                Ha1min2max=F0(h_test(j))+F1*alpha1_min+F2*alpha2_max;
                Ha1max2min=F0(h_test(j))+F1*alpha1_max+F2*alpha2_min;
                Ha1max2max=F0(h_test(j))+F1*alpha1_max+F2*alpha2_max;
                Ha1min2min=F0(h_test(j))+F1*alpha1_min+F2*alpha2_min;

                Ga1min2max=G0+G1*alpha1_min+G2*alpha2_max;
                Ga1max2min=G0+G1*alpha1_max+G2*alpha2_min;
                Ga1max2max=G0+G1*alpha1_max+G2*alpha2_max;
                Ga1min2min=G0+G1*alpha1_min+G2*alpha2_min;
                cvx_begin sdp        
                    variable P(3,3) symmetric semidefinite            
                    P>=eye(3)*0.00001;
                    (Ha1min2max-Ga1min2max*[K 0])'*P*(Ha1min2max-Ga1min2max*[K 0])-P <= -eye(3);
                    (Ha1max2min-Ga1max2min*[K 0])'*P*(Ha1max2min-Ga1max2min*[K 0])-P <= -eye(3);
                    (Ha1max2max-Ga1max2max*[K 0])'*P*(Ha1max2max-Ga1max2max*[K 0])-P <= -eye(3);
                    (Ha1min2min-Ga1min2min*[K 0])'*P*(Ha1min2min-Ga1min2min*[K 0])-P <= -eye(3);
                    cvx_quiet true
                cvx_end

                if strcmp(cvx_status,'Solved')
                    combis(i,j)=1;            
                end
            end
        end
    end
end

%% 8 vertices, this is going to be rough

taus=linspace(0.0001,0.15,60);
h_test=linspace(0.001,0.5,60);
combis2=zeros(60,60);

for j=1:length(taus)
    for i=1:length(h_test)
        if taus(i)<h_test(j)
            tau_mx=taus(i);
            tau_mn=taus(1);
            alpha1_min=alpha1(h_test(j),tau_mn);
            alpha1_max=alpha1(h_test(j),tau_mx);
            alpha12=alpha1(h_test(j),tau_mx/2);
            alpha2_min=alpha2(h_test(j),tau_mn);
            alpha2_max=alpha2(h_test(j),tau_mx);
            alpha22=alpha2(h_test(j),tau_mx/2);

            Ha1min2max=F0(h_test(j))+F1*alpha1_min+F2*alpha2_max;
            Ha1max2min=F0(h_test(j))+F1*alpha1_max+F2*alpha2_min;
            Ha1max2max=F0(h_test(j))+F1*alpha1_max+F2*alpha2_max;
            Ha1min2min=F0(h_test(j))+F1*alpha1_min+F2*alpha2_min;
            
            Ha1min22=F0(h_test(j))+F1*alpha1_min+F2*alpha22;
            Ha1max22=F0(h_test(j))+F1*alpha1_max+F2*alpha22;
            H12max=F0(h_test(j))+F1*alpha12+F2*alpha2_max;
            H12min=F0(h_test(j))+F1*alpha12+F2*alpha2_min;
            
            Ga1min2max=G0+G1*alpha1_min+G2*alpha2_max;
            Ga1max2min=G0+G1*alpha1_max+G2*alpha2_min;
            Ga1max2max=G0+G1*alpha1_max+G2*alpha2_max;
            Ga1min2min=G0+G1*alpha1_min+G2*alpha2_min;
            
            Ga1min22=G0+G1*alpha1_min+G2*alpha22;
            Ga1max22=G0+G1*alpha1_max+G2*alpha22;
            Ga12max=G0+G1*alpha12+G2*alpha2_max;
            Ga22min=G0+G1*alpha12+G2*alpha2_min;
            
            cvx_begin sdp        
                variable P(3,3) symmetric semidefinite            
                P>=eye(3)*0.00001;
                (Ha1min2max-Ga1min2max*[K 0])'*P*(Ha1min2max-Ga1min2max*[K 0])-P <= -eye(3);
                (Ha1max2min-Ga1max2min*[K 0])'*P*(Ha1max2min-Ga1max2min*[K 0])-P <= -eye(3);
                (Ha1max2max-Ga1max2max*[K 0])'*P*(Ha1max2max-Ga1max2max*[K 0])-P <= -eye(3);
                (Ha1min2min-Ga1min2min*[K 0])'*P*(Ha1min2min-Ga1min2min*[K 0])-P <= -eye(3);
                
                (Ha1min22-Ga1min22*[K 0])'*P*(Ha1min22-Ga1min22*[K 0])-P <= -eye(3);
                (Ha1max22-Ga1max22*[K 0])'*P*(Ha1max22-Ga1max22*[K 0])-P <= -eye(3);
                (H12max-Ga12max*[K 0])'*P*(H12max-Ga12max*[K 0])-P <= -eye(3);
                (H12min-Ga22min*[K 0])'*P*(H12min-Ga22min*[K 0])-P <= -eye(3);
                
                cvx_quiet true
            cvx_end

            if strcmp(cvx_status,'Solved')
                combis2(i,j)=1;            
            end
        end
    end
end
%%
combis2(combis2==0) = NaN;
surf(h_test,taus,combis2)
ylabel('$\tau$','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
title('Stable Combinations of $h$ and $\tau$','Interpreter','latex')
ylim([0 0.15])
xlim([0 0.5])
%%
for j=40:length(taus)
    for i=3:7
        for p=1:4
            if taus(i)<h_test(j)
                tau_mx=taus(i);
                tau_mn=taus(4-(p-1));
                alpha1_min=alpha1(h_test(j),tau_mn);
                alpha1_max=alpha1(h_test(j),tau_mx);
                alpha12=alpha1(h_test(j),tau_mx/2);
                alpha2_min=alpha2(h_test(j),tau_mn);
                alpha2_max=alpha2(h_test(j),tau_mx);
                alpha22=alpha2(h_test(j),tau_mx/2);

                Ha1min2max=F0(h_test(j))+F1*alpha1_min+F2*alpha2_max;
                Ha1max2min=F0(h_test(j))+F1*alpha1_max+F2*alpha2_min;
                Ha1max2max=F0(h_test(j))+F1*alpha1_max+F2*alpha2_max;
                Ha1min2min=F0(h_test(j))+F1*alpha1_min+F2*alpha2_min;

                Ha1min22=F0(h_test(j))+F1*alpha1_min+F2*alpha22;
                Ha1max22=F0(h_test(j))+F1*alpha1_max+F2*alpha22;
                H12max=F0(h_test(j))+F1*alpha12+F2*alpha2_max;
                H12min=F0(h_test(j))+F1*alpha12+F2*alpha2_min;

                Ga1min2max=G0+G1*alpha1_min+G2*alpha2_max;
                Ga1max2min=G0+G1*alpha1_max+G2*alpha2_min;
                Ga1max2max=G0+G1*alpha1_max+G2*alpha2_max;
                Ga1min2min=G0+G1*alpha1_min+G2*alpha2_min;

                Ga1min22=G0+G1*alpha1_min+G2*alpha22;
                Ga1max22=G0+G1*alpha1_max+G2*alpha22;
                Ga12max=G0+G1*alpha12+G2*alpha2_max;
                Ga22min=G0+G1*alpha12+G2*alpha2_min;

                cvx_begin sdp        
                    variable P(3,3) symmetric semidefinite            
                    P>=eye(3)*0.00001;
                    (Ha1min2max-Ga1min2max*[K 0])'*P*(Ha1min2max-Ga1min2max*[K 0])-P <= -eye(3);
                    (Ha1max2min-Ga1max2min*[K 0])'*P*(Ha1max2min-Ga1max2min*[K 0])-P <= -eye(3);
                    (Ha1max2max-Ga1max2max*[K 0])'*P*(Ha1max2max-Ga1max2max*[K 0])-P <= -eye(3);
                    (Ha1min2min-Ga1min2min*[K 0])'*P*(Ha1min2min-Ga1min2min*[K 0])-P <= -eye(3);

                    (Ha1min22-Ga1min22*[K 0])'*P*(Ha1min22-Ga1min22*[K 0])-P <= -eye(3);
                    (Ha1max22-Ga1max22*[K 0])'*P*(Ha1max22-Ga1max22*[K 0])-P <= -eye(3);
                    (H12max-Ga12max*[K 0])'*P*(H12max-Ga12max*[K 0])-P <= -eye(3);
                    (H12min-Ga22min*[K 0])'*P*(H12min-Ga22min*[K 0])-P <= -eye(3);

                    cvx_quiet true
                cvx_end

                if strcmp(cvx_status,'Solved')
                    combis2(i,j)=1;      
                end
            end
        end
    end
end

%% Question 5.1 & 5.2

% Define the static controller derived in ASS1 and define the system
% matrices
A= [ 5 8.5; 0 -3];
B= [ 0; 1];
K= place(A,B,[-2 -3]);

cvx_begin sdp        
    variable P(2,2) symmetric semidefinite
    variable Q(2,2) symmetric semidefinite
    subject to 
    P>=eye(2)*0.0001;
    Q>=eye(2)*0.0001;
    (A-B*K)'*P+P*(A-B*K)==-Q
    cvx_quiet true
cvx_end

events_amount = [];    

for n=1:50
	x_i = rand(2,1)*3;
	sigma=0.999;
    
    x_sk = x_i;
    t_start = 0;
    t_current = t_start;
    t_final = 3;
    phi = [A'*P+P*A+sigma*Q  -P*B*K; -(B*K)'*P  zeros(2,2)];
    tout = t_start;
    xout = x_i.';
    teout = [];
    xeout = [];
    ieout = [];

    for i = 1:500
        
        options = odeset('Events',@(t,x)odetrig(t,x,phi,x_sk));
        [t,x,te,xe,ie] = ode45(@(t,x) odefun(t,x,A,B,K,x_sk),[t_current t_final], x_sk, options);

        nt = length(t);
        tout = [tout; t(2:nt)];
        xout = [xout; x(2:nt,:)];
        teout = [teout; te];        
        xeout = [xeout; xe];
        ieout = [ieout; ie];

        x_sk = x(end,:)';
        t_current = t(nt);

        if t(end)==t_final
            break
        end
        
    end
 
	events_amount(n) = length(ieout); 
end  

plot([0;teout],[x_i(1) ;xeout(:,1)],'r*')
hold on
plot([0;teout],[x_i(2) ;xeout(:,2)],'k*')
hold on
plot([0;tout],[x_i(1) ;xout(:,1)],'r')
hold on
plot([0;tout],[x_i(2) ;xout(:,2)],'k')
xlabel('Time (s)','Interpreter','latex');
ylabel('$\xi$(s)','Interpreter','latex');
legend('$\xi(s_0)$','$\xi(s_1)$','Interpreter','latex')
title('Event-Triggering control with $\sigma=0.99$ ','Interpreter','latex')
hold off    



mean_events = mean(events_amount,2);

%% Question 5.3

events0_1 = []
h_inter=[]

for n=1:50
	x_i = randn(2,1);
	sigma=0.5;
    
    x_sk = x_i;
    t_start = 0;0.2
    t_current = t_start;
    t_final = 3;
    phi = [A'*P+P*A+sigma*Q  -P*B*K; -(B*K)'*P  zeros(2,2)];
    tout = t_start;
    xout = x_i.';
    teout = [];
    xeout = [];
    ieout = [];

    for i = 1:500
        
        options = odeset('Events',@(t,x)odetrig(t,x,phi,x_sk));
        [t,x,te,xe,ie] = ode45(@(t,x) odefun(t,x,A,B,K,x_sk),[t_current t_final], x_sk, options);

        nt = length(t);
        tout = [tout; t(2:nt)];
        xout = [xout; x(2:nt,:)];
        teout = [teout; te];        
        xeout = [xeout; xe];
        ieout = [ieout; ie];

        x_sk = x(end,:)';
        t_current = t(nt);

        if t(end)==t_final
            break
        end
        
        if x_sk'*P*x_sk<=0.1*x_i'*P*x_i
            break
        end      
    
        
    end
 
	events0_1(n) = length(ieout);
    h_inter(n)=t_current/events0_1(n);
end

plot([0;teout],[x_i(1) ;xeout(:,1)],'r*')
hold on
plot([0;teout],[x_i(2) ;xeout(:,2)],'k*')
hold on
plot([0;tout],[x_i(1) ;xout(:,1)],'r')
hold on
plot([0;tout],[x_i(2) ;xout(:,2)],'k')
xlabel('Time (s)','Interpreter','latex');
ylabel('$\xi$(s)','Interpreter','latex');
legend('$\xi(s_0)$','$\xi(s_1)$','Interpreter','latex')
hold off    

h_avg = mean(h_inter);


%% Question 5.4

sys = ss(A-B*K,[0;0],[0 0],0);
sysd=c2d(sys,h_avg);

x0=[1;-1];
time=linspace(0,3,3/h_avg)';
states=[];

for i=1:10   
    xk1=sysd.A*x0;
    states(i,1)=xk1(1);
    states(i,2)=xk1(2);
    x0=xk1;    
end


x_i = [1;-1];
sigma=0.5;

x_sk = x_i;
t_start = 0;
t_current = t_start;
t_final = 3;
phi = [A'*P+P*A+sigma*Q  -P*B*K; -(B*K)'*P  zeros(2,2)];
tout = t_start;
xout = x_i.';
teout = [];
xeout = [];
ieout = [];

for i = 1:500

    options = odeset('Events',@(t,x)odetrig(t,x,phi,x_sk));
    [t,x,te,xe,ie] = ode45(@(t,x) odefun(t,x,A,B,K,x_sk),[t_current t_final], x_sk, options);

    nt = length(t);
    tout = [tout; t(2:nt)];
    xout = [xout; x(2:nt,:)];
    teout = [teout; te];        
    xeout = [xeout; xe];
    ieout = [ieout; ie];

    x_sk = x(end,:)';
    t_current = t(nt);

    if t(end)==t_final
        break
    end

end

states=[ [1 -1];[states]];

plot([0;tout],[x_i(1) ;xout(:,1)],'r--')
hold on
plot([0;tout],[x_i(2) ;xout(:,2)],'r:')
hold on
plot(time,states(:,1),'k*')
hold on
plot(time,states(:,2),'k+')
xlabel('Time (s)','Interpreter','latex');
ylabel('$\xi$(s)','Interpreter','latex');
legend('$\xi(s_0)$ Event Triggered Control','$\xi(s_1)$ Event Triggered Control','$\xi(s_0)$ Time Triggered Control','$\xi(s_1)$ Time Triggered Control','Interpreter','latex')
hold off    


function dydt = odefun(t,x,A,B,K,x_sk)

dydt = A*x - B*K*x_sk;

end

function [condition,isterminal,direction] = odetrig(t,x,phi,x_sk)

condition = [x' x_sk']*phi*[x; x_sk];
condition = double(condition<=0);
isterminal = 1;
direction = -1;

end

