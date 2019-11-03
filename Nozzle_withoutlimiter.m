%% Project Shock tube 2nd order without limiter
clc
clear all
%% Defining Variables
L = 10;
N = 50;
y = 1.4;
R = 8.314;
dx = (L/(N-1));
%% Variables at left side
rho_l = 1.1308;
u_l = 80.054;
p_l = 9.6328*10^4;
M = 28.97 * 10^-3;

T_L = p_l*M/(rho_l*R);% Temprature required for speed of sound
a_L = sqrt(y*R*T_L);%Speed of sound
dt = 0.01*dx/a_L;% Using Courant number and CFL=0.01, its arbitrary
tf = 1;% Final time

Nt = fix(tf/dt);
%% Variables at right side
rho_r = 1.0275;
u_r = 150.535;
p_r = 8.4974*10^4;

E_l = ( p_l/(y-1) + 0.5*rho_l*(u_l^2) )/rho_l;
E_r = ( p_r/(y-1) + 0.5*rho_r*(u_r^2) )/rho_r;

H_l = E_l + p_l/rho_l;
H_r = E_r + p_r/rho_r;

u = zeros(N,Nt);
rho = zeros(N,Nt);
p = zeros(N,Nt);
E = zeros(N,Nt);
H = zeros(N,Nt);
Q1 = zeros(N,Nt);
Q2 = zeros(N,Nt);
Q3 = zeros(N,Nt);
F1 = zeros(N,Nt);
F2 = zeros(N,Nt);
F3 = zeros(N,Nt);
S=zeros(N,1);
dS=zeros(N,1);
g2=zeros(N,Nt);
%% Boundary conditions for primitive variables
x=0:dx:10;

for i=1:N/2
        rho(i,1) = rho_l;
        u(i,1) = u_l;
        p(i,1) = p_l;
        H(i,1) = E_l + p_l/rho_l;
        E(i,1) = E_l;       
        S(i)=1+1.5*(0.2*x(i)-1)^2;
        dS(i)=3*0.2*(0.2*x(i)-1);
        g2(i,1)=p(i)*dS(i);
end

for i=N/2+1:N
        rho(i,1) = rho_r;
        u(i,1) = u_r;
        p(i,1) = p_r;
        H(i,1) =  E_r + p_r/rho_r;
        E(i,1) = E_r;
        S(i)=1+0.5*(0.2*x(i)-1)^2;
        dS(i)=0.2*(0.2*x(i)-1);
        g2(i,1)=p(i)*dS(i);
end
%% Initial condition for primitive variables
u(1,:) = u_l;
u(N,:) = u_r;
rho(1,:) = rho_l;
rho(N,:) = rho_r;
p(1,:) = p_l;
p(N,:) = p_r;
E(1,:) = E_l;
E(N,:) = E_r;
H(1,:) = H_l;
H(N,:) = H_r;

%% IC and BC for Q
 Q1(:,1) = rho(:,1) .* S;
 Q2(:,1) = rho(:,1).*u(:,1).*S;
 Q3(:,1) = rho(:,1).*E(:,1).*S;
 
 Q1(1,:) = rho(1,:) .* S(1);
 Q2(1,:) = rho(1,:).*u(1,:).*S(1);
 Q3(1,:) = rho(1,:).*E(1,:).*S(1);
 
 Q1(N,:) = rho(N,:) .* S(N);
 Q2(N,:) = rho(N,:).*u(N,:).*S(N);
 Q3(N,:) = rho(N,:).*E(N,:).*S(N);
 
 %% IC and BC for F
 F1(:,1) = rho(:,1).*u(:,1).*S;
 F2(:,1) = ( rho(:,1).*(u(:,1).^2) + p(:,1) ) .* S;
 F3(:,1) = rho(:,1).*u(:,1).*H(:,1).*S;
 
 F1(1,:) = rho(1,:).*u(1,:)*S(1);
 F2(1,:) = ( rho(1,:).*(u(1,:).^2) + p(1,:) ) .* S(1);
 F3(1,:) = rho(1,:).*u(1,:).*H(1,:).*S(1);
 
 F1(N,:) = rho(N,:).*u(N,:).*S(N);
 F2(N,:) = ( rho(N,:).*(u(N,:).^2) + p(N,:) ) .* S(N);
 F3(N,:) = rho(N,:).*u(N,:).*H(N,:).*S(N);

 
 
 for k=1:Nt-1
     for i=2:N-1
         
         
         if i==2
        %% Roe averaging
        rho_plus2= sqrt(rho(i+2,k)*rho(i+1,k));%rho at i+3/2
        rho_plus1 = sqrt(rho(i,k)*rho(i+1,k));%rho at i+1/2
        rho_minus1 = sqrt(rho(i,k)*rho(i-1,k));%rho at i-1/2
        %rho_minus2 = sqrt(r*rho(i-1,k));%rho at i-3/2
        u_plus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i+1,k) )*u(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%Velocity at i+1/2
        u_plus2 = ( sqrt(rho(i+2,k))*u(i+2,k) + sqrt(rho(i+1,k))*u(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%Velocity at i+3/2
        u_minus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%Velocity at i-1/2
        %u_minus2 = ( sqrt(rho(i-2,k))*u(i-2,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%Velocity at i-/2
        H_plus2 = ( sqrt(rho(i+2,k))*H(i+2,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%H at i+3/2   
        H_plus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%H at i+1/2
        H_minus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%H at i-1/2
        %H_minus2 = ( sqrt(rho(i-2,k))*H(i-2,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%H at i-3/2
        a_plus2 = sqrt( (y-1)*(H_plus2 - 0.5*(u_plus2^2)) );
        a_plus1 = sqrt( (y-1)*(H_plus1 - 0.5*(u_plus1^2)) );
        a_minus1 = sqrt( (y-1)*(H_minus1 - 0.5*(u_minus1^2)) );
        %a_minus2 = sqrt( (y-1)*(H_minus2 - 0.5*(u_minus2^2)) );
        %% R inverse
        Rinv_plus = [1 (rho_plus1/(sqrt(2)*a_plus1)) (rho_plus1/(sqrt(2)*a_plus1));u_plus1 ((u_plus1 + a_plus1)*(rho_plus1/(sqrt(2)*a_plus1))) ((u_plus1 - a_plus1)*(rho_plus1/(sqrt(2)*a_plus1)));(0.5*(u_plus1^2)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) + a_plus1*u_plus1)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) - a_plus1*u_plus1))];                                     
        Rinv_minus = [1 (rho_minus1/(sqrt(2)*a_minus1)) (rho_minus1/(sqrt(2)*a_minus1));u_minus1 ((u_minus1 + a_minus1)*(rho_minus1/(sqrt(2)*a_minus1))) ((u_minus1 - a_minus1)*(rho_minus1/(sqrt(2)*a_minus1)));(0.5*(u_minus1^2)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) + a_minus1*u_minus1)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) - a_minus1*u_minus1))];
        %% R at i+3/2
        R_plus2 = [(1 - 0.5*(y-1)*(u_plus2^2)/(a_plus2^2)) ((y-1)*u_plus2/(a_plus2^2)) (1-y)/(a_plus2^2);((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) - a_plus2*u_plus2)) ((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 - u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1));((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) + a_plus2*u_plus2 )) -((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 + u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1))];
        %% R at i+1/2
        R_plus1 = [(1 - 0.5*(y-1)*(u_plus1^2)/(a_plus1^2)) ((y-1)*u_plus1/(a_plus1^2)) (1-y)/(a_plus1^2);((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) - a_plus1*u_plus1)) ((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 - u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1));((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) + a_plus1*u_plus1 )) -((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 + u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1))];
        %% R at i-1/2
        R_minus1 = [(1 - 0.5*(y-1)*(u_minus1^2)/(a_minus1^2)) ((y-1)*u_minus1/(a_minus1^2)) (1-y)/(a_minus1^2);((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) - a_minus1*u_minus1)) ((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 - u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1));((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) + a_minus1*u_minus1 )) -((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 + u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1))];
       
        %% R at i-3/2
        %R_minus2 = [(1 - 0.5*(y-1)*(u_minus2^2)/(a_minus2^2)) ((y-1)*u_minus2/(a_minus2^2)) (1-y)/(a_minus2^2);((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) - a_minus2*u_minus2)) ((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 - u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1));((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) + a_minus2*u_minus2 )) -((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 + u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1))];
        %% eigenmatrix
        L_plus = [abs(u_plus1) 0 0;0 abs(u_plus1 + a_plus1) 0;0 0 abs(u_plus1 - a_minus1)];
        L_minus = [abs(u_minus1) 0 0;0 abs(u_minus1 + a_minus1) 0;0 0 abs(u_minus1 - a_minus1)];
        
        Q_plus2 = [Q1(i+2,k)-Q1(i+1,k);Q2(i+2,k)-Q2(i+1,k);Q3(i+2,k)-Q3(i+1,k)];
        Q_plus1 = [Q1(i+1,k)-Q1(i,k);Q2(i+1,k)-Q2(i,k);Q3(i+1,k)-Q3(i,k)];
        Q_minus1 = [Q1(i,k)-Q1(i-1,k);Q2(i,k)-Q2(i-1,k);Q3(i,k)-Q3(i-1,k)];
        %Q_minus2 = [Q1(i-1,k)-Q1(i-2,k);Q2(i-1,k)-Q2(i-2,k);Q3(i-1,k)-Q3(i-2,k)];        
        
        M_plus2=R_plus2*Q_plus2;%M at i+3/2
        M_plus1=R_plus1*Q_plus1;%M at i+1/2
        M_minus1=R_minus1*Q_minus1;%M at i-1/2
        %M_minus2=R_minus2*Q_minus2;% M at i-3/2
        for v=1:3
            N_plus(v,1)=M_plus1(v)-minmod(M_minus1(v),M_plus1(v),M_plus2(v));%N i+1/2
            N_minus(v,1)=M_minus1(v)-minmod(M_minus1(v),M_plus1(v));%N at i-1/2, ignoring the M at the -1 node
        end
        AQ_plus = Rinv_plus*L_plus*N_plus;
        AQ_minus = Rinv_minus*L_minus*N_minus;
        Q1(i,k+1) = Q1(i,k) - dt * 0.5 * ( F1(i+1,k) - F1(i-1,k) - AQ_plus(1) + AQ_minus(1))/dx;
        Q2(i,k+1) = Q2(i,k) - dt *(-g2(i,k)+ 0.5 * ( F2(i+1,k) - F2(i-1,k) - AQ_plus(2) + AQ_minus(2))/dx);
        Q3(i,k+1) = Q3(i,k) - dt * 0.5 * ( F3(i+1,k) - F3(i-1,k) - AQ_plus(3) + AQ_minus(3))/dx;
        %% Extracting the primitive variables
        rho(i,k+1) =(Q1(i,k+1))./S(i);  
        u(i,k+1) = Q2(i,k+1)./(rho(i,k+1).*S(i));
        E(i,k+1) = Q3(i,k+1)./(rho(i,k+1).*S(i));
        p(i,k+1) = (y-1) * ( rho(i,k+1).*E(i,k+1) - 0.5*rho(i,k+1).*(u(i,k+1).^2) );
        H(i,k+1) = E(i,k+1) + p(i,k+1)./rho(i,k+1);
        g2(i,k+1)=p(i,k+1)*dS(i);
        F1(i,k+1) = rho(i,k+1).*u(i,k+1).*S(i);
        F2(i,k+1) = ( rho(i,k+1).*(u(i,k+1).^2) + p(i,k+1) ) .* S(i);
        F3(i,k+1) = rho(i,k+1).*u(i,k+1).*H(i,k+1).*S(i); 

         
         
         elseif i==N-1
        %% Roe averaging
        %rho_plus2= sqrt(rho(i+2,k)*rho(i+1,k));%rho at i+3/2
        rho_plus1 = sqrt(rho(i,k)*rho(i+1,k));%rho at i+1/2
        rho_minus1 = sqrt(rho(i,k)*rho(i-1,k));%rho at i-1/2
        rho_minus2 = sqrt(rho(i-2,k)*rho(i-1,k));%rho at i-3/2
        u_plus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i+1,k) )*u(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%Velocity at i+1/2
        %u_plus2 = ( sqrt(rho(i+2,k))*u(i+2,k) + sqrt(rho(i+1,k))*u(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%Velocity at i+3/2
        u_minus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%Velocity at i-1/2
        u_minus2 = ( sqrt(rho(i-2,k))*u(i-2,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%Velocity at i-/2
        %H_plus2 = ( sqrt(rho(i+2,k))*H(i+2,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%H at i+3/2   
        H_plus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%H at i+1/2
        H_minus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%H at i-1/2
        H_minus2 = ( sqrt(rho(i-2,k))*H(i-2,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%H at i-3/2
        %a_plus2 = sqrt( (y-1)*(H_plus2 - 0.5*(u_plus2^2)) );
        a_plus1 = sqrt( (y-1)*(H_plus1 - 0.5*(u_plus1^2)) );
        a_minus1 = sqrt( (y-1)*(H_minus1 - 0.5*(u_minus1^2)) );
        a_minus2 = sqrt( (y-1)*(H_minus2 - 0.5*(u_minus2^2)) );
        %% R inverse
        Rinv_plus = [1 (rho_plus1/(sqrt(2)*a_plus1)) (rho_plus1/(sqrt(2)*a_plus1));u_plus1 ((u_plus1 + a_plus1)*(rho_plus1/(sqrt(2)*a_plus1))) ((u_plus1 - a_plus1)*(rho_plus1/(sqrt(2)*a_plus1)));(0.5*(u_plus1^2)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) + a_plus1*u_plus1)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) - a_plus1*u_plus1))];                                     
        Rinv_minus = [1 (rho_minus1/(sqrt(2)*a_minus1)) (rho_minus1/(sqrt(2)*a_minus1));u_minus1 ((u_minus1 + a_minus1)*(rho_minus1/(sqrt(2)*a_minus1))) ((u_minus1 - a_minus1)*(rho_minus1/(sqrt(2)*a_minus1)));(0.5*(u_minus1^2)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) + a_minus1*u_minus1)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) - a_minus1*u_minus1))];
        %% R at i+3/2
        %R_plus2 = [(1 - 0.5*(y-1)*(u_plus2^2)/(a_plus2^2)) ((y-1)*u_plus2/(a_plus2^2)) (1-y)/(a_plus2^2);((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) - a_plus2*u_plus2)) ((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 - u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1));((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) + a_plus2*u_plus2 )) -((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 + u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1))];
        %% R at i+1/2
        R_plus1 = [(1 - 0.5*(y-1)*(u_plus1^2)/(a_plus1^2)) ((y-1)*u_plus1/(a_plus1^2)) (1-y)/(a_plus1^2);((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) - a_plus1*u_plus1)) ((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 - u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1));((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) + a_plus1*u_plus1 )) -((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 + u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1))];
        %% R at i-1/2
        R_minus1 = [(1 - 0.5*(y-1)*(u_minus1^2)/(a_minus1^2)) ((y-1)*u_minus1/(a_minus1^2)) (1-y)/(a_minus1^2);((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) - a_minus1*u_minus1)) ((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 - u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1));((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) + a_minus1*u_minus1 )) -((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 + u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1))];
       
        %% R at i-3/2
        R_minus2 = [(1 - 0.5*(y-1)*(u_minus2^2)/(a_minus2^2)) ((y-1)*u_minus2/(a_minus2^2)) (1-y)/(a_minus2^2);((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) - a_minus2*u_minus2)) ((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 - u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1));((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) + a_minus2*u_minus2 )) -((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 + u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1))];
        %% eigenmatrix
        L_plus = [abs(u_plus1) 0 0;0 abs(u_plus1 + a_plus1) 0;0 0 abs(u_plus1 - a_minus1)];
        L_minus = [abs(u_minus1) 0 0;0 abs(u_minus1 + a_minus1) 0;0 0 abs(u_minus1 - a_minus1)];
        
        %Q_plus2 = [Q1(i+2,k)-Q1(i+1,k);Q2(i+2,k)-Q2(i+1,k);Q3(i+2,k)-Q3(i+1,k)];
        Q_plus1 = [Q1(i+1,k)-Q1(i,k);Q2(i+1,k)-Q2(i,k);Q3(i+1,k)-Q3(i,k)];
        Q_minus1 = [Q1(i,k)-Q1(i-1,k);Q2(i,k)-Q2(i-1,k);Q3(i,k)-Q3(i-1,k)];
        Q_minus2 = [Q1(i-1,k)-Q1(i-2,k);Q2(i-1,k)-Q2(i-2,k);Q3(i-1,k)-Q3(i-2,k)];        
        
        %M_plus2=R_plus2*Q_plus2;%M at i+3/2
        M_plus1=R_plus1*Q_plus1;%M at i+1/2
        M_minus1=R_minus1*Q_minus1;%M at i-1/2
        M_minus2=R_minus2*Q_minus2;% M at i-3/2
        for v=1:3
            N_plus(v,1)=M_plus1(v)-M_minus1(v);%N i+1/2
            N_minus(v,1)=M_minus1(v)-M_minus2(v);%N at i-1/2, ignoring the M at the -1 node
        end
        AQ_plus = Rinv_plus*L_plus*N_plus;
        AQ_minus = Rinv_minus*L_minus*N_minus;
        Q1(i,k+1) = Q1(i,k) - dt * 0.5 * ( F1(i+1,k) - F1(i-1,k) - AQ_plus(1) + AQ_minus(1))/dx;
        Q2(i,k+1) = Q2(i,k) - dt *(-g2(i,k)+ 0.5 * ( F2(i+1,k) - F2(i-1,k) - AQ_plus(2) + AQ_minus(2))/dx);
        Q3(i,k+1) = Q3(i,k) - dt * 0.5 * ( F3(i+1,k) - F3(i-1,k) - AQ_plus(3) + AQ_minus(3))/dx;
        %% Extracting the primitive variables
        rho(i,k+1) =(Q1(i,k+1))./S(i);  
        u(i,k+1) = Q2(i,k+1)./(rho(i,k+1).*S(i));
        E(i,k+1) = Q3(i,k+1)./(rho(i,k+1).*S(i));
        p(i,k+1) = (y-1) * ( rho(i,k+1).*E(i,k+1) - 0.5*rho(i,k+1).*(u(i,k+1).^2) );
        H(i,k+1) = E(i,k+1) + p(i,k+1)./rho(i,k+1);
        g2(i,k+1)=p(i,k+1)*dS(i);
        F1(i,k+1) = rho(i,k+1).*u(i,k+1).*S(i);
        F2(i,k+1) = ( rho(i,k+1).*(u(i,k+1).^2) + p(i,k+1) ) .* S(i);
        F3(i,k+1) = rho(i,k+1).*u(i,k+1).*H(i,k+1).*S(i); 

         else
             %% Roe averaging
        rho_plus2= sqrt(rho(i+2,k)*rho(i+1,k));%rho at i+3/2
        rho_plus1 = sqrt(rho(i,k)*rho(i+1,k));%rho at i+1/2
        rho_minus1 = sqrt(rho(i,k)*rho(i-1,k));%rho at i-1/2
        rho_minus2 = sqrt(rho(i-2,k)*rho(i-1,k));%rho at i-3/2
        u_plus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i+1,k) )*u(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%Velocity at i+1/2
        u_plus2 = ( sqrt(rho(i+2,k))*u(i+2,k) + sqrt(rho(i+1,k))*u(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%Velocity at i+3/2
        u_minus1 = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%Velocity at i-1/2
        u_minus2 = ( sqrt(rho(i-2,k))*u(i-2,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%Velocity at i-/2
        H_plus2 = ( sqrt(rho(i+2,k))*H(i+2,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i+2,k)) + sqrt(rho(i+1,k)) );%H at i+3/2   
        H_plus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );%H at i+1/2
        H_minus1 = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );%H at i-1/2
        H_minus2 = ( sqrt(rho(i-2,k))*H(i-2,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i-2,k)) + sqrt(rho(i-1,k)) );%H at i-3/2
        a_plus2 = sqrt( (y-1)*(H_plus2 - 0.5*(u_plus2^2)) );
        a_plus1 = sqrt( (y-1)*(H_plus1 - 0.5*(u_plus1^2)) );
        a_minus1 = sqrt( (y-1)*(H_minus1 - 0.5*(u_minus1^2)) );
        a_minus2 = sqrt( (y-1)*(H_minus2 - 0.5*(u_minus2^2)) );
        %% R inverse
        Rinv_plus = [1 (rho_plus1/(sqrt(2)*a_plus1)) (rho_plus1/(sqrt(2)*a_plus1));u_plus1 ((u_plus1 + a_plus1)*(rho_plus1/(sqrt(2)*a_plus1))) ((u_plus1 - a_plus1)*(rho_plus1/(sqrt(2)*a_plus1)));(0.5*(u_plus1^2)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) + a_plus1*u_plus1)) ((rho_plus1/(sqrt(2)*a_plus1))*(0.5*(u_plus1^2) + (a_plus1^2)/(y-1) - a_plus1*u_plus1))];                                     
        Rinv_minus = [1 (rho_minus1/(sqrt(2)*a_minus1)) (rho_minus1/(sqrt(2)*a_minus1));u_minus1 ((u_minus1 + a_minus1)*(rho_minus1/(sqrt(2)*a_minus1))) ((u_minus1 - a_minus1)*(rho_minus1/(sqrt(2)*a_minus1)));(0.5*(u_minus1^2)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) + a_minus1*u_minus1)) ((rho_minus1/(sqrt(2)*a_minus1))*(0.5*(u_minus1^2) + (a_minus1^2)/(y-1) - a_minus1*u_minus1))];
        %% R at i+3/2
        R_plus2 = [(1 - 0.5*(y-1)*(u_plus2^2)/(a_plus2^2)) ((y-1)*u_plus2/(a_plus2^2)) (1-y)/(a_plus2^2);((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) - a_plus2*u_plus2)) ((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 - u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1));((1/(rho_plus2*a_plus2*sqrt(2)))*(0.5*(y-1)*(u_plus2^2) + a_plus2*u_plus2 )) -((1/(rho_plus2*a_plus2*sqrt(2)))*(a_plus2 + u_plus2*(y-1))) ((1/(rho_plus2*a_plus2*sqrt(2)))*(y-1))];
        %% R at i+1/2
        R_plus1 = [(1 - 0.5*(y-1)*(u_plus1^2)/(a_plus1^2)) ((y-1)*u_plus1/(a_plus1^2)) (1-y)/(a_plus1^2);((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) - a_plus1*u_plus1)) ((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 - u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1));((1/(rho_plus1*a_plus1*sqrt(2)))*(0.5*(y-1)*(u_plus1^2) + a_plus1*u_plus1 )) -((1/(rho_plus1*a_plus1*sqrt(2)))*(a_plus1 + u_plus1*(y-1))) ((1/(rho_plus1*a_plus1*sqrt(2)))*(y-1))];
        %% R at i-1/2
        R_minus1 = [(1 - 0.5*(y-1)*(u_minus1^2)/(a_minus1^2)) ((y-1)*u_minus1/(a_minus1^2)) (1-y)/(a_minus1^2);((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) - a_minus1*u_minus1)) ((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 - u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1));((1/(rho_minus1*a_minus1*sqrt(2)))*(0.5*(y-1)*(u_minus1^2) + a_minus1*u_minus1 )) -((1/(rho_minus1*a_minus1*sqrt(2)))*(a_minus1 + u_minus1*(y-1))) ((1/(rho_minus1*a_minus1*sqrt(2)))*(y-1))];
       
        %% R at i-3/2
        R_minus2 = [(1 - 0.5*(y-1)*(u_minus2^2)/(a_minus2^2)) ((y-1)*u_minus2/(a_minus2^2)) (1-y)/(a_minus2^2);((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) - a_minus2*u_minus2)) ((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 - u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1));((1/(rho_minus2*a_minus2*sqrt(2)))*(0.5*(y-1)*(u_minus2^2) + a_minus2*u_minus2 )) -((1/(rho_minus2*a_minus2*sqrt(2)))*(a_minus2 + u_minus2*(y-1))) ((1/(rho_minus2*a_minus2*sqrt(2)))*(y-1))];
        %% eigenmatrix
        L_plus = [abs(u_plus1) 0 0;0 abs(u_plus1 + a_plus1) 0;0 0 abs(u_plus1 - a_minus1)];
        L_minus = [abs(u_minus1) 0 0;0 abs(u_minus1 + a_minus1) 0;0 0 abs(u_minus1 - a_minus1)];
        
        Q_plus2 = [Q1(i+2,k)-Q1(i+1,k);Q2(i+2,k)-Q2(i+1,k);Q3(i+2,k)-Q3(i+1,k)];
        Q_plus1 = [Q1(i+1,k)-Q1(i,k);Q2(i+1,k)-Q2(i,k);Q3(i+1,k)-Q3(i,k)];
        Q_minus1 = [Q1(i,k)-Q1(i-1,k);Q2(i,k)-Q2(i-1,k);Q3(i,k)-Q3(i-1,k)];
        Q_minus2 = [Q1(i-1,k)-Q1(i-2,k);Q2(i-1,k)-Q2(i-2,k);Q3(i-1,k)-Q3(i-2,k)];        
        
        M_plus2=R_plus2*Q_plus2;%M at i+3/2
        M_plus1=R_plus1*Q_plus1;%M at i+1/2
        M_minus1=R_minus1*Q_minus1;%M at i-1/2
        M_minus2=R_minus2*Q_minus2;% M at i-3/2
        for v=1:3
            N_plus(v,1)=M_plus1(v)-M_minus1(v);%N i+1/2
            N_minus(v,1)=M_minus1(v)-M_minus2(v);%N at i-1/2
        end
        AQ_plus = Rinv_plus*L_plus*N_plus;
        AQ_minus = Rinv_minus*L_minus*N_minus;
        %% The scheme
        Q1(i,k+1) = Q1(i,k) - dt * 0.5 * ( F1(i+1,k) - F1(i-1,k) - AQ_plus(1) + AQ_minus(1))/dx;
        Q2(i,k+1) = Q2(i,k) - dt *(-g2(i,k)+ 0.5 * ( F2(i+1,k) - F2(i-1,k) - AQ_plus(2) + AQ_minus(2))/dx);
        Q3(i,k+1) = Q3(i,k) - dt * 0.5 * ( F3(i+1,k) - F3(i-1,k) - AQ_plus(3) + AQ_minus(3))/dx;
        %% Extracting the primitive variables
     rho(i,k+1) =(Q1(i,k+1))./S(i);  
     u(i,k+1) = Q2(i,k+1)./(rho(i,k+1).*S(i));
     E(i,k+1) = Q3(i,k+1)./(rho(i,k+1).*S(i));
     p(i,k+1) = (y-1) * ( rho(i,k+1).*E(i,k+1) - 0.5*rho(i,k+1).*(u(i,k+1).^2) );
     H(i,k+1) = E(i,k+1) + p(i,k+1)./rho(i,k+1);
     g2(i,k+1)=p(i,k+1)*dS(i);
     F1(i,k+1) = rho(i,k+1).*u(i,k+1).*S(i);
     F2(i,k+1) = ( rho(i,k+1).*(u(i,k+1).^2) + p(i,k+1) ) .* S(i);
     F3(i,k+1) = rho(i,k+1).*u(i,k+1).*H(i,k+1).*S(i); 
         end
     end
     %% Extracting the primitive variables
%      rho(3:N-2,k+1) =(Q1(3:N-2,k+1))./S(3:N-2);  
%      u(3:N-2,k+1) = Q2(3:N-2,k+1)./(rho(3:N-2,k+1).*S(3:N-2));
%      E(3:N-2,k+1) = Q3(2:N-2,k+1)./(rho(2:N-2,k+1).*S(3:N-2));
%      p(3:N-2,k+1) = (y-1) * ( rho(3:N-2,k+1).*E(3:N-2,k+1) - 0.5*rho(3:N-2,k+1).*(u(3:N-2,k+1).^2) );
%      H(3:N-2,k+1) = E(3:N-2,k+1) + p(3:N-2,k+1)./rho(3:N-2,k+1);
%      g2(3:N-2,k+1)=p(3:N-2,k+1)*dS(3:N-2);
%      F1(3:N-2,k+1) = rho(3:N-2,k+1).*u(3:N-2,k+1).*S(3:N-2);
%      F2(3:N-2,k+1) = ( rho(3:N-2,k+1).*(u(3:N-2,k+1).^2) + p(3:N-2,k+1) ) .* S(3:N-2);
%      F3(3:N-2,k+1) = rho(3:N-2,k+1).*u(3:N-2,k+1).*H(3:N-2,k+1).*S(3:N-2); 
 end
 %% Plotting
 plot(x,rho(:,Nt))
 figure
 plot(x,u(:,Nt))
 
 
 
 
 
 