%% Case 1 Laval Nozzle Euler
clc;
clear all;
%% Defining Variables
L = 10;
N = 100;
y = 1.4;
R = 8.314;
dx = (L/(N-1));
%% Variables at left side
rho_l = 1.1409;
u_l = 65.451;
p_l = 9.7534*10^4;
M = 28.97 * 10^-3;

T_L = p_l*M/(rho_l*R);% Temprature required for speed of sound
a_L = sqrt(y*R*T_L);%Speed of sound
dt = 0.01*dx/a_L;% Using Courant number and CFL=0.01, its arbitrary
tf = 0.00;% Final time

Nt = fix(tf/dt);
%% Variables at right side
rho_r = 1.1008;
u_r = 113.060;
p_r = 9.2772*10^4;

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
        S(i)=1+1.5*(1-0.2*x(i))^2;
        dS=3*(0.2*x(i)-1);
end

for i=N/2+1:N
        S(i)=1+0.5*(1-0.2*x(i))^2;
        dS(i)=(0.2*x(i)-1);
end

%% BCs for primitive variables
u(:,1)=u_l;
rho(:,1)=rho_l;
p(:,1)=p_l;
E(:,1)=E_l;
H(:,1)=H_l;
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
 
 F1(1,:) = rho(1,:).*u(1,:).*S(1);
 F2(1,:) = ( rho(1,:).*(u(1,:).^2) + p(1,:) ) .* S(1);
 F3(1,:) = rho(1,:).*u(1,:).*H(1,:).*S(1);
 
 F1(N,:) = rho(N,:).*u(N,:).*S(N);
 F2(N,:) = ( rho(N,:).*(u(N,:).^2) + p(N,:) ) .* S(N);
 F3(N,:) = rho(N,:).*u(N,:).*H(N,:).*S(N);
%% g2
g2(:,1)=p(:,1).*dS(:,1);
g2(1,:)=p(1,:).*S(1);
g2(N,:)=p(N,:).*S(N);
res=1;
k=1;
 while res>=10e-4
     for i=2:N-1
         %% Roe averaging
        rho_plus = sqrt(rho(i,k)*rho(i+1,k));
        rho_minus = sqrt(rho(i,k)*rho(i-1,k));
        u_plus = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i+1,k))*u(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );
        u_p2 = u_plus*u_plus;
        u_p3 = u_plus*u_plus*u_plus;
       
        u_minus = ( sqrt(rho(i,k))*u(i,k) + sqrt(rho(i-1,k))*u(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );
        u_m_sqr = u_minus*u_minus;
        u_m_cube = u_minus*u_minus*u_minus;
        
        H_plus = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i+1,k))*H(i+1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i+1,k)) );
        H_minus = ( sqrt(rho(i,k))*H(i,k) + sqrt(rho(i-1,k))*H(i-1,k) )/( sqrt(rho(i,k)) + sqrt(rho(i-1,k)) );
        
        a_plus = sqrt( (y-1)*(H_plus - 0.5*(u_plus^2)) );
        a_minus = sqrt( (y-1)*(H_minus - 0.5*(u_minus^2)) );
        %% R inverse
        Rinv_plus = [1 (rho_plus/(sqrt(2)*a_plus)) (rho_plus/(sqrt(2)*a_plus));u_plus ((u_plus + a_plus)*(rho_plus/(sqrt(2)*a_plus))) ((u_plus - a_plus)*(rho_plus/(sqrt(2)*a_plus)));(0.5*(u_plus^2)) ((rho_plus/(sqrt(2)*a_plus))*(0.5*(u_plus^2) + (a_plus^2)/(y-1) + a_plus*u_plus)) ((rho_plus/(sqrt(2)*a_plus))*(0.5*(u_plus^2) + (a_plus^2)/(y-1) - a_plus*u_plus))];                                     
        Rinv_minus = [1 (rho_minus/(sqrt(2)*a_minus)) (rho_minus/(sqrt(2)*a_minus));u_minus ((u_minus + a_minus)*(rho_minus/(sqrt(2)*a_minus))) ((u_minus - a_minus)*(rho_minus/(sqrt(2)*a_minus)));(0.5*(u_minus^2)) ((rho_minus/(sqrt(2)*a_minus))*(0.5*(u_minus^2) + (a_minus^2)/(y-1) + a_minus*u_minus)) ((rho_minus/(sqrt(2)*a_minus))*(0.5*(u_minus^2) + (a_minus^2)/(y-1) - a_minus*u_minus))];
        %% R
        R_plus = [(1 - 0.5*(y-1)*(u_plus^2)/(a_plus^2)) ((y-1)*u_plus/(a_plus^2)) (1-y)/(a_plus^2);((1/(rho_plus*a_plus*sqrt(2)))*(0.5*(y-1)*(u_plus^2) - a_plus*u_plus)) ((1/(rho_plus*a_plus*sqrt(2)))*(a_plus - u_plus*(y-1))) ((1/(rho_plus*a_plus*sqrt(2)))*(y-1));((1/(rho_plus*a_plus*sqrt(2)))*(0.5*(y-1)*(u_plus^2) + a_plus*u_plus )) -((1/(rho_plus*a_plus*sqrt(2)))*(a_plus + u_plus*(y-1))) ((1/(rho_plus*a_plus*sqrt(2)))*(y-1))];
        R_minus = [(1 - 0.5*(y-1)*(u_minus^2)/(a_minus^2)) ((y-1)*u_minus/(a_minus^2)) (1-y)/(a_minus^2);((1/(rho_minus*a_minus*sqrt(2)))*(0.5*(y-1)*(u_minus^2) - a_minus*u_minus)) ((1/(rho_minus*a_minus*sqrt(2)))*(a_minus - u_minus*(y-1))) ((1/(rho_minus*a_minus*sqrt(2)))*(y-1));((1/(rho_minus*a_minus*sqrt(2)))*(0.5*(y-1)*(u_minus^2) + a_minus*u_minus )) -((1/(rho_minus*a_minus*sqrt(2)))*(a_minus + u_minus*(y-1))) ((1/(rho_minus*a_minus*sqrt(2)))*(y-1))];
        %% eigenmatrix
        L_plus = [abs(u_plus) 0 0;0 abs(u_plus + a_plus) 0;0 0 abs(u_plus - a_minus)];
        L_minus = [abs(u_minus) 0 0;0 abs(u_minus + a_minus) 0;0 0 abs(u_minus - a_minus)];
        
        A_plus = Rinv_plus * L_plus * R_plus;
        A_minus = Rinv_minus * L_minus * R_minus;
        
        Q_plus = [Q1(i+1,k)-Q1(i,k);Q2(i+1,k)-Q2(i,k);Q3(i+1,k)-Q3(i,k)];
        Q_minus = [Q1(i,k)-Q1(i-1,k);Q2(i,k)-Q2(i-1,k);Q3(i,k)-Q3(i-1,k)];
        
        AQ_plus = A_plus * Q_plus;
        AQ_minus = A_minus * Q_minus;
        %% The scheme
        Q1(i,k+1) = Q1(i,k) - dt * 0.5 * ( F1(i+1,k) - F1(i-1,k) - AQ_plus(1) + AQ_minus(1))/dx;
        Q2(i,k+1) = Q2(i,k) - dt * 0.5 * ( F2(i+1,k) - F2(i-1,k) - AQ_plus(2) + AQ_minus(2)-g2(i,k))/dx;
        Q3(i,k+1) = Q3(i,k) - dt * 0.5 * ( F3(i+1,k) - F3(i-1,k) - AQ_plus(3) + AQ_minus(3))/dx;
        
     end
     %% Extracting the primitive variables
     rho(2:N-1,k+1) =(Q1(2:N-1,k+1))./S(2:N-1);
     u(2:N-1,k+1) = Q2(2:N-1,k+1)./(rho(2:N-1,k+1).*S(2:N-1));
     E(2:N-1,k+1) = Q3(2:N-1,k+1)./(rho(2:N-1,k+1).*S(2:N-1));
     p(2:N-1,k+1) = (y-1) * ( rho(2:N-1,k+1).*E(2:N-1,k+1) - 0.5*rho(2:N-1,k+1).*(u(2:N-1,k+1).^2) );
     H(2:N-1,k+1) = E(2:N-1,k+1) + p(2:N-1,k+1)./rho(2:N-1,k+1);
     
     F1(2:N-1,k+1) = rho(2:N-1,k+1).*u(2:N-1,k+1).*S(2:N-1);
     F2(2:N-1,k+1) = ( rho(2:N-1,k+1).*(u(2:N-1,k+1).^2) + p(2:N-1,k+1) ) .* S(2:N-1);
     F3(2:N-1,k+1) = rho(2:N-1,k+1).*u(2:N-1,k+1).*H(2:N-1,k+1).*S(2:N-1);
     g2(2:N-1,k+1)=p(2:N-1,k+1).*dS(2:N-1)';
     res=max((Q1(:,k+1)-Q1(:,k))./Q1(:,1));
     k=k+1;
 end
 %% Plotting
 plot(x,rho(:,Nt))
 figure
 
 
 
 
 
 
 
