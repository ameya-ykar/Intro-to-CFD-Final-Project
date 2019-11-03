clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%% INTITIALIZE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%



L = 10;
N =300;
gamma = 1.4;
R = 8.314;
delx = (L/(N-1));
x = 0:delx:10;
S = zeros(N,1);
dS = zeros(N,1);
    
for i=1:N
    
    if(i<=N/2)
        S(i,1) = 1 + 1.5*((1 - 0.2*x(i))^2);
        dS(i,1) = 0.6*(0.2*x(i) - 1);      
    else
        S(i,1) = 1 + 0.5*((1 - 0.2*x(i))^2);
        dS(i,1) = 0.2*(0.2*x(i) - 1);
    end
    
end


% rhoL = 1.1409;
% uL = 65.451;    % SUBSONIC
% pL = 9.7534 * 10^4;
% 
% rhoR = 1.1008;
% uR = 113.060;
% pR = 9.2772 * 10^4;

rhoL = 1.1308;
uL = 80.054;    % TRANSONIC
pL = 9.6328 * 10^4; 

rhoR = 1.0275;
uR = 150.535;
pR = 8.4974 * 10^4;

M = 28.97 * 10^-3; % Molar mass of air (Kg/mol)

TL = pL*M/(rhoL*R); % Temperature of left extreme
aL = sqrt(gamma*R*TL);
CFL=0.5;
%delt=10e-6;

delt = CFL*delx/688.2959;%Transonic
%delt=CFL*delx/523.1340;%Subsonic
tfinal = 1;

Nt = fix(tfinal/delt);

EL = ( pL/(gamma-1) + 0.5*rhoL*(uL^2) )/rhoL;
ER = ( pR/(gamma-1) + 0.5*rhoR*(uR^2) )/rhoR;

HL = EL + pL/rhoL;
HR = ER + pR/rhoR;

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
G2 = zeros(N,Nt);

%%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS FOR VARIABLES %%%%%%%%%%%%%%%%%%%%%
for i=1:N
    if (i <= N/2)  % 0 <= x <= 5
        rho(i,1) = rhoL;
        u(i,1) = uL; 
        p(i,1) = pL;
        H(i,1) = EL + pL/rhoL;
        E(i,1) = EL;       
    else 
        rho(i,1) = rhoR;  % 5 <= x <= 10
        u(i,1) = uR;
        p(i,1) = pR;
        H(i,1) =  ER + pR/rhoR;
        E(i,1) = ER;
    end
end
%%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS FOR VARIABLES %%%%%%%%%%%%%%%%%%%%
u(1,:) = uL;
u(N,:) = uR;
rho(1,:) = rhoL;
rho(N,:) = rhoR;
p(1,:) = pL;
p(N,:) = pR;
E(1,:) = EL;
E(N,:) = ER;
H(1,:) = HL;
H(N,:) = HR;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Q, F AND G MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Q1(:,1) = rho(:,1) .* S(:,1);
 Q2(:,1) = rho(:,1).*u(:,1).*S(:,1);  % INITIAL CONDITIONS
 Q3(:,1) = rho(:,1).*E(:,1).*S(:,1);

 Q1(1,:) = rho(1,:) .* S(1,1);
 Q2(1,:) = rho(1,:).*u(1,:).*S(1,1); % BOUNDARY CONDITIONS
 Q3(1,:) = rho(1,:).*E(1,:).*S(1,1);
 
 Q1(N,:) = rho(N,:) .* S(N,1);
 Q2(N,:) = rho(N,:).*u(N,:).*S(N,1); % BOUNDARY CONDITIONS
 Q3(N,:) = rho(N,:).*E(N,:).*S(N,1);
 
 
 F1(:,1) = rho(:,1).*u(:,1).*S(:,1);
 F2(:,1) = ( rho(:,1).*(u(:,1).^2) + p(:,1) ) .* S(:,1); % INITIAL CONDITIONS
 F3(:,1) = rho(:,1).*u(:,1).*H(:,1).*S(:,1);
 
 F1(1,:) = rho(1,:).*u(1,:).*S(1,1);
 F2(1,:) = ( rho(1,:).*(u(1,:).^2) + p(1,:) ) .* S(1,1); % BOUNDARY CONDITIONS
 F3(1,:) = rho(1,:).*u(1,:).*H(1,:).*S(1,1);
 
 F1(N,:) = rho(N,:).*u(N,:).*S(N,1);
 F2(N,:) = ( rho(N,:).*(u(N,:).^2) + p(N,:) ) .* S(N,1); % BOUNDARY CONDITIONS
 F3(N,:) = rho(N,:).*u(N,:).*H(N,:).*S(N,1);
 
 G2(:,1) = p(:,1).*dS(:,1);    % INITIAL CONDITIONS
 G2(1,:) = p(1,:).*dS(1,1);    % BOUNDARY CONDITIONS
 G2(N,:) = p(N,:).*dS(N,1);    % BOUNDARY CONDITIONS
 
 %%%%%%%%%%%%%%%%%%%%%%%%%% GHOST NODES VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
S_ghostL = 1 + 1.5*((1 - 0.2*(-delx))^2);    % LEFT
S_ghostR = 1 + 0.5*((1 - 0.2*(10+delx))^2);  % RIGHT
 
% LEFT

Q1_ghostL = rhoL * S_ghostL;
Q2_ghostL = rhoL * uL * S_ghostL;
Q3_ghostL = rhoL * EL * S_ghostL;

% RIGHT

Q1_ghostR = rhoR * S_ghostR;
Q2_ghostR = rhoR * uR * S_ghostR;
Q3_ghostR = rhoR * ER * S_ghostR;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Notation _m denotes i - 1/2
 % Notation _p denotes i + 1/2

 for n=1:Nt-1
     for i=2:N-1                
        
        rho_p = sqrt(rho(i,n)*rho(i+1,n));
        rho_m = sqrt(rho(i,n)*rho(i-1,n));
        
        u_p = ( sqrt(rho(i,n))*u(i,n) + sqrt(rho(i+1,n))*u(i+1,n) )/( sqrt(rho(i,n)) + sqrt(rho(i+1,n)) );
        u_p_sqr = u_p*u_p;
        u_p_cube = u_p*u_p*u_p;        
       
        u_m = ( sqrt(rho(i,n))*u(i,n) + sqrt(rho(i-1,n))*u(i-1,n) )/( sqrt(rho(i,n)) + sqrt(rho(i-1,n)) );        
        u_m_sqr = u_m*u_m;
        u_m_cube = u_m*u_m*u_m;        
        
        H_p = ( sqrt(rho(i,n))*H(i,n) + sqrt(rho(i+1,n))*H(i+1,n) )/( sqrt(rho(i,n)) + sqrt(rho(i+1,n)) );
        H_m = ( sqrt(rho(i,n))*H(i,n) + sqrt(rho(i-1,n))*H(i-1,n) )/( sqrt(rho(i,n)) + sqrt(rho(i-1,n)) );

        
        a_p = sqrt( (gamma-1)*(H_p - 0.5*(u_p^2)) );
        a_m = sqrt( (gamma-1)*(H_m - 0.5*(u_m^2)) );
        
        Rinv_p = [1 (rho_p/(sqrt(2)*a_p)) (rho_p/(sqrt(2)*a_p));u_p ((u_p + a_p)*(rho_p/(sqrt(2)*a_p))) ((u_p - a_p)*(rho_p/(sqrt(2)*a_p)));(0.5*(u_p^2)) ((rho_p/(sqrt(2)*a_p))*(0.5*(u_p^2) + (a_p^2)/(gamma-1) + a_p*u_p)) ((rho_p/(sqrt(2)*a_p))*(0.5*(u_p^2) + (a_p^2)/(gamma-1) - a_p*u_p))];                                     
        Rinv_m = [1 (rho_m/(sqrt(2)*a_m)) (rho_m/(sqrt(2)*a_m));u_m ((u_m + a_m)*(rho_m/(sqrt(2)*a_m))) ((u_m - a_m)*(rho_m/(sqrt(2)*a_m)));(0.5*(u_m^2)) ((rho_m/(sqrt(2)*a_m))*(0.5*(u_m^2) + (a_m^2)/(gamma-1) + a_m*u_m)) ((rho_m/(sqrt(2)*a_m))*(0.5*(u_m^2) + (a_m^2)/(gamma-1) - a_m*u_m))];
        
        R_p = [(1 - 0.5*(gamma-1)*(u_p^2)/(a_p^2)) ((gamma-1)*u_p/(a_p^2)) (1-gamma)/(a_p^2);((1/(rho_p*a_p*sqrt(2)))*(0.5*(gamma-1)*(u_p^2) - a_p*u_p)) ((1/(rho_p*a_p*sqrt(2)))*(a_p - u_p*(gamma-1))) ((1/(rho_p*a_p*sqrt(2)))*(gamma-1));((1/(rho_p*a_p*sqrt(2)))*(0.5*(gamma-1)*(u_p^2) + a_p*u_p )) (-(1/(rho_p*a_p*sqrt(2)))*(a_p + u_p*(gamma-1))) ((1/(rho_p*a_p*sqrt(2)))*(gamma-1))];
        R_m = [(1 - 0.5*(gamma-1)*(u_m^2)/(a_m^2)) ((gamma-1)*u_m/(a_m^2)) (1-gamma)/(a_m^2);((1/(rho_m*a_m*sqrt(2)))*(0.5*(gamma-1)*(u_m^2) - a_m*u_m)) ((1/(rho_m*a_m*sqrt(2)))*(a_m - u_m*(gamma-1))) ((1/(rho_m*a_m*sqrt(2)))*(gamma-1));((1/(rho_m*a_m*sqrt(2)))*(0.5*(gamma-1)*(u_m^2) + a_m*u_m )) (-(1/(rho_m*a_m*sqrt(2)))*(a_m + u_m*(gamma-1))) ((1/(rho_m*a_m*sqrt(2)))*(gamma-1))];
        
        LambdaA_p = [abs(u_p) 0 0;0 abs(u_p + a_p) 0;0 0 abs(u_p - a_p)];
        LambdaA_m = [abs(u_m) 0 0;0 abs(u_m + a_m) 0;0 0 abs(u_m - a_m)];    
        
        
        Q_p = [Q1(i+1,n)-Q1(i,n);Q2(i+1,n)-Q2(i,n);Q3(i+1,n)-Q3(i,n)];
        Q_m = [Q1(i,n)-Q1(i-1,n);Q2(i,n)-Q2(i-1,n);Q3(i,n)-Q3(i-1,n)];
        
        M_p = R_p*Q_p;
        M_m = R_m*Q_m;
        
        
        if(i == 2)
            
        rho_pp = sqrt(rho(i+1,n)*rho(i+2,n));
        rho_mm = rhoL;        

        u_pp = ( sqrt(rho(i+1,n))*u(i+1,n) + sqrt(rho(i+2,n))*u(i+2,n) )/( sqrt(rho(i+1,n)) + sqrt(rho(i+2,n)) );
        u_pp_sqr = u_pp*u_pp;
        u_pp_cube = u_pp*u_pp*u_pp;
        
        u_mm = uL;
        u_mm_sqr = u_mm*u_mm;
        u_mm_cube = u_mm*u_mm*u_mm;
        
        H_pp = ( sqrt(rho(i+1,n))*H(i+1,n) + sqrt(rho(i+2,n))*H(i+2,n) )/( sqrt(rho(i+1,n)) + sqrt(rho(i+2,n)) );
        H_mm = HL;        
        
        a_pp = sqrt( (gamma-1)*(H_pp - 0.5*(u_pp^2)) );
        a_mm = sqrt( (gamma-1)*(H_mm - 0.5*(u_mm^2)) );
        
        R_pp = [(1 - 0.5*(gamma-1)*(u_pp^2)/(a_pp^2)) ((gamma-1)*u_pp/(a_pp^2)) (1-gamma)/(a_pp^2);((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) - a_pp*u_pp)) ((1/(rho_pp*a_pp*sqrt(2)))*(a_pp - u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1));((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) + a_pp*u_pp )) (-(1/(rho_pp*a_pp*sqrt(2)))*(a_pp + u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1))];
        R_mm = [(1 - 0.5*(gamma-1)*(u_mm^2)/(a_mm^2)) ((gamma-1)*u_mm/(a_mm^2)) (1-gamma)/(a_mm^2);((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) - a_mm*u_mm)) ((1/(rho_mm*a_mm*sqrt(2)))*(a_mm - u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1));((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) + a_mm*u_mm )) (-(1/(rho_mm*a_mm*sqrt(2)))*(a_mm + u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1))];
        
        Q_pp = [Q1(i+2,n)-Q1(i+1,n);Q2(i+2,n)-Q2(i+1,n);Q3(i+2,n)-Q3(i+1,n)];               
        Q_mm = [Q1(i-1,n)-Q1_ghostL;Q2(i-1,n)-Q2_ghostL;Q3(i-1,n)-Q3_ghostL];  % CHANGED, USING GHOST Q AT LEFT EXTREME
        
        M_pp = R_pp*Q_pp;
        M_mm = R_mm*Q_mm;        
        
        [mod] = minmod(M_m(1),M_p(1),M_pp(1));
        modp1 = mod;
        [mod] = minmod(M_m(2),M_p(2),M_pp(2));
        modp2 = mod;
        [mod] = minmod(M_m(3),M_p(3),M_pp(3));
        modp3 = mod;
        
        [mod] = minmod(M_mm(1),M_m(1),M_p(1));
        modm1 = mod;
        [mod] = minmod(M_mm(2),M_m(2),M_p(2));
        modm2 = mod;
        [mod] = minmod(M_mm(3),M_m(3),M_p(3));
        modm3 = mod;


        Mod_p = [modp1;modp2;modp3];
        Mod_m = [modm1;modm2;modm3];

        
%          Mod_p = [M_m(1);M_m(2);M_m(3)];
%          Mod_m = [M_mm(1);M_mm(2);M_mm(3)];        
          
        
        elseif(i == (N-1))
            
        rho_pp = rhoR;
        rho_mm = sqrt(rho(i-1,n)*rho(i-2,n));        

        u_pp = uR;
        u_pp_sqr = u_pp*u_pp;
        u_pp_cube = u_pp*u_pp*u_pp;
        
        u_mm = ( sqrt(rho(i-1,n))*u(i-1,n) + sqrt(rho(i-2,n))*u(i-2,n) )/( sqrt(rho(i-1,n)) + sqrt(rho(i-2,n)) );
        u_mm_sqr = u_mm*u_mm;
        u_mm_cube = u_mm*u_mm*u_mm;
        
        H_pp = HR;
        H_mm = ( sqrt(rho(i-1,n))*H(i-1,n) + sqrt(rho(i-2,n))*H(i-2,n) )/( sqrt(rho(i-1,n)) + sqrt(rho(i-2,n)) );        
        
        a_pp = sqrt( (gamma-1)*(H_pp - 0.5*(u_pp^2)) );
        a_mm = sqrt( (gamma-1)*(H_mm - 0.5*(u_mm^2)) );
        
        R_pp = [(1 - 0.5*(gamma-1)*(u_pp^2)/(a_pp^2)) ((gamma-1)*u_pp/(a_pp^2)) (1-gamma)/(a_pp^2);((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) - a_pp*u_pp)) ((1/(rho_pp*a_pp*sqrt(2)))*(a_pp - u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1));((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) + a_pp*u_pp )) (-(1/(rho_pp*a_pp*sqrt(2)))*(a_pp + u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1))];
        R_mm = [(1 - 0.5*(gamma-1)*(u_mm^2)/(a_mm^2)) ((gamma-1)*u_mm/(a_mm^2)) (1-gamma)/(a_mm^2);((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) - a_mm*u_mm)) ((1/(rho_mm*a_mm*sqrt(2)))*(a_mm - u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1));((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) + a_mm*u_mm )) (-(1/(rho_mm*a_mm*sqrt(2)))*(a_mm + u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1))];
        
        Q_pp = [Q1_ghostR-Q1(i+1,n);Q2_ghostR-Q2(i+1,n);Q3_ghostR-Q3(i+1,n)];   % CHANGED
        Q_mm = [Q1(i-1,n)-Q1(i-2,n);Q2(i-1,n)-Q2(i-2,n);Q3(i-1,n)-Q3(i-2,n)];
        
        M_pp = R_pp*Q_pp;
        M_mm = R_mm*Q_mm;        
%         
        [mod] = minmod(M_m(1),M_p(1),M_pp(1));
        modp1 = mod;
        [mod] = minmod(M_m(2),M_p(2),M_pp(2));
        modp2 = mod;
        [mod] = minmod(M_m(3),M_p(3),M_pp(3));
        modp3 = mod;
        
        [mod] = minmod(M_mm(1),M_m(1),M_p(1));
        modm1 = mod;
        [mod] = minmod(M_mm(2),M_m(2),M_p(2));
        modm2 = mod;
        [mod] = minmod(M_mm(3),M_m(3),M_p(3));
        modm3 = mod;
% 
% 
        Mod_p = [modp1;modp2;modp3];
        Mod_m = [modm1;modm2;modm3];
        
%          Mod_p = [M_m(1);M_m(2);M_m(3)];
%          Mod_m = [M_mm(1);M_mm(2);M_mm(3)]; 


        else
            
        rho_pp = sqrt(rho(i+1,n)*rho(i+2,n));
        rho_mm = sqrt(rho(i-1,n)*rho(i-2,n));        

        u_pp = ( sqrt(rho(i+1,n))*u(i+1,n) + sqrt(rho(i+2,n))*u(i+2,n) )/( sqrt(rho(i+1,n)) + sqrt(rho(i+2,n)) );
        u_pp_sqr = u_pp*u_pp;
        u_pp_cube = u_pp*u_pp*u_pp;
        
        u_mm = ( sqrt(rho(i-1,n))*u(i-1,n) + sqrt(rho(i-2,n))*u(i-2,n) )/( sqrt(rho(i-1,n)) + sqrt(rho(i-2,n)) );
        u_mm_sqr = u_mm*u_mm;
        u_mm_cube = u_mm*u_mm*u_mm;
        
        H_pp = ( sqrt(rho(i+1,n))*H(i+1,n) + sqrt(rho(i+2,n))*H(i+2,n) )/( sqrt(rho(i+1,n)) + sqrt(rho(i+2,n)) );
        H_mm = ( sqrt(rho(i-1,n))*H(i-1,n) + sqrt(rho(i-2,n))*H(i-2,n) )/( sqrt(rho(i-1,n)) + sqrt(rho(i-2,n)) );        
        
        a_pp = sqrt( (gamma-1)*(H_pp - 0.5*(u_pp^2)) );
        a_mm = sqrt( (gamma-1)*(H_mm - 0.5*(u_mm^2)) );
        
        R_pp = [(1 - 0.5*(gamma-1)*(u_pp^2)/(a_pp^2)) ((gamma-1)*u_pp/(a_pp^2)) (1-gamma)/(a_pp^2);((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) - a_pp*u_pp)) ((1/(rho_pp*a_pp*sqrt(2)))*(a_pp - u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1));((1/(rho_pp*a_pp*sqrt(2)))*(0.5*(gamma-1)*(u_pp^2) + a_pp*u_pp )) (-(1/(rho_pp*a_pp*sqrt(2)))*(a_pp + u_pp*(gamma-1))) ((1/(rho_pp*a_pp*sqrt(2)))*(gamma-1))];
        R_mm = [(1 - 0.5*(gamma-1)*(u_mm^2)/(a_mm^2)) ((gamma-1)*u_mm/(a_mm^2)) (1-gamma)/(a_mm^2);((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) - a_mm*u_mm)) ((1/(rho_mm*a_mm*sqrt(2)))*(a_mm - u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1));((1/(rho_mm*a_mm*sqrt(2)))*(0.5*(gamma-1)*(u_mm^2) + a_mm*u_mm )) (-(1/(rho_mm*a_mm*sqrt(2)))*(a_mm + u_mm*(gamma-1))) ((1/(rho_mm*a_mm*sqrt(2)))*(gamma-1))];
        
        Q_pp = [Q1(i+2,n)-Q1(i+1,n);Q2(i+2,n)-Q2(i+1,n);Q3(i+2,n)-Q3(i+1,n)];
        Q_mm = [Q1(i-1,n)-Q1(i-2,n);Q2(i-1,n)-Q2(i-2,n);Q3(i-1,n)-Q3(i-2,n)];
        
        M_pp = R_pp*Q_pp;
        M_mm = R_mm*Q_mm;        
        
        [mod] = minmod(M_m(1),M_p(1),M_pp(1));
        modp1 = mod;
        [mod] = minmod(M_m(2),M_p(2),M_pp(2));
        modp2 = mod;
        [mod] = minmod(M_m(3),M_p(3),M_pp(3));
        modp3 = mod;
        
        [mod] = minmod(M_mm(1),M_m(1),M_p(1));
        modm1 = mod;
        [mod] = minmod(M_mm(2),M_m(2),M_p(2));
        modm2 = mod;
        [mod] = minmod(M_mm(3),M_m(3),M_p(3));
        modm3 = mod;


        Mod_p = [modp1;modp2;modp3];
        Mod_m = [modm1;modm2;modm3];
        
%         Mod_p = [M_m(1);M_m(2);M_m(3)];
%         Mod_m = [M_mm(1);M_mm(2);M_mm(3)];        
          
        end            
        
        N_p = M_p - Mod_p;
        N_m = M_m - Mod_m;
        
        A_Q_p = Rinv_p * LambdaA_p * N_p;
        A_Q_m = Rinv_m * LambdaA_m * N_m;
        
        Q1(i,n+1) = Q1(i,n) - delt * 0.5 * ( F1(i+1,n) - F1(i-1,n) - A_Q_p(1,1) + A_Q_m(1,1))/delx;
        Q2(i,n+1) = Q2(i,n) - delt * ( ( 0.5 * ( F2(i+1,n) - F2(i-1,n) - A_Q_p(2,1) + A_Q_m(2,1))/delx) - G2(i,n) );
        Q3(i,n+1) = Q3(i,n) - delt * 0.5 * ( F3(i+1,n) - F3(i-1,n) - A_Q_p(3,1) + A_Q_m(3,1))/delx;
        
     end
     % UPDATE VARIABLES
     rho(2:N-1,n+1) =(Q1(2:N-1,n+1))./S(2:N-1,1);
     u(2:N-1,n+1) = Q2(2:N-1,n+1)./(rho(2:N-1,n+1).*S(2:N-1,1));
     E(2:N-1,n+1) = Q3(2:N-1,n+1)./(rho(2:N-1,n+1).*S(2:N-1,1));
     p(2:N-1,n+1) = (gamma-1) * ( rho(2:N-1,n+1).*E(2:N-1,n+1) - 0.5*rho(2:N-1,n+1).*(u(2:N-1,n+1).^2) );
     H(2:N-1,n+1) = E(2:N-1,n+1) + p(2:N-1,n+1)./rho(2:N-1,n+1);
     
     F1(2:N-1,n+1) = rho(2:N-1,n+1).*u(2:N-1,n+1).*S(2:N-1,1);
     F2(2:N-1,n+1) = ( rho(2:N-1,n+1).*(u(2:N-1,n+1).^2) + p(2:N-1,n+1) ) .* S(2:N-1,1);
     F3(2:N-1,n+1) = rho(2:N-1,n+1).*u(2:N-1,n+1).*H(2:N-1,n+1).*S(2:N-1,1);
     
     G2(2:N-1,n+1) = p(2:N-1,n+1).*dS(2:N-1,1);
 end
% plot(x,u(:,Nt))
% hold on
    
% PLOTS
% Nt100=fix(Nt/100);
% Nt20=fix(Nt/20);
%  plot(x,u(:,Nt100),'b',x,u(:,Nt20),'g',x,u(:,Nt),'r');
%  xlabel('x')
%  ylabel('Velocity ')
%  title('Velocity VS Distance N = 300 at various times')
%  legend('t = 0.01','t=0.05','t=1')
%  figure;
%  plot(x,rho(:,Nt100),'b',x,rho(:,Nt20),'g',x,rho(:,Nt),'r');
%  xlabel('x')
%  ylabel('rho ')
%  title('Density VS Distance N = 300 at various times')
%  legend('t = 0.01','t=0.05','t=1')
%  figure;
%  plot(x,E(:,Nt100),'b',x,E(:,Nt20),'g',x,E(:,Nt),'r');
%  xlabel('x')
%  ylabel('Energy ')
%  title('Energy VS Distance N = 300 at various times')
%  legend('t = 0.01','t=0.05','t=1')
%  figure;
%  plot(x,p(:,Nt100),'b',x,p(:,Nt20),'g',x,p(:,Nt),'r');
%  xlabel('x')
%  ylabel('Pressure ')
%  title('Pressure VS Distance N = 300 at various times')
%  legend('t = 0.01','t=0.05','t=1')


function mod = minmod(x,y,z)

% THE MINMOD FUNCTION FOR EACH ROW OF THE 3 'M' MATRICES

if( x>0 && y>0 && z>0 )
    mod1 = min(x,y);
    mod = min(mod1,z);

elseif( x<0 && y<0 && z<0 )
    mod1 = max(x,y);
    mod = max(mod1,z);

else
    mod = 0;
end
end