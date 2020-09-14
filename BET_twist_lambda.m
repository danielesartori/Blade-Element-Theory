% Author: Daniele Sartori
% Date: October 2018
% References:
% - Hoffmann, Huang - Precision Flight Control for A Multi-Vehicle Quadrotor Helicopter Testbed
% - Hoffmann, Huang - Quadrotor Helicopter Flight Dynamics and Control, Theory and Experiment
% - Leishman - Principle of Helicopter Aerodynamics

clc
clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input vectors - Size = N x 1

%Rotor velocity vector module [m/s]
v=[2; 2.1; 2.2; 2.3; 2.4; 2.5];

%Rotor angle of attack [rad]
alpha=[0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

%Rotor angular speed [rpm]
Omega_rpm=[6200; 6220; 6240; 6260; 6280; 6300];

%Time [s]
t=[0; 0.1; 0.2; 0.3; 0.4; 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation properties

%Number of rotor radial elements [-]
Nr=100;

%Number of rotor angular elements [-]
Npsi=360;

%Iterative loop tolerance [-]
tol=1e-4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Propeller DJI 9450 properties

%Number of blades [-]
N_b=2;

%Rotor blade radius [m]
R_b=0.12;

%Rotor blade hub radius [m]
R_b_hub=0.007;

%Rotor blade average chord [m]
c_b=0.02074;

%Chord vector [m]
%Rotor blade poly fitting chord [m]
c_poly=[1.108365196825e-14	-6.286058435905e-12 1.483520999427e-09 -1.875804437693e-07 1.357510562304e-05 -0.0005515931891319 0.01126966744219 -0.09409957611780 0.9076873044609 8.571483392828];
r_v_mm=linspace(R_b_hub,R_b,Nr)*1000;
c_v=polyval(c_poly,r_v_mm)/1000;

%Blade pitch from thrust test experimental data fit [rad]
theta0=0.324;

%Theta twist [rad]
theta_tw=0.111; 

%Adimensional radius for reference theta
r_theta=0.25;

%Effective rotor hinge offset vector [0-1]
ef=0.65;

%Blade stiffness [N*m/rad]
k_b=0.2570;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ambient conditions

%Air density [kg/m3]
rho_air=1.225;

%Air temperature [°C]
Temp_air=15;

%Temperature at 0 deg C [K]
T0=273.15;

%Air dynamic viscosity (Sutherland formula) [kg/m·s]
k1=291.15;
k2=120;
visc_dyn=1.827e-5*(k1+k2)/(Temp_air+T0+k2)*((Temp_air+T0)/k1)^1.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis domain vectors

%Define rotor adimensional radius elements vector [-]
r_v=linspace(R_b_hub/R_b,1,Nr);
dr=(1-R_b_hub/R_b)/Nr;

%Define rotor angular position vector [rad]
psi=linspace(0,2*pi,Npsi);
dpsi=2*pi/Npsi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotor properties

%Rotor blade hub-tip radius [m]
R=R_b-R_b_hub;

%Rotor area [m^2]
Ar=pi*(R_b^2-R_b_hub^2);

%Blade aspect ratio [-]
AR=2*R_b/c_b;

%Blade Oswald factor (from Boehnke/Raymer) [-]
ew=1.78*(1-0.045*AR^0.68)-0.64;

%Induced drag constant [-]
Kar=1/(pi*ew*AR);

%Theoretical blade lift curve slope [1/rad]
Clalpha_b_theo=2*pi;

%Blade lift curve slope [1/rad]
Clalpha_b=Clalpha_b_theo/(1+Clalpha_b_theo*Kar);

%Rotor solidity [-]
sigma=N_b*c_b*R/Ar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forward flight conditions - General properties

%Conversion rotor rotational speed [rad/s]
Omega=2*pi*Omega_rpm/60;

%Rotor tip velocity [m/s]
v_tip=Omega*R_b;

%Rotor longitudinal advance ratio [-]
mu=v.*cos(alpha)./v_tip;

%Rotor vertical advance ratio [-]
lambda_c=-v.*sin(alpha)./v_tip;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotor flapping

%Effective moment of inertia of the blade about the hinge at ef [kg*m^2]
I_b_hinge=4/3*k_b./(N_b*Omega.^2*ef);

%Ratio flapping frequency/angular rate of rotor (also in Leishman pag 239 PDF) [-]
lambda_b=sqrt(1+1.5*ef+k_b./I_b_hinge./Omega.^2);

%Lock number (also in Leishman pag 217 PDF) [-]
lock=rho_air*Clalpha_b*c_b*(R_b^4-R_b_hub^4)./I_b_hinge;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forward flight conditions - Blade Element Theory for Rotating Blade (BETRB)

%Initialization
lambda_temp=0.05*ones(size(v,1),Nr);
lambda=0.05*ones(size(v,1),Nr);
dT=zeros(size(v,1),Npsi,Nr);
dHf=zeros(size(v,1),Npsi,Nr);
dYf=zeros(size(v,1),Npsi,Nr);
beta_dot=zeros(size(v,1),1);
err_check=1;

while err_check>tol
    
    
    %Evaluation of flapping angle    
    x=zeros(size(v,1),3);    
    theta_flap=mean(theta0+theta_tw*(r_v-r_theta));
    
    for i=1:size(v,1)    
        
        lambda_flap=mean(lambda(i,:),2);
        
        A=[lambda_b(i)^2       0                0;
            lock(i)/6*mu(i)        1-lambda_b(i)^2    -lock(i)/8;
            0                lock(i)/8           1-lambda_b(i)^2];
        
        B=[lock(i)/2*(theta_flap/4-lambda_flap/3); 0; lock(i)/3*mu(i)*theta_flap];

        x(i,:)=A\B;
        
    end
    
    %Flapping angles [rad]
    beta_con=x(:,1);
    beta_lon=x(:,2);
    beta_lat=x(:,3);
    
    
    for i=1:Npsi
        for j=1:Nr
            
            %Pitch angle [rad]
            theta=theta0+theta_tw*(r_v(j)-r_theta);
            
            %Inflow (Leishman pag 165 PDF) [-]         
            lambda(:,j)=sqrt( (sigma*Clalpha_b/16-lambda_c/2).^2 + sigma*Clalpha_b*theta*r_v(j)/8) -(sigma*Clalpha_b/16-lambda_c/2);
            
            
            %Flapping angle (sign change for beta_lon compared to Leishman, following Hoffmann convention)
            beta=beta_con-beta_lon*cos(psi(i))+beta_lat*sin(psi(i));
            
            %Flapping angle derivative [rad/s]
            beta_dot(1:end-1)=diff(beta)./diff(t);
            
            %Leishman pag 196 PDF
            ut_omegaR=r_v(j)+mu*sin(psi(i));
            up_omegaR=lambda(:,j)+mu.*sin(beta)*cos(psi(i))+r_v(j)*beta_dot./Omega;
            phi=atan(up_omegaR./ut_omegaR);
            U=sqrt( (up_omegaR*R_b.*Omega).^2 + (ut_omegaR*R_b.*Omega).^2 );
            
            %Reynolds number
            Re=rho_air*U*c_v(j)/visc_dyn;
            
            %Blade drag coefficient (Leishman pag 337 PDF) [-]
            Cd0_b=0.1166*Re.^-0.2;
            
            %Drag coefficient [-]
            Cd_b=Cd0_b+Kar*(Clalpha_b*(theta-phi)).^2;
            
            %Rotor tip loss coefficient (Leishman pag 114 PDF) [-]
            Bb=1-(1.386/N_b)*lambda(:,j);
            
            %Thrust
            dT(:,i,j)=N_b*c_v(j)*R*rho_air/(4*pi)*Bb.*U.^2.*(Clalpha_b*(theta-phi).*cos(phi)-Cd_b.*sin(phi)).*cos(beta)*dr*dpsi;
            
            %H-force
            dHf1=(Clalpha_b*(theta-phi).*sin(phi)+Cd_b.*cos(phi))*sin(psi(i));
            dHf2=(-Clalpha_b*(theta-phi).*cos(phi)+Cd_b.*sin(phi)).*cos(psi(i)).*sin(beta);
            dHf(:,i,j)=N_b*c_v(j)*R*rho_air/(4*pi)*Bb.*U.^2.*(dHf1+dHf2)*dr*dpsi;
            
            %Y-force
            dYf1=-(Clalpha_b*(theta-phi).*sin(phi)+Cd_b.*cos(phi))*cos(psi(i));
            dYf2=(-Clalpha_b*(theta-phi).*cos(phi)+Cd_b.*sin(phi)).*sin(psi(i)).*sin(beta);
            dYf(:,i,j)=N_b*c_v(j)*R*rho_air/(4*pi)*Bb.*U.^2.*(dYf1+dYf2)*dr*dpsi;
            
        end
    end
    
    
    %Loop convergence error
    err_check=max(abs(lambda-lambda_temp));
    
    %Assign found values to iteration parameters
    lambda_temp=lambda;
    
end


%Rotor thrust Blade Element Theory Rotating Blade [N]
T_BET=sum(sum(dT,3),2);

%Rotor H-force [N]
Hf=sum(sum(dHf,3),2);

%Rotor Y-force [N]
Yf=sum(sum(dYf,3),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure

figure(1)
plot(t,T_BET,'-o')
grid on
hold on
xlabel('t [s]')
ylabel('T [N]')

figure(2)
plot(t,Hf,'-o')
grid on
hold on
xlabel('t [s]')
ylabel('H [N]')

figure(3)
plot(t,Yf,'-o')
grid on
hold on
xlabel('t [s]')
ylabel('Y [N]')

