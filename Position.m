%%This code is for course AAE6102. 
%%Name: WANG Qipeng         Student ID:20031303R
clear;
load('rcvr.dat')
load('eph.dat')
%%% Two data sets should be matched
rcvr([1,2,3,4,5,6,7,8],:)=rcvr([8,1,2,3,4,5,6,7],:);
%%Step 1
c=299792458;%speed of light
Wedot=7.2921151467e-5;%rotation rate
mu=3.986005e14;%grativitation constant
F=-4.442807633e-10;%correction term constant
A=eph(:,10).^2;%Semi-major axis
n0=(mu./A.^3).^0.5;%mean motion
t=rcvr(:,1)-rcvr(:,3)/c;%calculate the GPS system time at time of transmition
tk=t-eph(:,4);   %Time from ephemeris reference epoch
n=n0+eph(:,11);%motion correction
Mk=eph(:,12)+n.*tk; %mean anomaly
%%%solve Kepler's equation
Ek=Mk;
for i =1:5 %%Actually, 2 iteration is enough
    Ek=Ek+(Mk-(Ek-eph(:,9).*sin(Ek)));
end
vk=atan2((1-eph(:,9).*eph(:,9)).^0.5.*sin(Ek),cos(Ek)-eph(:,9));
phi_k=vk+eph(:,13);
delta_mu_k=eph(:,18).*sin(2*phi_k)+eph(:,19).*cos(2*phi_k);
delta_r_k=eph(:,22).*sin(2*phi_k)+eph(:,23).*cos(2*phi_k);
delta_i_k=eph(:,20).*sin(2*phi_k)+eph(:,21).*cos(2*phi_k);
mu_k=phi_k+delta_mu_k;
rk=A.*(1-eph(:,9).*cos(Ek))+delta_r_k;
ik=eph(:,15)+delta_i_k+eph(:,17).*tk;
%%Positions in orbit plane
xk_p=rk.*cos(mu_k);
yk_p=rk.*sin(mu_k);
omega_k=eph(:,14)+(eph(:,16)-Wedot).*tk-Wedot.*eph(:,4);
x_k=xk_p.*cos(omega_k)-yk_p.*cos(ik).*sin(omega_k);
y_k=xk_p.*sin(omega_k)+yk_p.*cos(ik).*cos(omega_k);
z_k=yk_p.*sin(ik);
%%%%%Step 2: Determine the broadcast satellite clock error
Delta_t_r=F.*eph(:,9).*eph(:,10).*sin(Ek);
Delta_t_sv=eph(:,5)+eph(:,6).*(t-eph(:,3))+eph(:,7).*(t-eph(:,3)).*(t-eph(:,3))+Delta_t_r;
rho_corr=rcvr(:,3)+c.*Delta_t_sv;
%%%%%Step 3-6:
%Initialize algorithm
pos=[-2694685.473 -4293642.366 3857878.924 0]';
ini_bias=0;
Delta_x=1;
H=zeros(length(x_k),4);
while abs(Delta_x(1))>=1e-4
    %Linear GPS measurement equation
    Rho=sqrt((x_k-pos(1)).^2+(y_k-pos(2)).^2+(z_k-pos(3)).^2);
    a_x=(x_k-pos(1))./Rho;
    a_y=(y_k-pos(2))./Rho;
    a_z=(z_k-pos(3))./Rho;
    Delta_rho=Rho-rho_corr;
    H(:,1)=a_x;
    H(:,2)=a_y;
    H(:,3)=a_z;
    H(:,4)=-1*ones(length(x_k),1);
    Delta_x=pinv(H'*H)*H'*Delta_rho;
    pos=pos+Delta_x;
end
disp('Satellites Positions are:')
x_k
y_k
z_k
disp('User Position is:')
pos([1:3])