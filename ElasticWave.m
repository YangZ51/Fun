%% Elastic Wave Propagation in a Spherical Ball (Two Layers)

clear;
clf;
clc;

% Define the modeling frame
NT = 500;       % Number of Angular sectors
NR = 200;       % Number of Elements along Radius
nT = linspace(0,2*pi,NT);
nR = linspace(0,1,NR);
[T, R] = meshgrid(nT,nR);
[x,y] = pol2cart(T,R); % Convert grid to cartesian coordintes
% Define the steps
Dt = nT(2) - nT(1);
Dr = 0.005;

% Define the constants
% Pressure wave velocity
Vp1 = 8;
Vp2 = 10;
Vp3 = 12;
% Shear wave velocity
Vs1 = 1;
Vs2 = 3;
Vs3 = 5;
% Physical propertites
rho = 6;
dt = 100*Dr/((Vp1+Vp2+Vp3)/3);
dt1 = (0.25*Dr)/Vp1;
dt2 = (0.25*Dr)/Vp2;
% dts = (0.25*Dr)/Vps;
u1 = Vs1^2*rho;
u2 = Vs2^2*rho;
u3 = Vs3^2*rho;
a = (Vp1^2*rho)-(2*u1);
b = (Vp2^2*rho)-(2*u2);
c = (Vp3^2*rho)-(2*u3);

% Initialization
time = 0;
Ur = zeros(size(R)); Vr_new = Ur; Ur_old = Ur;
Ut = zeros(size(T)); Vt_new = Ut; Vt_old = Ut;
sigmarr = zeros(size(Ur));sigmarr_new = sigmarr; sigmarr_old = sigmarr;
sigmart = zeros(size(T));sigmart_new = sigmart; sigmart_old = sigmart;
sigmatt = zeros(size(T));sigmatt_new = sigmatt; sigmatt_old = sigmatt;

for i = 1:NR
    for j = 1:NT
     sigmarr(i,j)= 0.5*exp(-((i-180).^2+(j-200).^2)/8);
     sigmatt(i,j)= 0.5*exp(-((i-180).^2+(j-200).^2)/8);
    end
end

% Interations (stress & displacement)...
while time <= 1000
    for i = 2:NR-1
        for j = 1:NT
            if j == 1
                k1 = j+1; 
                k2 = NT;
            elseif j == NT
                k1 = 1;
                k2 = j-1;
            else
                k1 = j+1;
                k2 = j-1;
            end
            
            if nR(i) > 0.5 && nR(i) < 1
                sigmarr_new(i,j) = sigmarr(i,j) + (dt1)*(((a+(2*u1))/(2*Dr))*(Ur(i+1,j)-Ur(i-1,j))+(a/(nR(i)*2*Dt))*(Ut(i,k1)-Ut(i,k2))+(a/nR(i))*(Ur(i,j)));
                sigmatt_new(i,j) = sigmatt(i,j) + (dt1)*((a/(2*Dr))*(Ur(i+1,j)-Ur(i-1,j))+((a+2*u1)/(nR(i)*2*Dt))*(Ut(i,k1)-Ut(i,k2))+((a+2*u1)/nR(i))*(Ur(i,j)));
                sigmart_new(i,j) = sigmart(i,j) + (dt1)*((u1/(2*Dr))*(Ut(i+1,j)-Ut(i-1,j))+((u1/(nR(i)*2*Dt))*(Ur(i,k1)-Ur(i,k2)))-(u1/nR(i))*(Ut(i,j)));              
            end
            
            if nR(i) > 0.2 && nR(i) < 0.5
                sigmarr_new(i,j) = sigmarr(i,j) + (dt2)*(((b+(2*u2))/(2*Dr))*(Ur(i+1,j)-Ur(i-1,j))+(b/(nR(i)*2*Dt))*(Ut(i,k1)-Ut(i,k2))+(b/nR(i))*(Ur(i,j)));
                sigmatt_new(i,j) = sigmatt(i,j) + (dt2)*((b/(2*Dr))*(Ur(i+1,j)-Ur(i-1,j))+((b+2*u2)/(nR(i)*2*Dt))*(Ut(i,k1)-Ut(i,k2))+((b+2*u2)/nR(i))*(Ur(i,j)));
                sigmart_new(i,j) = sigmart(i,j) + (dt2)*((u2/(2*Dr))*(Ut(i+1,j)-Ut(i-1,j))+((u2/(nR(i)*2*Dt))*(Ur(i,k1)-Ur(i,k2)))-(u2/nR(i))*(Ut(i,j)));          
 
            elseif nR(i) == 1
                sigmart_new(i,j) = 0.5*(sigmart_old(i,j)-(rho*Vs1*Vt_old(i,j)));
                sigmarr_new(i,j) = 0.5*(sigmarr_old(i,j)-(rho*Vp1*Ur_old(i,j)));
                sigmatt_new(i,j) = sigmatt_old(i,j)+(sigmarr_new(i,j)-sigmarr_old(i,j))*(a/(a*2*u1));
 
            elseif nR(i) == 0.5
                sigmart_new(i,j) = 0.6*(sigmart_old(i,j)-(rho*Vs2*Vt_old(i,j)));
                sigmarr_new(i,j) = 0.6*(sigmarr_old(i,j)-(rho*Vp2*Ur_old(i,j)));
                sigmatt_new(i,j) = sigmatt_old(i,j)+(sigmarr_new(i,j)-sigmarr_old(i,j))*(b/(b*2*u2));
 
            elseif nR(i) == 0.2
                sigmatt_new(i,j) = sigmarr_old(i,j) +(c/(c+2*u3))*sigmarr_old(i,j);
                sigmarr_new(i,j) = 0;
                sigmart_new(i,j) = 0;       
                
            end
        end
    end
    for i = 2:NR-1
        for j = 1:NT
            if j == 1
                k1 = j+1;
                k2 = NT;
            elseif j == NT
                k1=1;
                k2=j-1;
            else
                k1=j+1;
                k2=j-1;
            end
            
            if nR(i) > 0.5 && nR(i) < 1
                Vr_new(i,j) = Ur(i,j) + ((dt1)*((1/(rho*2*Dr))*(sigmarr_new(i+1,j)-sigmarr_new(i-1,j))+(1/(rho*2*Dt*nR(i)))*((sigmart_new(i,k1))-(sigmart_new(i,k2)))+(1/(rho*nR(i)))*(sigmarr_new(i,j)-sigmatt_new(i,j))));
                Vt_new(i,j)= Ut(i,j) + ((dt1)*((1/(rho*2*Dr))*(sigmart_new(i+1,j)-sigmart_new(i-1,j))+(1/(rho*nR(i)*2*Dt))*(sigmatt_new(i,k1)- sigmatt_new(i,k2))+(2/(rho*nR(i)))*(sigmart_new(i,j))));
            end
            
            if nR(i) > 0.2 && nR(i) < 0.5
                Vr_new(i,j) = Ur(i,j) + ((dt2)*((1/(rho*2*Dr))*(sigmarr_new(i+1,j)-sigmarr_new(i-1,j))+(1/(rho*2*Dt*nR(i)))*((sigmart_new(i,k1))-(sigmart_new(i,k2)))+(1/(rho*nR(i)))*(sigmarr_new(i,j)-sigmatt_new(i,j))));
                Vt_new(i,j) = Ut(i,j) + ((dt2)*((1/(rho*2*Dr))*(sigmart_new(i+1,j)-sigmart_new(i-1,j))+(1/(rho*nR(i)*2*Dt))*(sigmatt_new(i,k1)- sigmatt_new(i,k2))+(2/(rho*nR(i)))*(sigmart_new(i,j))));
           
            elseif nR(i) == 1
                Vr_new(i,j) = 0.5*(Ur_old(i,j)-((1/(rho*Vp1))*sigmarr_old(i,j)));
                Vt_new(i,j) = 0.5*(Vt_old(i,j)-((1/(rho*Vs1))*simgart_old(i,j)));
          
            elseif nR(i) == 0.5
                Vr_new(i,j) = 0.6*(Ur_old(i,j)-((1/(rho*Vp2))*sigmarr_old(i,j)));
                Vt_new(i,j) = 0.6*(Vt_old(i,j)-((1/(rho*Vs2))*simgart_old(i,j)));
           
            elseif nR(i) == 0.2
                Vr_new(i,j) = Ur_old(i,j) + (1/(rho*Vp3))*sigmarr_old(i,j);
                Vt_new(i,j) = Vt_old(i,j) + (1/(rho*Vs3))*sigmart_old(i,j);
           
            end
        end
    end    
    
% Output update
time = time+10*dt;
Ur_old = Ur;
Ur = Vr_new;
Vt_old = Ut;
Ut = Vt_new;
sigmarr_old = sigmarr;
sigmarr = sigmarr_new;
sigmatt_old = sigmatt;
sigmatt = sigmatt_new;
sigmart_old = sigmart;
sigmart = sigmart_new;
[Xr,Yt]= pol2cart(Ut,Ur);
z = Xr./max(max(Xr));
figure(1), clf
pcolor(x,y,z);
shading interp
colorbar
hold on
xlabel('X [m]')
ylabel('Y [m]')
title(['Wave propagation after T = ',num2str(time),' sec'])
axis equal, axis tight
drawnow;
end
