% Acoustic Wave Modeling
% (NX,NZ,NT) - input- (Horizontal,Vertical) gridpt dimens. of vel
% model & # Time Steps
% FRE - input- Peak frequency of Ricker wavelet
% BVEL - input- NXxNZ matrix of background velocity model
% (dx,dt) - input- (space, time) sample intervals
% (xs,zs) - input- (x,z) coordinates of line source
% RICKER(NT) - input- NT vector of source time histories
% (p2,p1,p0) -calcul- (future,present,past) NXxNZ matrices of
% modeled pressure field
% (p0,p1) -output- Old and present pressure panels at time NT.
% REALDATA(NX,NT) -output- CSG seismograms at z=2

clear;
clf;
clc;

c=4.0;
% Peek frequency of the seismic source
FRE=20;
% Domain of line source (larger than that of background velocity model
NX=300;NZ=NX; 
% the grid size is dependent on the peek frequency of the source to
% stabilize the calculations
dx=c/FRE/20;
% time step (CFL condition)
dt=.5*dx/c;
%input- (x,z) coordinates of line source
xs=round(NX/2.3); 
zs=round(NX/2);
% Number of iterations in time
NT=900;
% time axis for source history and run time
t=[0:1:NT-1]*dt-0.95/FRE;
% Define the RICKER souce
RICKER=zeros(length(t));
RICKER= (1-t .*t * FRE^2 *pi^2 ) .*exp(- t.^2 * pi^2 * FRE^2 ) ;
% Plot the RICKER source
figure 
plot([0:NT-1]*dt,RICKER);
title('Ricker Wavelet');xlabel('Time (s)')
%input- NXxNZ matrix of background velocity model
BVEL=ones(NX,NZ)*c;
% Define two domains of different wave speed BVEL = 4 and BVEL = 4*1.2
BVEL(NX-round(NX/2):NX,:)= BVEL(NX-round(NX/2):NX,:)*1.2;
% Initialize seismogram
REALDATA=zeros(NX,NT);
% Initialize P_old = p0; P_current = p1; P_future=p2
p0=zeros(NX,NZ);p1=p0;p2=p0;
% constant from dericvation of finite difference scheme
cns=(dt/dx*BVEL).^2;
% Number of grid points in the domain
NX=200;NZ=NX;
figure
% loop over time (number of iteration, it) It is a vectorial calculation
for it=1:1:NT
p2 = 2*p1 - p0 + cns.*del2(p1);
% source added to the new pressure
p2(xs,zs) = p2(xs,zs) + RICKER(it);
% We store the seismogram through time at the source
REALDATA(:,it) = p2(xs,:)';
% update pressure arrays
p0=p1;
p1=p2;
% Plot the results
if mod(it,10)==0
    p00=p0/max(abs(p0(:))+.001);
    imagesc([1:NX]*dx,[1:NX]*dx,(p00+BVEL)); 
    colorbar;
    pause(.1);
end
end
p1=p0;p0=p2;
title('Snapshot of Acoustic Waves')
xlabel('X (km)')
ylabel('Z (km)')