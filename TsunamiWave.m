%Tsunami Wave Modeling

clear;
clf;
clc;

% Constants
g  = 9.81;
u0 = 0;                               
v0 = 0;
b  = 0;                               
h0 = 5030;                             

% Define the x domain
ni = 151;                               
xmax = 100000;                      
dx = xmax/(ni-1);
x  = [0:dx:xmax];

% Define the y domain
nj = 151;                              
ymax = 100000;                      
dy = ymax/(nj-1);
y  = [0:dy:ymax];

% Define the wavespeed
wavespeed = u0 + sqrt(g*(h0 - b));

% Define time domain
dt = 0.68*dx/wavespeed;
tmax = 1500;
t=[1:dt:tmax];
courant = wavespeed*dt/dx;

% Create initial u, v, b matrices
u=zeros(length(x), length(y), length(t));
v=zeros(length(x), length(y), length(t));
b=zeros(length(x), length(y));

% Define the depth of the water column (h)
h=zeros(length(x), length(y), length(t)); 
h(:,:,1) = 5000;                              
h((45000/100000*(length(x)-1)+1):floor(55000/100000*(length(x)-1)+1)...
    ,(45000/100000*(length(y)-1)+1):floor(55000/100000*(length(y)-1)+1),1) = 5030;

% Define the bathymetry surface (b)
for i = 1:length(x)
    if x(i) > 20001
        b(:,i) = 0;
    elseif x(i) < 20000
        b(:,i) = 5000/20000*(20000-x(i));
    end
end

% Interations...
for n=1:(length(t)-1)
    for i=2:(ni-1)
        for j=2:(nj-1)    
            u(i,j,n+1) = ((u(i+1,j,n) + u(i-1,j,n) + u(i,j+1,n) + u(i,j-1,n))/4)...
                - 0.5*(dt/dx)*((u(i+1,j,n)^2)/2 - (u(i-1,j,n)^2)/2)...
                - 0.5*(dt/dy)*(v(i,j,n))*(u(i,j+1,n) - u(i,j-1,n)) - 0.5*g*(dt/dx)*(h(i+1,j,n)-h(i-1,j,n));
            
         v(i,j,n+1) = ((v(i+1,j,n) + v(i-1,j,n) + v(i,j+1,n) + v(i,j-1,n))/4)...
                - 0.5*(dt/dy)*((v(i,j+1,n)^2)/2 - (v(i,j+1,n)^2)/2)...
                - 0.5*(dt/dx)*(u(i,j,n))*(v(i+1,j,n) - v(i-1,j,n)) - 0.5*g*(dt/dy)*(h(i,j+1,n)-h(i,j-1,n));
            
            h(i,j,n+1) = ((h(i+1,j,n) + h(i-1,j,n) + h(i,j+1,n) + h(i,j-1,n))/4)...
                - 0.5*(dt/dx)*(u(i,j,n))*((h(i+1,j,n)-b(i+1,j)) - (h(i-1,j,n)-b(i-1,j)))...
                - 0.5*(dt/dy)*(v(i,j,n))*((h(i,j+1,n)-b(i,j+1)) - (h(i,j-1,n)-b(i,j-1)))...
                - 0.5*(dt/dx)*(h(i,j,n)-b(i,j))*(u(i+1,j,n)- u(i-1,j,n))...
                - 0.5*(dt/dy)*(h(i,j,n)-b(i,j))*(v(i,j+1,n) - v(i,j-1,n));

        end
    end

    % Define Boundary Conditions
    u(1,:,n+1) = 2.5*u(2,:,n+1) - 2*u(3,:,n+1) + 0.5*u(4,:,n+1);
    u(length(x),:,n+1) = 2.5*u(ni-1,:,n+1) - 2*u(ni-2,:,n+1) + 0.5*u(ni-3,:,n+1);
    u(:,1,n+1) = 2.5*u(:,2,n+1) - 2*u(:,3,n+1) + 0.5*u(:,4,n+1);
    u(:,length(y),n+1) = 2.5*u(:,nj-1,n+1) - 2*u(:,nj-2,n+1) + 0.5*u(:,nj-3,n+1);

    v(1,:,n+1) = 2.5*v(2,:,n+1) - 2*v(3,:,n+1) + 0.5*v(4,:,n+1);
    v(length(x),:,n+1) = 2.5*v(ni-1,:,n+1) - 2*v(ni-2,:,n+1) + 0.5*v(ni-3,:,n+1);
    v(:,1,n+1) = 2.5*v(:,2,n+1) - 2*v(:,3,n+1) + 0.5*v(:,4,n+1);
    v(:,length(y),n+1) = 2.5*v(:,nj-1,n+1) - 2*v(:,nj-2,n+1) + 0.5*v(:,nj-3,n+1);

    h(1,:,n+1) = 2.5*h(2,:,n+1) - 2*h(3,:,n+1) + 0.5*h(4,:,n+1);
    h(length(x),:,n+1) = 2.5*h(ni-1,:,n+1) - 2*h(ni-2,:,n+1) + 0.5*h(ni-3,:,n+1);
    h(:,1,n+1) = 2.5*h(:,2,n+1) - 2*h(:,3,n+1) + 0.5*h(:,4,n+1);
    h(:,length(y),n+1) = 2.5*h(:,nj-1,n+1) - 2*h(:,nj-2,n+1) + 0.5*h(:,nj-3,n+1);

end

figure(1)

% Animation
for index=1:length(t)
    mesh(x,y,h(:,:,index))
    axis ([0 100000 0 100000 4990 5010])
    title ('AERSP 423 Computer Project Part II')
    xlabel('X Domain [m]')
    ylabel('Y Domain [m]') 
    zlabel('Height [m]')
    pause(0.02)
end
