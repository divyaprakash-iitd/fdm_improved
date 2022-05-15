clear; clc; close all;

% Constants
LEFT        = 1;
RIGHT       = 2;
TOP         = 3;
BOTTOM      = 4;

% Domain size
xmax = 15;
ymax = 15;
ar  = xmax/ymax;

% Mesh
nx      = 100;
ny      = nx/ar;
X       = linspace(0,xmax,nx);
Y       = linspace(0,ymax,ny);
[x,y]   = ndgrid(X,Y);

% Define solid(s) [Make sure that the solids do not overlap]
nsolid  = 2;
sloc    = nsolid+1;
dx      = xmax/sloc;
xp      = [];
yp      = [];
xip      = {};
yip      = {};
xg      = zeros(nx,ny,nsolid);
yg      = zeros(nx,ny,nsolid);
beta    = zeros(nx,ny,4,nsolid);

for isolid = 1:nsolid
    np = 60;
    cx = dx*(isolid);
    cy = ymax/2;

    [ixp,iyp]               = closed_curve('hypocycloid',np,2,[cx,cy]);
    [ixip,iyip,xg(:,:,isolid),yg(:,:,isolid),...
        beta(:,:,:,isolid)] = calculate_coefficients(ixp,iyp,x,y);
    
    % Collect the points
    xp = [xp;ixp];
    yp = [yp;iyp];
    xip = [xip;ixip];
    yip = [yip;iyip];
end

beta = sum(beta,4);
beta(beta<=1e-14) = 1;

[solid,liquid,boundary,inner] = generate_flags(xp,yp,x,y);  

% Initial and Boundary Conditions
T           = zeros(size(x));
T(boundary) = 0;

Tb = [600,600,300];
for isolid = 1:size(solid,3)
%     Tb = 500;%*isolid;%abs(500*sin(2*x(solid(:,:,isolid))));
    T(solid(:,:,isolid)) = Tb(isolid);
end


% Solve using finite difference method
for iter = 1:1000
    for i = 1:nx
        for j = 1:ny
            if inner(i,j)
                CC = 2*(1/beta(i,j,LEFT)/beta(i,j,RIGHT) + 1/beta(i,j,TOP)/beta(i,j,BOTTOM));
                CL = 2/beta(i,j,LEFT)*(1/(beta(i,j,LEFT) + beta(i,j,RIGHT)));
                CR = 2/beta(i,j,RIGHT)*(1/(beta(i,j,LEFT) + beta(i,j,RIGHT)));
                CT = 2/beta(i,j,TOP)*(1/(beta(i,j,TOP) + beta(i,j,BOTTOM)));
                CB = 2/beta(i,j,BOTTOM)*(1/(beta(i,j,TOP) + beta(i,j,BOTTOM)));
    
                T(i,j) = 1/CC * (T(i+1,j)*CR + T(i-1,j)*CL + T(i,j-1)*CB + T(i,j+1)*CT);
    
            end
        end
    end
end

% T(logical(sum(solid,3))) = nan;

% Plot figures
hold on
for isolid = 1:nsolid
    contourf(xg(:,:,isolid),yg(:,:,isolid),T,50,'edgecolor','none')
%     plot(xip{isolid},yip{isolid},'ko','MarkerSize',6)
end

plot(xp',yp','k-','LineWidth',2)
fill(xp',yp','k')

colormap("jet")
axis equal
% colorbar
mesh(x,y,ones(size(x)),'facealpha',0,'edgealpha',0.3)
axis off

%% Shapes examples
% [xp,yp] = closed_curve('cardoid',np,2,[cx,cy]);
% [xp,yp] = closed_curve('cardoid',np,2,[cx,cy]);
% [xp,yp] = closed_curve('hypocycloid',np,3,[cx,cy]);
% [xp,yp] = closed_curve('tear',np,3,[cx,cy]);
% [xp,yp] = closed_curve('ellipse',np,3,[cx,cy]);
