function timeEvolution2D


%     freeParticlePipeline
%      doubleSlitPipeline;
    tunnellingPipeline

% ------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATIONS EXPERIMENTS  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freeParticlePipeline()
% FREEPARTICLEPIPELINE

    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 1;
    
    % Momentum
    px = 3.5;
    py = 0;
    
    % Gaussian wave packet width
    sigmax = 3.5;
    sigmay = 3.5;
    
    % Gaussian wave packet starting position
    xo = -40;
    yo = 0;
    
    % Space lattice size
    L = 50;
    N = 250;
    
    % Time lattice size
    nt = 200;
    dt = 0.05;

    % Space-time discretization init
    t =timeGridInit(nt);
    [x,y,dx,dy] = spaceGridInit(L,N);
    [k2,~,~] = momentumGridInit(N,dx,dy);

    % Gaussian wave packet
    psi = psiInit(x,y,xo,yo,px,py,sigmax,sigmay);

    % Potential grid
    U = potentialInit(N);

    % Surf Plot
    axis = [-L , L , -10 , 10 , 0 , 0.1];
    surfPlot = plotInit(x,y,psi,axis);
    view(0,90)

    % Run simulation
    core(psi,U,k2,t,dt,nplot,surfPlot,recordMovie,writerObj)

% ------------------------------------

function doubleSlitPipeline()
% DOUBLESLITPIPELINE
    type = 'doubleSlit';

    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 3;
    
    % Momentum
    px = 3.5;
    py = 0;
    
    % Gaussian wave packet width
    sigmax = 1;
    sigmay = 1;
    
    % Gaussian wave packet starting position
    xo = -10;
    yo = 0;
    
    % Space lattice size
    L = 50;
    N = 301;
    
    % Time lattice size
    nt = 1000;
    dt = 0.01;

    % Space-time discretization init
    t =timeGridInit(nt);
    [x,y,dx,dy] = spaceGridInit(L,N);
    [k2,~,~] = momentumGridInit(N,dx,dy);
    
    yo = 0.5*dy;
    
    % Gaussian wave packet
    psi = psiInit(x,y,xo,yo,px,py,sigmax,sigmay);
    
    % Potential grid
    U = potentialInit(N,type);
    
    % Surf Plot
    axis = [-L+10 , L , -25 , 25 , 0 , 0.5];
    surfPlot = plotInit(x,y,psi,axis);
    view(0,90)

    % Run simulation
    core(psi,U,k2,t,dt,nplot,surfPlot,recordMovie,writerObj);
% ------------------------------------

function tunnellingPipeline()
% TUNNELLINGPIPELINE
    type = 'tunnelling';
    
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 1;
    
    % Momentum
    px = 3;
    py = 0;
 
    % Gaussian wave packet width
    sigmax = 4;
    sigmay = 4;
  
    % Gaussian wave packet starting position
    xo = -10;
    yo = 0;

    % Space lattice size
    L = 50;
    N = 250;
    
    % Time lattice size
    nt = 200;
    dt = 0.05;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,y,dx,dy] = spaceGridInit(L,N);
    [k2,~,~] = momentumGridInit(N,dx,dy);

    % Gaussian wave packet
    psi = psiInit(x,y,xo,yo,px,py,sigmax,sigmay);

    % Potential grid
    U = potentialInit(N,type);

    % Surf Plot
    axis = [-L+10 , L , -10 , 10 , 0 , 0.1];
    surfPlot = plotInit(x,y,psi,axis);
    view(0,90)

    % Run simulation
    core(psi,U,k2,t,dt,nplot,surfPlot,recordMovie,writerObj)

% ------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATION INIT UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function core(psi,U,k2,t,dt,nplot,h,flag,writer)
% CORE run simulation
    pause(2)
    upsi = psi;
    for j = t
       % 1. 
       upsi = exp(-1i*dt*U/2).*upsi;

       % 2. Compute the Fourier transform of upsi to go to the momentum space
       upsiN1 = fft2(upsi);

       % 3. Multiply upsiN1 by the kinetic part of the evolution equation
       upsiN2 = exp(-1i*dt*k2/2).*upsiN1;

       % 4. Go back to x space with the inverse Fourier transform
       upsi2 = ifft2(upsiN2);

       % 5. Finally
       upsi = exp(-1i*dt*U/2).*upsi2;

       if mod(j,nplot) == 0
           set(h,'Zdata',abs(upsi).^2);
           captureFrame(flag,writer,getframe(gcf))
%            pause(0.001)
           drawnow
        end
    end
    pause(3)
    close(gcf)
% ------------------------------------

function [t]=timeGridInit(nt)
% TIMEGRIDINIT
    t = 1:nt;
% ------------------------------------

function [x,y,dx,dy] = spaceGridInit(l,n)
% SPACEGRIDINIT
    [x,y] = meshgrid(linspace(-l,l,n),linspace(-l,l,n));
    dx = x(1,2) - x(1,1);
    dy = y(2,1) - y(1,1);
% ------------------------------------

function [k2,kx,ky] = momentumGridInit(N,dx,dy)
% MOMENTUMGRIDINIT
    n = floor(N/2);
    nn = floor((N-1)/2);
    kx = (2*pi/(dx*N))*(-n:nn)';
    ky = (2*pi/(dy*N))*(-n:nn)';
    kx = fftshift(kx);
    ky = fftshift(ky);
    [kx , ky] = meshgrid(kx,ky);
    k2 = kx.^2+ky.^2;
% ------------------------------------

function psi = psiInit(x,y,xo,yo,px,py,sigmax,sigmay)
% PSIINIT
% Wawe packet init
    psi = exp(-(x-xo).^2/2/sigmax -(y-yo).^2/2/sigmay +1i*px*x+1i*py*y )/(pi*sigmax*pi*sigmay).^(1/4);
% ------------------------------------

function U = potentialInit(n,type)
% POTENTIALINIT
    if(nargin < 2)
        type = 'freeParticle';
    end
    switch type
        case 'freeParticle'
            U = zeros(n,n);
        case 'doubleSlit'
            barrierHeight = 800;
            barrierWidth  = 5;
            slitsDistance = 4;
            slitsWidth    = 8;
            slitsDistance = slitsDistance/2;
            U = zeros(n,n);
            l = (n-1)/2 + 1;
            ll = (n -1)/2 + 1 + barrierWidth;
            U(1:n,l:ll) = barrierHeight;
            U(l-slitsDistance-slitsWidth+1:l-slitsDistance+1,l:ll) = 0;
            U(l+slitsDistance+1:l+slitsWidth+slitsDistance+1,l:ll) = 0;
        case 'tunnelling'
            barrierHeight = 5;
            barrierWidth  = 5;
            U = zeros(n,n);
            l = n/2;
            ll = n/2 + barrierWidth;
            U(1:n,l:ll) = barrierHeight;
    end
% ------------------------------------

function [surfPlot] = plotInit(x,y,psi,ax)
% PLOTINIT
   clf
   shg
%    set(gcf,'menubar','none');
   set(gcf,'Color',[0.6 0.6 0.6]);
   surfPlot = surf(x,y,abs(psi).^2);
   surfPlot.FaceColor = 'interp';
   surfPlot.EdgeColor = 'none';
   if(nargin == 4)
       axis([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6)])
       box on
   end
   grid
% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   VIDEO RECORDING UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function writerObj = initVidRecording(flag,path,frameRate)
% INITVIDRECORDING video recording 
% properties initialization.
    if(flag)
        if(nargin < 2)
            path = 'out.avi';
        end
        if(nargin < 3)
            frameRate = 25;
        end
        writerObj = VideoWriter(path); % Name it.
        writerObj.FrameRate = frameRate; % How many frames per second.
        open(writerObj);
    else
        writerObj = false;
    end
% ------------------------------------

function captureFrame(flag,writer,frame)
% CAPTUREFRAME writeVideo wrapper:
% Write the image frame to the video file
    if(flag)
        writeVideo(writer, frame)
    end
% ------------------------------------
