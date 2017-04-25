function timeEvolution1D

%     freeParticlePipeline
%    harmonicOscillatorPipeline
%     potentialStepPipeline
%     potentialBarrierPipeline
    potentialWellPipeline

% ------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         EXPERIMENTS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freeParticlePipeline()
% FREEPARTICLEPIPELINE
    type = 'freeParticle';
    
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 5;
    
    % Momentum
    px = 4;
    
    % Gaussian wave packet width
    sigmax = 4;
    
    % Gaussian wave packet starting position
    xo = -40;
    
    % Space lattice size
    L = 50;
    N = 1024;
    
    % Time lattice size
    nt = 1300;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(N,x,L,type);

    % Plot
    axis = [-L L -0.6 0.60];
    plot = plotInit(x,psi,U,axis);

    % Run simulation
    core(psi,U,k2,t,dt,nplot,plot,recordMovie,writerObj)
% ------------------------------------

function harmonicOscillatorPipeline()
% HARMONICOSCILLATORPIPELINE
    type = 'harmonicOscillator';
    
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 5;
    
    % Momentum
    px = 0.1;
    
    % Gaussian wave packet width
    sigmax = 1;
    
    % Gaussian wave packet starting position
    xo = -6;
    
    % Space lattice size
    L = 10;
    N = 1024;
    
    % Time lattice size
    nt = 2000;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(N,x,L,type);

    % Surf Plot
    axis = [-L L -0.6 0.60];
    plot = plotInit(x,psi,U,axis);
    legend('V = Vo*0.5*x^2','real(\Psi)','imag(\Psi)')

    % Run simulation
    core(psi,U,k2,t,dt,nplot,plot,recordMovie,writerObj)
% ------------------------------------

function potentialStepPipeline()
% POTENTIALSTEPPIPELINE
    type = 'potentialStep';
    
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 5;
    
    % Momentum
    px = 3;
    
    % Gaussian wave packet width
    sigmax = 4;
    
    % Gaussian wave packet starting position
    xo = -10;
    
    % Space lattice size
    L = 50;
    N = 1024;
    
    % Time lattice size
    nt = 1100;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(N,x,L,type);

    % Plot
    axis = [-L+10 L-10 -0.5 0.50];
    plot = plotInit(x,psi,U,axis);
    legend('V = Vo*heaviside(x)','real(\Psi)','imag(\Psi)')

    % Run simulation
    core(psi,U,k2,t,dt,nplot,plot,recordMovie,writerObj)
% ------------------------------------

function potentialBarrierPipeline()
% POTENTIALBARRIEPIPELINE
    type = 'potentialBarrier';
    
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 5;
    
    % Momentum
    px = 3;
    
    % Gaussian wave packet width
    sigmax = 4;
    
    % Gaussian wave packet starting position
    xo = -10;
    
    % Space lattice size
    L = 50;
    N = 1024;
    
    % Time lattice size
    nt = 1000;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(N,x,L,type);

    % Plot
    axis = [-L+10 L-10 -0.5 0.50];
    plot = plotInit(x,psi,U,axis);
    legend('V = vo*(heaviside(x+w)-heaviside(x-w))','real(\Psi)','imag(\Psi)')

    % Run simulation
    core(psi,U,k2,t,dt,nplot,plot,recordMovie,writerObj)
% ------------------------------------

function potentialWellPipeline()
% POTENTIALWELLPIPELINE
    type = 'potentialWell';
    
    % Plot & Movie
    recordMovie = true;
    writerObj = initVidRecording(recordMovie);
    nplot = 5;
    
    % Momentum
    px = 3;
    
    % Gaussian wave packet width
    sigmax = 4;
    
    % Gaussian wave packet starting position
    xo = -13;
    
    % Space lattice size
    L = 50;
    N = 1024;
    
    % Time lattice size
    nt = 1000;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(N,x,L,type);

    % Plot
    axis = [-L+10 L-10 -0.5 0.50];
    plot = plotInit(x,psi,U,axis);
    legend('V = -vo*(heaviside(x+w)-heaviside(x-w))','real(\Psi)','imag(\Psi)')

    % Run simulation
    core(psi,U,k2,t,dt,nplot,plot,recordMovie,writerObj)
% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATION INIT UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function core(psi,U,k2,t,dt,nplot,plt,flag,writer)
% CORE run simulation
    pause(2)
    upsi = psi;
    for j = t
       % 1. 
       upsi = exp(-1i*dt*U/2).*upsi;

       % 2. Compute the Fourier transform of upsi to go to the momentum space
       upsiN1 = fft(upsi);

       % 3. Multiply upsiN1 by the kinetic part of the evolution equation
       upsiN2 = exp(-1i*dt*k2/2).*upsiN1;

       % 4. Go back to x space with the inverse Fourier transform
       upsi2 = ifft(upsiN2);

       % 5. Finally
       upsi = exp(-1i*dt*U/2).*upsi2;
       
       if mod(j,nplot) == 0
           set(plt(1),'Ydata',abs(upsi).^2);
           set(plt(2),'Ydata',real(upsi));
           set(plt(3),'Ydata',imag(upsi));
           captureFrame(flag,writer,getframe(gcf))
           drawnow
        end
    end
% ------------------------------------

function [t]=timeGridInit(nt)
% TIMEGRIDINIT
    t = 1:nt;
% ------------------------------------

function [x,dx] = spaceGridInit(l,n)
% SPACEGRIDINIT
    x = linspace(-l,l,n);
    dx = x(2) - x(1);
% ------------------------------------

function [k2,kx] = momentumGridInit(N,dx)
% MOMENTUMGRIDINIT
    n = floor(N/2);
    nn = floor((N-1)/2);
    kx = (2*pi/(dx*N))*(-n:nn)';
    kx = fftshift(kx);
    k2 = (kx.^2)';
% ------------------------------------

function psi = psiInit(x,xo,px,sigmax)
% PSIINIT
% Wawe packet init
    psi = exp(-(x-xo).^2/2/sigmax+1i*px*x )/(2*pi*sigmax).^(1/4);
% ------------------------------------

function U = potentialInit(n,x,l,type)
% POTENTIALINIT
    if(nargin < 2)
        type = 'freeParticle';
    end
    switch type
        case 'freeParticle'
            U = zeros(1,n);
        case 'harmonicOscillator'
            Uo = 0.4;
            U = Uo*x.^2/2;
        case 'potentialStep'
            Uo = 5.4;
            U = Uo*heaviside(x);
        case 'potentialBarrier'
            Uo = 5.4;
            w = 1;
            U = Uo*(heaviside(x+w)-heaviside(x-w));
        case 'potentialWell'
            Uo = -5.4;
            w = 2;
            U = Uo*(heaviside(x+w)-heaviside(x-w));

    end
% ------------------------------------

function [plt] = plotInit(x,psi,U,ax)
% PLOTINIT
   clf
   shg
%    set(gcf,'menubar','none');
   set(gcf,'Color',[0.6 0.6 0.6]);
   psiPlot = abs(psi).^2;
   Usc = U*max(psiPlot)/max(abs(U));
   plot(x,Usc,'--k');
   hold on
   box on
   hh = plot(x,real(psi),'--y');
   hhh = plot(x ,imag(psi),'--b');
   h = fill(x,psiPlot,'r');
   grid
   ylabel('|\Psi|^2','FontSize' ,22)
   xlabel('X','FontSize' ,20)
   if(nargin == 4)
       axis([ax(1) ax(2) ax(3) ax(4)])
       box on
   end
   grid
   set( gca, 'Color', [0.5,0.5,0.5] );
   plt = [h hh hhh];
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