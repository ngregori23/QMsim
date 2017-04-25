function imaginaryTimePropagation
% IMAGINARYTIMEPROPAGATION
    harmonicOscillatorGSPipeline
% ------------------------------------

function harmonicOscillatorGSPipeline()
% HARMONICOSCILLATORPIPELINE
    % Plot & Movie
    recordMovie = false;
    writerObj = initVidRecording(recordMovie);
    nplot = 1;
    
    % Momentum
    px = 22;
    
    % Gaussian wave packet width
    sigmax = 1;
    
    % Gaussian wave packet starting position
    xo = 0;
    
    % Space lattice size
    L = 10;
    N = 1024;
    
    % Time lattice size
    nt = 600;
    dt = 0.01;

    % Space-time discretization init
    t = timeGridInit(nt);
    [x,dx] = spaceGridInit(L,N);
    [k2,~] = momentumGridInit(N,dx);

    % Gaussian wave packet
    psi = psiInit(x,xo,px,sigmax);

    % Potential grid
    U = potentialInit(x);

    % Surf Plot
    plot = plotInit(x,psi,U,t);
    
    % Run simulation
    core(psi,dx,U,k2,t,dt,nplot,plot,recordMovie,writerObj);
% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATION INIT UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function upsi = core(upsi,dx,U,k2,t,dt,nplot,plt,flag,writer)
% CORE run simulation
    pause(2)
    for j = t
       % 1. 
       upsi = exp(-dt*U/2).*upsi;

       % 2. Compute the Fourier transform of upsi to go to the momentum space
       phiN1 = fft(upsi);

       % 3. Multiply upsiN1 by the kinetic part of the evolution equation
       phiN2 = exp(-dt*k2/2).*phiN1;

       % 4. Go back to x space with the inverse Fourier transform
       upsi2 = ifft(phiN2);

       % 5. Finally
       upsi = exp(-dt*U/2).*upsi2;
       if mod(j,nplot) == 0
           Eo = totalEnergy(upsi,k2,U)
           upsiPlot=upsi/sqrt(sum(upsi.*conj(upsi))*dx);
           set(plt(1),'Ydata',real(upsiPlot));
           set(plt(2),'Ydata',imag(upsiPlot));
           addpoints(plt(3),j,Eo)
           captureFrame(flag,writer,getframe(gcf))
           str = sprintf('Eo = %8.1f',Eo);
           set(plt(4),'string',str);
           drawnow
           pause(0.01)
        end
    end
    pause(3)
% ------------------------------------

function E = totalEnergy(psi,k2,U,mass)
% TOTALENERGY
    if(nargin < 4)
        mass = 1;
    end
    E = kEnergy(psi,k2,mass) + pEnergy(psi,U);
% ------------------------------------

function Ek = kEnergy(psi,k2,mass)
% KENERGY
    if(nargin < 3)
        mass = 1;
    end
    phi=fft(psi);
    Ek=sum(k2.*abs(phi).^2/(2*mass))/sum(abs(phi).^2);
% ------------------------------------

function Ep = pEnergy(psi,U)
% PENERGY
    Ep=sum(U.*abs(psi).^2)/sum(abs(psi).^2);
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
    psi = exp(-(x-xo).^2/2/sigmax + 1i*px*x )/(2*pi*sigmax).^(1/4);
% ------------------------------------

function U = potentialInit(x)
% POTENTIALINIT
	Uo = 1;
	U = Uo*x.^2/2;
% ------------------------------------

function [plt] = plotInit(x,psi,U,t)
% PLOTINIT
   clf
   shg
   set(gcf,'menubar','none');
   set(gcf,'Color',[0.6 0.6 0.6]);
   psiPlot = abs(psi).^2;
   Usc = U*max(psiPlot)/max(abs(U));
   subplot(1,2,1)
   plot(x,Usc,'--k');
   hold on
   box on
   grid
   realP = plot(x,real(psi),'--y');
   imagP = plot(x ,imag(psi),'--b');
   legend('V = Vo*0.5*x^2','real(\Psi)','imag(\Psi)')
   ylabel('|\Psi|^2','FontSize' ,22)
   xlabel('X','FontSize' ,20)
   set( gca, 'Color', [0.5,0.5,0.5] );
   axis([-10,10,-0.8,0.8]);
   
   subplot(1,2,2)
   Eo = ones(size(t))*0.5;
   plot(t,Eo,'y','LineWidth',1);
   set( gca, 'Color', [0.5,0.5,0.5] );
   xlabel('Imaginary time','FontSize' ,20)
   ylabel('Ground state Energy','FontSize' ,22)
   legend('Theoretical Eo : 0.5')
   energyP = animatedline('LineWidth',2,'Color', [ 0.9, 0.9, 0.9 ] );
   titl2 = title(sprintf('Eo = %8.1f',0));
   ylim([0.49,0.51])
   
   plt = [realP imagP energyP titl2];
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