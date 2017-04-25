function schrodinger1D
% SCHRODINGER1D
    particleInABoxPipeline
%     figure
%     harmonicOscillatorPipeline

% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         EXPERIMENTS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function particleInABoxPipeline
% PAETICLEINABOX
% Solve eigenvalue problem for 1D Schrodinger equation.
% In this experiment U(x) is a Infinite square well.
% 
% the technique is extendable to arbitrary potentials.
% 
% Then the numerical solution is compared with the analitical solution for the
% infinite square well problem.
    type = 'particleInABox';

    nEig = 100;
    nEig2Plot = 3;
    
    % Problem properties
    hbar = 1;
    mass = 1;
    
    % Space lattice size & discratization 
    L = 1;
    N = 1*1000;
    x = linspace(-L,L,N);
    
     % Hemiltonian Matrix
    [H,~,U] = hemiltonianMatrixInit(type,x,N,L,hbar,mass);
    
    % Eigenvalue & eigenstates solver
    opts = eigsOptionsInit();
    [V,En] = eigs(H, nEig ,'SA',opts);
    
    % Plot U & first nEig2Plot eigensolutions
    plotter(type,x,U,V,En,nEig2Plot,L);
% ------------------------------------

function harmonicOscillatorPipeline
% Solve eigenvalue problem for 1D Schrodinger equation.
% In this experiment U(x) is a quadratic harmonic oscillator,
%
% The problem is solved by tridiagonal laplacian matrix operator.
% the technique is extendable to arbitrary potentials.
% 
% Then the numerical solution is compared with the analitical solution for the
% harmonic oscillator problem.
    type = 'harmonicOscillator';

    nEig = 100;
    nEig2Plot = 3;

    % Problem properties
    hbar = 1;
    mass = 1;
    
    % Space lattice size & discratization
    L = 10;
    N = 1000;
    x = linspace(-L,L,N);
    
    % Hemiltonian Matrix
    [H,~,U] = hemiltonianMatrixInit(type,x,N,L,hbar,mass);
    
    % Eigenvalue & eigenstates solver
    [V,En] = eigs(H, nEig ,'SA');
    
    % Plot U & first nEig2Plot eigensolutions
    plotter(type,x,U,V,En,nEig2Plot,L);
% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATION INIT UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,K,U] = hemiltonianMatrixInit(type,x,N,L,hbar,mass)
% HEMILTONIANMATRIXINIT
    % Potential definition
    U = potentialInit(type,x);
    
    % Kinetic part, momentum operator with tridiagonal laplacian matrix
    % operator
    lap = laplacianOperator(N,-L,L);
    K = -1/2*(hbar^2/mass)*lap;
    H = K + diag(U)';
% ------------------------------------

function U = potentialInit(type,x)
% POTENTIALINIT
	switch type
        case 'particleInABox'
            U = x.*0;
        case 'harmonicOscillator'
            Uo = 1;
            U = Uo*(0.5*x.^2);
	end
% ------------------------------------

function lap = laplacianOperator( N , mL , L )
% LAPLACIANOPERATOR
	x = linspace(mL,L,N);
	dx = x(2) - x(1);
    lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) + diag(ones((N-1),1),-1))/(dx^2);
% ------------------------------------

function opts = eigsOptionsInit()
% EIGSOPTIONINIT
	opts.tol = 1e-9;
	opts.disp = 0;
	opts.isreal = true;
	opts.issym  = true;
	opts.p = 130 ; 
% ------------------------------------

function plotter(type,x,U,V,E,n,L)
	switch type
        case 'harmonicOscillator'
        	clf
            shg
            [E,ind] = sort(diag(E));
            V = V(:,ind);
            vPlot = V(:,1:n);
            ePlot = E(1:n);
            Usc = U*max(abs(vPlot(:)))/max(abs(U));
            subplot(2,1,1)
            plot(x,vPlot,x,Usc,'--k','LineWidth',1);
            title('Harmonic Oscillator')
            leg = [repmat('E = ',n,1),num2str(ePlot)];
            set(gcf,'menubar','none');
            set(gcf,'Color',[0.6 0.6 0.6]);
            set(gca, 'Color', [0.5,0.5,0.5] );
            legend(leg)
            subplot(2,1,2)
            nPoints = length(E);
            set(gcf,'menubar','none');
            set(gcf,'Color',[0.6 0.6 0.6]);
            set(gca, 'Color', [0.5,0.5,0.5] );
            axis([0 nPoints  0  max(E)])
            exactSolution = animatedline('LineWidth',1,'Color', 'y' );
            numericalSolution = animatedline('LineWidth',1,'Color', 'c' );
            absDiff = animatedline('LineWidth',1,'Color', 'r' );
            legend('Exact solution','Numerical solution','Absolute difference','Location','northwest');
            x = linspace(0,nPoints,nPoints);
            for n=1:nPoints
                addpoints(exactSolution,x(n),harmonicOscillatorExactSolution(n,1,1))
                addpoints(numericalSolution,x(n),(E(n)))
                addpoints(absDiff,x(n),abs(harmonicOscillatorExactSolution(n,1,1) - E(n)))
                pause(0.02)
            end
        case 'particleInABox'
        	clf
            shg
            [E,ind] = sort(diag(E));
            V = V(:,ind);
            vPlot = V(:,1:n);
            ePlot = E(1:n);
            Usc = U*max(abs(vPlot(:)))/max(abs(U));
            subplot(2,1,1)
            plot(x,vPlot,x,Usc,'--k','LineWidth',1);
            title('Particle in a box')
            leg = [repmat('E = ',n,1),num2str(ePlot)];
            set(gcf,'menubar','none');
            set(gcf,'Color',[0.6 0.6 0.6]);
            set(gca, 'Color', [0.5,0.5,0.5] );
            legend(leg)
            subplot(2,1,2)
            nPoints = length(E);
            set(gcf,'menubar','none');
            set(gcf,'Color',[0.6 0.6 0.6]);
            set(gca, 'Color', [0.5,0.5,0.5] );
            axis([0 nPoints  0  max(E)])
            exactSolution = animatedline('LineWidth',1,'Color', 'y' );
            numericalSolution = animatedline('LineWidth',1,'Color', 'c' );
            absDiff = animatedline('LineWidth',1,'Color', 'r' );
            legend('Exact solution','Numerical solution','Absolute difference','Location','northwest');
            x = linspace(0,nPoints,nPoints);
            for n=1:nPoints
                addpoints(exactSolution,x(n),particleInABoxExactSolution(n,L,1,1))
                addpoints(numericalSolution,x(n),(E(n)))
                addpoints(absDiff,x(n),abs(particleInABoxExactSolution(n,L,1,1) - E(n)))
                pause(0.02)
            end
	end
% ------------------------------------

function En = harmonicOscillatorExactSolution(n,hbar,omega)
% HARMONICOSCILLATOREXACTSOLUTION
    En = hbar*omega*(n + 0.5);
% ------------------------------------

function En = particleInABoxExactSolution(n,l,hbar,m)
% PARTICLEINABOXEXACTSOLUTION
    En = n.^2*(hbar^2*pi^2)/(2*m*(2*l)^2);
% ------------------------------------

