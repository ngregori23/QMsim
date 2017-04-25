function schrodinger2D
% SCHRODINGER2D
particleInABoxPipeline

% ------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         EXPERIMENTS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function particleInABoxPipeline
% PARTICLEINABOXPIPELINE
    type = 'particleInABox';
    
    n = 2;
    m = 3;
    
    nmodes = max(n,m);
    
    % Problem properties
    hbar = 1;
    mass = 1;
    
    % Space lattice size & discratization 
    L = 1;
    N = 512;
    X = linspace( -L, L , N);
    Y = linspace( -L ,L , N);
    
	% Hemiltonian Matrix
	[Hx,Hy] = hemiltonianMatrixInit(N,L,hbar,mass);
    
	% Eigenvalue & eigenstates solver
    [Psix,Ex] = eigs(Hx,nmodes,'sa');
    evalx = diag(Ex);
    [Ex, permx] = sort(evalx);
    Psix = Psix(:,permx);

    [Psiy,Ey] = eigs(Hy,nmodes,'sa');
    evaly = diag(Ey); 
    [Ey, permy] = sort(evaly);
    Psiy = Psiy(:,permy);
    
    px = Psix(:,n);
    py = Psiy(:,m);
    
    
    Vnm = px*py';
    
	% Plot the eigensolutions
    plotter(X,Y,Vnm,n,m);
  % ------------------------------------  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIMULATION INIT UTILS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Hx,Hy] = hemiltonianMatrixInit(N,L,hbar,mass)
% HEMILTONIANMATRIXINIT
    % Kinetic part, momentum operator with tridiagonal laplacian matrix
    % operator
    lap = laplacianOperator(N,-L,L);
    Hx = -1/2*(hbar^2/mass)*lap ;
    Hy = -1/2*(hbar^2/mass)*lap ;
% ------------------------------------

function lap = laplacianOperator( N , mL , L )
% LAPLACIANOPERATOR
    x = linspace(mL,L,N);
    dx = x(2) - x(1);
    lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) + diag(ones((N-1),1),-1))/(dx^2);
% ------------------------------------

function plotter(x,y,v,n,m)
% PLOTTER
	clf
	shg
    v = v/max(max(v));
    subplot(1,2,1)
    s = surf(x,y,v);
	set(gcf,'Color',[0.6 0.6 0.6]);
	set(gca, 'Color', [0.5,0.5,0.5] );
    xlabel('X','FontSize',22)
    ylabel('Y','FontSize',22)
    zlabel(['\Psi',num2str(n),num2str(m)],'FontSize',22)
	s.EdgeColor = 'interp';
    subplot(1,2,2)
    s = surf(x,y,v);
	set(gcf,'Color',[0.6 0.6 0.6]);
	set(gca, 'Color', [0.5,0.5,0.5] );
    xlabel('X','FontSize',22)
    ylabel('Y','FontSize',22)
    zlabel(['\Psi',num2str(n),num2str(m)],'FontSize',22)
	s.EdgeColor = 'interp';
    view(0,90)
% ------------------------------------
