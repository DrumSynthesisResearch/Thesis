%%%%%%%%%%%%%%%%%%%%%%FDTD Code example%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
%Grid dimensions
Nx=40;
Ny=40;

%Sample rate
Fs=44100;

%Basic Operators
Ix=eye(Nx+1);
Iy=eye(Ny+1);
h=1;
Dxx=(-2.*Ix+diag(ones(Nx,1),1)+diag(ones(Nx,1),-1))./(h^2);
Dyy=(-2.*Ix+diag(ones(Ny,1),1)+diag(ones(Ny,1),-1))./(h^2);
Dxx2=kron(Dxx,Iy);
Dyy2=kron(Ix,Dyy);
Lap=Dxx2+Dyy2;
Bih=Lap*Lap;

%Excite with...
option=1
%option=2

if option==1
% ...modal surface of (m,n) where m is # modal diameters 
% and n is # nodal circles
k=1/Fs*5000; %Time step

    m=3;
    n=2;
    Rstrike=0;
    X=(-Nx:2:Nx);
    Y=(-Ny:2:Ny);
    [x,y] = meshgrid(X,Y);
    [THETA,R]=cart2pol( x,y);
    RR=R/Nx;
    Xmn=max(besselzero(m,n)); %nth non-trivial zero of Jm
    Amp=5 %5cm Amplitude 
    W=(Amp.*cos(m*THETA).*besselj(m,RR*Xmn)).*(R<Nx);
    surf(W)
end

if option==2
%... Central strike (raised cosine)
k=1/Fs*1000; %Time step

    X=(-Nx:2:Nx);
    Y=(-Ny:2:Ny);
    [x,y] = meshgrid(X,Y);
    [THETA,R]=cart2pol( x,y);
    RR=R/Nx;
    W=5*cos(RR*20).*(RR<pi/(2*20));
    surf(W)
end

%Matrix to concatenated vector
w=reshape(W,[(Nx+1)*(Ny+1) 1]); %Grid vector
BC=reshape(R<=Nx,[(Nx+1)*(Ny+1) 1]); %Circular boundary condition vector

%Constants
sig0=1; %Damping 1
sig1=1; %Damping 2
r=0.2; %radius
H=1.75E-4; %Thickness
v=0.25; %Poisson's ratio (0-0.5)
A=pi*r^2; %Membrane area
E=3.5E9; %Young's Modulus
D=(E*(H^3))/(12*(1-(v^2))); %Bending stiffness
xi=6*(D/(A*H^2)); %Berger Plate coefficient
rho=1400; %Surface density
T=1140; %Tension
kap=sqrt(D/(rho*H)); %Kappa
I=eye((Nx+1)*(Ny+1)); %Identity matrix
q=(h*k*sqrt(xi/2)).*Lap*w;


%Set number of iterations and preallocate movie frames
loops = 0.5*Fs;
clear M
M(ceil(loops/100)) = struct('cdata',[],'colormap',[]);

%Update scheme
B=2*I-(k.^2)*(kap.^2).*Bih;
wMinus1=w;
wCurrent=w;
j=1
for i=1:loops
wPlus1=B*wCurrent-wMinus1;
wMinus1=wCurrent.*BC;
wCurrent=wPlus1.*BC;
if floor(i/100)==i/100
j=j+1
MCurrent=reshape(wCurrent,[(Nx+1) (Ny+1)]);
surf(MCurrent);
xlim([0 Nx]);
ylim([0 Ny]);
zlim([-5 5]);
M(j) = getframe;
end
end
