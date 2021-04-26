
clear all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            Unified earthquake cycles of               %
%          fault slip and viscoelastic strain           %
%                                                       %
%                  Valere Lambert, 2016                 %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Simulates unified earthquake cycles including slip on 
% two parallel strike-slip faults and viscoelastic strain
% beneath the brittle-ductile transition in 2D antiplane.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%            P H Y S I C A L   M O D E L                %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% rigidity (MPa)
G=30e3;

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% Shear Stress kernels for distributed deformation
s1312=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    +log((x2-L/2).^2+(x3+D+W).^2) - log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    -log((x2-L/2).^2+(x3+D).^2) + log((x2+L/2).^2+(x3+D).^2));

s1212=@(D,L,W,x2,x3) G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2)) ...
    -atan((x3+D+W)./(x2-L/2))-atan((x3+D)./(x2+L/2)) ...
    +atan((x3+D)./(x2-L/2))+atan((x3+D+W)./(x2+L/2)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

s1213=@(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3+D+W).^2) + log((x2+L/2).^2+(x3+D+W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2) ...
    +log((x2-L/2).^2+(x3+D).^2) - log((x2+L/2).^2+(x3+D).^2));

s1313=@(D,L,W,x2,x3) G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W)) ...
    +atan((x2+L/2)./(x3+D))  -atan((x2-L/2)./(x3+D)) ...
    -atan((x2+L/2)./(x3+D+W))+atan((x2-L/2)./(x3+D+W)))...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

% Stress kernels for fault slip
s12h=@(x2,x3,y2,y3,Wf) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-Wf)./((x2-y2).^2+(x3-y3-Wf).^2)-(x3+y3+Wf)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,Wf) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-Wf).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+Wf).^2) ...
    )/2/pi;

% Displacement Kernels for distributed deformation
uk12=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x3-D-W).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x3-D-W).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x3-D).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x3-D).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x2-L/2).*(atan((x3-D-W)./(x2-L/2))-atan((x3-D)./(x2-L/2))) ...
    +2*(x2+L/2).*(atan((x3-D)./(x2+L/2))-atan((x3-D-W)./(x2+L/2))) ...
    +(x3+D+W).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x3+D+W).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x3+D).*log((x2+L/2).^2+(x3+D).^2) ...
    +(x3+D).*log((x2-L/2).^2+(x3+D).^2) ...
    +2*(x2+L/2).*(atan((x3+D+W)./(x2+L/2))-atan((x3+D)./(x2+L/2))) ...
    +2*(x2-L/2).*(atan((x3+D)./(x2-L/2))-atan((x3+D+W)./(x2-L/2))) ...
    );

uk13=@(D,L,W,x2,x3) 1/(2*pi)*( ...
    (x2-L/2).*log((x2-L/2).^2+(x3-D-W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3-D-W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3-D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3-D).^2) ...
    +2*(x3-W-D).*(atan((x2-L/2)./(x3-D-W))-atan((x2+L/2)./(x3-D-W))) ...
    +2*(x3-D).*(atan((x2+L/2)./(x3-D))-atan((x2-L/2)./(x3-D))) ...
    +(x2-L/2).*log((x2-L/2).^2+(x3+D+W).^2) ...
    -(x2+L/2).*log((x2+L/2).^2+(x3+D+W).^2) ...
    -(x2-L/2).*log((x2-L/2).^2+(x3+D).^2) ...
    +(x2+L/2).*log((x2+L/2).^2+(x3+D).^2) ...
    +2*(x3+W+D).*(atan((x2-L/2)./(x3+D+W))-atan((x2+L/2)./(x3+D+W))) ...
    +2*(x3+D).*(atan((x2+L/2)./(x3+D))-atan((x2-L/2)./(x3+D))) ...
    );

% Displacement kernels for fault slip
u1h=@(x2,x3,y2,y3,W) ...
    (+atan((x3-y3)./(x2-y2))-atan((x3+y3)./(x2-y2)) ...
     -atan((x3-y3-W)./(x2-y2))+atan((x3+y3+W)./(x2-y2)) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % %  %

%            -25 km       0       25 km
% _______________________________________________ 0 km    --------> x2
%               |                   |                    |
%               |                   |                    |
%               |      Brittle      |                    |
%               |                   |                    V
%               |                   |                    x3
%          West | fault        East | fault
% ----------------------------------------------- 35 km
%
%              Viscoelastic Shear Zones
%

% Fault Meshes
% West Fault (m)
y2_W=-25e3;
y3=0e3;

% East Fault (m)
y2_E=25e3;

% Brittle-Ductile Tranisition Depth (m)
Transition=35e3;

% Number of patches
ss.M=120;

dz=Transition/ss.M;
fpoles=y3+(0:ss.M)'*dz;
% Tops of fault patches
ss.y3f=fpoles(1:end-1);
% Width of fault patches
Wf=ones(ss.M,1)*dz;
%%
% Shear Zone Mesh
ss.Nx=50;
ss.Nz=50;
eps=1e-12;

% Shear zone Grid
% edges along x3
ss.polesz=Transition+tan((0:ss.Nz)'*pi/(2.2*(ss.Nz+eps)))*Transition;

% center of shear zone (x3)
ss.x3c=(ss.polesz(2:end)+ss.polesz(1:end-1))/2;

% shear zone width 
W=ss.polesz(2:end)-ss.polesz(1:end-1);

xx3t=repmat(ss.polesz(1:end-1)',ss.Nx,1);  % tops
xx3b=repmat(ss.polesz(2:end)',  ss.Nx,1);  % bottoms
xx3c=repmat(ss.x3c',ss.Nx,1);              % centers

% edges along x2
ss.polesxc=(2*y2_W/1e3:(2*y2_E-2*y2_W)/(26e3):2*y2_E/1e3)'*1e3;
edges= floor(ss.Nx-length(ss.polesxc)+1)/2; 
ss.polesxl=flipud(min(ss.polesxc)-tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition);
ss.polesxr=max(ss.polesxc)+tan((0:edges)'*pi/(2.2*(edges)+eps))*Transition;
ss.polesx=[ss.polesxl(1:end-1);ss.polesxc;ss.polesxr(2:end)];

% center of shear zone (x2)
ss.x2c=(ss.polesx(2:end)+ss.polesx(1:end-1))/2;

% shear zone length 
L=ss.polesx(2:end)-ss.polesx(1:end-1);
xx2l=repmat(ss.polesx(1:end-1),1,length(ss.x3c));  % left edge
xx2r=repmat(ss.polesx(2:end),  1,length(ss.x3c));  % right edge
xx2c=repmat(ss.x2c,1,length(ss.x3c));              % center

%%

% Grid points
xpointsFW=[repmat(y2_W-0.15e3,length(ss.y3f),1),...
          repmat(y2_W+0.15e3,length(ss.y3f),1),...
          repmat(y2_W+0.15e3,length(ss.y3f),1),...
          repmat(y2_W-0.15e3,length(ss.y3f),1)]';
      
xpointsFEV=[repmat(y2_E-0.15e3,length(ss.y3f),1),...
  repmat(y2_E+0.15e3,length(ss.y3f),1),...
  repmat(y2_E+0.15e3,length(ss.y3f),1),...
  repmat(y2_E-0.15e3,length(ss.y3f),1)]';
      
zpointsF=[repmat(fpoles(1:end-1),1,1),...
          repmat(fpoles(1:end-1),1,1),...
          repmat(fpoles(2:end),1,1),...
          repmat(fpoles(2:end),1,1)]';

xpoints=[xx2l(:)';xx2r(:)';xx2r(:)';xx2l(:)'];
zpoints=[xx3t(:)';xx3t(:)';xx3b(:)';xx3b(:)'];

%%     
% Surface points
x2GPS=(-400:10:399)'*1e3;
%%
% Stress kernels from fault
ss.k12W=zeros(length(xx2c(:)),ss.M);
ss.k13W=zeros(length(xx2c(:)),ss.M);
ss.k12E=zeros(length(xx2c(:)),ss.M);
ss.k13E=zeros(length(xx2c(:)),ss.M);

ss.KWE=zeros(ss.M,ss.M);
ss.KWW=zeros(ss.M,ss.M);
ss.KEW=zeros(ss.M,ss.M);
ss.KEE=zeros(ss.M,ss.M);

ss.ku1W=zeros(length(x2GPS),ss.M);
ss.ku1E=zeros(length(x2GPS),ss.M);
% Fields from Faults
for k=1:ss.M
    % we evaluate the stress at the center of the shear zones
    ss.k12W(:,k)=s12h(xx2c(:),xx3c(:),y2_W,ss.y3f(k),Wf(k));  
    ss.k13W(:,k)=s13h(xx2c(:),xx3c(:),y2_W,ss.y3f(k),Wf(k));
    
    ss.k12E(:,k)=s12h(xx2c(:),xx3c(:),y2_E,ss.y3f(k),Wf(k));  
    ss.k13E(:,k)=s13h(xx2c(:),xx3c(:),y2_E,ss.y3f(k),Wf(k));
    
    % we evaluate the stress at the center of the fault patches
    ss.KWE(:,k)=s12h(y2_E,ss.y3f+dz/2,y2_W,ss.y3f(k),Wf(k));
    ss.KWW(:,k)=s12h(y2_W,ss.y3f+dz/2,y2_W,ss.y3f(k),Wf(k));
    ss.KEW(:,k)=s12h(y2_W,ss.y3f+dz/2,y2_E,ss.y3f(k),Wf(k));
    ss.KEE(:,k)=s12h(y2_E,ss.y3f+dz/2,y2_E,ss.y3f(k),Wf(k));
    
    % Displacement kernels
    ss.ku1W(:,k)=u1h(x2GPS,0,y2_W,ss.y3f(k),Wf(k));
    ss.ku1E(:,k)=u1h(x2GPS,0,y2_E,ss.y3f(k),Wf(k));
end

% Stress kernels from shear zones
ss.k1312=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1313=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1212=zeros(length(ss.x2c)*length(ss.x3c));
ss.k1213=zeros(length(ss.x2c)*length(ss.x3c));

ss.k1212fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1312fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1213fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1313fW=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));

ss.k1212fE=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1312fE=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1213fE=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));
ss.k1313fE=zeros(length(ss.y3f),length(ss.x2c)*length(ss.x3c));

% Establish displacement kernels
ss.ku12=zeros(length(x2GPS),length(ss.x2c)*length(ss.x3c));
ss.ku13=zeros(length(x2GPS),length(ss.x2c)*length(ss.x3c));

% Fields from Shear zones
for kx=1:length(ss.x2c)
    for kz=1:length(ss.x3c)
        % we evaluate the stress at the center of the shear zones
        ss.k1212(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1312(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1213(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));
        ss.k1313(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),xx2c(:)-ss.x2c(kx),xx3c(:));

        % we evaluate stress at the center of the fault patches 
        ss.k1212fW(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),y2_W-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1312fW(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),y2_W-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1213fW(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),y2_W-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1313fW(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),y2_W-ss.x2c(kx),ss.y3f(:)+dz/2);
        
        ss.k1212fE(:,(kz-1)*ss.Nx+kx)=s1212(ss.polesz(kz),L(kx),W(kz),y2_E-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1312fE(:,(kz-1)*ss.Nx+kx)=s1312(ss.polesz(kz),L(kx),W(kz),y2_E-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1213fE(:,(kz-1)*ss.Nx+kx)=s1213(ss.polesz(kz),L(kx),W(kz),y2_E-ss.x2c(kx),ss.y3f(:)+dz/2);
        ss.k1313fE(:,(kz-1)*ss.Nx+kx)=s1313(ss.polesz(kz),L(kx),W(kz),y2_E-ss.x2c(kx),ss.y3f(:)+dz/2);
        
        % we evaluate the displacement at the top of the patches
        ss.ku12(:,(kz-1)*ss.Nx+kx)=uk12(ss.polesz(kz),L(kx),W(kz),x2GPS(:)-ss.x2c(kx),0*x2GPS(:));
        ss.ku13(:,(kz-1)*ss.Nx+kx)=uk13(ss.polesz(kz),L(kx),W(kz),x2GPS(:)-ss.x2c(kx),0*x2GPS(:));
    end
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Confining pressure (MPa) and Temperature (K)
k  = 3.138; % thermal conductivity (W/m/K)
Cp = 1171 ; % specific heat (J/kg/K)
Rm = 3330 ; % mantle density (kg/m^3)

Pconf       = Rm*9.8*ss.x3c/1e6;  % Shear zones
Pconf_fault = Rm*9.8*(ss.y3f+dz); % Faults

Kappa     = k / (Rm * Cp); % Thermal diffusivity (m^2/s)
Age_plate = 2e15; % seconds
ss.Tprof  = 300+1380*erf(ss.x3c/(sqrt(4* Kappa * Age_plate)));  % Kelvin

% default friction properties (velocity-weakening friction)
% effective confining pressure (MPa)
ss.sigmab = 1000;

% frictional parameters
ss.aW = 1e-3*ones(size(ss.y3f));
ss.bW = ss.aW+2.1e-4*ones(size(ss.y3f));

ss.aE = 1e-3*ones(size(ss.y3f));
ss.bE = ss.aE+2.1e-4*ones(size(ss.y3f));

% static friction coefficient
ss.mu0 = 0.2*ones(size(ss.y3f));

% characteristic weakening distance (m)
ss.L = 0.012*ones(size(ss.y3f));

% plate velocity (m/s)
ss.V_plate = 1e-9*ones(size(ss.y3f));

% reference slip rate (m/s)
ss.Vo = 1e-6*ones(size(ss.y3f));

% shear wave speed (m/s)
ss.Vs = 3e3*ones(size(ss.y3f));

% Velocity-strengthening at edges 
% West fault: 5-15km velocity-weakening
topW    = floor(5e3/(Transition/ss.M));
bottomW = ceil(15e3/(Transition/ss.M));
ss.bW(1:topW)      = ss.aW(1:topW)-2.1e-4*ones(topW,1);
ss.bW(bottomW:end) = ss.aW(bottomW:end)-2.1e-4*ones(length(ss.aW(bottomW:end)),1);

% East fault: 3-20km velocity-weakening
topE    = floor(3e3/(Transition/ss.M));
bottomE = ceil(20e3/(Transition/ss.M));
ss.bE(1:topE)      = ss.aE(1:topE)-2.1e-4*ones(topE,1);
ss.bE(bottomE:end) = ss.aE(bottomE:end)-2.1e-4*ones(length(ss.aE(bottomE:end)),1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Values taken for wet olivine - Hirth, G. and D. Kohlstedt (2003)

% Driving strain rate (1/s)
ss.e12p_plate = 1e-14*ones(length(ss.x2c)*length(ss.x3c),1);
ss.e13p_plate =      zeros(length(ss.x2c)*length(ss.x3c),1);

% Rheological Parameters
% Reference Strain Rate (for stress in MPa)
ss.Adif = 1e6*ones(length(ss.x3c)*length(ss.x2c),1);
ss.Adis = 90 *ones(length(ss.x3c)*length(ss.x2c),1);

% Power-Law Exponent
ss.n = 3.5*ones(length(ss.x3c)*length(ss.x2c),1);

% Activation Energy Wet Oliving (J/mol)
ss.Qdif = 335e3*ones(length(ss.x3c)*length(ss.x2c),1);
ss.Qdis = 480e3*ones(length(ss.x3c)*length(ss.x2c),1);

% Activation Volume (m^3/mol)
ss.Voldif = 4e-6*ones(length(ss.x3c)*length(ss.x2c),1);  
ss.Voldis = 11e-6*ones(length(ss.x3c)*length(ss.x2c),1);  

% Grain size (m)
ss.d    = 1e-2*ones(length(ss.x3c)*length(ss.x2c),1);
ss.pexp = 3*ones(length(ss.x3c)*length(ss.x2c),1);

% Water fugacity (H/10^6 Si)
ss.COH = 1000*ones(length(ss.x3c)*length(ss.x2c),1);
ss.r   = 1.2*ones(length(ss.x3c)*length(ss.x2c),1);

% Pressure (Pa)
ss.P = repmat(1e6*Pconf',length(ss.x2c),1);
ss.P = reshape(ss.P,[length(ss.x2c)*length(ss.x3c),1]);

% Temperature (K)
Te0 = repmat(ss.Tprof',length(ss.x2c),1);
Te0 = reshape(Te0,[length(ss.x2c)*length(ss.x3c),1]);

% Coefficients for dislocation and diffusion creep
ss.Const_dis = ss.Adis.*exp(-(ss.Qdis+ss.P.*ss.Voldis)./(8.314.*Te0)).*ss.COH.^(ss.r);
ss.Const_diff = ss.Adif.*exp(-(ss.Qdif+ss.P.*ss.Voldif)./(8.314.*Te0)).*ss.COH.^(ss.r).*ss.d.^(-ss.pexp);

% Strengh profile
s120 = (ss.e12p_plate./ss.Const_dis).^(1./ss.n);
s130 = zeros(size(s120));
e120 = zeros(size(s120));
e130 = zeros(size(s120));

% Fault Strength
ss.strength_E = ss.sigmab.*(ss.mu0+(ss.aE-ss.bE).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);
ss.strength_W = ss.sigmab.*(ss.mu0+(ss.aW-ss.bW).*log(ss.V_plate./ss.Vo))+G*ss.V_plate./(2*ss.Vs);

% Plot strength profiles
figure(1);clf;
subplot(2,1,1);
plot(ss.y3f/1e3,ss.strength_E,ss.y3f/1e3,ss.strength_W)
xlabel('Depth (km)')
ylabel('Strength (MPa)');
subplot(2,1,2);
plot(ss.x3c/1e3,log10(s120(1:length(ss.x2c):end)))
xlim([ss.x3c(1)/1e3,ss.x3c(end)/1e3]);
xlabel('Depth (km)')
ylabel('Strength (MPa) log10');


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
%         N U M E R I C A L   S O L U T I O N           %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% state parameters
ss.dgfF=3;
ss.dgfS=4;
%% Initialize State Vector
Y0=zeros(ss.M*ss.dgfF+ss.M*ss.dgfF+length(ss.x2c)*length(ss.x3c)*ss.dgfS,1);

% Fault patches
Y0(1:ss.dgfF:ss.M*ss.dgfF)=zeros(size(ss.y3f));
Y0(2:ss.dgfF:ss.M*ss.dgfF)=ss.strength_W;
Y0(3:ss.dgfF:ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

Y0(ss.M*ss.dgfF+1:ss.dgfF:2*ss.M*ss.dgfF)=zeros(size(ss.y3f));
Y0(ss.M*ss.dgfF+2:ss.dgfF:2*ss.M*ss.dgfF)=ss.strength_E;
Y0(ss.M*ss.dgfF+3:ss.dgfF:2*ss.M*ss.dgfF)=log(ss.Vo./ss.V_plate);

% Shear zones
Y0(2*ss.M*ss.dgfF+1:ss.dgfS:end)=s120;
Y0(2*ss.M*ss.dgfF+2:ss.dgfS:end)=s130;
Y0(2*ss.M*ss.dgfF+3:ss.dgfS:end)=e120;
Y0(2*ss.M*ss.dgfF+4:ss.dgfS:end)=e130;

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) odeViscoelastic(t,y,ss);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6); 
[t,Y]=ode45(yp,[0 1e10],Y0,options);
toc
%%
% Compute the instantaneous derivative
Yp=zeros(length(t)-1,size(Y,2));
for k=1:length(t)-1
    Yp(k,:)=(Y(k+1,:)-Y(k,:))/(t(k+1)-t(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Strain rate at center
Ep=sqrt(Yp(:,2*ss.M*ss.dgfF+floor(length(ss.x2c)/2)*ss.dgfS+3:ss.dgfS*length(ss.x2c):end)'.^2 +...
        Yp(:,2*ss.M*ss.dgfF+floor(length(ss.x2c)/2)*ss.dgfS+4:ss.dgfS*length(ss.x2c):end)'.^2);
    
% Velocity
V_W=Yp(:,1:ss.dgfF:ss.M*ss.dgfF);
V_E=Yp(:,ss.M*ss.dgfF+1:ss.dgfF:2*ss.M*ss.dgfF);

% Maximum Velocity
VWmax=zeros(length(t)-1,1);
VEmax=zeros(length(t)-1,1);
for ts=1:length(t)-1
    VWmax(ts)=max(V_W(ts,:));
    VEmax(ts)=max(V_E(ts,:));
end

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
figure(1);clf;set(gcf,'name','Time Evolution')
f1=subplot(5,1,1);cla;
pcolor(t(1:end-1)/3.15e7,ss.y3f/1e3,log10(V_W')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V_W))) max(max(log10(V_W)))]);
colormap(f1,parula);
title(h,'Slip Rate West (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

f2=subplot(5,1,2);cla;
pcolor(t(1:end-1)/3.15e7,ss.y3f/1e3,log10(V_E')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V_W))) max(max(log10(V_W)))]);
colormap(f2,parula);
title(h,'Slip Rate East (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

f3=subplot(5,1,3);cla;
pcolor(t(1:end-1)/3.15e7,ss.x3c/1e3,log10(Ep)), shading flat
set(gca,'YDir','reverse');

caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
h1=colorbar('Location','NorthOutside');
colormap(f3,hot);
title(h1,'Strain Rate (1/s)')
xlabel('Time (Yr)')
ylabel('Depth (km)');

subplot(5,1,4);cla;
plot(t(1:end-1)/3.15e7,log10(VWmax),t(1:end-1)/3.15e7,log10(VEmax))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Maximum slip rates on faults')

subplot(5,1,5);cla;
plot(t(1:end-1)/3.15e7,log10(V_W(:,floor((topW+bottomW)/2))),...
     t(1:end-1)/3.15e7,log10(V_E(:,floor((topE+bottomE)/2))))
xlabel('Time (Yr)')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zones')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f1=figure(2);clf;set(gcf,'name','Time Step Evolution')

subplot(5,1,1);cla;
pcolor(1:length(t)-1,ss.y3f/1e3,log10(V_W')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V_W))) max(max(log10(V_W)))]);
colormap(f1,parula);
title(h,'Slip Rate West (m/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

f2=subplot(5,1,2);cla;
pcolor(1:length(t)-1,ss.y3f/1e3,log10(V_E')), shading flat
set(gca,'YDir','reverse');

h=colorbar('Location','NorthOutside');
caxis([min(min(log10(V_W))) max(max(log10(V_W)))]);
colormap(f2,parula);
title(h,'Slip Rate East (m/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

f3=subplot(5,1,3);cla;
pcolor(1:length(t)-1,ss.x3c/1e3,log10(Ep)), shading flat
set(gca,'YDir','reverse');

caxis([log10(min(min(Ep))) log10(max(max(Ep)))]);
h1=colorbar('Location','NorthOutside');
colormap(f3,hot);
title(h1,'Strain Rate (1/s)')
xlabel('Time Steps')
ylabel('Depth (km)');

subplot(5,1,4);cla;
plot(1:length(t)-1,log10(VWmax),1:length(t)-1,log10(VEmax))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Maximum slip rates on faults')

subplot(5,1,5);cla;
plot(1:length(t)-1,log10(V_W(:,floor((topW+bottomW)/2))),...
     1:length(t)-1,log10(V_E(:,floor((topE+bottomE)/2))))
xlabel('Time Steps')
ylabel('Velocity (m/s) log10')
title('Time series at center of seismogenic zones')
