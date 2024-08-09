%--------------------------------------------------------------------------
% ceren.m
%
% (Particle Physics and Optics)
%
% This program simulates generation of Cherenkov photons when a particle
% is passing through an optical material (Radiator) having rectangular
% prism shape. These photons are generated inside the radiator when 
% particle's speed is greater then that of speed of light in radiator.
% Cherenkov photons are created while the following inequality is satisfied:
%     Beta*n > 1
% (where Beta = v_particle/c and n = c/v_light, and c = speed of light in vacuum.
% Production points of photons and their directions are written to a file
% called "cherenkov.dat" which can be used in NSC mode of Zemax OpticStudio.
% Note that in order to use these Cherenkov Photons, you need to fill
% "Inside Of" column in Zemax NSC Mode. Then, Zemax can do Ray Tracing
% and model absorption and reflection effects (not scattering) for each photon.
%
% In the program, we assume that
%   All distance units are in mm
%   All photon wavelength units are in um (micrometer)
%   All particle units are in GeV
%
% Inputs:
%   p        = momentum of the incident particle
%   m        = mass of the indident particle
%   theta    = particle direction w.r.t. z-axis
%   phi      = particle direction w.r.t. x-axis
%   Radiator = optical material name (to get index data)
%   Width    = width of the optical material in mm
%   Heigth   = height of the optical material in mm
%   Thick    = thickness of the optical material in mm
%
% Output file format:
%
%   Number_of_lines Dimension_Flag
%   x1 y1 z1 l1 m1 n1 intensity1 wavelength1
%   x2 y2 z2 l2 m2 n2 intensity2 wavelength2
%   ...
%
%   where
%   Number_of_Lines are number of photons
%   Dimension_Flag is 4 for mm (lens unit in Zemax)
%   x1 y1 z1       are the generation point of photon_1 in 3D space
%   l1 m1 n1       are the direction cosines of photon 1
%   wavelength1    is the wavelength of the photon 1 in um
%   intensity1 = intensity2 = ... = 1
%
%   First Developed : Jul 2023
%   Last Modified   : Feb 2024 (wavelength2color function is added).
%   Author          : Ahmet.Bingul@cern.ch
%--------------------------------------------------------------------------
clear; clc;

% We have some global variables
global Radiator
global Beta
global SF11 SF10 BK7 PMMA AEROGEL WATER H2O PS G15

% Output file name
filename = 'cherenkov.dat';
fileID = fopen(filename,'w');
printToFile = 1; % 0=no or 1=yes

% Particle info
p     = 1.0;              % particle momentum in GeV
m     = 0.1056583755;     % particle mass in GeV
theta = (10*pi/180)*rand; % particle (random) direction w.r.t z-axis
phi   = 2*pi*rand;        % particle (random) direction w.r.t x-axis

% Setup radiator material type and geometry
SetRadiatorNames();
Radiator = PS;
Width = 50;
Height= 50;
Thick = 50;

%  -----                        -----
% --- simulations start from here ---
%  -----                        -----

% Particle's 4-vector and related calculations
E     = sqrt(p*p+m*m); % paricle energy
px    = p*sin(theta)*cos(phi);
py    = p*sin(theta)*sin(phi);
pz    = p*cos(theta);
p4    = [px, py, pz, E];
Beta  = p/E;
gamma = 1/sqrt(1-Beta^2);
% unit vector in particle direction
u     = p4(1:3)/norm(p4(1:3));

% Particle initial and final hit positions on the optical material
% (Particle enters to glass at z = 0 and leaves it at z = Thick)
x0 = -Width/2  + Width*rand;
y0 = -Height/2 + Height*rand;
z0 = 0;
xf = x0 + u(1)*(Thick-z0)/u(3);
yf = y0 + u(2)*(Thick-z0)/u(3);
zf = z0 + u(3)*(Thick-z0)/u(3);

% Total track length of the particle in the material
TOTL = sqrt((xf-x0)^2+(yf-y0)^2+(zf-z0)^2);

% Rotation matrix and its transpose (= inverse)
R  = Rotate(u);
RT = R';

% Draw glass and particle trajectory (optional)
figure(1)
plot3(0,0,0);
hold on
DrawGlass(Width,Height,Thick)
DrawParticle(x0,y0,z0,xf,yf,zf)
axis equal
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)'); grid on
title('Only 2% of photons are shown.')

% Here we calculate number of Cherenkov photons generated in the 
% given wavelength range. As described in Geat4 Manual, the final number 
% of photons produced is calculated from a Poisson distribution with a mean 
% of <n> = StepLength * dN/dx
WLmin = 0.4; % wavelength
WLmax = 0.7; % range
Ncpmm = integrate(WLmin, WLmax); % number of photons per mm
Nph   = Ncpmm * TOTL;            % total (mean) number of photons
Nc    = round(poissrnd(Nph));    % final number of photons generated

% Print the header info to the top of the output file
if printToFile ~= 0
    fprintf(fileID,'%d 4\n',Nc);
end

% Construct Cherenkov photon data
lambda   = zeros(Nc,1); % wavelength from wavelength distribution
theta_ch = zeros(Nc,1); % Cherenckov angle for each photon
theta_cr = zeros(Nc,1); % Critical angle for each photon

% Photon loop to generate location and direction of Cherenkov lights
for i = 1:Nc
   lambda(i)   = GetWavelength(WLmin,WLmax);
   theta_ch(i) = acos(1/(index(lambda(i))*Beta));
   theta_cr(i) = asin(1/index(lambda(i)));
   % local coordinates of photon
   lx0 = 0.0;
   ly0 = 0.0;
   lz0 = rand * TOTL;
   phi = 2*pi*rand;
   lx  = sin(theta_ch(i))*cos(phi)+lx0;
   ly  = sin(theta_ch(i))*sin(phi)+ly0;
   lz  = cos(theta_ch(i))+lz0;
   % transform local coordinates to global
   r0  = RT * [lx0 ly0 lz0]' + [x0 y0 z0]';
   r1  = RT * [lx  ly  lz ]' + [x0 y0 z0]';
   % unit vector in the direction of photon
   U = (r1-r0) / norm(r1-r0);
   U = U';
   % direction cosine of photons
   cosTx = dot(U,[1 0 0]);
   cosTy = dot(U,[0 1 0]);
   cosTz = dot(U,[0 0 1]);
   % draw photons
   if rand<0.02
     s = (Thick-r0(3))/U(3);
     rph = r0' + U*s;
     colorCode = wavelength2color(lambda(i)*1000, 'gammaVal', 1, ...
                 'maxIntensity', 255, 'colorSpace', 'rgb');
     colorCode = colorCode / 255;
     line([r0(1) rph(1)],[r0(2) rph(2)],[r0(3) rph(3)],'color',colorCode)
   end
   % print photon's phusical data to the file
   if printToFile ~= 0
     fprintf(fileID,'%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %d %6.4f\n', ...
                    r0(1),r0(2),r0(3), cosTx,cosTy,cosTz, 1, lambda(i));
   end
end

if printToFile ~= 0
  fclose(fileID);
end
hold off

fprintf('Final Report:\n');
fprintf(' Particle Momentum      : %f GeV\n',p);
fprintf(' Particle Mass          : %f GeV\n',m);
fprintf(' Particle Energy        : %f GeV\n',E);
fprintf(' Particle Beta (=v/c)   : %f \n',Beta);
fprintf(' Radiator material      : %s \n',GetRadiatorName());
fprintf(' TOTL in radiator       : %f mm\n',TOTL);
fprintf(' Wavelength range       : [%5.2f, %5.2f] um\n',WLmin,WLmax);
fprintf(' # of photons/mm        : %f \n',Ncpmm);
fprintf(' Total # of photons     : %d \n',Nc);
fprintf(' Output source filename : %s \n',filename);

% plot the wavelength, cherenkov angle, etc. distributions
figure(2)
histogram(lambda,20);
xlabel('wavelength (um)'); 
ylabel('Entries');

figure(3)
histogram(theta_ch*180/pi,20);
xlabel('Cherenkov angle (deg)');
ylabel('Entries');

figure(4)
histogram(theta_cr*180/pi,20);
xlabel('Critical angle (deg)');
ylabel('Entries');

figure(5)
plot(lambda,theta_cr*180/pi,'.');
xlabel('Wavelength');
ylabel('cherenkov Angle');
grid on


% End of main part of the program
%**************************************************************************

% -----                         -----
%--- Sub programs starts from here ---
% -----                         -----

%--------------------------------------------------------------------------
% Returns index of refraction for the given wavelenth x in um.
%--------------------------------------------------------------------------
function n = index(x)
  global SF11 SF10 BK7 PMMA AEROGEL WATER H2O PS G15
  global Radiator
  n = 1.0;
  if BK7~=2
      disp('You should first call SetRadiatorNames().');
      return
  end

      if Radiator == SF11 && x>=0.37 && x<=2.5
         n=sqrt(1+1.73759695./(1-0.013188707./x.^2)+0.313747346./(1-0.0623068142./x.^2)+1.89878101./(1-155.23629./x.^2));
  elseif Radiator == SF10 && x>=0.38 && x<=2.5
         n=sqrt(1+1.62153902./(1-0.0122241457./x.^2)+0.256287842./(1-0.0595736775./x.^2)+1.64447552./(1-147.468793./x.^2));
  elseif Radiator == BK7 && x>=0.3 && x<=2.5
         n=sqrt(1+1.03961212./(1-0.00600069867./x.^2)+0.231792344./(1-0.0200179144./x.^2)+1.01046945./(1-103.560653./x.^2));
  elseif Radiator == PMMA && x>=0.4047 && x <=1.083
          n=sqrt(1+0.99654./(1-0.00787./x.^2)+0.18964./(1-0.02191./x.^2)+0.00411./(1-3.85727./x.^2));
  elseif Radiator == AEROGEL && x>=0.35 && x <=0.70
         % see https://cds.cern.ch/record/1431905/files/LHCb-TALK-2007-072.pdf
         x = x*1000; % convert to nm
         n=sqrt(1+0.05639*x.^2./(x.^2-6.9256e+03));
  elseif (Radiator == WATER || Radiator == H2O) && x>=0.2 && x <=1.0
         n=0.937*sqrt(1+9.311*x.^2./(9.311*x.^2-0.06303));
  elseif Radiator == PS && x>=0.37 && x<=1.0
         n=1.576-0.01434./x+0.0144./x.^2;
  elseif Radiator == G15
         n = 1.5;
  end
end
%--------------------------------------------------------------------------
% integrant for the calcuation of number of Cherenkov photons/unit length
%--------------------------------------------------------------------------
function y = f(x)
    global Beta
    alpha = 1/137.035999084;
    n = index(x);
    y = (2*pi*alpha/(x*x)) * (1-1/(Beta*Beta * n*n));
end
%--------------------------------------------------------------------------
% Calculate number of cherenkov photons in the wavelength range [a,b]
% wavelenth must be in um.
%--------------------------------------------------------------------------
function Nph = integrate(a,b)
    n = 1000;
    s = 0.5*(f(a)+f(b));
    h = abs(b-a)/n;
    factor = 1e+3;
    for i =1:n-1
        s = s + f(a+i*h);
    end
  Nph = s * h * factor;
end
%--------------------------------------------------------------------------
% Returns a random wavelength in um between [Lmin,Lmax] from the Cherenkov
% distibution using rejection method (from the integrant function, f(x)).
%--------------------------------------------------------------------------
function x = GetWavelength(Lmin,Lmax)
   fmax = f(Lmin);
   while 1
       x = Lmin + (Lmax-Lmin)*rand;
       y = f(x);
       ptest = fmax*rand;
       if y > ptest
           return
       end
   end
end
%--------------------------------------------------------------------------
% Setup radiator names and codes.
%--------------------------------------------------------------------------
function SetRadiatorNames()
   global SF11 SF10 BK7 PMMA AEROGEL WATER H2O PS G15
   SF11=0; SF10=1; BK7=2; PMMA=3; AEROGEL=4; WATER=5; H2O=5; PS=6; G15=7;
end
%--------------------------------------------------------------------------
% Setup radiator names and codes.
%--------------------------------------------------------------------------
function name = GetRadiatorName()
   name = 'none';
   global Radiator
   if Radiator == 0
       name = 'SF11';
   elseif Radiator == 1
       name = 'SF10';
   elseif Radiator == 2
       name = 'BK7';
   elseif Radiator == 3
       name = 'PMMA';
   elseif Radiator == 4
       name = 'AEROGEL';
   elseif Radiator == 5
       name = 'WATER (H2O)';
   elseif Radiator == 6
       name = 'POLYSTYRENE';
   elseif Radiator == 7
       name = 'G15 (constant index n = 1.5)';
   end
end
%--------------------------------------------------------------------------
% Setup a special rotation matrix using zyz convention.
% In new rotated coordinate, z'-axis is selected in the direction of input
% vector u = (ux,uy,uz) (parallel to incident particle direction).
% So, after rotation, new coordinates of u becomes u' = (0,0,uz).
%--------------------------------------------------------------------------
function R = Rotate(u)
   R = zeros(3,3);
   lm = length(u);
   if lm ~= 3
       disp('Input vector must be 3x1 (or 1x3) vector! zeros(3,3) is returned.')
       return
   end
   u  = u / norm(u); % normalize input vector
   A  = atan(u(2)/u(1));
   if (u(2)<0 && u(1)<0) || (u(2)>0 && u(1)<0)
      A = A + pi;
   end
   B  = acos(u(3));
   C  = 0;
   Rz = [cos(A) sin(A) 0 ; -sin(A) cos(A) 0; 0      0 1];
   Ry = [cos(B) 0 -sin(B);  0      1      0; sin(B) 0 cos(B)];
   Rzz= [cos(C) sin(C) 0 ; -sin(C) cos(C) 0; 0      0 1];
   R  = Rzz*Ry*Rz;
   return
   % matrix form for C/C++ or Python
   %ca = cos(A); sa = sin(A);
   %cb = cos(B); sb = sin(B);
   %cg = cos(C); sg = sin(C);
   %R[0][0] =  cg*cb*ca-sg*sa;
   %R[0][1] =  cg*cb*sa+sg*ca;
   %R[0][2] = -cg*sb;
   %R[1][0] = -sg*cb*ca-cg*sa;
   %R[1][1] = -sg*cb*sa+cg*ca;
   %R[1][2] =  sg*sb;
   %R[2][0] =  sb*ca;
   %R[2][1] =  sb*sa;
   %R[2][2] =  cb;
end

%--------------------------------------------------------------------------
% Draws the optical material.
%--------------------------------------------------------------------------
function DrawGlass(w,h,t)
  c = 'blue';
  line([+w -w], [+h +h], [t t],'color',c,'lineWidth',2)
  line([-w -w], [+h -h], [t t],'color',c,'lineWidth',2)
  line([-w +w], [-h -h], [t t],'color',c,'lineWidth',2)
  line([+w +w], [-h +h], [t t],'color',c,'lineWidth',2)
  line([+w -w], [+h +h], [0 0],'color',c,'lineWidth',2)
  line([-w -w], [+h -h], [0 0],'color',c,'lineWidth',2)
  line([-w +w], [-h -h], [0 0],'color',c,'lineWidth',2)
  line([+w +w], [-h +h], [0 0],'color',c,'lineWidth',2)
  line([+w +w], [+h +h], [0 t],'color',c,'lineWidth',2)
  line([-w -w], [+h +h], [0 t],'color',c,'lineWidth',2)
  line([-w -w], [-h -h], [0 t],'color',c,'lineWidth',2)
  line([+w +w], [-h -h], [0 t],'color',c,'lineWidth',2)
  %s = 1.2;
  %line(s*[+w -w], s*[+h +h], s*[t t],'color','w','lineWidth',0.5)
  %line(s*[-w -w], s*[+h -h], s*[t t],'color','w','lineWidth',0.5)
  %line(s*[-w +w], s*[-h -h], s*[t t],'color','w','lineWidth',0.5)
  %line(s*[+w +w], s*[-h +h], s*[t t],'color','w','lineWidth',0.5)
  %line(s*[+w -w], s*[+h +h], s*[-t -t],'color','w','lineWidth',0.5)
  %line(s*[-w -w], s*[+h -h], s*[-t -t],'color','w','lineWidth',0.5)
  %line(s*[-w +w], s*[-h -h], s*[-t -t],'color','w','lineWidth',0.5)
  %line(s*[+w +w], s*[-h +h], s*[-t -t],'color','w','lineWidth',0.5)
end

%--------------------------------------------------------------------------
% Draws the particle trajectory in the glass.
%--------------------------------------------------------------------------
function DrawParticle(x0,y0,z0,xf,yf,zf)
  line([x0 xf], [y0 yf], [z0 zf],'color','k','lineWidth',2);
end

%--------------------------------------------------------------------------
% Function: wavelength2color
% Author: Urs Hofmann
% Mail: hofmannu@biomed.ee.ethz.ch
% Date: 08.07.2020
% Version: 1.0
% Description: converts a wavelength in nm into a color
% input arguments:
%   - maxIntensity: maximum intensity of colorspace
%   - gammaVal 
%   - used colorSpace (either 'rgb', or 'hsv')
% example usage
%  wavelength2color(532, 'gammaVal', 1, 'maxIntensity', 255, 'colorSpace', 'rgb')
% principle stolen from:
%    https://academo.org/demos/wavelength-to-colour-relationship/
function colorCode = wavelength2color(wavelength, varargin)
  % default arguments
  maxIntensity = 1;
  gammaVal = 0.8;
  colorSpace = 'rgb';
  for iargin=1:2:(nargin-1)
    switch varargin{iargin}
      case 'maxIntensity' 
        maxIntensity = varargin{iargin + 1};
      case 'gammaVal'
        gammaVal = varargin{iargin + 1};
      case 'colorSpace'
        switch varargin{iargin + 1}
          case 'rgb'
            colorSpace = 'rgb';
          case 'hsv'
            colorSpace = 'hsv';
          otherwise
            error('Invalid colorspace defined');
        end
      otherwise
        error('Invalid argument passed');
    end
  end
	function outputVal = adjust(inputVal, factor)
		if (inputVal == 0)
	  	outputVal = 0;
	  else
			outputVal = (inputVal * factor)^gammaVal;
	  end
	end
	if (wavelength >= 380) && (wavelength < 440)
		r = -(wavelength - 440) / (440 - 380);
    g = 0;
    b = 1;
	elseif (wavelength >= 440) && (wavelength < 490)
		r = 0;
 		g = (wavelength - 440) / (490 - 440);
    b = 1;
  elseif (wavelength >= 490) && (wavelength < 510)
  	r = 0;
    g = 1;
    b = -(wavelength - 510) / (510 - 490);
  elseif (wavelength >= 510) && (wavelength < 580)
  	r = (wavelength - 510) / (580 - 510);
    g = 1;
    b = 0;
  elseif (wavelength >= 580) && (wavelength < 645)
    r = 1;
    g = -(wavelength - 645) / (645 - 580);
    b = 0;
  elseif (wavelength >= 645) && (wavelength < 780)
    r = 1;
    g = 0;
    b = 0;
  else
  	r = 0;
    g = 0;
    b = 0;
  end
    
  if (wavelength >= 380) && (wavelength < 420)
  	factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);
  elseif (wavelength >=  420) && (wavelength < 700)
    factor = 1;
  elseif (wavelength >= 700) && (wavelength < 780)
  	factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
  else
    factor = 0;
  end
  r = adjust(r, factor);
  g = adjust(g, factor);
  b = adjust(b, factor);
  rgbCode = [r, g, b];
  switch colorSpace
    case 'rgb'
      colorCode = rgbCode;
    case 'hsv'
      colorCode = rgb2hsv(rgbCode);
  end
  colorCode = colorCode * maxIntensity;
end
