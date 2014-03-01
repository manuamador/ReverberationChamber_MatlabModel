%ExampleFDvect.m sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%   	            Frequency Domain (vect)                     %
%                    			                        %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

global Lt c Rf POS f

tic
l=8.7;
p=3.7;
h=2.9;


c = 299792458;%
Lt = 1e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)

f=1e6:1e6:150e6; %frequency range

%%Q factor simulation
% mur=150000;
% mu0=4*pi*10^-7;
%mu=1%mur*mu0;
%sigma=0.05%10^7;

%delta=sqrt(2./(2*pi*f*mu*sigma));

%V=l*p*h;
%S=2*(l*p+p*h+l*h);
%L=4*V/S;

%Q=1./(1./(3/2*V/S./delta)+1./(16*pi^2*V*(f/c).^3));

%tau=Q./(2*pi*f);

%%Q facor from measurements Qm(fq)
% f=fq(1):1e6:fq(end);
% Q=interp1(fq,Qm,f);
%tau=Q./(2*pi*f);


%Loss coefficient 
%Rf=exp(-L./(2*c*tau));
Rf = 0.998*ones(1,length(f));




%Reception point rectangular coordinates
X_1 = 4.5;
Y_1 = 3;
Z_1 = 1.5;

[Ex,Ey,Ez] = FRvect(X_1,Y_1,Z_1);
toc
%CIR
figure(1)
subplot(3,1,1)
plot(f,20*log10(abs(Ex)))
title('E_x')
xlabel('f [Hz]')
ylabel('[dB/Hz]')

subplot(3,1,2)
plot(f,20*log10(abs(Ey)))
title('E_y')
xlabel('f [Hz]')
ylabel('[dB/Hz]')

subplot(3,1,3)
plot(f,20*log10(abs(Ez)))
title('E_z')
xlabel('f [Hz]')
ylabel('[dB/Hz]')

figure(2)
subplot(1,3,1)
plot(real(Ez),imag(Ex),'.')
title('E_x')
grid on
axis equal

subplot(1,3,2)
plot(real(Ez),imag(Ey),'.')
title('E_y')
grid on
axis equal

subplot(1,3,3)
plot(real(Ez),imag(Ez),'.')
title('E_z')
grid on
axis equal
