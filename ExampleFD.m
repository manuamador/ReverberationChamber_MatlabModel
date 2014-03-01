%ExampleFD.m sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%                     FREQUENCY RESPONSE                        %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

global f c Rf POS va

%dimensions
l=8.7;
p=3.7;
h=2.9;

tic
c = 299792458;%
Lt = 1e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)

f=1e6:1e6:150e6;

%Loss coefficient 
Rf = 0.998*ones(1,length(f)); %the loss coefficient Rf can be a function of the frequency

% %if a measurement of the quality factor is available Q(fq), use:
% f=fq(1):1e6:fq(end);
% Qi=interp1(fq,Q,f);
% L=2*(l*p*h)/(lp+lh+ph); %reverberation distance
% Rf=exp(-L*f.*2*pi./(2*c*Qf));

va=0*ones(1,length(f)); %volumetric absorption due to air absorption


%Reception point rectangular coordinates
X_1 = 4.5;
Y_1 = 3;
Z_1 = 1.5;
tic
[Ex,Ey,Ez] = FR(X_1,Y_1,Z_1);
toc


%Frequency response figure
figure(1)
subplot(3,1,1)
plot(f/1e6,20*log10(abs(Ex)))
title('FFT_x')
xlim([0 max(f)/1e6])
grid on
xlabel('frequency in MHz')
ylabel('dB')

subplot(3,1,2)
plot(f/1e6,20*log10(abs(Ey)))
title('FFT_y')
xlim([0 max(f)/1e6])
grid on
xlabel('frequency in MHz')
ylabel('dB')

subplot(3,1,3)
plot(f/1e6,20*log10(abs(Ez)))
title('FFT_z')
xlim([0 max(f)/1e6])
grid on
xlabel('frequency in MHz')
ylabel('dB')
