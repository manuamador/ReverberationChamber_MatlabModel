%ExampleTDvect.m sample program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%          CHANNEL IMPULSE RESPONSE, ANGLES of ARRIVAL          %
%                          (vect)                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

global Lt c R POS

tic
c = 299792458;%
Lt = 1e-6; %Time-window length in seconds
nbre_elements = 1; %number of radiating elements
dmax = c*Lt; %maximal distance

%loading the Position matrix from the image generator
filename = sprintf('%delem_%dns1s8.mat',nbre_elements,round(Lt/(1e-9))); 
load(filename)


%Loss coefficient 
R = 0.998;

%Reception point rectangular coordinates
X_1 = 4.5;
Y_1 = 3;
Z_1 = 1.5;

[Sx,Sy,Sz,t,azim,elev] = CIRvect(X_1,Y_1,Z_1);

%CIR
figure(1)
subplot(3,1,1)
plot(t,Sx,'.','MarkerSize',2)
xlim([0 Lt])
title('E_x')
xlabel('time [s]')
ylabel('[V/m]')

subplot(3,1,2)
plot(t,Sy,'.','MarkerSize',2)
xlim([0 Lt])
title('E_y')
xlabel('time [s]')
ylabel('[V/m]')
subplot(3,1,3)
plot(t,Sz,'.','MarkerSize',2)
xlim([0 Lt])
title('E_z')
xlabel('time [s]')
ylabel('[V/m]')


figure(2)
subplot(2,1,1)
plot(t,azim,'.','MarkerSize',2)
xlim([0 Lt])
title('azimuth')
xlabel('time [s]')
ylabel('[rad]')

subplot(2,1,2)
plot(t,elev,'.','MarkerSize',2)
xlim([0 Lt])
title('elevation')
xlabel('time [s]')
ylabel('[rad]')

