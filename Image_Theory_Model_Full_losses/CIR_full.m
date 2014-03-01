%CIR_full.m Computes the channel impulse response at the rectangular coordinates X_1,Y_1,Z_1 for a given sources matrix POS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%            CHANNEL IMPULSE RESPONSE FUNCTION V 2.0            %
%                                                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Sx,Sy,Sz]=CIR(X_1,Y_1,Z_1)

global Lt c Rx Ry Rz N POS

%Channel Impulse Response computation, we compute 8 CIR corresponding to
%each eighth of the system 

Sx = zeros(1,N);
Sy = zeros(1,N);
Sz = zeros(1,N);

Sx1 = zeros(1,N);
Sy1 = zeros(1,N);
Sz1 = zeros(1,N);
Sx2 = zeros(1,N);
Sy2 = zeros(1,N);
Sz2 = zeros(1,N);
Sx3 = zeros(1,N);
Sy3 = zeros(1,N);
Sz3 = zeros(1,N);
Sx4 = zeros(1,N);
Sy4 = zeros(1,N);
Sz4 = zeros(1,N);
Sx5 = zeros(1,N);
Sy5 = zeros(1,N);
Sz5 = zeros(1,N);
Sx6 = zeros(1,N);
Sy6 = zeros(1,N);
Sz6 = zeros(1,N);
Sx7 = zeros(1,N);
Sy7 = zeros(1,N);
Sz7 = zeros(1,N);
Sx8 = zeros(1,N);
Sy8 = zeros(1,N);
Sz8 = zeros(1,N);



%1/4
%CIR 1/8, images in the space defined by x>=0, y>=0, z>=0
[Sx1,Sy1,Sz1] = CIR8th_full(X_1,Y_1,Z_1);

%CIR 2/8, images in the space defined by x>=0, y>=0, z<=0, we need to add
%one order to the images previously used, use negative values for their z
%coordinates and add pi to their azimut angle

POS(:,4) = POS(:,4)+1; %order+1 
POS(:,3) = -POS(:,3); 
POS(:,6) = POS(:,6)+pi;
[Sx2,Sy2,Sz2] = CIR8th_full(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3); %change again the signs of the z coordinates because the next eighth is in the half-space z>0


%2/4
%CIR 3/8, images in the space defined by x>=0, y<=0, z>=0
POS(:,2) = -POS(:,2); %change the signs of the y coordinates
POS(:,5) = pi-POS(:,5); %change the tilt angle
[Sx3,Sy3,Sz3] = CIR8th_full(X_1,Y_1,Z_1);

%CIR 4/8, images in the space defined by x>=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3); 
POS(:,6) = POS(:,6)+pi;
[Sx4,Sy4,Sz4] = CIR8th_full(X_1,Y_1,Z_1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3); 

%3/4
%CIR 5/8, images in the space defined by x<=0, y>=0, z>=0
POS(:,4) = POS(:,4)-1; %order-1
POS(:,1) = -POS(:,1);
[Sx5,Sy5,Sz5] = CIR8th_full(X_1,Y_1,Z_1);

%CIR 6/8, images in the space defined by x<=0, y>=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3); 
POS(:,6) = POS(:,6)+pi;
[Sx6,Sy6,Sz6] = CIR8th_full(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3); 
POS(:,5) = pi-POS(:,5);

%4/4
%CIR 7/8, images in the space defined by x<=0, y<=0, z>=0
POS(:,2) = -POS(:,2);
POS(:,6) = -POS(:,6);
[Sx7,Sy7,Sz7] = CIR8th_full(X_1,Y_1,Z_1);

%CIR 8/8, images in the space defined by x<=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3); 
POS(:,6) = POS(:,6)+pi;
[Sx8,Sy8,Sz8] = CIR8th_full(X_1,Y_1,Z_1);
POS(:,1) = -POS(:,1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3); 
POS(:,4) = POS(:,4)-3;
POS(:,6) = mod(POS(:,6),2*pi);

Sx = Sx1+Sx2+Sx3+Sx4+Sx5+Sx6+Sx7+Sx8;
Sy = Sy1+Sy2+Sy3+Sy4+Sy5+Sy6+Sy7+Sy8;
Sz = Sz1+Sz2+Sz3+Sz4+Sz5+Sz6+Sz7+Sz8;
