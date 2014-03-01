%CIRvect.m Computes the channel impulse response at the rectangular coordinates X_1,Y_1,Z_1 for a given sources matrix POS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%           CHANNEL IMPULSE RESPONSE FUNCTION V 2.0             %
%                          (vect)                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sx,Sy,Sz,t,azim,elev]=CIRvect(X_1,Y_1,Z_1)

global Lt c R POS

%Channel Impulse Response computation, we compute 8 CIR corresponding to
%each eighth of the system


%1/4
%CIR 1/8, images in the space defined by x>=0, y>=0, z>=0
[Sx1,Sy1,Sz1,dist1,azim1,elev1] = CIR8thvect(X_1,Y_1,Z_1);

%CIR 2/8, images in the space defined by x>=0, y>=0, z<=0, we need to add
%one order to the images previously used, use negative values for their z
%coordinates and add pi to their azimut angle

POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Sx2,Sy2,Sz2,dist2,azim2,elev2] = CIR8thvect(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3); %change again the signs of the z coordinates because the next eighth is in the half-space z>0

%2/4
%CIR 3/8, images in the space defined by x>=0, y<=0, z>=0
POS(:,2) = -POS(:,2); %change the signs of the y coordinates
POS(:,5) = pi-POS(:,5); %change the tilt angle
[Sx3,Sy3,Sz3,dist3,azim3,elev3] = CIR8thvect(X_1,Y_1,Z_1);

%CIR 4/8, images in the space defined by x>=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Sx4,Sy4,Sz4,dist4,azim4,elev4] = CIR8thvect(X_1,Y_1,Z_1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3);

%3/4
%CIR 5/8, images in the space defined by x<=0, y>=0, z>=0
POS(:,4) = POS(:,4)-1; %order-1
POS(:,1) = -POS(:,1);
[Sx5,Sy5,Sz5,dist5,azim5,elev5] = CIR8thvect(X_1,Y_1,Z_1);

%CIR 6/8, images in the space defined by x<=0, y>=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Sx6,Sy6,Sz6,dist6,azim6,elev6] = CIR8thvect(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3);
POS(:,5) = pi-POS(:,5);

%4/4
%CIR 7/8, images in the space defined by x<=0, y<=0, z>=0
POS(:,2) = -POS(:,2);
POS(:,6) = -POS(:,6);
[Sx7,Sy7,Sz7,dist7,azim7,elev7] = CIR8thvect(X_1,Y_1,Z_1);

%CIR 8/8, images in the space defined by x<=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Sx8,Sy8,Sz8,dist8,azim8,elev8] = CIR8thvect(X_1,Y_1,Z_1);
POS(:,1) = -POS(:,1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3);
POS(:,4) = POS(:,4)-3;
POS(:,6) = mod(POS(:,6),2*pi);

Sx =[Sx1; Sx2; Sx3; Sx4; Sx5; Sx6; Sx7; Sx8];
Sy =[Sy1; Sy2; Sy3; Sy4; Sy5; Sy6; Sy7; Sy8];
Sz =[Sz1; Sz2; Sz3; Sz4; Sz5; Sz6; Sz7; Sz8];
t =[dist1; dist2; dist3; dist4; dist5; dist6; dist7; dist8]/c;
azim=[azim1; azim2; azim3; azim4; azim5; azim6; azim7; azim8];
elev=[elev1; elev2; elev3; elev4; elev5; elev6; elev7; elev8];

