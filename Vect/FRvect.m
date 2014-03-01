%FRvect.m Computes the channel impulse response at the rectangular coordinates X_1,Y_1,Z_1 for a given sources matrix POS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%             FREQUENCY RESPONSE FUNCTION V 2.0                 %
%                          (vect)                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ex,Ey,Ez]=FRvect(X_1,Y_1,Z_1)

global Lt c Rf POS f

%Frequency Response computation, we compute 8 FR corresponding to
%each eighth of the system


%1/4
%FR 1/8, images in the space defined by x>=0, y>=0, z>=0
[Ex1,Ey1,Ez1] = FR8thvect(X_1,Y_1,Z_1);

%FR 2/8, images in the space defined by x>=0, y>=0, z<=0, we need to add
%one order to the images previously used, use negative values for their z
%coordinates and add pi to their azimut angle

POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Ex2,Ey2,Ez2] = FR8thvect(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3); %change again the signs of the z coordinates because the next eighth is in the half-space z>0

%2/4
%FR 3/8, images in the space defined by x>=0, y<=0, z>=0
POS(:,2) = -POS(:,2); %change the signs of the y coordinates
POS(:,5) = pi-POS(:,5); %change the tilt angle
[Ex3,Ey3,Ez3] = FR8thvect(X_1,Y_1,Z_1);

%FR 4/8, images in the space defined by x>=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Ex4,Ey4,Ez4] = FR8thvect(X_1,Y_1,Z_1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3);

%3/4
%FR 5/8, images in the space defined by x<=0, y>=0, z>=0
POS(:,4) = POS(:,4)-1; %order-1
POS(:,1) = -POS(:,1);
[Ex5,Ey5,Ez5] = FR8thvect(X_1,Y_1,Z_1);

%FR 6/8, images in the space defined by x<=0, y>=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Ex6,Ey6,Ez6] = FR8thvect(X_1,Y_1,Z_1);
POS(:,3) = -POS(:,3);
POS(:,5) = pi-POS(:,5);

%4/4
%FR 7/8, images in the space defined by x<=0, y<=0, z>=0
POS(:,2) = -POS(:,2);
POS(:,6) = -POS(:,6);
[Ex7,Ey7,Ez7] = FR8thvect(X_1,Y_1,Z_1);

%FR 8/8, images in the space defined by x<=0, y<=0, z<=0
POS(:,4) = POS(:,4)+1; %order+1
POS(:,3) = -POS(:,3);
POS(:,6) = POS(:,6)+pi;
[Ex8,Ey8,Ez8] = FR8thvect(X_1,Y_1,Z_1);
POS(:,1) = -POS(:,1);
POS(:,2) = -POS(:,2);
POS(:,3) = -POS(:,3);
POS(:,4) = POS(:,4)-3;
POS(:,6) = mod(POS(:,6),2*pi);

Ex =[Ex1+Ex2+Ex3+Ex4+Ex5+Ex6+Ex7+Ex8];
Ey =[Ey1+Ey2+Ey3+Ey4+Ey5+Ey6+Ey7+Ey8];
Ez =[Ez1+Ez2+Ez3+Ez4+Ez5+Ez6+Ez7+Ez8];

