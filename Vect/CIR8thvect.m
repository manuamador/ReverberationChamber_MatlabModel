
%CIR8thvect.m computes the CIR of a given 8th of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%        CHANNEL IMPULSE RESPONSE Eighth FUNCTION V 2.0         %
%                          (vect)                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sx8th,Sy8th,Sz8th,dist,azim,elev]=CIR8thvect(X_1,Y_1,Z_1)

global Lt c R POS

DX = X_1-POS(:,1);
DY = Y_1-POS(:,2);
DZ = Z_1-POS(:,3);

dist = sqrt(DX.^2+DY.^2+DZ.^2);

ca    = cos(POS(:,5));
sa    = sin(POS(:,5));
cb    = cos(POS(:,6));
sb    = sin(POS(:,6));

distx = ((-sb).^2+(1-(-sb).^2).*ca).*DX+(-sb.*cb.*(1-ca)).*DY+(cb.*sa).*DZ;
disty = (-sb.*cb.*(1-ca)).*DX+((cb).^2+(1-cb.^2).*ca).*DY+(sb.*sa).*DZ;
distz = (-cb.*sa).*DX+(-sb.*sa).*DY+ca.*DZ;
DXY=sqrt(DX.^2+DY.^2);
azim=(acos((POS(:,1)-X_1)./DXY)).*sign(asin((POS(:,2)-Y_1)./DXY));
elev=asin((POS(:,3)-Z_1)./dist);
clear DX DY DZ DXY

% distxy = sqrt(distx.^2+disty.^2);
% azim=(acos((POS(:,1)-X_1)./distxy)).*sign(asin((POS(:,2)-Y_1)./distxy));
% elev=asin((POS(:,3)-Z_1)./dist);

distxy = sqrt(distx.^2+disty.^2);

costheta = distz./dist;
sintheta = distxy./dist;
cosphi   = distx./distxy;
sinphi   = disty./distxy;
L = (R.^POS(:,4))./dist; %Attenuation
clear distx disty distz distxy

%Projection in the usual rectangular coordinates
Sx8th =L.*((((-sb).^2+(1-(-sb).^2).*ca).*(-sintheta.*costheta.*cosphi)+(-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*sinphi)+(-cb.*sa).*(-sintheta.*(-sintheta))));
Sy8th = L.*(((-sb.*cb.*(1-ca)).*(-sintheta.*costheta.*cosphi)+((cb).^2+(1-(cb).^2).*ca).*(-sintheta.*costheta.*sinphi)+(-sb.*sa).*(-sintheta.*(-sintheta))));
Sz8th = L.*(((cb.*sa).*(-sintheta.*costheta.*cosphi)+(sb.*sa).*(-sintheta.*costheta.*sinphi)+ca.*(-sintheta.*(-sintheta))));

clear L ca cb sa sb costheta sintheta cosphi sinphi 

