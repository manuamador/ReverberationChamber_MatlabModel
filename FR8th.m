
%FR8th.m computes the FR of a given 8th of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%           Frequency response Eighth FUNCTION V 2.0            %
%                                                               %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ex8th,Ey8th,Ez8th]=FR8th(X_1,Y_1,Z_1)

global f c Rf POS va

Ex8th = zeros(1,length(f));
Ey8th = zeros(1,length(f));
Ez8th = zeros(1,length(f));

for j=1:1:length(POS) %loop over the image-currents... a vector version exists and will be released soon
    
    DX = X_1-POS(j,1);
    DY = Y_1-POS(j,2);
    DZ = Z_1-POS(j,3);
    
    dist = sqrt(DX^2+DY^2+DZ^2);
    phase = 2*pi*dist*f/c;
   
  
        alpha = POS(j,5);
        beta  = POS(j,6);
        ca    = cos(alpha);
        sa    = sin(alpha);
        cb    = cos(beta);
        sb    = sin(beta);
        
        %Ralpha/beta, coordinates transformation matrix:
        %                           Ralpbeta=[(-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha) -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
        %                                     -sin(beta)*cos(beta)*(1-cos(alpha))   (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
        %                                     -cos(beta)*sin(alpha)                 	-sin(beta)*sin(alpha)                        	cos(alpha)];
        
        %rectangular coordinates calculation in the local system attached to the considered current (developped expressions, a matrix product is slower)
        
        distx = ((-sb)^2+(1-(-sb)^2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ;
        disty = (-sb*cb*(1-ca))*DX+((cb)^2+(1-cb^2)*ca)*DY+(sb*sa)*DZ;
        distz = (-cb*sa)*DX+(-sb*sa)*DY+(ca)*DZ;
        
        distxy = sqrt(distx^2+disty^2);
        
        LL = Rf.^POS(j,4); %Attenuation
        E = 1/dist;%*dist.^va; %Free-space attenuation add *dist.^va to include volumetric absorption by the air
        
        costheta = distz/dist;
        sintheta = distxy/dist;
        cosphi   = distx/distxy;
        sinphi   = disty/distxy;
        Antth    = -sintheta; %dipole radiation pattern
        
        %Reverse transformation matrix
        %                         Ralpbetainv=[((-sb)^2+(1-(-sb)^2)*ca)	(-sb*cb*(1-ca))			(-cb*sa);
        %                                     (-sb*cb*(1-ca))			((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
        %                                     (cb*sa) 					(sb*sa) 					 ca];
        
      
        %Projection in the usual rectangular coordinates
        
        Vx = (((-sb)^2+(1-(-sb)^2)*ca)*(Antth*costheta*cosphi)+(-sb*cb*(1-ca))*(Antth*costheta*sinphi)+(-cb*sa)*(-sintheta*Antth))*exp(1i*phase);
        Vy = ((-sb*cb*(1-ca))*(Antth*costheta*cosphi)+((cb)^2+(1-(cb)^2)*ca)*(Antth*costheta*sinphi)+(-sb*sa)*(-sintheta*Antth))*exp(1i*phase);
        Vz = ((cb*sa)*(Antth*costheta*cosphi)+(sb*sa)*(Antth*costheta*sinphi)+ca*(-sintheta*Antth))*exp(1i*phase);
        
        %Three axis channel impulse response construction
        Ex8th = Ex8th+LL.*E.*Vx;
        Ey8th = Ey8th+LL.*E.*Vy;
        Ez8th = Ez8th+LL.*E.*Vz;
 
end
