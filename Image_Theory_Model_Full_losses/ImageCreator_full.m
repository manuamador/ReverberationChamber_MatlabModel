%ImageCreator_full.m standalone script that generates the images for a l,p,h
%rectangular cavity use this script if you need to discrimates losses from
%the direction x, direction y and direction z (a door can be opened in the cavity and
%losses in a direction can be greater).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%                    IMAGE CREATOR V 2.0                        %
%        by E. Amador (emmanuel.amador@insa-rennes.fr)          %
%                         IETR/DGA                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
tic
Lt = 1e-6 %Time-window length in seconds

%Physical dimensions of the cavity (length, width, heigth) in meters
l = 8.7;
p = 3.7;
h = 2.9;

c = 299792458; %speed of light
dmax = Lt*c %maximal distance for the choses time-window
Ltp = Lt+3*max([l p h])/c; %security margin in order to not truncate the sources too much

%Radiating element(s), use columns here for more than one radiating elements
% rectanguar coordinates
X = [1];
Y = [2];
Z = [1];

%angular orientation
tilt = 0%[pi/2-acos(sqrt(2/3))];
azimut = 0%=[pi/4];

%Source creation main parameters
ordre = round(dmax/min([l;p;h]))+1; %Maximum order

Memo = 4/3/8*pi*(Ltp*c)^3/l/p/h*8*8; %memory usage
disp([num2str(round(Memo/1e6)),' MB needed'])


disp('Horizontal plane generation')

POS = [X Y Z zeros(length(X),1) zeros(length(X),1) zeros(length(X),1) tilt azimut];
POSP = POS;

for z=1:1:length(X)
    for i=0
        for j=1:ordre
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP = [POSP; X(z)	2*j*p+Y(z)	Z(z)	0 abs(2*j)      0   tilt(z)     mod(azimut(z),2*pi);
                             +X(z)	2*j*p-Y(z)	Z(z)	0 abs(2*j)-1    0   pi-tilt(z)  mod(pi-azimut(z),2*pi)]; 
            end
        end
    end
    
end


for z=1:1:length(X)
    for i=1:ordre
        for j=0
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP = [POSP; 2*i*l-X(z)	Y(z)	Z(z)    abs(2*i)-1  0   0		pi-tilt(z)  mod(2*pi-azimut(z),2*pi);
                              2*i*l+X(z)	Y(z)	Z(z)	abs(2*i)    0   0      	tilt(z)     mod(azimut(z),2*pi)];
            end
        end
    end
    
end


for z=1:1:length(X)
    
    
    
    for i=1:ordre
        for j=1:ordre
            if ((i*l)^2+(j*p)^2)<(dmax+l+p)^2
                
                POSP = [POSP;   2*i*l-X(z)	2*j*p+Y(z)	Z(z)    abs(2*i)-1  abs(2*j)    0     pi-tilt(z)  mod(2*pi-azimut(z),2*pi);
                                2*i*l+X(z)	2*j*p+Y(z)	Z(z)	abs(2*i)    abs(2*j)    0     tilt(z)     mod(azimut(z),2*pi);
                                2*i*l-X(z)	2*j*p-Y(z)	Z(z)    abs(2*i)-1	abs(2*j)-1  0     tilt(z)     mod(pi+azimut(z),2*pi);
                                2*i*l+X(z)	2*j*p-Y(z)	Z(z)	abs(2*i)    abs(2*j)-1  0     pi-tilt(z)  mod(pi-azimut(z),2*pi)];
            end
        end
    end
    
end


%odd horizontal plane generation
POSI = POSP;
POSI(:,3) = h-POSI(:,3);
POSI(:,6) = mod(POSI(:,6)+pi,2*pi);%azimuth flipping

disp('Vertical duplication')

for k=0:ordre
    if mod(k,2)==0
        POS1 = [POSP+[zeros(length(POSP),1) zeros(length(POSP),1) (k*h)*ones(length(POSP),1)    zeros(length(POSP),1)   zeros(length(POSP),1)    abs(k)*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
    else
        POS1 = [POSI+[zeros(length(POSI),1) zeros(length(POSI),1) (k*h)*ones(length(POSI),1)    zeros(length(POSP),1)   zeros(length(POSP),1)    abs(k)*ones(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1)]];
    end
    
    POS1(:,9) = sqrt(POS1(:,1).^2+POS1(:,2).^2+POS1(:,3).^2)./c; %compute the distance/c of each source
    
    U = find(POS1(:,9)<Ltp); %find the sources within the chosen time window
    POS1(:,9) = [];
    POS1 = POS1(U,:); %keeps the sources in the time window
    POS = [POS; POS1];
    clear POS1
end


POS(1,:)=[];
disp('Saving...') %saving the POSITION matrix
filename = sprintf('%delem_%dns1s8_full.mat',length(X),round(Lt/(1e-9)));
save(filename,'POS')
toc
disp('Done.')
