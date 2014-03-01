%ImageCreator.m standalone script that generates the images for a l,p,h
%rectangular cavity

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
tilt = [pi/2-acos(sqrt(2/3))];
azimut = [pi/4];

%Source creation main parameters
ordre = round(dmax/min([l;p;h]))+1; %Maximum order

Memo = 1/6*pi*(Ltp*c)^3/l/p/h*6*8; %memory usage
disp([num2str(round(Memo/1e6)),' MB needed'])

disp('Horizontal plane generation')

POS = [X Y Z zeros(length(X),1) tilt azimut];
POSp = POS;

for z=1:1:length(X)
    for i=1:ordre
        POSp = [POSp; 	2*i*l-X(z)	Y(z)	Z(z)    abs(2*i)-1		pi-tilt(z)  mod(2*pi-azimut(z),2*pi);
            			2*i*l+X(z)	Y(z)	Z(z)	abs(2*i)      	tilt(z)     mod(azimut(z),2*pi)];
    end
    
end

POSi = POSp;
POSi(:,2) = p-POSi(:,2);
POSi(:,5) = pi-POSi(:,5);
POSi(:,6) = pi-POSi(:,6);
POSP=[];

for j=0:2:ordre
    
    POS1 = [		POSp+[zeros(length(POSp),1) (j*p)*ones(length(POSp),1)  	zeros(length(POSp),1)       abs(j)*ones(length(POSp),1) zeros(length(POSp),1) zeros(length(POSp),1)]];
    POS1 = [POS1; 	POSi+[zeros(length(POSi),1) ((j+1)*p)*ones(length(POSi),1)  zeros(length(POSi),1)    	abs(j+1)*ones(length(POSi),1) zeros(length(POSi),1) zeros(length(POSi),1)]];
    
    POS1(:,7) = sqrt(POS1(:,1).^2+POS1(:,2).^2)./c; %compute the distance/c of each source
    U = find(POS1(:,7)<Ltp); %find the sources within the chosen time window
    POS1(:,7) = [];
    POS1 = POS1(U,:); %keeps the sources in the time window
    POSP = [POSP; POS1];
    clear POS1
end
toc

%odd horizontal plane generation
POSI = POSP;
POSI(:,3) = h-POSI(:,3);
POSI(:,6) = mod(POSI(:,6)+pi,2*pi);%azimuth flipping

disp('Vertical duplication')

for k=0:2:ordre    
    POS1 = [		POSP+[zeros(length(POSP),1) zeros(length(POSP),1) (k*h)*ones(length(POSP),1)        abs(k)*ones(length(POSP),1) zeros(length(POSP),1) zeros(length(POSP),1)]];
    POS1 = [POS1;	POSI+[zeros(length(POSI),1) zeros(length(POSI),1) ((k+1)*h)*ones(length(POSI),1)    abs(k+1)*ones(length(POSI),1) zeros(length(POSI),1) zeros(length(POSI),1)]];
    
    POS1(:,7) = sqrt(POS1(:,1).^2+POS1(:,2).^2+POS1(:,3).^2)./c; %compute the distance/c of each source
    U = find(POS1(:,7)<Ltp); %find the sources within the chosen time window
    POS1(:,7) = [];
    POS1 = POS1(U,:); %keeps the sources in the time window
    POS = [POS; POS1];
    clear POS1
end
POS(1,:)=[];
disp('Saving...') %saving the POSITION matrix
filename = sprintf('%delem_%dns1s8.mat',length(X),round(Lt/(1e-9)));
save(filename,'POS')
toc
disp('Done.')
