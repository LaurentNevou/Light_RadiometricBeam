%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% last update 15 March 2021, LNEV %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Power=1;                % light power [Watt]

xmax = 5;      Nx = 50; % xmax is the last point of the x vector and Nx is number of point of this vector [mm]
ymax = 4;      Ny = 40; % ymax is the last point of the y vector and Nx is number of point of this vector [mm]
zmax = 2*xmax; Nz = 81; % same for z [mm]

FOI  = 10;              % angle of the beam FOI
w0=0.1;                 % waist of the Gaussian beam
z0=2;                   % z position of the waist

zcut = 5;               % 2D xy-plot at z=zcut
xscale=[-1 1]*xmax;     % ploting scale
tilt=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = logspace(-3,log10(xmax),Nx);
x = sort([-x 0 x]);

y = logspace(-3,log10(ymax),Ny);
y = sort([-y 0 y]);

z = logspace(-2,log10(zmax),Nz);
if isempty(find(z==zcut))
  z=sort([z zcut ]);
end

[Y,XX,ZZ]=meshgrid(y,x,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tilt=tilt*pi/180;
X =  XX*cos(tilt) + ZZ*sin(tilt);
Z = -XX*sin(tilt) + ZZ*cos(tilt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Gaussian beam calculus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll=pi*w0*tan(FOI*pi/180);
zr=pi*w0^2./ll;

w=w0*sqrt(1+((Z-z0)/zr).^2);
L = Power*2/pi/w0^2 * (w0./w).^2 .* exp(-2*(X.^2 + Y.^2)./w.^2) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% indexes and values for plotting %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idz=find(z==zcut);
idy=find(y==0);

idxL=find(abs(x-xscale(1))==min(abs(x-xscale(1))) ) ;
idxR=find(abs(x-xscale(2))==min(abs(x-xscale(2))) ) ;

M=max(  [ max(L(:,idy,end))  max(L(idxL,idy,:))  max(L(idxR,idy,:)) ] ) ;
m=max(L(:,idy,idz));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('position',[50 50 1400 1000],'color','w')
%figure('position',[-3500 300 1400 1000],'color','w')

FS=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1,'fontsize',FS)
hold on; grid on; box on;

pcolor(x,z,squeeze(L(:,idy,:))' )
plot(xscale,[1 1]*z(idz),'r--')

v=sort(logspace(log10(M), log10(m),3));
contour(x,z,squeeze(L(:,idy,:))',v,'showtext','off','linecolor','m','linewidth',1)
plot(+squeeze(w(1,1,:))',z,'w','linewidth',2)
plot(-squeeze(w(1,1,:))',z,'w','linewidth',2)

colormap(jet)
colorbar
shading flat
xlim(xscale)
ylim([0 2*xscale(end)])
caxis([0 5*M])

xlabel('x (mm)')
ylabel('z (mm)')

title(strcat('FOI = ',num2str(FOI,'%.0f'),'deg @y=',num2str(y(idy)),'mm' ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2,'fontsize',FS)
hold on; grid on; box on;

pcolor(x,y,squeeze(L(:,:,idz))')

v=sort(logspace(log10(M), log10(m),3)) ;
contour(x,y,squeeze(L(:,:,idz))',v,'showtext','off','linecolor','w','linewidth',2)


colormap(jet)
colorbar
shading flat
xlim(xscale)
ylim(xscale)
xlabel('x (mm)')
ylabel('y (mm)')
title(strcat('FOI = ',num2str(FOI,'%.0f'),'deg @z=',num2str(z(idz)),'mm' ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3,'fontsize',FS)
hold on; grid on; box on;

plot(x,L(:,idy,idz)/max(L(:,idy,idz)) , 'r.-' )

xlim(xscale)
xlabel('x (mm)')
ylabel('Amplitude (norm. u.)')
title(strcat('@z=',num2str(z(idz)), 'mm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I just check that the integral (the power) is constant over z

s=squeeze( trapz(y,trapz(x,L,1),2) );

subplot(2,2,4,'fontsize',FS)
hold on; grid on; box on;

plot(z,s,'ro-')

xlabel('z (mm)')
ylabel('Power (W)')
title(strcat('Total integral, P=', num2str(Power),'W'))
