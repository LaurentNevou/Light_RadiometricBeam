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

rmax = 5;      Nr = 50; % rmax is the last point of the r vector and Nr is number of point of this vector [mm]
zmax = 2*rmax; Nz = 81; % same for z [mm]

FOI  = 10;              % angle of the beam FOI
w0=0.1;                 % waist of the Gaussian beam
z0=2;                   % z position of the waist

zcut = 5;               % 2D xy-plot at z=zcut
xscale=[-1 1]*rmax;     % ploting scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = logspace(-3,log10(rmax),Nr);
r = sort([0 r]);

z = logspace(-2,log10(zmax),Nz);
if isempty(find(z==zcut))
  z=sort([z zcut ]);
end

[Z,R] = meshgrid(z,r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Gaussian beam calculus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll=pi*w0*tan(FOI*pi/180);
zr=pi*w0^2./ll;

w=w0*sqrt(1+((Z-z0)/zr).^2);
L = Power*2/pi/w0^2 * (w0./w).^2 .* exp(-2*R.^2./w.^2) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% indexes and values for plotting %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idz=find(z==zcut);

idxL=find(abs(r-xscale(1))==min(abs(r-xscale(1))) ) ;
idxR=find(abs(r-xscale(2))==min(abs(r-xscale(2))) ) ;

M=max(  [ max(L(:,end))  max(L(idxL,:))  max(L(idxR,:)) ] ) ;
m=max(L(:,idz));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[50 50 1400 1000],'color','w')
%figure('position',[-3500 300 1400 1000],'color','w')

FS=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1,'fontsize',FS)
hold on; grid on; box on;

pcolor(+r,z,squeeze(L)' )
pcolor(-r,z,squeeze(L)' )
plot(xscale,[1 1]*z(idz),'r--')

v=sort(logspace(log10(M), log10(m),3));
contour(+r,z,squeeze(L)',v,'showtext','off','linecolor','m','linewidth',1)
contour(-r,z,squeeze(L)',v,'showtext','off','linecolor','m','linewidth',1)
plot(+squeeze(w(1,:))',z,'w','linewidth',2)
plot(-squeeze(w(1,:))',z,'w','linewidth',2)

colormap(jet)
colorbar
shading flat
xlim(xscale)
ylim([0 2*xscale(end)])
m=150*Power/max(z)^2;
caxis([0 m])

xlabel('r (mm)')
ylabel('z (mm)')

title(strcat('FOI = ',num2str(FOI,'%.0f'),'deg @y=0mm' ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphi=50;
phi=linspace(0,2*pi,Nphi);
[PHI,RR] = meshgrid(phi,r);
Lp= repmat(L(:,:,1),[1 1 Nphi]);

X=RR.*cos(PHI);
Y=RR.*sin(PHI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2,'fontsize',FS)
hold on; grid on; box on;

pcolor(X,Y,squeeze(Lp(:,idz,:)  ))

colormap(jet); colorbar; shading flat
xlim(xscale)
ylim(xscale)
axis equal
xlabel('r (mm)')
ylabel('r (mm)')
title(strcat('@z=',num2str(z(idz)),'mm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3,'fontsize',FS)
hold on; grid on; box on;

plot(+r,L(:,idz)/max(L(:,idz)) , 'r.-' )
plot(-r,L(:,idz)/max(L(:,idz)) , 'r.-' )

xlim(xscale)
xlabel('r (mm)')
ylabel('Amplitude (norm. u.)')
title(strcat('@z=',num2str(z(idz)), 'mm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I just check that the integral (the power) is constant over z

s = squeeze( trapz(r,2*pi*R.*squeeze(L),1) ) ;

subplot(2,2,4,'fontsize',FS)
hold on; grid on; box on;

plot(z,s,'ro-')

xlabel('z (mm)')
ylabel('Power (W)')
title(strcat('Total integral, P=', num2str(Power),'W'))
