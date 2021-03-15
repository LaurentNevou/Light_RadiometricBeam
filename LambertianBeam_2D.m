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

rmax = 2;      Nr = 90; % xmax is the last point of the x vector and Nx is number of point of this vector [mm]
zmax = 2*rmax; Nz = 80; % same for z [mm]
Nphi = 35;

FOI  = 30;              % angle of the beam FOI

zcut = 1;               % 2D xy-plot at z=zcut
xscale=[-1 1]*rmax;     % ploting scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = logspace(-2,log10(zmax),Nz);
if isempty(find(z==zcut))
  z=sort([z zcut ]);
end

r = logspace(-3,log10(rmax),Nr);
r = [0 r];
[Z,R] = meshgrid(z,r);

O = linspace(-pi/2 , pi/2, 200);
p = log(0.5)./log(cos(FOI*pi/180));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Lambertian beam calculus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cosO = Z ./ sqrt( R.^2  + Z.^2 )    ;
cosO(cosO<0)=0;

trans = cosO.^3 ./ (Z.^2) ; %% it is a function that transform the Intensity (W/sr) in Irradiance (W/m2)

L = Power / pi *(p+1)/2 * cosO.^p .* trans ;

% the pi(p+1)/2 is a normalization constant so that the integral over the solide angle (half sphere) = the power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% indexes and values for plotting %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idz=find(z==zcut);

idxL=find(abs(-r-xscale(1))==min(abs(-r-xscale(1))) ) ;
idxR=find(abs(r-xscale(2))==min(abs(r-xscale(2))) ) ;

M=max(  [ max(L(:,end))  max(L(idxL,:))  max(L(idxR,:)) ] ) ;
m=max(L(:,idz));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[50 50 1400 1000],'color','w')
%figure('position',[-3500 300 1400 1000],'color','w')

FS=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1,'fontsize',FS)
hold on; grid on; box on;

plot( O*180/pi , cos(O).^p , 'r.-')

xlim([-1 1]*90)
set (gca , 'xtick', [-90:15:90]);
set (gca , 'ytick', 0:0.25:1);

xlabel('Angle (deg)')
ylabel('Amplitude (norm. u.)')
title(strcat('FOI=',num2str(FOI),'deg <=> cosO\^',num2str(p,'%.1f')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2,'fontsize',FS)
hold on; grid on; box on;

plot(r,L(:,idz)/max(L(:,idz)) , 'r.-' )
plot(-r,L(:,idz)/max(L(:,idz)) , 'r.-' )

xlim(xscale)
xlabel('r (mm)')
ylabel('Amplitude (norm. u.)')
title(strcat('@z=',num2str(z(idz)), 'mm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3,'fontsize',FS)
hold on; grid on; box on;

pcolor(+r,z,L')
pcolor(-r,z,L')
plot(xscale,[1 1]*z(idz),'r--')

v=sort(logspace(log10(M), log10(m),6));
contour(+r,z,squeeze(L)',v,'showtext','off','linecolor','w','linewidth',2)
contour(-r,z,squeeze(L)',v,'showtext','off','linecolor','w','linewidth',2)

colormap(jet)
colorbar
shading flat
xlim(xscale)
ylim([0 2*xscale(end)])
axis equal
m=10*Power*sqrt(p)/max(z)^2;
caxis([0 m])

xlabel('r (mm)')
ylabel('z (mm)')

title(strcat('FOI = ',num2str(FOI,'%.0f'),'deg'  ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=linspace(0,2*pi,Nphi);
[PHI,RR] = meshgrid(phi,r);
Lp= repmat(L(:,:,1),[1 1 Nphi]);

X=RR.*cos(PHI);
Y=RR.*sin(PHI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4,'fontsize',FS)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I just check that the integral (the power) is constant over z
% just it is in polar...

figure('position',[1400 100 450 500],'color','w')

subplot(1,1,1,'fontsize',FS)
hold on; grid on; box on;

s = squeeze( trapz(r,2*pi*R.*squeeze(L),1) ) ;
plot(z,s,'ro-')

xlabel('z (mm)')
ylabel('Power (W)')
title(strcat('Total integral, P=', num2str(Power),'W'))
