close all
[NZcoastLat,NZcoastLong]=NZCoastlineData

% plot(NZcoastLong,NZcoastLat)
% xlim([171.5 173.5])
% ylim([-44.25 -43])
% hold on;
% Srn = readtable('SOURCE_checker.txt');
% Stn = readtable('STATION_checker_25stat.txt');

Srn = readtable('SOURCE_SRF.txt');
% Srn = readtable('SOURCE_SRF.27s.txt');
Stn = readtable('STATION_dh_4km.txt');
% Stn = readtable('STATION_dh_4km.extended.txt');
% Stn = readtable('STATION_dh_4km_strong_motion.txt');
% load Cij.mat;

Rec = 4e3*[Stn.Var1 Stn.Var2 -Stn.Var3];
Src = 4e3*[Srn.Var1 Srn.Var2 -Srn.Var3];

nx=88;ny=88;nz=60;dx=4;
dxx=(175.167000-170.92780)/((nx-1)*4e3);
dyy=(-40.509468-(-43.651422))/((ny-1)*4e3);

xx = 0:4e3:(nx-1)*4e3;
yy = 0:4e3:(ny-1)*4e3;

xx=xx*(175.167000-170.92780)/((nx-1)*4e3)+170.92780;
yy=yy*(-40.509468-(-43.651422))/((ny-1)*4e3)+(-43.651422);

Rec(:,1) = Rec(:,1)*(175.167000-170.92780)/((nx-1)*4e3)+170.92780;
Rec(:,2) = Rec(:,2)*(-40.509468-(-43.651422))/((ny-1)*4e3)+(-43.651422);

Src(:,1) = Src(:,1)*(175.167000-170.92780)/((nx-1)*4e3)+170.92780;
Src(:,2) = Src(:,2)*(-40.509468-(-43.651422))/((ny-1)*4e3)+(-43.651422);

n1 = length(xx);    %X resolution
n2 = length(yy);   %Y resolution

% Plot Density
[X,Y] = meshgrid(xx,yy);

Rec(:,2) = Y(end)-Rec(:,2)+Y(1);
Src(:,2) = Y(end)-Src(:,2)+Y(1);
%Velocity profile
str2='z0=2';eval(str2);
str1='x0=44';
str2='z0=2';
str3='y0=44';
% % % 
% str1='x0=100';
% str2='z0=5';text(Rec(4,1)+4e3,Rec(4,2),'MQZ')
% str3='y0=100';

eval(str1);
eval(str2);
eval(str3);

str1=strcat('x0=',num2str(x0*dx),'[km]');
str3=strcat('y0=',num2str(y0*dx),'[km]');
str2=strcat('z0=',num2str(z0*dx),'[km]');

% snap_file1='NewVM_20200207_4KM/3366146/vs3dfile.s';
% snap_file2='NewVM_20200207_4KM/3366146/vp3dfile.p';
% snap_file1='NewVM_20200207_4KM/vs3dfile_smooth_4km.s';
% snap_file2='NewVM_20200207_4KM/vp3dfile_smooth_4km.p';
snap_file1='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/NewVM_20200207_4KM/NZVM_2020/vs3dfile_2020.s';
snap_file2='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/NewVM_20200207_4KM/NZVM_2020/vp3dfile_2020.p';
Vs0= Vsp_read_2(nx,ny,nz,snap_file1);
Vp0= Vsp_read_2(nx,ny,nz,snap_file2);

% snap_file1='NewVM_20200207_4KM/vs3dfile_true_square_2x_G2x4_no_smooth.s';
% snap_file2='NewVM_20200207_4KM/vp3dfile_true_square_2x_G2x4_no_smooth.p';
% snap_file1='NewVM_20200207_4KM/vs3dfile_true_square_4x_G5x4.s';
% snap_file2='NewVM_20200207_4KM/vp3dfile_true_square_4x_G5x4.p';

% snap_file1='NewVM_20200207_4KM/vs3dfile_iter_5.s';
% snap_file2='NewVM_20200207_4KM/vp3dfile_iter_5.p';

% filepath='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD/Dump_25_06_20_15iters_16s_checker_4x_G515_cc09_bfgs_dt_adj_default_continue_iter5/';
% filepath='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/NewVM_20200207_4KM/'
% filepath='INVERSION_GOOD/Dump_28_05_20_5iters_27s_BB_data_G212_bfgs_dt_vs_am_adj_ntape5_dtts5_cc09_Tmax_5s_kn_new_2nd_run_good/';
% filepath='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/NewVM_20200207_4KM/NZVM_2020/2nd_RUN/'
filepath='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/NewVM_20200207_4KM/NZVM_2020/'
% snap_file1=strcat(filepath,'vs3dfile_iter_10.s');
% snap_file2=strcat(filepath,'vp3dfile_iter_10.p');
snap_file1=strcat(filepath,'vs3dfile_iter_9.s');
snap_file2=strcat(filepath,'vp3dfile_iter_9.p');

% snap_file1='NewVM_20200207_4KM/NZVM_2020/vs3dfile_2020.s';
% snap_file2='NewVM_20200207_4KM/NZVM_2020/vp3dfile_2020.p';

Vs= Vsp_read_2(nx,ny,nz,snap_file1);
Vp= Vsp_read_2(nx,ny,nz,snap_file2);
% Vsc_min=2.5; Vsc_max=3.5;
% Vsc_min=2; Vsc_max=4.5;
% Vsmean_arr =[3.1863, 3.5786, 3.8046] - mean Vs_initial at 4km, 12km, 20km
% Vsmean_arr =[3.1863, 3.5786, 3.8046];
% ss = [15, 32, 16, 29, 14];
ss = 1:13;
% ss = 1:27;

LL = squeeze( log(Vs(:,z0,:)./Vs0(:,z0,:)));
LL = flip(LL,2);
% Vsc_mean=mean(LL(:));
% Vsc_mean = Vsmean_arr(3);
% disp(Vsc_mean)% 
% Vsc_min=Vsc_mean*0.8; Vsc_max=Vsc_mean*1.2;

fig = figure;
% set(gcf,'Position',[100 100 600 500])
set(gcf,'Position',[100 100 500 600])
% ColorMap  = loadcmap('WhBlGrYeRe.c3g');
% colormap(polarmap(64));
colormap(flipud(polarmap(64)));
% colormap(jet(64));
% colormap(flipud(jet(64)));
% surf(X,Y,LL'), shading interp , view (2), hold on
surf(X,Y,LL'), shading interp , view (2), hold on;

plot3(NZcoastLong,NZcoastLat,1000*ones(size(NZcoastLat)),'k','linewidth',1);
plot3(Rec(:,1),Rec(:,2),-Rec(:,3),'v','markersize',8,'MarkerEdgeColor','k',...
      'MarkerFaceColor',[0, 0.4470, 0.7410],'linewidth',1), hold on
plot3(Src(ss,1),Src(ss,2),-Src(ss,3),'p','markersize',10,'MarkerEdgeColor','k',...
      'MarkerFaceColor',[0.6350 0.0780 0.1840],'linewidth',1), %axis equal
% plot3(Rec(3,1)*ones(size(Y(1:end))),Y(1:end),1000*ones(size(Y(1:end))),'--r','linewidth',1.5);  
% text(Rec(4,1),Rec(4,2),'MQZ')
ylim([Y(1) Y(end)]), xlim([X(1) X(end)])
% colorbar('vertical');
C = colorbar('location','SouthOutside');
% set(get(C,'XLabel'),'String',strcat('Vs [km/s] (',num2str(double(Vsc_mean), 2),' \pm 20%)'),'fontsize',14)
% set(gca,'CLim',[Vsc_min Vsc_max])
set(gca,'CLim',[-0.2 0.2])

% xlabel('Longitude ','fontsize',14)
% ylabel('Latitude','fontsize',14)

text(173.3,-43.5,1000,strcat('z = ',{' ' }, num2str((z0-1)*dx),' km'),'fontweight','bold','fontsize',14)
% text(173.3,-40.6,1000,strcat('x=',num2str((x0)*dx),'km'), 'Color','r','fontweight','bold','fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{true}/m_{00})'),'fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{00:NZVM2020}/m_{00:NZVM2010})'),'fontsize',14)
set(get(C,'XLabel'),'String',strcat('ln(m_{09}/m_{00})'),'fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{10}/m_{00})'),'fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{12}/m_{00})'),'fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{15}/m_{00})'),'fontsize',14)

% set(get(C,'XLabel'),'String',strcat('ln(m_true/m0)'),'fontsize',14)
% set(get(C,'XLabel'),'String',strcat('ln(m_{true}/m0)'),'fontsize',14)
set(C,'FontSize',14);

text(170.8,-40,'(a)','fontsize',12);

% xlabel(strcat('Vs [km/s] ),'km'),'fontsize',14)
% title(strcat('m10 at z=',num2str((z0-1)*dx)),'fontsize',14)

% hold on
% for i=1:length(ss) 
%     text(Src(ss(i),1)+2*dxx,Src(ss(i),2)+2*dyy,num2str(ss(i)));
% end

% for i=1:size(Stn,1)   
%     text(Rec(i,1)+2*dxx,Rec(i,2)+2*dyy,Stn.Var4(i));
% end

saveas(fig,'V4_real.jpg');