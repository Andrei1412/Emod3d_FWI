clear all
close all
str1='shot=6';
eval(str1);
% matrix_file1=strcat('GS_shot',num2str(shot),'.txt');
% matrix_file2=strcat('GP_shot',num2str(shot),'.txt');

% matrix_file1='KS.txt';
% matrix_file2='KP.txt';

% matrix_file1='Tp_rec.txt';
% matrix_file2='Tp_xyz.txt';
% 
matrix_file1='Grads_iter_1.txt';
matrix_file2='Gradp_iter_1.txt';
% % % % % % %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
% matrix_file1='Gra_S.txt';
% matrix_file2='Gra_P.txt';
% % % % % 
% matrix_file1='Results/Syn_VM_B_2layers/CASE1/Dump_26_02/Gra_K.txt.it1';
% matrix_file2='Results/Syn_VM_B_2layers/CASE1/Dump_26_02/Gra_M.txt.it1';

fid1 = fopen(matrix_file1,'r');
matrix_dummy1 = fscanf(fid1,'%f');

fid2 = fopen(matrix_file2,'r');
matrix_dummy2 = fscanf(fid2,'%f');

% % % % % max1=0.01*max(abs(matrix_dummy1));
% max1=0.05*abs(min(matrix_dummy1));
% for i=1:length(matrix_dummy1)
%     if matrix_dummy1(i)>max1
%         matrix_dummy1(i)=max1;
%     end
%     if matrix_dummy1(i)<-max1
%         matrix_dummy1(i)=-max1;
%     end    
% end
% 
% % max2=0.01*max(abs(matrix_dummy2));
% max2=0.05*abs(min(matrix_dummy2));
% for i=1:length(matrix_dummy2)
%     if matrix_dummy2(i)>max2
%         matrix_dummy2(i)=max2;
%     end
%     if matrix_dummy2(i)<-max2
%         matrix_dummy2(i)=-max2;
%     end    
% end
%  nx=269;ny=269;nz=68;
nx=88;ny=88;nz=60;
% nx=120;ny=120;nz=60;
%  nx=90;ny=120;nz=60;


nts=5000;
% nt=floor(nts/20)+1;
nt=250;

dx=2.0;dy=dx;dz=dx;
x=0:dx:(nx-1)*dx;
y=0:dy:(ny-1)*dy;
z=0:dz:(nz-1)*dz;
% snap_file1='synth_homogeneous/Velocity_Model/vs3dfile_h.s';
% snap_file2='synth_homogeneous/Velocity_Model/vp3dfile_h.p';
% snap_file3='synth_homogeneous/Velocity_Model/rho3dfile.d';
% % 
% Vs= Vsp_read_2(nx,ny,nz,snap_file1);
% Vp= Vsp_read_2(nx,ny,nz,snap_file2);
% rho= Vsp_read_2(nx,ny,nz,snap_file3);
Vs=ones(nx,nz,ny)*2.0;
Vp=ones(nx,nz,ny)*4.0;
rho=2.7*ones(nx,nz,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_Vs=reshape(matrix_dummy1,[nx nz ny]);
G_Vp=reshape(matrix_dummy2,[nx nz ny]);

% G_Vs(:,1:5,:)=0;
% G_Vp(:,1:5,:)=0;

G_Vp_max=max((G_Vp(:)));
G_Vs_max=max((G_Vs(:)));

% % 
% str1='x0=50';
% str2='z0=10';
% str3='y0=50';

str1='x0=44';
str2='z0=2';
str3='y0=44';

% str1='x0=200';
% str2='z0=25';
% str3='y0=200';
% 
% str1='x0=100';
% str2='z0=25';
% str3='y0=100';

eval(str1);
eval(str2);
eval(str3);

Vsc_min=-1e-11 ;
Vsc_max=-Vsc_min;

% figure(5)
% subplot(211)
% % imagesc(squeeze(Vs (:,:,c1)))
% imagesc(y,z,squeeze(G_Vp(x0,:,:)));axis xy;zoom on;
% title('Gra__Vs')
% colormap(jet(64));
% % colormap(flipud(jet(64)));
% colorbar('vertical');
% set(gca,'Ydir','reverse')
% % set(gca,'CLim',[Vsc_min Vsc_max])
% xlabel('y-axis [km]')
% ylabel('z-axis [km]')
% subplot(212)
% % imagesc(squeeze(Vp (:,:,c1)))  
% imagesc(y,z,squeeze(G_Vs(x0,:,:)));axis xy;zoom on;
% title('Gra__Vs')
% colormap(jet(64));
% % colormap(flipud(jet(64)));
% colorbar('vertical');
% set(gca,'Ydir','reverse')
% % set(gca,'CLim',[Vpc_min Vpc_max])
% xlabel('y-axis [km]')
% ylabel('z-axis [km]')
% xlabel(str1)

G_Vp_z=squeeze(G_Vp(:,z0,:));
G_Vs_z=squeeze(G_Vs(:,z0,:));
fig = figure(6);
set(gcf,'Position',[100 100 500 900])
subplot(211)
% imagesc(squeeze(Vs (:,:,c1)))
imagesc(x,y,G_Vp_z');axis xy;zoom on;
% imagesc(squeeze(G_Vp(:,z0,:)))
title('Gra__Vp')
colormap(jet(64));
% caxis([-1e-12, 1e-12 ] )
% colormap(flipud(jet(64)));
% colormap(flipud(polarmap(64)));
colorbar('vertical');
% set(gca,'Ydir','reverse')
% set(gca,'CLim',[Vsc_min Vsc_max])
xlabel('x-axis [km]')
ylabel('y-axis [km]')
subplot(212)
% imagesc(squeeze(Vp (:,:,c1)))  
imagesc(x,y,G_Vs_z');axis xy;zoom on;
% imagesc(squeeze(G_Vs(:,z0,:)))
title('Gra__Vs')
colormap(jet(64));
% caxis([-1e-12, 1e-12 ] )
% colormap(flipud(jet(64)));
% colormap(flipud(polarmap(64)));
colorbar('vertical');
% set(gca,'Ydir','reverse')
% set(gca,'CLim',[Vsc_min Vsc_max])
xlabel('x-axis [km]')
ylabel('y-axis [km]')
xlabel(str2)
% saveas(fig,strcat('kernel_shot_',num2str(shot),'.jpg'));
saveas(fig,strcat('kernel_shot_i.jpg'));

G_Vp_y=squeeze(G_Vp(:,:,y0));
G_Vs_y=squeeze(G_Vs(:,:,y0));
figure(7)
subplot(211)
% imagesc(squeeze(Vs (:,:,c1)))
imagesc(x,z,G_Vp_y');axis xy;zoom on;
title('Gra__Vp')
colormap(jet(64));
% colormap(flipud(jet(64)));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
% colormap(flipud(polarmap(64)));
colorbar('vertical');
set(gca,'Ydir','reverse')
% set(gca,'CLim',[Vsc_min Vsc_max])
xlabel('x-axis [km]')
ylabel('z-axis [km]')
subplot(212)
% imagesc(squeeze(Vp (:,:,c1)))  
imagesc(x,z,G_Vs_y');axis xy;zoom on;
title('Gra__Vs')
colormap(jet(64));
% colormap(flipud(jet(64)));
% colormap(flipud(polarmap(64)));
colorbar('vertical');
set(gca,'Ydir','reverse')
% set(gca,'CLim',[Vsc_min Vsc_max])
caxis([-1e-12, 1e-12 ] )
xlabel('x-axis [km]')
ylabel('z-axis [km]')
xlabel(str3)

% % return
% Plot_3D_nofixcolor_Gra
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            