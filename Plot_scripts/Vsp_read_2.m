function Vp= Vsp_read_2(nx,ny,nz,snap_file3)

%nx=140;ny=120;nz=47;%ny=120; nx*nz=6580

%snap_file3='vp3dfile.p';
fid3=fopen(snap_file3,'r','ieee-le');
[sek3,count3]=fread(fid3,'float');

fclose(fid3);

%nv=count3/(nx*ny*nz); 

num=0;
Vp=zeros(nx,nz,ny);
%for kk=1:nv
for i=1:ny
    for j=1:nz
        for k=1:nx
            num=num+1;
            Vp(k,j,i)=sek3(num);
        end
    end
end

%end