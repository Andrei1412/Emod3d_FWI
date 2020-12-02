import numpy as np
import os

def read_srf_source(source_file):
    #read source file in cartesian coordinates
    with open(source_file, 'r') as f:
        lines = f.readlines()
    line0 = lines[0].split()
    nShot = int(line0[0])
    S = np.zeros((nShot, 3))
    sNames = []
    for i in range(1, nShot + 1):
        line_i = lines[i].split()
        S[i - 1, 0] = line_i[0]
        S[i - 1, 1] = line_i[1]
        S[i - 1, 2] = line_i[2]
        sNames.append(line_i[3])

    return nShot, S, sNames

def read_source_new_utm(source_file):
    #read source file in lat/lon/depth (km)
    with open(source_file, 'r') as f:
        lines = f.readlines()
    line0 = lines[0].split()
    nShot = int(line0[0])
    S = np.zeros((nShot, 3))
    sNames = []
    for i in range(1, nShot + 1):
        line_i = lines[i].split()
        S[i - 1, 0] = float(line_i[4])
        S[i - 1, 1] = float(line_i[5])
        S[i - 1, 2] = float(line_i[6])
        sNames.append(line_i[3])

    return nShot, S, sNames

def read_stat_name(station_file):
    #read station file in cartesian coordinates
    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nRec=int(line0[0])
    R=np.zeros((nRec,3))
    statnames = []
    for i in range(1,nRec+1):
        line_i=lines[i].split()
        R[i-1,0]=int(line_i[0])
        R[i-1,1]=int(line_i[1])
        R[i-1,2]=int(line_i[2])
        statnames.append(line_i[3])
    return nRec, R, statnames

def read_stat_name_utm(station_file):
    # read station file in lat/lon/depth (km)
    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0 = lines[0].split()
    nRec = int(line0[0])
    R = np.zeros((nRec, 2))
    statnames = []
    for i in range(1, nRec + 1):
        line_i = lines[i].split()
        R[i - 1, 0] = float(line_i[3])
        R[i - 1, 1] = float(line_i[4])
        statnames.append(line_i[2])
    return nRec, R, statnames

def Vsp_read(nx, ny, nz, snap_file):
    #Read binary files for model: Vs,Vp,rho etc.
    fid = open(snap_file, 'r')
    sek = np.fromfile(fid, dtype='<f4')
    Vp = np.reshape(sek, [ny, nz, nx])
    return Vp


def write_model_Vsp(rho, Vs, Vp):
    #write updated Vs,Vp,rho
    [ny, nz, nx] = Vs.shape

    model_file0 = 'rho3dfile1.d'
    model_file1 = 'vs3dfile1.s'
    model_file2 = 'vp3dfile1.p'

    fid = open(model_file0, 'wb')
    sek1 = np.array(rho, dtype=np.float32)
    sek1.astype('float32').tofile(fid)

    fid00 = open(model_file1, 'wb')
    sek2 = np.array(Vs, dtype=np.float32)
    sek2.astype('float32').tofile(fid00)

    fid01 = open(model_file2, 'wb')
    sek3 = np.array(Vp, dtype=np.float32)
    sek3.astype('float32').tofile(fid01)


def write_par_opt_source_i(Si):
    #write e3d.par file from template for Emod3d: V3.0.7-dc_stf_file (forward simulation with waveform output only)
    #only source location in cartesian coordinated changed, other parameter kept the same.
    #e3d.par for each source was prepared in davanced and store at FWD/All_srf/
    fname = 'FWD/e3d_mysource_i_opt.par'
    # file_default='e3d_mysource.par'
    os.system('cp FWD/e3d_opt_default.par FWD/e3d_mysource_i_opt.par')
    fid = open(fname, 'a')
    fid.write("%s\n" % ('xsrc=' + str(Si[0])))
    fid.write("%s\n" % ('ysrc=' + str(Si[1])))
    fid.write("%s\n" % ('zsrc=' + str(Si[2])))

def write_par_source_i(Si):
    #write e3d.par file from template for Emod3d: V3.0.7-xyz (forward simulation with waveform and forward strain tensor output)
    sname = 'FWD/e3d_mysource_i.par'

    os.system('cp FWD/e3d_mysource_xyz_default.par FWD/e3d_mysource_i.par')
    fid = open(sname, 'a')
    fid.write("%s\n" % ('xsrc=' + str(int(Si[0]))))
    fid.write("%s\n" % ('ysrc=' + str(int(Si[1]))))
    fid.write("%s\n" % ('zsrc=' + str(int(Si[2]))))