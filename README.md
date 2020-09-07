==========================================================================
                        multiecho_B0_recon
==========================================================================

prepare_data.m
==============
The Matlab script 'prepare_data(id)' prepares data for multi-echo 
reconstruction including B0 estimation based on
- raw data for reconstruction stored in 'data/raw/[id].mat' and
- k-space data stored in '../13C_MRI_simulation/results/[id].mat',
and stores the processed data in 'data/proc/[id].mat'.
    Example: prepare_data('invivo1')

Raw data 'data/raw/[id].mat':
       - coil_map:    sensitivity map      4-D array: [Nx,Ny,Nz,Nc]
       - alpha_map:   flip angle map [rad] 7-D array: [Nx,Ny,Nz,1,1,1,Nm]
       - T2s_map:     T2* map [s]          7-D array: [Nx,Ny,Nz,1,1,1,Nm]
       - B0_map:      initial estimate [Hz]3-D array: [Nx,Ny,Nz]
                      here: B0_map=0
       - k:           sample points [1/m]  2-D array: [Ns,dim]
       - ts:          sampling time [s]    1-D array: [Ns,1]
       - fm:          chemical shifts [Hz] 1-D array: [1,Nm]
       - TE:          echo times [s]       1-D array: [1,Ne]
       - FOV:         field-of-view [m]    1-D array: [1,dim]

K-space data '../13C_MRI_simulation/results/[id].mat':
       - s:           k-space              7-D array: [Ns,1,1,Nc,Nd,1,Ne]
	
Image-domain data is given as 7-D array [x,y,z,c,d,hp,m] of 
size [Nx,Ny,Nz,Nc,Nd,Np,Nm] with attributes:
       - x,y,z        voxel
       - c            coils
       - d            dynamics
       - hp           heart phases
       - m            metabolites

K-space data is given as 7-D array [ks,1,1,c,d,hp,e] of 
size [Ns,1,1,Nc,Nd,Np,Nm] with attributes:
       - ks           samples
       - c            coils
       - d            dynamics
       - hp           heart phases
       - e            echoes

Processed data 'data/proc/[id].mat':
       - E:           encoding operator    2-D array: [Ns*Ne,Nv*Nm]
       - E_nrm:       normalization factor of E
       - s:           k-space              3-D array: [Ns*Ne,Nc,Nd]
       - coil_map:    coil sensitivity map 2-D array: [Nv,Nc]
       - t:           time vector [s]      1-D array: [Ns*Ne,1]
       - B0_map:      initial B0 [Hz]      2-D array: [Nx,Ny]
                      here: B0_map=0

Dimensions:
       - Nx,Ny,Nz     number of voxels per dimension
       - Nv           total number of voxels Nv=Nx*Ny*Nz
       - Nc           number of coils
       - Nd           number of dynamics
       - Nm           number of metabolites
       - Ns           number of k-space samples
       - Ne           number of echoes
       - dim          spatial encoding dimensionality (2D or 3D)

Exemplary raw data for reconstruction based on in vivo measurements using 
hyperpolarized [1-13C]pyruvate in the pig heart is available in 
'data/raw/invivo1-4.mat'.
Corresponding synthetic k-space data can be found in the repository 
13C_MRI_simulation (https://github.com/jtraechtler/13C_MRI_simulation.git)
in '13C_MRI_simulation/results/invivo1-4.mat'.


run_recon.sh
============
The shell script 'run_recon id' specifies reconstruction parameters 
(e.g. regularization parameters, dynamic time frame XX) and runs the 
reconstruction with B0 estimation for the multi-echo data stored 
in 'data/proc/[id].mat'.
Reconstructed data is stored in 'recon/[id]_dyn_XX.mat'.
    Example: ./run_recon.sh "data/proc/invivo1.mat"

Reconstruced data 'recon/[id]_dyn_XX.mat':
       - img:         reconstruced image   7-D array: [Nx,Ny,1,1,1,1,Nm]
       - B0_map:      estimated B0 map [Hz]2-D array: [Nx,Ny]
       - f_val:       error
       - lambda_rho:  image regularization parameter
       - lambda_b0:   B0 regularization parameter
