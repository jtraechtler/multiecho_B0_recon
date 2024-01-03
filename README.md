# multiecho_B0_recon

## UPDATE

The repository was moved to https://gitlab.ethz.ch/ibt-cmr-public/multiecho_b0_recon.git.

## Description

Python code to reconstruct hyperpolarized 13C metabolic MRI data with B0 correction from echo shift-encoded (multi-echo) k-space data without initial B0 estimate. Metabolite images and B0 map are jointly estimated.\
Reference: https://doi.org/10.1002/mrm.28710

## Requirements

Python\
MATLAB\
Repository: https://gitlab.ethz.ch/ibt-cmr-public/13c_mri_simulation.git

## Usage

1. Run Matlab script ```prepare_data(id)``` to prepare data for multi-echo reconstruction with B0 estimation based on raw data for reconstruction stored in *'data/raw/[id].mat'* and based on k-space data stored in *'../13C_MRI_simulation/results/[id].mat'*. The processed data is stored in *'data/proc/[id].mat'*.
2. Run ```run_recon id``` to reconstruct data using multi-echo reconstruction with B0 estimation. The shell script ```run_recon.sh``` specifies reconstruction parameters (e.g. regularization parameters, dynamic time frame *XX*) and calls the python script ```recon_B0.py``` to run the reconstruction with B0 estimation for the multi-echo data stored in *'data/proc/[id].mat'*. Reconstructed data is stored in *'recon/[id]\_dyn\_[XX].mat'*.

Example: 
```
prepare_data('invivo1')
```
```
./run_recon.sh "data/proc/invivo1.mat"
```

## Data

### Raw data *'data/raw/[id].mat'*:

| parameter   | description | dimension     |
| :---        |    :---     |   :---        |
| *coil_map* |    sensitivity map |      4-D array: [*Nx*,*Ny*,*Nz*,*Nc*] |
| *alpha_map* |   flip angle map [rad] | 7-D array: [*Nx*,*Ny*,*Nz*,1,1,1,*Nm*] |
| *T2s_map* |     T2* map [s] |          7-D array: [*Nx*,*Ny*,*Nz*,1,1,1,*Nm*] |
| *B0_map* |      initial B0 estimate [Hz], here: *B0_map*=0 |       3-D array: [*Nx*,*Ny*,*Nz*] |
| *k* |           sample points [1/m] |  2-D array: [*Ns*,*dim*] |
| *ts* |          sampling time [s] |    1-D array: [*Ns*,1] |
| *fm* |          chemical shifts [Hz] | 1-D array: [1,*Nm*] |
| *TE* |          echo times [s] |       1-D array: [1,*Ne*] |
| *FOV* |         field-of-view [m] |    1-D array: [1,*dim*] |
                     
### K-space data *'../13C_MRI_simulation/results/[id].mat'*:

| parameter   | description | dimension     |
| :---        |    :---     |   :---        |
| *s* |           k-space |              7-D array: [*Ns*,1,1,*Nc*,*Nd*,1,*Ne*] |

Exemplary raw data for reconstruction based on in vivo measurements using hyperpolarized [1-13C]pyruvate in the pig heart is available in *'data/raw/invivo1-4.mat'*.
Corresponding synthetic k-space data can be found in the repository 13C_MRI_simulation (https://gitlab.ethz.ch/ibt-cmr-public/13c_mri_simulation.git) in *'13C_MRI_simulation/results/invivo1-4.mat'.*

### Processed data *'data/proc/[id].mat'*:
| parameter   | description | dimension     |
| :---        |    :---     |   :---        |
| *E* |         encoding operator |    2-D array: [*NsNe*,*NvNm*] |
| *E_nrm* |     normalization factor of *E* | scalar |
| *s* |         k-space |              3-D array: [*NsNe*,*Nc*,*Nd*] |
| *coil_map* |  coil sensitivity map | 2-D array: [*Nv*,*Nc*] |
| *t* |         time vector [s] |      1-D array: [*NsNe*,1] |
| *B0_map* |    initial B0 estimate [Hz], here: *B0_map*=0 |      2-D array: [*Nx*,*Ny*] |
                      
### Reconstruced data *'recon/[id]\_dyn\_[XX].mat'*:
| parameter   | description | dimension     |
| :---        |    :---     |   :---        |
| *img* |       reconstructed image |   7-D array: [*Nx*,*Ny*,1,1,1,1,*Nm*] |
| *B0_map* |    estimated B0 map [Hz] | 2-D array: [*Nx*,*Ny*] |
| *f_val* |     error | scalar |
| *lambda_rho* |image regularization parameter | scalar |
| *lambda_b0* | B0 regularization parameter | scalar |
		      
### Dimensions:

| parameter   | description |
| :---        |    :---     |
| *Nx*,*Ny*,*Nz* |     number of voxels per dimension |
| *Nv* |           total number of voxels *Nv*=*NxNyNz* |
| *Nc* |           number of coils |
| *Nd* |           number of dynamics |
| *Nm* |           number of metabolites |
| *Ns* |           number of k-space samples |
| *Ne* |           number of echoes |
| *dim* |          spatial encoding dimensionality (2D or 3D) |

**Image-domain** data is given as 7-D array [*x*,*y*,*z*,*c*,*d*,*hp*,*m*] of size [*Nx*,*Ny*,*Nz*,*Nc*,*Nd*,*Np*,*Nm*] with attributes:
| parameter   | description |
| :---        |    :---     |
| *x,y,z* |        voxels |
| *c* |            coils |
| *d* |            dynamics |
| *hp* |           heart phases |
| *m* |            metabolites |

**K-space** data is given as 7-D array [*ks*,1,1,*c*,*d*,*hp*,*e*] of size [*Ns*,1,1,*Nc*,*Nd*,*Np*,*Nm*] with attributes:
| parameter   | description |
| :---        |    :---     |
| *ks* |           samples |
| *c* |            coils |
| *d* |            dynamics |
| *hp* |           heart phases |
| *e* |            echoes |
