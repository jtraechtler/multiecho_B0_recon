function prepare_data(id)
% function prepare_data(id)
%=========================================================================
%
%	TITLE:
%       prepare_data.m
%
%	DESCRIPTION:
%       Example script to prepare data for reconstruction with joint
%       estimation of metabolite images and field map based on
%       pre-processed raw data stored in 'data/raw/[id].mat' and  
%       k-space data stored in '../13C_MRI_simulation/recon/[id].mat'.
%        - coil_map:    sensitivity map      4-D array: [Nx,Ny,Nz,Nc]
%        - alpha_map:   flip angle map [rad] 7-D array: [Nx,Ny,Nz,1,1,1,Nm]
%        - T2s_map:     T2* map [s]          7-D array: [Nx,Ny,Nz,1,1,1,Nm]
%        - B0_map:      initial estimate [Hz]3-D array: [Nx,Ny,Nz]
%                       here: B0_map=0
%        - k:           sample points [1/m]  2-D array: [Ns,dim]
%        - ts:          sampling time [s]    1-D array: [Ns,1]
%        - fm:          chemical shifts [Hz] 1-D array: [1,Nm]
%        - TE:          echo times [s]       1-D array: [1,Ne]
%        - FOV:         field-of-view [m]    1-D array: [1,dim]
%        - s:           k-space              7-D array: [Ns,1,1,Nc,Nd,1,Ne]
%
%	INPUT:
%       id:             data set ID
%
%	SAVED FILES:							
%       Simulation results are stored in 'results/[id].mat'
%        - E:           encoding operator    2-D array: [Ns*Ne,Nv*Nm]
%        - E_nrm:       normalization factor of E
%        - s:           k-space              3-D array: [Ns*Ne,Nc,Nd]
%        - coil_map:    coil sensitivity map 2-D array: [Nv,Nc]
%        - t:           time vector [s]      1-D array: [Ns*Ne,1]
%        - B0_map:      initial B0 [Hz]      2-D array: [Nx,Ny]
%                       here: B0_map=0
%
%	VERSION HISTORY:
%       200821JT Initial version for release
%
%	    JULIA TRAECHTLER (TRAECHTLER@BIOMED.EE.ETHZ.CH)
%
%=========================================================================

%% input
if nargin < 1
    id = 'invivo1';
end

%% add path
addpath('../13C_MRI_simulation/code/')

%% get_data
load(['data/raw/',id,'.mat'],'coil_map','alpha_map','T2s_map','B0_map',...
    'k','ts','fm','TE','FOV');
load(['../13C_MRI_simulation/results/',id,'.mat'],'s');

%% dimensions
dim = size(k,2);        % spatial encoding dimensionality (2D)
N = size(B0_map);       % number of pixels per dimension [Nx,Ny]
Nv = prod(N);           % total number of pixels Nx*Ny
Ns = length(ts);        % number of k-space samples
Ne = length(TE);        % number of echoes
Nc = size(coil_map,4);  % number of coils
Nd = size(s,5);         % number of dynamics

%% t: [Ns,Ne] time vector including sampling time & shifted echo times [s]
t = ts+TE;

%% r: [Nv,dim] spatial coordinates [m]
[r] = generate_spatial_grid(N,FOV,dim);

%% F: [Ns,1,1,Nv,1] spatial encoding (Fourier transform)
[F] = generate_F_operator(k,r);

%% M: [Ns,Ne,1,1,Nm] chemical shift encoding
[M] = generate_M_operator(fm,t);

%% A: [1,Ne,1,Nv*Nz,Nm] depolarization
[A] = generate_A_operator(alpha_map,Ne);

%% T2s: [Ns,Ne,1,Nv*Nz,Nm] transverse relaxation
[T2s] = generate_T2s_operator(T2s_map,t);

%% E: [Ns*Ne*Nc,Nv*Nm] extended encoding operator without B0 and without W
[E]   = generate_E_operator(F,M,1,A,T2s,1);
E_nrm = svds(double(E),1);       % normalization factor
E     = E/E_nrm;                 % normalize E

%% coil_map: [Nv,Nc] coil sensitiviy weighting
coil_map = reshape(coil_map,[Nv,Nc]);

%% reshape s: [Ns*Ne,Nc,Nd] multi-echo, multi-coil k-space
s = reshape(permute(s,[1,7,4,5,2,3,6]),[Ns*Ne,Nc,Nd]);

%% reshape t: [Ns*Ne,1]
t = reshape(t,[Ns*Ne,1]);

%% save data
save(['data/proc/',id,'.mat'],'E','E_nrm','s','coil_map','t','B0_map')

end