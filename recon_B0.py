import random
random.seed(5)
import tensorflow as tf
import scipy.io as sio
import numpy as np
from scipy.optimize import minimize
import argparse

#--------------------------------------------------------------------------------------
#--------------------------------- define signal model --------------------------------
#--------------------------------------------------------------------------------------
# signal model for the multi-echo, multi-coil k-space signal
def signal_model(rho, B0_map, t, E, Nv, Nm, Nc, Nse):
    # reshape input
    rho = tf.reshape(rho, [Nc, -1])
    t 	= tf.reshape(tf.constant(t, dtype = tf.float32), [Nse, 1])

    # B0-induced phase-offset
    exp = 2 * np.pi * t * tf.reshape(B0_map, [1, Nv]) 
    B0 	= tf.reshape(tf.complex(tf.cos(exp), tf.sin(exp)), [Nse, 1, Nv])

    # extended encoding matrix
    E = tf.reshape(tf.constant(E, dtype = tf.complex64), [Nse, Nm, Nv])

    # extended encoding matrix including B0-induced phase-offset
    EB0 = tf.reshape(E * B0, [Nse, Nv*Nm])

    # k-space signal
    s_model = tf.matmul(EB0, rho, transpose_b = True) # Nse - Nc
    return s_model

#--------------------------------------------------------------------------------------
#-------------------------------- get and prepare data --------------------------------
#--------------------------------------------------------------------------------------
# parser
parser = argparse.ArgumentParser('''
Reconstruct multi-echo k-space data jointly estimating for metabolite images and field map. Solves the following optimization problem:\n
min_{rho,b0}  || E(b0) * rho - s||_2^2 + lambda_rho * TV(rho) + lambda_b0 * TV(b0)

Format of input mat file:
E        -- Ns*Ne x Nv*Nm complex encoding matrix that encapsulates Fourier encoding and weighting, 
            expected to have unit norm >> E = E/E_nrm
E_nrm    -- normalization factor of encoding operator: E_nrm = svds(double(E),1)
s        -- Ns*Ne x Nc x Nd complex k-space measurement vector
coil_map -- Nv x Nc complex coil sensitivity map
t        -- Ns*Ne time vector including sample times for all echoes [s]
B0_map 	 -- Nx x Ny initial B0 map [Hz]

''')
parser.add_argument('--input', type=str, help = 'path to input mat file', required = True)
parser.add_argument('--output', type=str, help = 'path to output mat file', required = True)
parser.add_argument('--lambda_rho', type=float, help = 'image regularization weight', required = True)
parser.add_argument('--lambda_b0', type=float, help = 'b0 regularization weight', required = True)
parser.add_argument('--max_iter', type=int, default = 3000, help = 'max number of optimizer iterations')
parser.add_argument('--tune_b0', type=int, default = 1, help = '1 - to optimize wrt b0, 0 - just use initialization for b0')
parser.add_argument('--dyn', type=int, default = 1, help = 'k-space dynamic')
args = parser.parse_args()

# arguments
lambda_rho 	= args.lambda_rho
lambda_b0  	= args.lambda_b0
tune_b0 	= bool(args.tune_b0)
dyn 		= args.dyn

# avoid division by 0, approximates norm |x| by sqrt(x^2+epsilon)
epsilon = 1e-5

# read data
mat 	 = sio.loadmat(args.input)
E 		 = mat['E'].astype(np.complex64)
s 		 = mat['s'].astype(np.complex64)
t 		 = mat['t'].astype(np.float32).ravel()
B0_init  = mat['B0_map'].astype(np.float32).transpose()
coil_map = mat['coil_map'].astype(np.complex64)
E_nrm 	 = mat['E_nrm'].astype(np.float32)

# dimensions
N = B0_init.shape		# number of pixels per dimension
Nv = N[0]*N[1] 			# total number of pixels
Nse = E.shape[0] 		# number of k-space samples * number of echoes
Nm = int(E.shape[1]/Nv)	# number of metabolites
Nc = s.shape[1]			# number of coils

# measured k-space signal for current dynamic
s = s[:,:,dyn-1]

# normalize k-space
s_nrm = np.sqrt(s.size) / np.linalg.norm(s)
s = s * s_nrm

# configure tensorflow
ftype = tf.float32
ntype = np.float32
config = tf.compat.v1.ConfigProto()
config.intra_op_parallelism_threads = 12
config.inter_op_parallelism_threads = 6
sess = tf.compat.v1.InteractiveSession(config=config)

#--------------------------------------------------------------------------------------
#------------------------- components for minimization problem ------------------------
#--------------------------------------------------------------------------------------
# regularization parameters
ph_lambda_rho = tf.compat.v1.placeholder(tf.float32, shape = [])
ph_lambda_b0 = tf.compat.v1.placeholder(tf.float32, shape = [])

# image estimate
rho_re = tf.Variable(np.random.rand(1, Nm, N[0], N[1])*0.00001, dtype = tf.float32, trainable = True)
rho_im = tf.Variable(np.random.rand(1, Nm, N[0], N[1])*0.00001, dtype = tf.float32, trainable = True)
rho = tf.complex(rho_re, rho_im)

# field map estimate
B0_map = tf.Variable(B0_init, dtype = tf.float32, trainable = tune_b0)

# coil sensitivity weighting
coil_map = tf.reshape(tf.constant( np.transpose(coil_map), dtype = tf.complex64), [Nc, 1, N[0], N[1]])
coil_nrm = tf.sqrt(tf.reduce_sum(tf.abs(coil_map)**2, axis = 0, keepdims = True )) + 1e-6
coil_map = coil_map / tf.cast(coil_nrm, tf.complex64 )

# coil sensitivity weighted image
rhoW = coil_map * rho

# k-space signal model
s_model = signal_model(rhoW, B0_map, t, E, Nv, Nm, Nc, Nse)

#--------------------------------------------------------------------------------------
#-------------------------------- minimization problem --------------------------------
#--------------------------------------------------------------------------------------
# residual
residual = (s_model) - tf.constant(s, dtype = tf.complex64)

# data error
data_error = tf.reduce_sum(tf.sqrt( tf.math.real(residual)**2 + tf.math.imag(residual)**2 + epsilon )**2)

# discrete spatial derivatives of image
drho_1 = rho - tf.roll(rho, 1, 2)
drho_2 = rho - tf.roll(rho, 1, 3)

# total variation of image
TV_rho = tf.reduce_sum( tf.sqrt( tf.math.real(drho_1)**2 + tf.math.real(drho_2)**2 + tf.math.imag(drho_1)**2 + tf.math.imag(drho_2)**2 + epsilon) )

# discrete spatial derivatives of field map
dB0_1 = B0_map - tf.roll(B0_map, 1, 0)
dB0_2 = B0_map - tf.roll(B0_map, 1, 1)

# total variation of field map
TV_b0 = tf.reduce_sum( tf.sqrt( tf.math.real(dB0_1)**2 + tf.math.real(dB0_2)**2 + tf.math.imag(dB0_1)**2 + tf.math.imag(dB0_2)**2 + epsilon) )

# objective function
f = data_error + ph_lambda_rho * TV_rho + ph_lambda_b0 * TV_b0

#--------------------------------------------------------------------------------------
#----------------------------- solve minimization problem -----------------------------
#--------------------------------------------------------------------------------------
learning_rate = tf.compat.v1.placeholder(tf.float32, shape=[])
train_op = tf.contrib.opt.ScipyOptimizerInterface( \
				f, \
				method='L-BFGS-B', \
				options={'maxiter': args.max_iter, 'maxcor' : 30, 'ftol':1e-10, 'eps':1e-7, 'iprint':1, 'disp':False})

# build graph
sess.run(tf.global_variables_initializer())

# optimize
ferrs = []
feed = {ph_lambda_b0 : lambda_b0, ph_lambda_rho : lambda_rho}
FVAL_TRACE = []
def print_loss(loss):
	global FVAL_TRACE
	FVAL_TRACE.append(loss)
train_op.minimize(sess,
				 feed_dict = feed,
				 loss_callback=print_loss,
				 fetches=[f])
ferrs = FVAL_TRACE

#--------------------------------------------------------------------------------------
#-------------------------- solution of minimization problem --------------------------
#--------------------------------------------------------------------------------------
rho = rho.eval().transpose() / s_nrm / E_nrm # unnormalize
rho = np.reshape(rho,[N[1], N[0], 1, 1, 1, 1, Nm])
B0_map = B0_map.eval().transpose()

# save results
sio.savemat(args.output, { 'img' : rho, \
                            'B0_map' : B0_map, \
                            'f_val' : ferrs, \
                            'lambda_rho' : lambda_rho, \
                            'lambda_b0' : lambda_b0 } )