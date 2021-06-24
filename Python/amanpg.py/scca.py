import numpy as np
import time
from numpy import linalg as LA

def scca_amanpg(a, b, n, p, q, b1, b2,
                maxiter, tol, inner_tol, 
                sp_type):
    
    start = time.perf_counter()
    
    # initial value assignment
    m = np.shape(a)[0]
    tau1 = b1 * np.sqrt((n + np.log(p)) / m)
    tau2 = b2 * np.sqrt((n + np.log(p)) / m)
    
    
    if sp_type == 'l1':
        h1 = lambda x: tau1 * np.sum(np.abs(x), axis = 0)
        h2 = lambda x: tau2 * np.sum(np.abs(x), axis = 0)
        prox.func = lambda b, lamb, r: prox_l1(b, lamb, r)
    elif sp_type == 'l21':
        h1 = lambda x: tau1 * np.sum(LA.norm(x,2,1))
        h2 = lambda x: tau2 * np.sum(LA.norm(x,2,1))
        prox.func = lambda b, lamb, r: prox_l21(b, lamb, r)
         
    inner_flag1 = 0
    #Dn = sparse(DuplicationM(n))
    pDn = LA.solve(Dn.T @ Dn, Dn.T) 
    
    #center data
    a -= np.mean(a, axis=1, keepdims=True)
    b -= np.mean(b, axis=1, keepdims=True)
    t_min = 1e-4 # minimum stepsize
    atb = A.T @ B / (m-1)
    bta = atb.T
    ata = a.T @ a / (m-1)
    btb = b.T @ b / (m-1)
    
    #set type
    gamma = 1e-4 / (m-1)
    u, s, v = LA.svd(a, full_matrices=False)
    if s[-1] < 1e-4:
        M1, pda = (1 - gamma)*ata + gamma*np.identity(p), 0
    else:
        M1, pda = ata, 1
        
    u, s, v = LA.svd(b, full_matrices=False)
    if s[-1] < 1e-4:
        M2, pdb = (1 - gamma)*btb + gamma*np.identity(q), 0
    else:
        M2, pda = btb, 1 
    
    typea = 1 if m > p/2 else 0
    typeb = 1 if m > q/2 else 0
