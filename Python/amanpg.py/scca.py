{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import time\n",
    "from numpy import linalg as LA\n",
    "\n",
    "def scca_amanpg(a, b, n, p, q, b1, b2,\n",
    "                maxiter, tol, inner_tol, \n",
    "                sp_type):\n",
    "    \n",
    "    start = time.perf_counter()\n",
    "    \n",
    "    # initial value assignment\n",
    "    m = np.shape(a)[0]\n",
    "    tau1 = b1 * np.sqrt((n + np.log(p)) / m)\n",
    "    tau2 = b2 * np.sqrt((n + np.log(p)) / m)\n",
    "    \n",
    "    \n",
    "    if sp_type == 'l1':\n",
    "        h1 = lambda x: tau1 * np.sum(np.abs(x), axis = 0)\n",
    "        h2 = lambda x: tau2 * np.sum(np.abs(x), axis = 0)\n",
    "        prox.func = lambda b, lamb, r: prox_l1(b, lamb, r)\n",
    "    elif sp_type == 'l21':\n",
    "        h1 = lambda x: tau1 * np.sum(LA.norm(x,2,1))\n",
    "        h2 = lambda x: tau2 * np.sum(LA.norm(x,2,1))\n",
    "        prox.func = lambda b, lamb, r: prox_l21(b, lamb, r)\n",
    "         \n",
    "    inner_flag1 = 0\n",
    "    #Dn = sparse(DuplicationM(n))\n",
    "    pDn = LA.solve(Dn.T @ Dn, Dn.T) \n",
    "    \n",
    "    #center data\n",
    "    a -= np.mean(a, axis=1, keepdims=True)\n",
    "    b -= np.mean(b, axis=1, keepdims=True)\n",
    "    t_min = 1e-4 # minimum stepsize\n",
    "    atb = A.T @ B / (m-1)\n",
    "    bta = atb.T\n",
    "    ata = a.T @ a / (m-1)\n",
    "    btb = b.T @ b / (m-1)\n",
    "    \n",
    "    #set type\n",
    "    gamma = 1e-4 / (m-1)\n",
    "    u, s, v = LA.svd(a, full_matrices=False)\n",
    "    if s[-1] < 1e-4:\n",
    "        M1, pda = (1 - gamma)*ata + gamma*np.identity(p), 0\n",
    "    else:\n",
    "        M1, pda = ata, 1\n",
    "        \n",
    "    u, s, v = LA.svd(b, full_matrices=False)\n",
    "    if s[-1] < 1e-4:\n",
    "        M2, pdb = (1 - gamma)*btb + gamma*np.identity(q), 0\n",
    "    else:\n",
    "        M2, pda = btb, 1 \n",
    "    \n",
    "    typea = 1 if m > p/2 else 0\n",
    "    typeb = 1 if m > q/2 else 0"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
