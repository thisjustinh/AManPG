import numpy as np
import time
from numpy import linalg as LA
from numpy.core.numeric import full

def spca_amanpg(b, mu, lamb, n, x0, y0, f_palm, 
                gamma=0.5, type=0, maxiter=1e4, tol=1e-5, 
                verbose=False):
    
    start = time.perf_counter()
    
    m, d = np.shape(b)
    
    # anonymous function for sum of matrix mu times colsums of x
    h = lambda x: np.sum(mu.T * np.sum(np.abs(x), axis=0, keepdims=True))

    if d < m * 2:
        b = b.T @ b
        type = 1

    u, s, v = LA.svd(b, full_matrices=False)
    
    ly = 2 * s[0] ** 2 + 2 * lamb if not type else 2 * s[0] + 2 * lamb
    
    ### initialization ###
    x, y = x0, y0
    total_linesearch = 0
    linesearch_flag = 1
    min_step = 0
    linesearch_flag_y = 1
    t = 1 / ly
    tau = 100 / d

    if not type:
        ay = b.T @ (b @ y)
        ax = b.T @ (b @ x)
    else:
        ay = b @ y
        ax = b @ x

    ### Main Loop ###
    if not lamb == np.inf:
        fx = -2 * np.sum(x * ay)
        fy = np.sum(y * ay) + lamb * LA.norm(y, 'fro') ** 2 + h(y)
        f_rgd = [fx + fy]

        for i in range(1, int(maxiter)):
            if verbose:
                loop_start = time.perf_counter()
                print("=========================")
                print("On iteration", i)
            
            ### update y ###
            t = t * 1.01 if not linesearch_flag_y else max(1 / ly, t / 1.01)
            linesearch_flag_y = 0

            y_t = prox_l1(y - t * 2 * (ay - ax + lamb * y), mu * t, n)
            ayt = b.T @ (b @ y_t) if not type else b @ y_t
            f_ytrial = -2 * np.sum(x * ayt) + np.sum(y_t * ayt) + lamb * LA.norm(y_t, 'fro') ** 2 + h(y_t)
            normpg = LA.norm(y_t - y, 'fro') ** 2 / t ** 2

            # adjust step size to be appropriate and recalculate values
            while (f_ytrial > f_rgd[i - 1] - 1e-3 * t * normpg):
                t *= gamma
                if t < 1e-5 / d:
                    break

                y_t = prox_l1(y - t * 2 * (ay - ax + lamb * y), mu * t, n)
                ayt = b.T @ (b @ y_t) if not type else b @ y_t
                f_ytrial = -2 * np.sum(x * ayt) + np.sum(y_t * ayt) + lamb * LA.norm(y_t, 'fro') ** 2 + h(y_t)
                linesearch_flag_y = 1

            y, ay = y_t, ayt  # assign updated values from loop

            ### update x ###
            if not linesearch_flag:
                tau *= 1.1
            if min_step == 1:
                tau = 1 / d
            linesearch_flag = 0
            min_step = 0
            gx = -2 * ay
            xgx = gx.T @ x
            rgx = gx - 0.5 * x @ (xgx + xgx.T)  # Projected gradient
            
            tx = x - tau * rgx
            sigma, u = LA.eig(tx.T @ tx)
            j = u @ np.diag(np.sqrt(1 / sigma)) @ u.T
            x_trial = tx @ j
            f_xtrial = -2 * np.sum(x_trial * ay)
            fxval = -2 * np.sum(x * ay)
            normpg = LA.norm(rgx, 'fro') ** 2

            while (f_xtrial > fxval - 1e-3 * tau * normpg):
                tau *= gamma
                if tau < 1e-5 / d:
                    min_step = 1
                    break

                tx = x - tau * rgx
                sigma, u = LA.eig(tx.T @ tx)
                j = u @ np.diag(np.sqrt(1 / sigma)) @ u.T
                x_trial = tx @ j
                total_linesearch += 1
                linesearch_flag = 1
                f_xtrial = -2 * np.sum(x_trial * ay)

            x, fx = x_trial, f_xtrial  # assign updated values from loop
            ax = b.T @ (b @ x) if not type else b @ x
            fy = np.sum(y * ay) + lamb * LA.norm(y, 'fro') ** 2 + h(y)
            f_rgd.append(fx + fy)

            if verbose:
                print("fx: ", fx)
                print("fy: ", fy)
                print("Finished with value", f_rgd[i], "and difference", abs(f_rgd[i]-f_rgd[i - 1]))
                print("Done in time", time.perf_counter() - loop_start)
            
            if abs(f_rgd[i] - f_rgd[i - 1]) < tol and f_rgd[i] < f_palm or abs(f_rgd[i] - f_rgd[i - 1]) < 1e-12:
                if verbose:
                    print("Final difference of", abs(f_rgd[i] - f_rgd[i - 1]))
                break
    else:  # lamb is inf
        fx = -2 * np.sum(x * ay)
        fy = LA.norm(y, 'fro') ** 2 + h(y)
        f_rgd = [fx + fy]

        for i in range(1, int(maxiter)):
            if verbose:
                loop_start = time.perf_counter()
                print("=========================")
                print("On iteration", i)

            if not linesearch_flag:
                tau *= 1.1
            if min_step == 1:
                tau = 1 / d

            min_step = 0
            linesearch_flag = 0

            ### update y ###
            t = 1/2
            y = prox_l1(y - 2 * t * (-ax + y), mu * t, n)
            ay = b.T @ (b @ y) if not type else b @ y

            ### update x ###
            gx = -2 * ay
            xgx = gx.T @ x
            rgx = gx - x @ xgx  # Canonical Riemannian gradient
            
            tx = x - tau * rgx
            u, s, v = LA.svd(tx, full_matrices=False)
            x_trial = u @ v  # note that v is already transposed in numpy
            f_xtrial = -2 * np.sum(x_trial * ay)
            fxval = -2 * np.sum(x * ay)
            normpg = LA.norm(rgx, 'fro') ** 2

            while f_xtrial > fxval - 1e-3 * tau * normpg:
                tau *= gamma
                if tau < 1e-3 / d:
                    min_step = 1
                    break

                tx = x - tau * rgx
                u, s, v = LA.svd(tx, full_matrices=False)
                x_trial = u @ v  # note that v is already transposed in numpy
                f_xtrial = -2 * np.sum(x_trial * ay)
                total_linesearch += 1
                linesearch_flag = 1
                f_xtrial = -2 * np.sum(x_trial * ay)

            x, fx = x_trial, f_xtrial  # assign updated values from loop
            ax = b.T @ (b @ x) if not type else b @ x
            fy = LA.norm(y, 'fro') ** 2 + h(y)
            f_rgd.append(fx + fy)

            if verbose:
                print("fx: ", fx)
                print("fy: ", fy)
                print("Finished with value", f_rgd[i], "and difference", abs(f_rgd[i]-f_rgd[i - 1]))
                print("Done in time", time.perf_counter() - loop_start)
            
            if abs(f_rgd[i] - f_rgd[i - 1]) < tol:
                if verbose:
                    print("Final difference of", abs(f_rgd[i] - f_rgd[i - 1]))
                break
    
    y_norm = np.sqrt(np.sum(y ** 2, axis=0)).reshape(4,1)
    y_norm[y_norm == 0] = 1

    sparsity = np.sum(y == 0) / (d * n)
    y_man = np.divide(y, np.ones((d, 1)) @ y_norm.T)

    return i, f_rgd[i], sparsity, time.perf_counter() - start, x, y_man


def prox_l1(b, lamb, r):
    a = np.abs(b) - lamb.T

    act_set = (a > 0).astype(float) if r < 15 else a > 0
    x_prox = act_set * np.sign(b) * a
    inact_set = a <= 0  # if required
    return x_prox


def normalize(x):
    x -= np.mean(x, axis=1, keepdims=True)  # center

    # normalize rows using l2 norm
    x /= LA.norm(x, axis=1, keepdims=True)
    return x


if __name__ == '__main__':
    maxiter = 1e4
    tol = 1e-5
    n = 4  # columns
    d = 500  # dimensions
    m = 1000  # sample size
    mu = 0.1 * np.ones((n, 1))
    type = 0
    lamb = np.inf
    f_palm = 1e5
    
    test_finite = False

    if not test_finite:
        # testing for lambda = inf
        for i in range(1, 11):
            np.random.seed(i)
            a = np.random.rand(m, d)
            a = normalize(a)
            _, _, v = LA.svd(a)
            x0 = v[:, 0:n]
            # print(np.mean(a, axis=1)[0], LA.norm(a, axis=1)[0])
            iter, f_amanpg, sparsity, timediff, x, y_man = spca_amanpg(a, mu, lamb, n, x0, x0, f_palm, verbose=False)
            print(f"{iter} iterations with final value {f_amanpg}, sparsity {sparsity}, timediff {timediff}.")
    else:
        a = np.loadtxt(open('A.csv', 'rb'), delimiter=",")
        _, _, v = LA.svd(a, full_matrices=True)
        x0 = v.T[:, 0:n]
        iter, f_amanpg, sparsity, timediff, x, y_man = spca_amanpg(a, mu, 1, n, x0, x0, f_palm, verbose=False)
        print(f"{iter} iterations with final value {f_amanpg}, sparsity {sparsity}, timediff {timediff}.")
