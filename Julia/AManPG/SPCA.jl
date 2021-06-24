module SparsePCA

using Base: Bool
using LinearAlgebra
using Printf: @printf

export spca_amanpg, prox_l1

function spca_amanpg(a::Matrix{Float64}, mu::Matrix{Float64}, lambda::Number, n::Integer, 
                     x0::Matrix{Float64}, y0::Matrix{Float64}, f_palm::Number;
                     type::Integer=0, maxiter::Integer=1e4, tol::Number=1e-5, gamma::Number=0.5, verbose::Bool=false)
    start = time()

    b = copy(a)
    m, d = size(b)

    # anonymous function to get l1 norm per column times mu
    h = x -> sum(mu * sum(abs(x), 1))

    if d < m * 2
        b = b' * b
        type = 1
    end

    _, s, _ = svd(b)  # used to get largest singular value
    ly = type == 0 ? 2 * s[1] ^ 2 + 2 * lambda : 2 * s[1] + 2 * lambda

    ### initial point ###
    x, y = x0, y0
    total_linesearch = 1
    linesearch_flag = 1
    min_step = 0
    linesearch_flag_y = 1
    t, tau = 1 / ly, 100 / d
    if type == 0
        ay = b' * (b * y)
        ax = b' * (b * x)
    else
        ay = b * y
        ax = b * x
    end

    ### main loop ###
    if lambda != Inf
        fx = -2 * sum(x .* ay)
        fy = -2 * sum(y .* ay) + lambda * norm(y) ^ 2 + h(y)
        f_rgd = [fx + fy]

        for i = 2:maxiter
            if verbose
                i_start = time()
                println("=========================")
                @printf("On iteration %i\n", i)
            end

            ### update y ###
            t = linesearch_flag_y == 0 ? t * 1.01 : max(1/ly, t/1.01)
            linesearch_flag_y = 0
            y_t = prox_l1(y - t * 2 * (ay - ax + lambda * y), mu * t, n)
            ayt = type == 0 ? b' * (b * y_t) : b * y_t
            f_ytrial = -2 * sum(x .* ayt) + sum(y_t .* ayt) + lambda * norm(y_t) ^ 2 + h(y_t)
            normpg = norm(y_t - y) ^ 2 / t ^ 2

            while f_ytrial > f_rgd[i - 1] - 1e-3 * t * normpg
                t *= gamma
                if t < 1e-5 / d
                    break
                end

                y_t = prox_l1(y - t * 2 * (ay - ax + lambda * y), mu * t, n)
                ayt = type == 0 ? b' * (b * y_t) : b * y_t
                f_ytrial = -2 * sum(x .* ayt) + sum(y_t .* ayt) + lambda * norm(y_t) ^ 2 + h(y_t)
                linesearch_flag_y = 1
            end

            y, ay = y_t, ayt  # assign updated values from loop

            ### update x ###
            if linesearch_flag == 0
                tau *= 1.1
            end
            if min_step == 1
                tau = 1/d
            end

            linesearch_flag = 0
            min_step = 0
            gx = -2 * ay
            xgx = gx' * x
            rgx = gx - 0.5 * x * (xgx + xgx')  # projected gradient
            
            tx = x - tau * rgx
            vals, vecs = eigen(tx' * tx)
            j = vecs * Diagonal(sqrt(1 ./ vals)) * vecs'
            x_trial = tx * j
            fx_trial = -2 * sum(x_trial .* ay)
            fxval = -2 * sum(x .* ay)
            normpg = norm(rgx) ^ 2

            while fx_trial > fxval - 1e-3 * tau * normpg
                tau *= gamma
                if tau < 1e-5 / d
                    min_step = 1
                    break
                end

                tx = x - tau * rgx
                vals, vecs = eigen(tx' * tx)
                j = vecs * Diagonal(sqrt(1 ./ vals)) * vecs'
                x_trial = tx * j
                fx_trial = -2 * sum(x_trial .* ay)
                total_linesearch += 1
                linesearch_flag = 1
            end

            x, fx = x_trial, fx_trial  # assign updated values from loop
            ax = type == 0 ? b' * (b * x) : b * x
            fy = sum(y .* ay) + lambda * norm(y)^2 + h(y)
            push!(f_rgd, fx + fy)

            if verbose
                @printf("fx: %.4f, fy: %.4f\n", fx, fy)
                @printf("Finished with difference %g\n", abs(f_rgd[i] - f_rgd[i - 1]))
                @printf("Done in time %f\n", time() - i_start)
            end

            # if normDsquared < tol^2
            if (abs(f_rgd[i] - f_rgd[i -1]) < tol && f_rgd[i] < f_palm) || abs(f_rgd[i] - f_rgd[i - 1]) < 1e-12
                if verbose
                    printf("Final difference of %g\n", abs(f_rgd[i] - f_rgd[i - 1]))
                end
                break
            end
        end
    else  # elastic net parameter is Inf
        fx = -2 * sum(x .* ay)
        fy = norm(y) ^ 2 + h(y)
        f_rgd = [fx + fy]

        for i = 2:maxiter
            if verbose
                i_start = time()
                println("=========================")
                @printf("On iteration %i\n", i)
            end

            if linesearch_flag == 0
                tau *= 1.1
            end
            if min_step == 1
                tau = 1 / d
            end

            min_step = 0
            linesearch_flag = 0

            ### update y ###
            t = 1/2
            y = prox_l1(y - 2 * t * (-ax + y), mu * t, n)
            ay = type == 0 ? b' * (b * y) : b * y

            ### update x ###
            gx = -2 * ay
            xgx = gx' * x
            rgx = gx - x * xgx  # canonical Riemannian gradient
            
            tx = x - tau * rgx
            u, _, v = svd(tx)
            x_trial = u * v'
            fx_trial = -2 * sum(x_trial .* ay)
            fxval = -2 * sum(x .* ay)
            normpg = norm(rgx)^2

            while fx_trial > fxval - 1e-3 * tau * normpg
                tau *= gamma
                if tau < 1e-3 / d
                    min_step = 1
                    break
                end

                tx = x - tau * rgx
                u, _, v = svd(tx)
                x_trial = u * v'
                fx_trial = -2 * sum(x_trial .* ay)
                total_linesearch += 1
                linesearch_flag = 1
            end

            x, fx = x_trial, fx_trial  # assign updated values from loop
            ax = type == 0 ? b' * (b * x) : b * x
            fy = norm(y)^2 + h(y)
            push!(f_rgd, fx + fy)

            if verbose
                @printf("fx: %.4f, fy: %.4f\n", fx, fy)
                @printf("Finished with difference %g\n", abs(f_rgd[i] - f_rgd[i - 1]))
                @printf("Done in time %f\n", time() - i_start)
            end

            if abs(f_rgd[i] - f_rgd[i - 1]) < tol
                if verbose
                    printf("Final difference of %g\n", abs(f_rgd[i] - f_rgd[i - 1]))
                end
                break
            end
        end
    end

    y_norm = sqrt(sum(y .^ 2, axis=1))
    y_norm[y_norm == 0] = 1
    y_man = y ./ (ones(d, 1) * y_norm)
    sparsity = sum(y == 0) / (d * n)

    return i, f_rgd[i], sparsity, time() - start, x, y_man
end

function prox_l1(b, lambda, r)
    a = abs(b) - lambda
    act_set = r < 15 ? float(a > 0) : a > 0
    return act_set .* sign(a) .* a, act_set
end

end  # end module