function [xc, nsol, max_iter] = ELLIPTICFixedPoint1d(tau, eps, a_, nxdx, bcond, ifexact, ifplot, ifdemo)

    if nargin < 5, ifdemo = 0; end
    if nargin < 4, ifplot = 1; end

    % Setting up grid and initial values
    if length(nxdx) > 1
        x = nxdx';
        nx = length(x) - 1;
        dx = diff(x);
    else
        nx = nxdx;
        dx = 1 / nx * ones(nx, 1);
        a = 0; b = 1;
        x = 0 * dx;
        x(1) = a;
        for j = 2:nx + 1
            x(j) = x(j - 1) + dx(j - 1);
        end
    end

    xc = x(1:nx) + dx / 2;
    h = max(dx);

    % Set up boundary conditions
    bflag_left = bcond(1);
    bval_left = bcond(2);
    bflag_right = bcond(3);
    bval_right = bcond(4);

    % Transmissibilities for finite volume
    kcof = ones(nx, 1);
    tx = zeros(nx + 1, 1);
    for j = 2:nx
        tx(j) = 2 / (dx(j - 1) / kcof(j - 1) + dx(j) / kcof(j));
    end

    if bflag_left == 0
        j = 1;
        tx(j) = 2 / (dx(j) / kcof(j));
    end
    if bflag_right == 0
        j = nx + 1;
        tx(j) = 2 / (dx(j - 1) / kcof(j - 1));
    end

    % Diffusion matrix
    diffmat = sparse(nx, nx);
    for j = 2:nx
        gl = j - 1; gr = j;
        diffmat(gl, gl) = diffmat(gl, gl) + tx(j);
        diffmat(gl, gr) = diffmat(gl, gr) - tx(j);
        diffmat(gr, gl) = diffmat(gr, gl) - tx(j);
        diffmat(gr, gr) = diffmat(gr, gr) + tx(j);
    end

    if bflag_left == 0
        j = 1;
        gr = 1;
        diffmat(gr, gr) = diffmat(gr, gr) + tx(j);
    end
    if bflag_right == 0
        j = nx + 1;
        gl = nx;
        diffmat(gl, gl) = diffmat(gl, gl) + tx(j);
    end

    tmat = eye(nx, nx);
    allmat = (dx .* tmat) + (tau * eps * diffmat);

    % Initial conditions
    nsol = sin(4 * pi * (xc + sin(xc))) + pi * xc .^ 4;
    tol = 1e-12;     
    t = 0; T = 0.1; 
    max_iter = 0; % Variable to store maximum iterations used

    while t < T - tau
        t = t + tau;
        oldnsol = nsol;
        
        % Fixed-point iteration using specified while loop structure
        iter = 0;
        while 1
            rhs = allmat \ (dx .* oldnsol) + dx .* tau * a_ .* (nsol - nsol.^3);

            % Apply boundary conditions
            if bflag_left == 0
                if ifexact
                    [dval_left, ~] = exfun(x(1));
                else
                    dval_left = bval_left;
                end
                qleft = tx(1) * dval_left;
            else
                if ifexact
                    [~, nval_left] = exfun(x(1));
                    qleft = -kcof(1) * nval_left;
                else
                    qleft = bval_left;
                end
            end

            if bflag_right == 0
                if ifexact
                    [dval_right, ~] = exfun(x(nx + 1));
                else
                    dval_right = bval_right;
                end
                qright = -tx(nx + 1) * dval_right;
            else
                if ifexact
                    [~, nval_right] = exfun(x(nx + 1));
                    qright = -kcof(nx) * nval_right;
                else
                    qright = bval_right;
                end
            end

            rhs(1) = rhs(1) + qleft;
            rhs(nx) = rhs(nx) - qright;

            % Calculate new guess
            nsol_new = rhs;
            resnorm = norm(nsol_new - nsol);
            
            % Check for convergence
            if resnorm < tol
                nsol = nsol_new;
                break;
            else
                iter = iter + 1;
                nsol = nsol_new;
            end
        end
        
        % Update max_iter_used if the current iteration count is greater
        max_iter = max(max_iter, iter);
    end

    % Plotting the final solution
    if ifplot
        plot(xc, nsol, 'r-', 'LineWidth', 1.5);
        title(sprintf('Solution at a = %g, \\epsilon = %g, t = %g, dt = %g, h = %g', a_, eps, T, tau, h));
        axis([0 1 0 2]);
        xlabel('x');
        ylabel('u(x,t)');
        grid on;
    end

    % Calculate error if exact solution is provided
    if ifexact
        exsol = exfun(xc);
        if ifplot
            plot(xc, exsol, 'k', xc, nsol, 'r*-');
            legend('exact', 'numerical', 'location', 'best');
            title(sprintf('Solution '));
        end
        l2err = sqrt(sum((nsol - exsol) .^ 2 .* dx));
        fprintf('Error is err_l2=%g\n', l2err);
    end

    % Print maximum number of iterations used
    fprintf('Maximum number of iterations used: %d\n', max_iter);
end

% function [v, dv] = exfun(x)
%     v = (exp(10 .* x .* x)) / (exp(10));
%     dv = (20 .* x .* exp(10 .* x .* x)) / (exp(10));
% end
% 
% function v = rhsfun(x)
%     v = x - x .^ 3;
% end