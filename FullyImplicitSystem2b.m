function FullyImplicitSystem2b(alpha, c, tau)
    if nargin < 3
        alpha = 1;
        c = 15;
        tau = 0.025;
    end

    % Discretization parameters
    T = 1;
    t = 0;
    max_tau = tau; % Maximum allowed tau
    min_tau = 0.001; % Minimum allowed tau for stability

    % Initial conditions
    uold = 1;
    vold = 0.1; 

    % Initialize storage for solutions
    t_values = t;
    u_values = uold;
    v_values = vold;

    % Define the function g(t)
    g = @(t) -sin(4 * t);

    % Newton's method parameters
    tolerance = 1e-8;
    maxiter = 100;

    % Main loop for time stepping with variable tau
    while t < T
        unext = uold; % Initial guess for u
        iter = 0;

        % Newton's method for solving u
        while iter < maxiter
            F = unext - uold - tau * g(t + tau) + tau * alpha * unext^3 + c * tau * (uold - vold);
            F_prime = 1 + 3 * tau * alpha * unext^2;

            % Newton's update
            unextnew = unext - F / F_prime;

            % Check for convergence
            if abs((unextnew - unext) / (1e-5 + abs(unext))) < tolerance
                break; % Convergence achieved
            end  

            % Update guess for next iteration
            unext = unextnew;
            iter = iter + 1;
        end

        % Adaptive time-stepping: Adjust tau based on Newton solver performance
        if iter >= maxiter / 2  % If Newton took many iterations, reduce tau
            tau = max(tau * 0.5, min_tau);
        elseif iter < maxiter / 4 && tau < max_tau % If convergence was fast, increase tau
            tau = min(tau * 1.2, max_tau);
        end

        % Store new value of u after convergence
        unew = unextnew;

        % Update V
        vnew = (1 / (1 + c * tau)) * (c * tau * unew + vold);

        % Update time and append results
        t = t + tau;
        t_values = [t_values, t];
        u_values = [u_values, unew];
        v_values = [v_values, vnew];

        % Update old values for next time step
        uold = unew;
        vold = vnew;

        % Check if reached end of T
        if t >= T
            break;
        end
    end    

    % Plot the results
    figure;
    plot(t_values, u_values, 'b', 'LineWidth', 1.5, 'DisplayName', 'u(t)');
    hold on;
    plot(t_values, v_values, 'r--', 'LineWidth', 1.5, 'DisplayName', 'v(t)');
    xlabel('Time');
    ylabel('Solution');
    legend('u(t)', 'v(t)');
    title(sprintf('Case alpha = %g, c = %g, maxdt = %.3f', alpha, c, max_tau));
    % xlim([0 1]);
    % ylim([0 1]);
   
end
