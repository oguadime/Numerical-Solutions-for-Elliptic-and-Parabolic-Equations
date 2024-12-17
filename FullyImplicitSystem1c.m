function FullyImplicitSystem1c(alpha, c, tau)
    if nargin < 3
        alpha = 4;
        c = 10;
        tau = 0.025;
    end

    % Discretization parameters
    T = 1;
    dt = tau;
    t = 0:dt:T;
    Nsteps = length(t);

    % Initial conditions
    uold = 1;
    vold = 0.1; 

    unew = zeros(1, Nsteps);
    vnew = zeros(1, Nsteps);
    newton_iters = zeros(1, Nsteps); % To store Newton iterations per step

    unew(1) = uold;
    vnew(1) = vold;

    % Define the function g(t)
    g = @(t) -sin(4 * t);

    % Newton's method parameters
    tolerance = 1e-15;
    maxiter = 100;

    % Main loop for time stepping
    for i = 2:Nsteps
        unext = uold; % Initial guess
        iter = 0;

        % Update U (using Newton's method) 
        while iter < maxiter
            F = unext - uold - dt*g(t(i)) + dt*alpha*unext^3 + c*dt*(unext - ((c*dt*unext) + vold)/(1 + c*dt));
            F_prime = 3*dt*alpha*unext^2 + (1+2*c*dt)/(1+c*dt);

            % Newton's update
            unextnew = unext - F/F_prime;

            % Check for convergence
            if abs((unextnew-unext)/(1e-5 + abs(unext))) < tolerance
                break; % convergence achieved
            end  

            % Update guess for next iteration
            unext  = unextnew;
            iter = iter + 1;
        end

        % Store the number of iterations for reporting
        newton_iters(i) = iter;

        % Store new value of u after convergence
        unew(i) = unextnew;

        % Update V
        vnew(i) = (1/(1+c*dt))*(c*dt*unew(i) + vold);

        % Update old values for next time step
        uold = unew(i);
        vold = vnew(i);
    end    

    % Plot the results
    % figure;
    % plot(t, unew, 'b', 'LineWidth', 1.5, 'DisplayName', 'u(t)');
    % hold on;
    % plot(t, vnew, 'r--', 'LineWidth', 1.5, 'DisplayName', 'v(t)');
    % xlabel('Time');
    % ylabel('Solution');
    % legend('u(t)', 'v(t)');
    % title(sprintf('Case alpha = %g, c = %g, maxdt = %.3f', alpha, c, tau));
    % xlim([0 1]);
    % ylim([0 1]);

    % Reporting on the solver performance
    avg_iters = mean(newton_iters(2:end)); % Average iterations per time step
    max_iters = max(newton_iters);         % Maximum iterations in any time step
    min_iters = min(newton_iters(2:end));  % Minimum iterations, excluding the first step
    fprintf('Average Newton iterations per step: %.2f\n', avg_iters);
    fprintf('Maximum Newton iterations in any step: %d\n', max_iters);
    fprintf('Minimum Newton iterations (excluding first step): %d\n', min_iters);
end
