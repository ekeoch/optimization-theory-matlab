x = sym('x', [1 1]);
x_val = 1; %random starting point x =/= 0

y = (-3*x^2) + (12/(x^2)) + 2*exp(x^2);

gradient = diff(y);
hessian = diff(gradient);

g_val = double(subs(gradient, x, x_val));
h_val = double(subs(hessian, x, x_val));
step = NaN(1, 1);

e = 1e-2; 
iter = 0;
err = inf;
max_iter = 50;

while ((e < err) & (iter <= max_iter))
    rnk = rank(h_val);
    if (rnk == 1)
        step = h_val\g_val;
    else
        hess_p = h_val(:,1:rnk);
        step(1:rnk) = hess_p\g_val;
        step((rnk+1):1) = zeros(1-rnk,1);
    end
    
    ind = find(isnan(step) | isinf(step));
    
    if ~isempty(ind)
        break;
    end
    
    x_val = x_val - step; % new value of X
    g_val = double(subs(gradient, x, x_val));  % compute gradient
    h_val= double(subs(hessian, x, x_val)); % compute hessian
    err = ((g_val')*g_val)^0.5; % compute norm of gradient
    iter = iter + 1;
end

fprintf('optimal point by newton: %s\n\n', x_val);

opt_solution = double(subs(y, x, x_val));

fprintf('optimal solution by newton: %s\n\n', opt_solution);





