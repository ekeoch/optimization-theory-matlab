x = sym('x', [1 1]);
l_val = 0;
r_val = 10;

eps = 1e-3;
lam = 1e-2;

e_check = inf;
abs_grad = inf;

y = (-3*x^2) + (12/(x^2)) + 2*exp(x^2);

grad = diff(y);

while ( (e_check > eps) & (abs_grad > lam))
    x_val = (l_val + r_val)/2;
    grad_x = double(subs(grad, x, x_val));
    e_check = r_val - l_val;
    abs_grad = abs(grad_x);
    
    if(grad_x >= 0)
        r_val = x_val;
    else
        l_val = x_val;
    end
    
end

fprintf('optimal point by bisection: %s\n\n', x_val);

opt_solution = double(subs(y, x, x_val));

fprintf('optimal solution by bisection: %s\n\n', opt_solution);


