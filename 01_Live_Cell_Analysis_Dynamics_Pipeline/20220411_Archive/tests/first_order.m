function dxdt = first_order(x)

    y = x(2);
    u = x(1);
    z = x(3);
    
    tau = 5;
    K = 2;
    
    dydt = (-y + K*u)/tau;
    dzdt = (-z + y)/tau;
    
    dxdt(1) = dydt;
    dxdt(2) = dzdt;
    
end
