function impedance_delta = impedance_difference(omega)
    R = 525;  
    C = 7e-5; 
    L = 3;    
    M = 75;

    if omega <= 0
        error('Omega musi być większa od zera.');
    end

    impedance = 1 / sqrt(1/(R^2) + (omega*C - 1/(omega*L))^2);

    impedance_delta = impedance - M;
end
