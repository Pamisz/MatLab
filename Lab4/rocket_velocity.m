function velocity_delta = rocket_velocity(t)
% velocity_delta - różnica pomiędzy prędkością rakiety w czasie t oraz zadaną prędkością M
% t - czas od rozpoczęcia lotu rakiety dla którego ma być wyznaczona prędkość rakiety
M = 750; % [m/s]
g = 1.622; %[m/s^2]
m0 = 150000; %[kg]
u = 2000; %[m/s]
q = 2700; %[kg/s]

if t <= 0
    error("Nieprawidłowa wartość czasu.")
end

velocity_delta = u * log(m0/(m0-q*t))-g*t;
velocity_delta = velocity_delta - M;

end