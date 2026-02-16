clear
% Fluid Properties (Gnielinski)
T_g = 1282.32;
Pr = 0.6133;
mu = 3.5192e-5;
k_g = 0.4615;
Pc = 2.068e6;
R = 8314/4.032;
rhog = Pc/(R*T_g);
gam = 1.336;
D = 0.0087;
u = 0.018/(rhog*pi*(D/2)^2);
Ma = u/1879.8;

% Wall & Boundary Properties 
L = 0.0443738;
n = 1000;
dy = L / n;
h_nat = 5; 
T_amb = 300.0;
k_w = 328.84;
rho = 8885;
cp = 389.37;
tf = 5;
dt = 0.001;
phi = 300 * ones(n, 1);

% Convection Coefficients
Re = rhog * u * D / mu;
f = (0.79 * log(Re) - 1.64)^(-2);
Nu = ((f / 8) * (Re - 1000) * Pr) / (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1));
hG = Nu * k_g / D;
Taw = T_g * (1 + (gam - 1) / 2 * Ma^2 * Pr^(1/3));

% Tridiagonal Matrix Setup
DD = k_w / dy;
ap0 = rho * cp * dy / dt;
d = (ap0 + 2 * DD) * ones(n, 1);
d(1) = ap0 + DD + hG;
d(n) = ap0 + DD + h_nat;
a = -DD * ones(n - 1, 1);
b = -DD * ones(n - 1, 1);

% Video Recording Setup
writerObj = VideoWriter('Heat_Transfer_Simulation.mp4', 'MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

% Plotting Setup
fig = figure('Color', 'w');
x_axis = linspace(0, L, n);
hP = plot(x_axis, phi, 'r', 'LineWidth', 2);
grid on;
axis([0 L 280 Taw]);
xlabel('Depth (m)');
ylabel('Temperature (K)');

% Transient Loop 
for t = 0:dt:tf
    % Source term construction
    c = ap0 * phi;
    c(1) = c(1) + hG * Taw;
    c(n) = c(n) + h_nat * T_amb;
    
    % Thomas Algorithm Solver
    dp = d;
    cp = c;
    for i = 2:n
        m = b(i - 1) / dp(i - 1);
        dp(i) = dp(i) - m * a(i - 1);
        cp(i) = cp(i) - m * cp(i - 1);
    end
    
    phi(n) = cp(n) / dp(n);
    for i = n - 1:-1:1
        phi(i) = (cp(i) - a(i) * phi(i + 1)) / dp(i);
    end
    
    % Update plot and title
    set(hP, 'YData', phi);
    title(sprintf('Time: %.3f s | Wall Temp: %.2f K', t, phi(1)));
    drawnow expose;
    
    % Store Video Frames
    if mod(round(t/dt), 10) == 0
        frame = getframe(fig);
        writeVideo(writerObj, frame);
    end
end

% Finalize Video 
close(writerObj);
disp('Simulation complete. Video saved as Heat_Transfer_Simulation.mp4');