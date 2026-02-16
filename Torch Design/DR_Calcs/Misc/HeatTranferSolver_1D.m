%% 1-D Convection Heat Tranfer Solver %%
% Author: Rafael Macia Titos %

%%
%contour = importdata("contour.csv");        
%x = contour(:,1)./1000;                            %Distance From Injector Plate [m]
%r = contour(:,2)./1000;     
clear

INtoM = 0.0254;

% Fluid Properties 
T_g = 1093.65;                  %Gas Temperature (k)
Pr = 0.6064;                    %Prandlt Number
mu = 3.0790e-5;                 %Gas Dynamic Viscosity (Pa*s)
k_g = 0.403589;                 %Gas Thermal Conductivity (W/(m*k)
Pc = 1.109812725e6;             %Combustion Chamber Pressure (Pa)
MW = 4.032;                     %Molecular Weight   (kg/kmol)
R = 2062.003968253968;          %Gas Constant   (J/(k*kg))
cp_g = 7.9487632e3;             %Gas Specific Heat   (J/(k*kg))
rhog = 0.49206;                 %Gas Density    (kg/m^3)
gam = 1.3503;                   %Specific heat ratio
D = 0.296 * INtoM;              %Throat Diameter (m)
A = pi*(D/2)^2;                 %Throat Area (m^2)
Ma = 1.0;                       %Mach Number
u = 1745.02176071493;           %Gas Velocity (m/s)
Re = rhog * u * D / mu;         %Reynolds Number


% Extra Properties
Throat_D = 0.296 * INtoM;           %Throat Diameter (m)
Throat_A = pi*(Throat_D/2)^2;       %Throat Area (m)
Cstar = 1353.4928267626089564;      %Characteristic Velocity (m/s)

% Wall & Boundary Properties 
L = 0.323 * INtoM;                  %Wall Thickness (m)
n = 1000;                           %Wall Section
dy = L / n;                         %wall Section Size
h_nat = 25;                          %Natural Convection Heat Transfer Coefficient
T_amb = 300.0;                      %Ambient Temperature
k_w = 389.37;                         %Wall Thermal Conductivity  (W/(m*k)
rho = 8885;                         %Wall Density    (kg/m^3)
cp = 420;                           %Wall Specific Heat (J/(k*kg))
tf = 5;                             %Simulation Time    (s)
dt = 0.0005;                        %Time Step (s)
phi = 298 * ones(n, 1);             %Initial Wall Temperature   (k)
T_melt = 1673;                      %Wall Melting Temperature   (k)

% Convection Coefficients (Gnielinski) [Only use under certain Conditions]
%f = (0.79 * log(Re) - 1.64)^(-2);   %Friction Factor
%Nu = ((f / 8) * (Re - 1000) * Pr) / (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1));  %Nusselt Number
% hG = Nu * k_g / D;          %Forced Convection Heat Transfer Coeffienct

% Convection Coefficients (Bartz)
hG_base = (0.026 / Throat_D^0.2) * ...
              ((mu^0.2 * cp_g) / (Pr^0.6)) * ...
              (Pc / Cstar)^0.8*(Throat_A / A)^0.9; 
Taw = T_g * (1 + (gam - 1) / 2 * Ma^2 * Pr^(1/3));  %Adiabatic Wall/Recovery Temperature


% Plotting & Video Setup
writerObj = VideoWriter('Heat_Transfer_Simulation.avi','Motion JPEG AVI');
writerObj.FrameRate = 30;
open(writerObj);

fig = figure('Color', 'w');
x_axis = linspace(0, L, n);
hP = plot(x_axis, phi, 'r', 'LineWidth', 2); 
grid on;
axis([0 L 280 3200]); 
yline(1673, '--k', ['Melting Point (Steel 304)'], 'LabelVerticalAlignment', 'bottom');
xlabel('Depth (m)');
ylabel('Temperature (K)');

% Matrix Pre-calculations
DD = k_w / dy;
ap0 = rho * cp * dy / dt;

%Melting Flag
Melt = 0;

% Transient Loop
for t = 0:dt:tf
    Tw = phi(1);
    term_Ma = (1 + (gam-1)/2 * Ma^2);
    
    % Full Bartz Correction
    sigma = 1 / ( (0.5 * (Tw/T_g) * term_Ma + 0.5)^0.68 * (term_Ma)^0.12 );
    hG = hG_base * sigma;
    
    % Matrix update 
    d = (ap0 + 2 * DD) * ones(n, 1);
    d(1) = (ap0 / 2) + DD + hG;
    d(n) = (ap0 / 2) + DD + h_nat;
    a = -DD * ones(n - 1, 1);
    b = -DD * ones(n - 1, 1);
    
    % Source term
    c = ap0 * phi;
    c(1) = (ap0 / 2) * phi(1) + hG * Taw;
    c(n) = (ap0 / 2) * phi(n) + h_nat * T_amb;
    
    % Thomas Algorithm
    dp = d;
    cp_vec = c; 
    for i = 2:n
        m = b(i - 1) / dp(i - 1);
        dp(i) = dp(i) - m * a(i - 1);
        cp_vec(i) = cp_vec(i) - m * cp_vec(i - 1);
    end
    
    phi(n) = cp_vec(n) / dp(n);
    for i = n - 1:-1:1
        phi(i) = (cp_vec(i) - a(i) * phi(i + 1)) / dp(i);
    end
    
    % Visualization and Video Recording
    if mod(round(t/dt), 5) == 0
        set(hP, 'YData', phi);
        title(sprintf('Time: %.2fs | h_G: %.1f | Wall: %.1fK', t, hG, phi(1)));
        drawnow expose;
        writeVideo(writerObj, getframe(fig));
        
        % Check for melting
        if (phi(1) > T_melt) && (Melt == 0);
            fprintf('WARNING: Wall surface melting at t = %.3f s\n', t);
            Melt = 1;
        end
    end
end

close(writerObj);
disp('Simulation complete.');