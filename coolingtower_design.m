% ================================================================================
% ==                                                                            ==
% == Title:             Cooling towers design                                   ==
% == Problem Author:    Prof. Alessandro Di Pretoro (SuPER,POLIMI)              ==
% ==                    alessandro.dipretoro@polimi.it                          ==
% ==                    https://super.chem.polimi.it/ing-alessandro-di-pretoro  ==
% == Solution:          @ahmedsabriz                                            ==
% ==                                                                            ==
% ================================================================================

clc
close all
clear vars

% Physical Properties
% L refers to liquid water
% G refers to evaporated water
L_in=18000;             % [kg/h]
G_in=13000;             % [kg/h]
TL_in=60;               % [°C]
TL_out=25;              % [°C]
TG_in=22;               % [°C]
P=760;      		 	% [mmHg]
Z=0.6;      			% [-] relative humidity
Cu=0.26;        		% [kcal/kg*K]
CpL=1;          		% [kcal/kg*K]
DHev=580;       		% [kcal/kg]

% Heat transfer coefficients per unit of length
hga=5000;       		% [kcal/m/h/K] 
hla=28000;              % [kcal/m/h/K]

% Vapor Pressure Correlation
A=-49.705;              % Pev 1st coeff	mmHg vs °C
B=2.71;                 % Pev 2nd coeff mmHg vs °C

% Humidity
Ps=A+B*TG_in;           % [mmHg]
Us=0.62*Ps/P;           % [mol_w / mol_air]
U_in=Z*Us;              % [mol_w / mol_air]
Gdry=G_in/(1+U_in);     % [kg/h]

% Mass transfer Coefficient from Lewis equation
Kua=hga/Cu;             % [kg/m/h]

% Countercurrent iterations
max_iterations = 100;
step = 0.005;
height_guess = 10;
L_out_guess = L_in;     % First guess value for water outlet flowrate
tollerance = 0.5;

% Integration
for i = 1:max_iterations
    zspan=0:step:height_guess;
    y0 = [U_in L_out_guess G_in TL_out TG_in];
    [z, y] = ode15s(@(z,y) odefun(z, y, hla, hga, Kua, DHev, A, P, B, Gdry, CpL, Cu), zspan, y0);
    
    % Assign output to vectors
    U   = y(:,1);
    L   = y(:,2);
    G   = y(:,3);
    TL  = y(:,4);
    TG  = y(:,5);
    
    % Determine the index I of zspan where the TL_in condition is met
    [M, I] = min(abs(TL_in - TL)); 
    height_calculated = I * step;
    L_in_calculated = L(I);
    error = L_in - L_in_calculated;

    % Print Iterations
    fprintf('Iteration: %03i - L_out_guess: %.0f - L_in_calculated: %.0f - Error: %.1f \n', ...
                 i, L_out_guess, L_in_calculated, error);
    
    % Monitor convergence
    if abs(error) < tollerance
       break
    end

    L_out_guess = L_out_guess + 0.5 * error;

end

fprintf('Packing Height: %.1f\n', height_calculated);

% Trimming result vectors:
z   = z(1:I);
U   = U(1:I);
L   = L(1:I);
G   = G(1:I);
TL  = TL(1:I);
TG  = TG(1:I);

% Plots
figure(1)
plot(z,TL,'-b', z,TG,'-r','LineWidth',1)
title("Temperature Profile Along the Column")
legend("Water Temperature", "Air Temperature")
xlabel("Height [m]")
ylabel('Temperature [°C]')

figure(2)
plot(z,G,'-b', z,L,'-r','LineWidth',1)
title("Flowrates Along the Column")
legend("Air Flow", "Water Flow")
xlabel("Height [m]")
ylabel('Flowrate [kg/h]')
ytickformat('%d')
ax = gca;
ax.YAxis.Exponent = 0;

function dy = odefun(z, y, hla, hga, Kua, DHev, A, P, B, Gdry, CpL, Cu)
% Assigning initial condition parameters
U   = y(1);
L   = y(2);
G   = y(3);
TL  = y(4);
TG  = y(5);

% Temperature and Humidity at interface
% From heat balance on the adiabatic system
Ti=(hla*TL+hga*TG+Kua*DHev*(U-0.62*A/P))/(hla+hga+Kua*DHev*0.62*B/P);
Ui=0.62*(A+B*Ti)/P;

% ODEs
dU  = Kua*(Ui-U)/Gdry;
dL  = Kua*(Ui-U);
dG  = Kua*(Ui-U);
dTL = hla*(TL-Ti)/L/CpL;
dTG = hga*(Ti-TG)/Gdry/Cu;

dy = [dU; dL; dG; dTL; dTG];
end
