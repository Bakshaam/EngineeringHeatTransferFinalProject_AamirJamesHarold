%% Aamir Baksh, James Shaw, Harold Griffin                  Heat Transfer Final Project
clc; clear;

%% Question 1

%% Initializing Variables
m = 10;                                  %kg/s           Flow rate of Silica Particles
D_ext = 16;                              %m              External diameter
alpha = .95;                             %Ratio          Absorbtivity
h = 90;                                  %m              Tower height
Cp = 1130;
U_inf = 1.5;                             %m/s            Wind speed
L = 9;             %VARIABLE             %m              Reciever Height
H = 0.254;                                      % height of the channel in (m)
w = (3/10)*H;                                   % width of the channel in (m)
A_channel = (0.5*w*H);                   %m^2            Cross-sectional area of a single channel
P_channel = 23.22 / 39.37;               %m              Perimeter for channel walls
D_eff = 4 * (A_channel/P_channel);       %m              Effective diameter
A_in = pi * D_eff * L;                   %m^2
h_in = 913.6;                            %W/m^2*K        Silica particles heat transfer coeff (from linear interpolation)
A_s = 2*(10.11/39.37)*L;                 %m^2            Exposed surface area

T_air = 300;                             %Kelvin         Assumed Air Temperature
T_in = 750 + 273.15;                     %Kelvin         Inlet temperature
T_out = 1300;                            %Kelvin         Assumed outlet temperature
T_f = (T_out + T_in)/2;                  %Kelvin         Average fluid temperature

epsilon = 0.58;
sigma = 5.67*10^-8;                      %W/m^2*K^4

kine_visc = 1.57*10^-5;                  %m^2/s          Kinematic Viscocity of air at room temp
Pr = 0.728;
g = 9.8;

Q_solar_recieved = alpha * 1000000 * epsilon;
error = 100000;
errorPlot = [];
count = 1;

numPanels = round((D_ext*pi)/w);
totalSurfaceArea = numPanels * A_s;
Q_solar_recieved = Q_solar_recieved * totalSurfaceArea;
Q_solar_recieved = Q_solar_recieved / numPanels;


%% Loop to find actual T_out
while (abs(error) > 1)
    % Finding wall temperatures
    Q_absorb_in = m * Cp * (T_out - T_in);
    T_wall_in = (Q_absorb_in / (A_s * h_in)) + T_f;
    T_wall_out = T_wall_in;             %For this case, our T_wall_out is equal to our T_wall_in since we ignore conduction effects

    %Finding Q_loss_radiation
    Q_loss_rad = epsilon * sigma * A_s * ((T_wall_out)^4 - (T_air)^4);

    %Finding Q_loss_convection (forced)
    ReL = (U_inf * L) / kine_visc;
    if ReL<500000
        NuL1 = 0.664*(ReL^0.5) * (Pr^(1/3)); %Laminar flow over flat plate correlation
    else
        xc = 50000*kine_visc/U_inf;
        if xc/L>0.95
            NuL1 = 0.664*ReL^0.5*Pr^(1/3);
        else
            NuL1 = ((0.037*ReL^0.8)-871)*Pr^(1/3); %Mixed flow over flat plate correlation
        end
    end
    kair = 0.02638;
    h_forced = NuL1*kair/L;

    %Finding Q_loss_convection (natural)
    gPr = 0.75*Pr^0.5/(0.609+(1.221*Pr^0.5)+(1.238*Pr))^0.25;
    B = 1 / ((T_wall_out + T_air)/2);
    GrL = (g * B * (T_wall_out - T_air) * L^3) / kine_visc^2;
    NuL2 = (4/3) * (GrL/4)^0.25*gPr;
    h_natural = NuL2*kair/L;

    %Find hcombined and Qloss-conv
    h_combined = h_forced + h_natural;
    Q_loss_convection = A_s * h_combined * (T_wall_out - T_air);

    %Calculating energy balance
    Q_solar_calc = Q_loss_convection + Q_loss_rad + Q_absorb_in;

    error = (Q_solar_recieved - Q_solar_calc);
    errorPlot(count) = error;

    %Adjusting compared to expected
    if (abs(error) > 1)
        if (error > 0)
            T_out = T_out + rand(1);
        else
            T_out = T_out - rand(1);
        end
        T_f = (T_out + T_in)/2;
    end
    error
    count = count + 1;

end

T_out

