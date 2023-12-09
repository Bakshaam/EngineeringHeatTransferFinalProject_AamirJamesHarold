%% Aamir Baksh, James Shaw, Harold Griffin                  Heat Transfer Final Project
clc; clear;

%% Question 4

%% Initializing Variables
m = 10;                                  %kg/s           Flow rate of Silica Particles
D_ext = 16;                              %m              External diameter
alpha = .95;                             %Ratio          Absorbtivity
h = 90;                                  %m              Tower height
Cp = 1130;
U_inf = 1.5;                             %m/s            Wind speed
L = 6;             %VARIABLE             %m              Reciever Height
t = 0.001;                               %m              Wall thickness

T_air = 300;                             %Kelvin         Assumed Air Temperature
T_in = 750 + 273.15;                     %Kelvin         Inlet temperature
T_out = 900 + 273.15;                    %Kelvin         Assumed outlet temperature
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

%% Bending tip calcuations:
H = 0.254;                                      % height of the channel in (m)
w = (3/10)*H;                                     % width of the channel in (m)
% prompt = "What is the bend radius in (m)? "
rb = 0.002;                                     % radius of the bend in (m)
w1 = 2*rb;                                      % width of the subtracted tip in (m)
H1 = rb*(10/1.5);                               % height of the subtracted tip in (m)
Ac_tip1 = (0.5)*(w1*H1);                        % area of the subtracted tip, including the round in (m^2)
delta_Ac_tip1 = (0.5*pi*rb^2);                  % area of the rounded section maintained in channel area in (m^2)
Ac_tip = Ac_tip1 - delta_Ac_tip1;               % Total cross sectional area removed after bend radius rb is applied in (m^2)


A_channel = ((0.5*H*w))-Ac_tip;       %m^2            Cross-sectional area of a single channel

S2 = sqrt(H1^2 + (rb)^2)- rb*tan(0.0872665);     %m              Side length that is lost when rounding the corner - rb*tan(0.0872665)
S1 = sqrt(H^2 + (w/2)^2);                   %m              Original side length
arclength = 2*pi*rb*(170/360);              %m              Arc length of the rounded corner
P_channel = 2*(S1-S2) + w + arclength;      %m              Perimeter for channel walls
D_eff = 4 * (A_channel/P_channel);          %m              Effective diameter
A_in = pi * D_eff * L;                      %m^2
h_in = 913.6;                               %W/m^2*K        Silica particles heat transfer coeff (from linear interpolation)

%% CALCULATE NEW SURFACE AREA

A_removed = 2*sqrt((H1)^2+rb^2);
A_add = pi*rb;

A_s1 = (S1*2)*L;                              %m^2            Exposed surface area
A_s = A_s1 - A_removed*L + A_add*L;

numPanels = round((D_ext*pi)/w);
totalSurfaceArea = numPanels * A_s;
Q_solar_recieved = Q_solar_recieved * totalSurfaceArea;
Q_solar_recieved = Q_solar_recieved / numPanels;


while (abs(error) > 1)

    %% Loop to find mass flow rate needed to maintain 900 degrees T_out

    % Finding wall temperatures
    Q_absorb_in = m * Cp * (T_out - T_in);
    T_wall_in = (Q_absorb_in / (A_s * h_in)) + T_f;

    %Conduction effects based on wall thickness 1 mm
    k = 22.525;                 %Interpolated Thermal conductivity for T_f for Inconel
    T_wall_out = (t*Q_absorb_in)/(k*A_s) + T_wall_in;

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

    %Adjusting compared to expected
    error = (Q_solar_recieved - Q_solar_calc);
    if (abs(error) > 1)
        if (error > 0)
            m = m + rand(1)*(Q_solar_recieved / Q_solar_calc);
        else
            m = m - rand(1)*(Q_solar_recieved / Q_solar_calc);
        end
        %T_f = (T_out + T_in)/2;
    end
    error;
    count = count + 1;

end
m

