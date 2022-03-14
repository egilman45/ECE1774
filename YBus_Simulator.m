%------------------------------------------------------%
%ECE1774 Advanced Power Systems Analysis
%Spring 2022
%------------------------------------------------------%

clc

%---Definitions---%

j = 1i;
w = 377;

%Defining conductors data
%Bundles of 2 Subconductors per phase 1.5' spacing
partridge_cond.r = 0.385; %ohm/mi
partridge_cond.GMR = 0.0217; %ft
partridge_cond.diam = 0.642/12; %ft 
partridge_cond.d = 1.5; %ft
partridge_cond.nb = 2; 

%Define Line Data in feet
line_data.Dab = 19.5;
line_data.Dbc = 19.5;
line_data.Dca = 39;

%Define Line Lengths
line_length.L1 = 10;
line_length.L2 = 25;
line_length.L3 = 20;
line_length.L4 = 10;
line_length.L5 = 35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Find Sequence Impedance---%

[Z1, Z2, Z3, Z4, Z5] = find_sequence_impedance(line_data, partridge_cond, line_length, w, j);
Y1 = 1/Z1;
Y2 = 1/Z2;
Y3 = 1/Z3;
Y4 = 1/Z4; 
Y5 = 1/Z5;
disp(Y1)
disp(Y2)
disp(Y3)
disp(Y4)
disp(Y5)


%---Find Shunt Admittance---%
%Even of medium length lines

[Y1_shunt, Y2_shunt, Y3_shunt, Y4_shunt, Y5_shunt] = find_shunt_admittances(line_data, partridge_cond, line_length, w, j);

disp(Y1_shunt)
disp(Y2_shunt)
disp(Y3_shunt)
disp(Y4_shunt)
disp(Y5_shunt)

G1_y = 0;
G2_y = 0;
T1_y = 1.46-14.63j;
T2_y = 1.582+18.98j;

Y_12 = -1 * T1_y;
Y_11 = G1_y + T1_y;
Y_34 = -1 * T2_y;
Y_22 = (Y1_shunt/2) + Y1 + T1_y;
Y_33 = (Y2_shunt/2) + Y2 + (Y3_shunt/2) + Y3 + T2_y;
Y_44 = T2_y + G2_y;
Y_55 = (Y1_shunt/2) + Y1 + (Y2_shunt/2) + Y2 + (Y4_shunt/2) + Y4;
Y_25 = -1 * Y1;
Y_35 = -1 * Y2;
Y_36 = -1 * Y3;
Y_57 = -1 * Y4;
Y_67 = -1 * Y5;
Y_66 = (Y3_shunt/2) + Y3 + (Y5_shunt/2) + Y5;
Y_77 = (Y4_shunt/2) + Y4 + (Y5_shunt/2) + Y5;



%Defining Ybus matrix
Ybus = [Y_11 Y_12 0 0 0 0 0;
        Y_12 Y_22 0 0 Y_25 0 0;
        0 0 Y_33 Y_34 Y_35 Y_36 0;
        0 0 Y_34 Y_44 0 0 0;
        0 Y_25 Y_35 0 Y_55 0 Y_57;
        0 0 Y_36 0 0 Y_66 Y_67;
        0 0 0 0 Y_57 Y_67 Y_77];
    
Ybus

%---Functions---%


function Deq = find_Deq(data)
    Deq = (data.Dab*data.Dbc*data.Dca)^(1/3);
end

function D_SL = find_D_SL(cond) 
    if (cond.nb ==1)
       D_SL = cond.GMR; 
    end
    if (cond.nb == 2) 
        D_SL = sqrt(cond.GMR * cond.d);
    end 
    if(cond.nb == 3)
       D_SL = (cond.GMR * (cond.d)^2)^(1/3);
    end
    if (cond.nb == 4)
       D_SL = (1.091 *((cond.GMR * (cond.d)^3)^(1/4)));
    end
end 

function D_SC = find_D_SC(cond)
    D_SC = sqrt((cond.diam/2) * cond.d);
end

function R = find_resistance(cond, length) 
    %length should be in miles
    R = (cond.r / cond.nb) * length;
end

function X = find_reactance(w, Deq, D_SL,length)
    %Length should be in meters
    new_length = length * 1609.34;
    X = w * (2e-7) * log(Deq/D_SL) * new_length;
end

function Y = calc_shunt_line(length, D_SC, Deq, j, w)
    
    Y = ((j * w * 2 * pi * (8.8541878e-12)) / (log(Deq/D_SC))) * 1609.34 * length;

end 

function [Z1pu, Z2pu, Z3pu, Z4pu, Z5pu] = find_sequence_impedance(line_data, partridge_cond, line_length, w, j)
    
    %Find Deq
    Deq = find_Deq(line_data);
    
    %Find Equivalent GMR
    D_SL = find_D_SL(partridge_cond);

    %Find resistance
    R1 = find_resistance(partridge_cond, line_length.L1);
    R2 = find_resistance(partridge_cond, line_length.L2);
    R3 = find_resistance(partridge_cond, line_length.L3);
    R4 = find_resistance(partridge_cond, line_length.L4);
    R5 = find_resistance(partridge_cond, line_length.L5);

    %printR_results(R1, R2, R3, R4, R5);

    %Find Reactances
    X1 = find_reactance(w, Deq, D_SL, line_length.L1);
    X2 = find_reactance(w, Deq, D_SL, line_length.L2);
    X3 = find_reactance(w, Deq, D_SL, line_length.L3);
    X4 = find_reactance(w, Deq, D_SL, line_length.L4);
    X5 = find_reactance(w, Deq, D_SL, line_length.L5);

    %printX_results(X1, X2, X3, X4, X5);
    
    %Build Impedances
    Z1 = R1 + (j * X1);
    Z2 = R2 + (j * X2);
    Z3 = R3 + (j * X3);
    Z4 = R4 + (j * X4);
    Z5 = R5 + (j * X5);
    
    %Zbase Calculation
    Zbase = (230000)^2 / (100e6);
    
    %Per Unit of Impedance
    Z1pu = Z1 / Zbase;
    Z2pu = Z2 / Zbase;
    Z3pu = Z3 / Zbase; 
    Z4pu = Z4 / Zbase;
    Z5pu = Z5 / Zbase;
    
end

function [Y1, Y2, Y3, Y4, Y5] = find_shunt_admittances(line_data, partridge_cond, line_length, w, j)
    
    Deq = find_Deq(line_data);
    D_SC = find_D_SC(partridge_cond);
    
    Y1 = calc_shunt_line(line_length.L1, D_SC, Deq, j, w);
    Y2 = calc_shunt_line(line_length.L2, D_SC, Deq, j, w);
    Y3 = calc_shunt_line(line_length.L3, D_SC, Deq, j, w);
    Y4 = calc_shunt_line(line_length.L4, D_SC, Deq, j, w);
    Y5 = calc_shunt_line(line_length.L5, D_SC, Deq, j, w);
    
    
end 

