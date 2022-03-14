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

%Bus Power Information
P = [0 0 0 0 1 1.1 1]; 
Q = [0 0 0 0 .65 .50 .70]; 

%Flat Start Information
V_o = [1.0 1.0 1.0 1.0 1.0 1.0 1.0];
delta = [0 0 0 0 0 0 0];

%Number of Buses
N = 7; 

%The amount of iterations to do
number = 4;

%The convergence criteria in pu
criteria = 0.0001;

%Set iteration counter
counter = 0;

%----------Execution----------%

%Find the Ybus
Ybus = find_Ybus(line_data, partridge_cond, line_length, w, j);

%First Iteration from Flat Start
[newV, newDelta, newcounter, met] = iteration(Ybus, V_o, delta, P, Q, N, criteria, counter);
disp("Converge?: " + met)

%Continue until convergance or max iterations
while(newcounter < number && met == false)
    [newV, newDelta, newcounter, met] = iteration(Ybus, newV, newDelta, P, Q, N, criteria, newcounter);
end

%Find Current in the Lines
Y_lines = find_line_admittances(line_data, partridge_cond, line_length, w, j);
current_lines = find_current_lines(Y_lines, newV, newDelta);
disp("Current in Lines")
disp(current_lines)

%Find Power in Lines
power_lines = find_power_lines(current_lines, newV);

%Find Power Loss
power_loss = find_power_loss(Ybus, V_o, delta, P, Q, N);
disp("Power Loss:")
disp(power_loss)


%---Functions---%

function Deq = find_Deq(data)
    Deq = (data.Dab*data.Dbc*data.Dca)^(1/3);
end

function D_SL = find_D_SL(cond) 
    if (cond.nb == 2) 
        D_SL = sqrt(cond.GMR * cond.d);
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

function Y_lines = find_line_admittances(line_data, partridge_cond, line_length, w, j)
    [Z1, Z2, Z3, Z4, Z5] = find_sequence_impedance(line_data, partridge_cond, line_length, w, j);
    Y_lines(1) = 1/Z1;
    Y_lines(2) = 1/Z2;
    Y_lines(3) = 1/Z3;
    Y_lines(4) = 1/Z4; 
    Y_lines(5) = 1/Z5;
end

function Ybus = find_Ybus(line_data, partridge_cond, line_length, w, j)
    %---Find Sequence Impedance---%

    [Z1, Z2, Z3, Z4, Z5] = find_sequence_impedance(line_data, partridge_cond, line_length, w, j);
    Y1 = 1/Z1;
    Y2 = 1/Z2;
    Y3 = 1/Z3;
    Y4 = 1/Z4; 
    Y5 = 1/Z5;

    %---Find Shunt Admittance---%
    %Even of medium length lines

    [Y1_shunt, Y2_shunt, Y3_shunt, Y4_shunt, Y5_shunt] = find_shunt_admittances(line_data, partridge_cond, line_length, w, j);

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

end

function J1_kn = find_J1_kn_element(V_o, Ybus, delta, k, n)
    
    J1_kn = V_o(k) * abs(Ybus(k,n)) * V_o(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
    
end 

function J2_kn = find_J2_kn_element(V_o, Ybus, delta, k, n)
    
    J2_kn = V_o(k) * abs(Ybus(k,n)) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
    
end 

function J3_kn = find_J3_kn_element(V_o, Ybus, delta, k, n)
    
    J3_kn = -1 * V_o(k) * abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
    
end 

function J4_kn = find_J4_kn_element(V_o, Ybus, delta, k, n)
    
    J4_kn = V_o(k) * abs(Ybus(k,n)) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
    
end 

function J1_kk = find_J1_kk_element(V_o, Ybus, delta, k, N)
    result = 0;
    for n = 1:N
        
        num = abs(Ybus(k,n)) * V_o(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
        if n == k
            num = 0;
        end
        result = result + num;
    end
    
    J1_kk = -1 * V_o(k) * result;
    
end 

function J2_kk = find_J2_kk_element(V_o, Ybus, delta, k, N)
    result = 0;
    for n = 1:N
        num = abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
        result = result + num;
    end
    
    J2_kk = V_o(k) * abs(Ybus(k,k)) * cos(angle(Ybus(k,k))) + result;
    
end 

function J3_kk = find_J3_kk_element(V_o, Ybus, delta, k, N)
    result = 0;
    for n = 1:N
        num = abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
        if n == k
            num = 0;
        end
        result = result + num;
    end
    
    J3_kk = V_o(k) * result;
    
end 

function J4_kk = find_J4_kk_element(V_o, Ybus, delta, k, N)
    result = 0;
    for n = 1:N
        num = abs(Ybus(k,n)) * V_o(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
        result = result + num;
    end
    
    J4_kk = -1 * V_o(k) * abs(Ybus(k,k)) * sin(angle(Ybus(k,k))) + result;
    
end 

function J1 = find_J1(V_o, Ybus, delta, N)
    
    J1 = zeros(N,N);
    
    for k = 1:N
        for n = 1:N
           if k == n 
               J1(k,k) = find_J1_kk_element(V_o, Ybus, delta, k, N);
           else
               J1(k,n) = find_J1_kn_element(V_o, Ybus, delta, k, n);
           end
        end
    end 
    
    
end

function J2 = find_J2(V_o, Ybus, delta, N)
    
    J2 = zeros(N,N);
    
    for k = 1:N
        for n = 1:N
           if k == n 
               J2(k,k) = find_J2_kk_element(V_o, Ybus, delta, k, N);
           else
               J2(k,n) = find_J2_kn_element(V_o, Ybus, delta, k, n);
           end
        end
    end 
    
    
end

function J3 = find_J3(V_o, Ybus, delta, N)
    
    J3 = zeros(N,N);
    
    for k = 1:N
        for n = 1:N
           if k == n 
               J3(k,k) = find_J3_kk_element(V_o, Ybus, delta, k, N);
           else
               J3(k,n) = find_J3_kn_element(V_o, Ybus, delta, k, n);
           end
        end
    end 
    
    
end

function J4 = find_J4(V_o, Ybus, delta, N)
    
    J4 = zeros(N,N);
    
    for k = 1:N
        for n = 1:N
           if k == n 
               J4(k,k) = find_J4_kk_element(V_o, Ybus, delta, k, N);
           else
               J4(k,n) = find_J4_kn_element(V_o, Ybus, delta, k, n);
           end
        end
    end 
    
    
end

function J = find_jacobian(V_o, Ybus, delta, N) 


    J1 = find_J1(V_o, Ybus, delta, N);
    J2 = find_J2(V_o, Ybus, delta, N);
    J3 = find_J3(V_o, Ybus, delta, N);
    J4 = find_J4(V_o, Ybus, delta, N);
    
    J = [J1 J2; J3 J4];
    
    %Delete rows and columns
    J(1, :) = [];
    J(7, :) = [];
    J(9, :) = [];
    J(:, 1) = [];
    J(:, 7) = [];
    J(:, 9) = [];

end

function mismatch = find_mismatch_matrix(Ybus, V_o, delta, P, Q, N)
    
    result = 0;
    mismatch_P = zeros(7, 1);
    
    for k = 1:N
        for n = 1:N
            num = abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
            result = result + num; 
        end
            mismatch_P(k,1) = V_o(k) * result;
    end
   
    
    result = 0;
    mismatch_Q = zeros(7, 1);
    
    for k = 1:N
        for n = 1:N
            num = abs(Ybus(k,n)) * V_o(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
            result = result + num;
        end
            mismatch_Q(k,1) = V_o(k) * result;
    end
    
    mismatch_x = [mismatch_P ; mismatch_Q];
    
    expected = zeros(14,1);
    for i = 1:7
        expected(i, 1) = P(i);
    end
    for y = 8:14
        expected(y, 1) = Q(y-7);
    end
    
    m = expected - mismatch_x;
    m(1, :) = [];
    m(7, :) = [];
    m(9, :) = [];
    
    mismatch = m;
    
end

function [V_delta_mismatch, newV, newDelta] = find_V_delta_mismatch(mismatch, J, N, V_o, delta)
    
    newV = zeros(1,N);
    newDelta = zeros(1,N);
    mismatch_new = zeros(14, 1);
    
    mismatch_new(1) = 0;
    mismatch_new(2) = mismatch(1,1);
    mismatch_new(3) = mismatch(2,1);
    mismatch_new(4) = mismatch(3,1);
    mismatch_new(5) = mismatch(4,1);
    mismatch_new(6) = mismatch(5,1); 
    mismatch_new(7) = mismatch(6,1);
    mismatch_new(8) = mismatch(7,1);
    mismatch_new(9) = 0;
    mismatch_new(10) = mismatch(8,1);
    mismatch_new(11) = mismatch(9,1);
    mismatch_new(12) = mismatch(10,1);
    mismatch_new(13) = 0;
    mismatch_new(14) = mismatch(11,1);

    V_delta_mismatch = inv(J) * mismatch;
    
    for i = 1:N
        newV(i) = mismatch_new(i,1) + V_o(i);
    end
    for k = 1:N
       newDelta(i) = mismatch_new(k+N,1) + delta(k);
    end
    
end 

function [newV, newDelta, newcounter, met] = iteration(Ybus, V_o, delta, P, Q, N, criteria, counter)
        
        mismatch = find_mismatch_matrix(Ybus, V_o, delta, P, Q, N);

        J = find_jacobian(V_o, Ybus, delta, N);

        [V_delta_mismatch, newV, newDelta] = find_V_delta_mismatch(mismatch, J, N, V_o, delta);
        
        met = find_convergance(V_delta_mismatch, criteria, N);
        
        print_results(V_delta_mismatch, mismatch, counter);
        
        newcounter = counter+1;
end

function met = find_convergance(mismatch, criteria, N)
    mettemp = true;
    disp("Mismatch")
    for i = 1:N
        if(-1 * mismatch(i,1) <= criteria)
            
            disp(mismatch(i,1))
        else
            mettemp = false;
        end
    end
    met = mettemp;
end

function print_results(V_delta_mismatch, mismatch, counter)
    disp("For Iteration " + counter)
    disp("The Voltage and Angle Mismatch is: ")
    disp(V_delta_mismatch)
    disp("The Power Mismatch is: ")
    disp(mismatch)
end

function current_lines = find_current_lines(Y_lines, newV, newDelta)
    current_lines = zeros(5,1);
    current_lines(1) = Y_lines(1) * (newV(2) - newV(5)); 
    
    %What is up with this angle
    if((angle(newDelta(2) - newDelta(5))-180) < 0)
        disp("Current Flowing from 5 to 2");
    end
    current_lines(2) = Y_lines(2) * (newV(3) - newV(5));
    if((angle(newDelta(3) - newDelta(5))-180) < 0)
        disp("Current Flowing from 3 to 5");
    end
    current_lines(3) = Y_lines(3) * (newV(3) - newV(6));
    if((angle(newDelta(3) - newDelta(6))-180) < 0)
        disp("Current Flowing from 6 to 3");
    end
    current_lines(4) = Y_lines(4) * (newV(5) - newV(7));
    if((angle(newDelta(3) - newDelta(6))-180) < 0)
        disp("Current Flowing from 7 to 5");
    end
    current_lines(5) = Y_lines(5) * (newV(6) - newV(7));
    if((angle(newDelta(3) - newDelta(6))-180) < 0)
        disp("Current Flowing from 7 to 6");
    end
    
    
    
end

function power_lines = find_power_lines(current_lines, newV)
    power_lines = zeros(5,1);
    power_lines(1) = current_lines(1) * (newV(2) - newV(5)); 
    power_lines(2) = current_lines(2) * (newV(3) - newV(5));
    power_lines(3) = current_lines(3) * (newV(3) - newV(6));
    power_lines(4) = current_lines(4) * (newV(5) - newV(7));
    power_lines(5) = current_lines(5) * (newV(6) - newV(7));    
    
end

function power_loss = find_power_loss(Ybus, V_o, delta, P, Q, N)
    result = 0;
    P_found = zeros(7, 1);
    
    for k = 1:N
        for n = 1:N
            num = abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
            result = result + num; 
        end
            P_found(k,1) = V_o(k) * result;
    end
   
    
    result = 0;
    Q_found = zeros(7, 1);
    
    for k = 1:N
        for n = 1:N
            num = abs(Ybus(k,n)) * V_o(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
            result = result + num;
        end
            Q_found(k,1) = V_o(k) * result;
    end
    
    mismatch_x = [P_found ; Q_found];
    
    expected = zeros(14,1);
    for i = 1:14
       expected(i) = P(1); 
    end
    
    m = expected - mismatch_x;
    
    power_loss = m;
end

