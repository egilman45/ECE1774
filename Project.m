

clear 
clc

%------DEFINITIONS AND DATA COLLECTION--------%

w = 377;
j = 1i;


%Find the Number of Nodes and Lines in the System
N = 7;
numLines = 7;

%Find Bus Connections
bus_connections = find_bus_connections(N);

%Have User Define the Connections
line_bus_connections = find_line_connections(N, numLines);

%Find the Bus Types
bus_type = find_bus_types();

%Line Data
conductor.r = 0.385; %ohm/mi
conductor.GMR = 0.0217; %ft
conductor.diam = 0.642/12; %ft 
conductor.d = 1.5; %ft
conductor.nb = 2; 

%Find More Line Data
[line_lengths, phase_distance] = find_line_lengths(numLines);

%Bus Power Information
P = [0 0 0 2 -1 -1.1 -1]; 
Q = [0 0 0 0 -.65 -.50 -.70]; 

%Flat Start
V = [1 1 1 1 1 1 1];
delta = [0 0 0 0 0 0 0];

N = 7;
flow_type = 1;
criteria = 0.0001;
counter = 0;
numIterations = 10;
met = false;

%Milestone 1----------------------------------------------------------

%--------SOLVING FOR YBUS---------%
%Find Transformer Information
[TX, TX_connections, numTX, TX_lines] = find_transformers();

%Find Admittance of the Lines
Y_lines = find_line_admittance(w, j, numLines, phase_distance, conductor, line_lengths, TX, numTX, TX_lines);

%Find the Shunt Admittance of the Lines
shunt = find_shunt_admittance(numLines, phase_distance, conductor, line_lengths, w, j);

%Find the Diagonal Elements of the Ybus
diagMatrix = find_diag_elements(line_bus_connections, Y_lines, shunt, N, numLines, TX, TX_connections, numTX);

%Find Outside Diaglonal Elements of Ybus
off_diag = find_off_diag_elements(line_bus_connections, Y_lines, numLines, N, TX, TX_connections, numTX);

%Find the Ybus
Ybus = find_Ybus_matrix(diagMatrix, off_diag, bus_connections, N);

% Ybus = [1.46883-j*14.6235 -1.46883+j*14.62351 0 0 0 0 0;
%        -1.46883+j*14.62351 27.06428-j*94.94624 0 0 -25.59545+j*80.34127 0 0;
%     	0 0 25.00760-j*91.23075	-1.59730+j*18.98608	-10.36387+j*32.11662 -13.04643+j*40.21160 0;	
% 		0 0 -1.59730+j*18.98608	1.59730-j*18.98608 0 0 0;		
% 		0 -25.59545+j*80.34127 -10.36387+j*32.11662 0 61.55477-j*192.71566 0 -25.59545+j*80.34127;
% 		0 0 -13.04643+j*40.21160 0 0 20.42355-j*63.05407 -7.37712+j*22.94457;
% 		0 0 0 0 -25.59545+j*80.34127 -7.37712+j*22.94457 32.97256-j*103.20234];

% Ybus_str = string(Ybus);
% Ybus_table = array2table(Ybus_str);
% % disp(Ybus_table)
% filename_output = 'Project_Outputs.xlsx';
% writetable(Ybus_table, filename_output, 'Sheet', 1, 'Range', 'B3');

%Newton Rasphon Power Flow
cap_set = 0;
testing_var_control = 0;
print = 1;
newton_rasp_power_flow(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, Y_lines, line_bus_connections, numLines, cap_set, numIterations, met, testing_var_control, print)


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

function bus_type = find_bus_types()
%     prompt = {'Label the Bus Types (1) Slack (2) PV (0) Other: : '};
%     dlgtitle = 'General System Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     bus_type = str2num(answer{1});

    bus_type = [1 0 0 2 0 0 0];
end

function [line_lengths, phase_distance] = find_line_lengths(numLines)
%     prompt = {'Enter Each Line Length (in miles) Seperated by a Space: '};
%     dlgtitle = 'Line Length Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     string = answer;
%     line_lengths = str2num(string{1});
%     
%     prompt = {'Enter the distance from Phase A to B (in feet): ', 'Enter the distance from Phase B to C (in feet): ', 'Enter the distance from Phase C to A (in feet): '};
%     dlgtitle = 'Phase Distance Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     phase_distance.Dab = str2num(answer{1});
%     phase_distance.Dbc = str2num(answer{2});
%     phase_distance.Dca = str2num(answer{3});

    line_lengths = zeros(1,numLines);
    line_lengths(1) = 10;
    line_lengths(2) = 25;
    line_lengths(3) = 20;
    line_lengths(4) = 10;
    line_lengths(5) = 35;
    phase_distance.Dab = 19.5;
    phase_distance.Dbc = 19.5;
    phase_distance.Dca = 39;
    
    
end

function conductor = find_phase_conductor_info(input)
   
    conductor.r = input(20,1); %ohm/mi
    conductor.GMR = input(20,2); %ft
    conductor.diam = input(20,3); %ft 
    conductor.d = input(20,4); %ft
    conductor.nb = input(20,5); 
end

function [X, R] = get_X_R_pu()
    prompt = {'Enter Each Line Resistance (R in pu) Seperated by a Space: ', 'Enter Each Line Reactance (X in pu) Seperated by a Space: '};
    dlgtitle = 'Per Unit Resistance Collection';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle,dims);
    string = answer;
    R = str2num(string{1});
    X = str2num(string{2});

end 

function Y_lines = find_line_admittance(w, j, numLines, phase_distance, conductor, line_lengths, TX, numTX, TX_lines)
    
%     prompt = {'Type 1 if you have the per unit X and R already, Type 2 if you have conductor information to input: '};
%     dlgtitle = 'Line Admittance Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     type = str2num(answer{1});
% 
%     if (type == 1)
%         [X, R] = get_X_R_pu();
%         Y_lines = find_line_admittances_with_pu(X, R, numLines);
%         line_lengths = 0;
%         phase_distance = 0; 
%         conductor = 0;
%     else 
%         %Find Length of Lines
%         [line_lengths, phase_distance] = find_line_lengths();
% 
%         %Find Conductor Information
%         conductor = find_conductor_info();
%         
%         Y_lines = find_line_admittances_with_line_data(phase_distance, conductor, line_lengths, w, j, numLines);
%     end

    Y_lines = find_line_admittances_with_line_data(phase_distance, conductor, line_lengths, w, j, numLines);
    for i = 1:numTX
        index = TX_lines(i);
        Y_lines(index) = TX(i);
    end
end

function Y_lines = find_line_admittances_with_pu(X, R, numLines) 

    Y_lines = zeros(1, numLines);
    for i = 1:numLines
       Y_lines(i) = 1/(R(i) + j*X(i));    
    end

end

function Ypu = find_line_admittances_with_line_data(phase_distance, conductor, line_lengths, w, j, numLines)
    %Find Deq
    Deq = find_Deq(phase_distance);
    
    %Find Equivalent GMR
    D_SL = find_D_SL(conductor);

    %Find resistance
    R = zeros(1,numLines);
    for i = 1:numLines
        R(i) = find_resistance(conductor, line_lengths(i));
    end

    %printR_results(R1, R2, R3, R4, R5);

    %Find Reactances
    X = zeros(1,numLines);
    for k = 1:numLines
        X(k) = find_reactance(w, Deq, D_SL, line_lengths(k));
    end

    %Build Impedances
    Z = zeros(1,numLines);
    for q = 1:numLines
        Z(q) = R(q) + (j * X(q));
    end

    
    %Get System Information for Zbase Calculation
%     prompt = {'What is the Voltage Base (kV) of the System? : ', 'What is the Apparent Power (S MVA) the System? : '};
%     dlgtitle = 'General System Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     V_sys = str2num(answer{1});
%     S_sys = str2num(answer{2});
    
    V_sys = 230;
    S_sys = 100;
    
    Zbase = (V_sys * 1000)^2 / (S_sys * 1e6);
    
    %Per Unit of Impedance
    Zpu = zeros(1,numLines);
    for p = 1:numLines
        Zpu(p) = Z(p) / Zbase;
    end
    
    %Convert to Admittance
    Ypu = zeros(1,numLines);
    for t = 1:numLines
        if(Zpu(t) ~= 0)
            Ypu(t) = (1/Zpu(t)); 
        else
            Ypu(t) = 0 + i*0;
        end
    end    
end


function bus_connections = find_bus_connections(N)
    
%     bus_connections = zeros(N,N);
%     string1 = 'Which Buses are Connected to Bus ';
%     string2 = ' Put in Order 1 for Connected and 0 for Not, Space Inbetween Each: ';
%     for i = 1:N
%         string3 = num2str(i);
%         string = strcat(string1, string2, string3);
%         prompt = {string};
%         dlgtitle = 'Bus Connection Data Collection';
%         dims = [1 50];
%         answer = inputdlg(prompt, dlgtitle,dims);
%         bus_connections(i,:) = str2num(answer{1});
%     end

    bus_connections = [1 1 0 0 0 0 0;
                       1 1 0 0 1 0 0;
                       0 0 1 1 1 1 0;
                       0 0 1 1 0 0 0;
                       0 1 1 0 1 0 1;
                       0 0 1 0 0 1 1;
                       0 0 0 0 1 1 1];

end

function line_bus_connections = find_line_connections(N, numLines)
%Fill In Matrix With Which Lines are Connected to Which Buses

    %Define Matrix of Line to Bus connections
%     line_bus_connections = zeros(N, numLines);
% 
%     for i = 1:N 
%        for k = 1:numLines
%           ansCon = input(sprintf("Is Line %d Connected to Bus %d ? y/n : ", k, i), 's'); 
%           if (ansCon == 'y') 
%              line_bus_connections(i,k) = 1; 
%           end
%        end
%     end

    line_bus_connections = [0 0 0 0 0 1 0;
                            1 0 0 0 0 1 0;
                            0 1 1 0 0 0 1;
                            0 0 0 0 0 0 1;
                            1 1 0 1 0 0 0;
                            0 0 1 0 1 0 0;
                            0 0 0 1 1 0 0];
end


function B = get_B_pu(numLines)
    prompt = {'Enter Each Line Shunt Admittance (B in pu) Seperated by a Space: '};
    dlgtitle = 'Shunt Admittance Collection';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle,dims);
    string = answer;
    B_not = str2num(string{1});
    B = zeros(1,numLines);
    for i = 1:numLines
        B(i) = B_not(i) * 0.5;
    end
end

function D_SC = find_D_SC(cond)
    D_SC = sqrt((cond.diam/2) * cond.d);
end

function Y_shunt_lines = find_shunt_admittances_with_line_data(phase_distance, conductor, line_lengths, j, w, numLines)
    Deq = find_Deq(phase_distance);
    D_SC = find_D_SC(conductor);
    
    Y_shunt_lines = zeros(1,numLines);
    for i = 1:numLines
        Y_shunt_lines(i) = calc_shunt_line(line_lengths(i), D_SC, Deq, j, w);
    end
    
    Z_shunt_lines = zeros(1,numLines);
    for j = 1:numLines
       Z_shunt_lines(j) = (1/Y_shunt_lines(j));
    end
    
    V_sys = 230;
    S_sys = 100;
    
    Zbase = (V_sys * 1000)^2 / (S_sys * 1e6);
    
    %Per Unit of Impedance
    Zpu = zeros(1,numLines);
    for p = 1:numLines
        Zpu(p) = Z_shunt_lines(p) / Zbase;
    end
    
    for k = 1:numLines
       Y_shunt_lines(k) = (1/Zpu(k));
    end
    
end

function Y_shunt = find_shunt_admittance(numLines, phase_distance, conductor, line_lengths, w, j)
%     prompt = {'Type 1 if you have the per unit B already, Type 2 if you want to use conductor information: '};
%     dlgtitle = 'Shunt Admittance Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     type = str2num(answer{1});
% 
%     if (type == 1)
%         Y_shunt = get_B_pu(numLines);
%     else 
%         Y_shunt = find_shunt_admittances_with_line_data(phase_distance, conductor, line_lengths, w, j, numLines);
%     end

    Y_shunt = find_shunt_admittances_with_line_data(phase_distance, conductor, line_lengths, w, j, numLines);
end

function [TX, TX_connections, numTX, TX_lines] = find_transformers()

%     prompt = {'Type How Many Transformers There? '};
%     dlgtitle = 'Number of Transformers';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     numTX = str2num(answer{1});
%     
%     TX = [1, numTX];
%     TX_connections = [2, numTX];
%     
%     if (numTX > 0) 
%         for i = 1:numTX
%             prompt = {'Type 1 if You Know the Impedance Already: '};
%             dlgtitle = 'Number of Transformers';
%             dims = [1 35];
%             answer = inputdlg(prompt, dlgtitle,dims);
%             ansTX = str2num(answer{1});
%             
%             if(ansTX == 1)
%                 prompt = {'Type Real Part then Space then Complex: '};
%                 dlgtitle = 'Transformer Input';
%                 dims = [1 35];
%                 answer = inputdlg(prompt, dlgtitle,dims);
%                 realTX = str2num(answer{1});
%                 complexTX = str2num(answer{2});
%             else
%                 prompt = {'What is MVA base of TX: ', 'What is MVA base of System: ', 'What is current pu: ', 'What is x/r Ratio: '};
%                 dlgtitle = 'Transformer Calculation';
%                 dims = [1 35];
%                 answer = inputdlg(prompt, dlgtitle,dims);
%                 MVA_TX = str2num(answer{1});
%                 MVA_base = str2num(answer{2});
%                 pu_old = str2num(answer{3});
%                 xr = str2num(answer{4});
%                 
%                 realTX = (pu_old) * (MVA_base / MVA_TX);
%                 complexTX = atan(xr);
%             end
%             
%             prompt = {'Type First Bus Then Space Then Second Bus: '};
%             dlgtitle = 'Transformer Bus Connections';
%             dims = [1 35];
%             answer = inputdlg(prompt, dlgtitle,dims);
%             bus1 = str2num(answer{1});
%             bus2 = str2num(answer{2});
%             
%             TX(i) = realTX + complexTX;
%             TX_connections(1,i) = bus1;
%             TX_connections(2,i) = bus2;
%             
%         end
%     end
    j = 1i;

    TX = [(1.463-j*14.63) (1.582-j*18.98)];
    TX_connections = [1 3; 2 4];
    TX_lines = [6 7];
    numTX = 2;
    
end

function diagMatrix = find_diag_elements(connection, Y_lines, shunt, N, numLines, TX, TX_connections, numTX)

    j = 1i;

    matrix = zeros(N, 1);
    sum = 0;
    for i = 1:N
        bus_info = connection(i,:);
        element1 = (bus_info .* Y_lines);
        element2 = (bus_info .* (j*shunt));
        for k = 1:numLines
            sum = sum + element1(k) + element2(k);
        end
        matrix(i) = sum;
        sum = 0;
    end
    
%     if(numTX > 0)        
%         for q = 1:numTX
%             bus = TX_connections(1,q);
%             matrix(bus) = matrix(bus) + TX(q);
%         end
%     end    
    diagMatrix = matrix;

end

function off_diag = find_off_diag_elements(line_bus_connections, Y_lines, numLines, N, TX, TX_connections, numTX)
    
    off_diag = zeros(N,N);

    for i = 1:numLines
        count = 0;
        while count < 1
            for k = 1:N
                if(line_bus_connections(k,i)==1)
                    bus1 = k;
                    count = count + 1;
                end
            end
        end
        
        count = 0;
        while count < 1
            for h = 1:N
                if(line_bus_connections(h,i)==1 && h ~= bus1)
                    bus2 = h;
                    count = count + 1;
                end
            end
        end
        
        off_diag(bus1, bus2) = -1 * Y_lines(i);
        off_diag(bus2, bus1) = -1 * Y_lines(i);
        
        if(numTX > 0)        
            for q = 1:numTX
                if(TX_connections(1,q) == bus1 && TX_connections(2,q) == bus2)
                    off_diag(bus1, bus2) = off_diag(bus1, bus2) + TX(q);
                    off_diag(bus2, bus1) = off_diag(bus2, bus1) + TX(q);
                end
            end
        end
        
    end
    
end

function Ybus = find_Ybus_matrix(diagMatrix, off_diag, bus_connections, N)
    Ybus = zeros(N,N);
    for i = 1:N
       for k = 1:N
          if(k==i)
             Ybus(k,k) = diagMatrix(k); 
          else 
             Ybus(i,k) = off_diag(i,k); 
          end
       end
    end
    
    Ybus = Ybus .* bus_connections;
    
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
    
    J2_kk = (V_o(k) * abs(Ybus(k,k)) * cos(angle(Ybus(k,k)))) + result;
    
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
    
    J4_kk = (-1 * V_o(k) * abs(Ybus(k,k)) * sin(angle(Ybus(k,k)))) + result;
    
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

function [J, delta_y] = find_jacobian(V, Ybus, delta, N, bus_type, flow_type, delta_y) 
    
    if(flow_type == 1)
        J1 = find_J1(V, Ybus, delta, N);
        J2 = find_J2(V, Ybus, delta, N);
        J3 = find_J3(V, Ybus, delta, N);
        J4 = find_J4(V, Ybus, delta, N);

        J = [J1 J2; J3 J4];
        
%         J_t = array2table(J);
%         disp(J_t)
        
        count = 0;
        %Find slack bus 
        for i = 1:N
           if (bus_type(i) == 1)
              %Delete row and column for slack bus
               J(i, :) = [];
               delta_y(i) = [];
               J(:, i) = [];
               count = count + 1;
               J((i+N-count), :) = [];
               delta_y((i+N-count)) = [];
               J(:, (i+N-count)) = [];
               count = count + 1;
           elseif(bus_type(i) == 2)
                J((i+N-count), :) = [];
                delta_y(i+N-count) = [];
                J(:, (i+N-count)) = [];
           end
        end
    end
end

function delta_y = find_delta_y(Ybus, V, delta, N, P, Q, bus_type)
    Px = zeros(1,N);
    Qx = zeros(1,N);
    
    for k = 1:N
        if(bus_type(k) == 1)
          Px(k) = 0;
          Qx(k) = 0;
        elseif(bus_type(k) == 2)
          Px(k) = find_Pk(Ybus, V, delta, N, k); 
          Qx(k) = 0;
        else
        Px(k) = find_Pk(Ybus, V, delta, N, k);
        Qx(k) = find_Qk(Ybus, V, delta, N, k);
        end
    end
    

    
    Px = Px.';
    Qx = Qx.';
    m = [Px; Qx];
    n = [P.'; Q.'];
    
    delta_y = n- m;

end

function a = find_Pk(Ybus, V, delta, N, k)
    result = 0;
    for n = 1:N
       num = abs(Ybus(k,n)) * V(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
       result = result + num;
    end
    a = V(k) * result;
end

function a = find_Qk(Ybus, V, delta, N, k)
    result = 0;
    for n = 1:N
       num = abs(Ybus(k,n)) * V(n) * sin(delta(k) - delta(n) - angle(Ybus(k,n)));
       result = result + num;
    end
    a = V(k) * result;
end

function [newV, newDelta] = find_V_delta_mismatch(mismatch, J, N, V, delta, bus_type)
    newV = zeros(1,N);
    newDelta = zeros(1,N);
    bus_type_d = zeros(1, N);
    
    for n = 1:N
       if(bus_type(n) == 2)
          bus_type_d(n) = 0; 
       else
           bus_type_d(n) = bus_type(n);
       end
    end
    
    bus_type_v = bus_type;
    
    count = 1;
    V_delta_mismatch = (J)^(-1) * mismatch;
    
    for j = 1:N
        if(bus_type_d(j) == 1)
            newDelta(j) = 0;
        elseif(bus_type_d(j) == 2)
            newDelta(j) = 1;
        else
            newDelta(j) = V_delta_mismatch(count) + delta(j);
            count = count + 1;
        end
    end        
    
    for i = 1:N
       if(bus_type_v(i) == 1)
          newV(i) = 1;
       elseif(bus_type_v(i) == 2)
          newV(i) = 1;
       else
          newV(i) = V_delta_mismatch(count) + V(i);
          count = count + 1;
       end
    end                 
    
end 

function met = find_convergance(power, criteria)
    met = true;
    if(max(abs(power)) > criteria)
    met = false;
    end
end

function [newV, newDelta, newcounter, met, delta_y_org] = iteration(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, print)
        
        %Compute delta(y(i)) = [delta(P(i)) delta(Q(i))]       
        delta_y = find_delta_y(Ybus, V, delta, N, P, Q, bus_type);
        delta_y_org = delta_y;
        
        %Find the Jacobian
        [J, delta_y] = find_jacobian(V, Ybus, delta, N, bus_type, 1, delta_y);

        %Find delta and voltage mismatch
        [newV, newDelta] = find_V_delta_mismatch(delta_y, J, N, V, delta, bus_type);
        
        %Check for convergance
        met = find_convergance(delta_y, criteria);
        
        %Print Results of Iteration
        if(print == 1)
            print_results(J, delta_y_org, met, newV, newDelta, counter);
        end
        
        %Set Count to Next Iteration
        newcounter = counter+1;
        
        delta_y = delta_y_org;
end

function print_results(J, delta_y, met, newV, newDelta, counter)
    
    delta_y = delta_y.';

    disp("For Iteration " + counter)
    disp('')
    
    disp("Jacobian: ")
    J_table = array2table(J);
    disp(J_table)
   
    disp("Delta Y Power Mismatches")
    delta_y_table = array2table(delta_y);
    disp(delta_y_table)
    
%     disp("The Real Power Mismatch is: ")
%     realP = delta_y(1:N);
%     real_p_table = array2table(realP);
%     disp(real_p_table)
%     disp("The Reactive Power Mismatch is: ")
%     reactiveQ = delta_y(N+1:N*2);
%     reactive_p_table = array2table(reactiveQ);
%     disp(reactive_p_table)
    disp('Voltage for Next Iteration:')
    newV_table = array2table(newV);
    disp(newV_table)
    disp('Angles for Next Iteration:')
    newDelta_table = array2table(newDelta);
    disp(newDelta_table)
    
    disp(' ')
    if(met == true)
        disp(["Convergance Found On Iteration " + num2str(counter)])
    else
       disp("Convergance Not Found on This Iteration.") 
    end
    
end

function Q = add_cap(Q_start, var)
    Q = Q_start + var;
end

function newton_rasp_power_flow(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, Y_lines, line_bus_connections, numLines, cap_set, numIterations, met, testing_var_control, print)
    
    %Define Data
    Q_limit = -1.75;
    bus_add_cap = 7;
    raise_var_bus = 6;
    
    %Add capactior if needed
    if(cap_set == 1)
        disp("Turning on Capacitor")
        var = -.1;
        Q(bus_add_cap) = add_cap(Q(bus_add_cap), var);
        cap_set = 0;
    end
    
    if(testing_var_control == 1)
        print = 0;
        while (cap_set ~= 1)
        %Raise the vars on Bus 6 until the var control is switched on for
        %the generator on PV bus 4
            %Raise VARs
            disp(Q)
            Q(raise_var_bus) = Q(raise_var_bus) + (-.1);
            disp('Q')
            disp(Q)
            
            while (met ~= true && counter < numIterations)
                [V, delta, counter, met] = iteration(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, print);
            end
            
            %Check to see if PV Bus is over limit
            for y = 1:N
                if(bus_type(y) == 2)
                    if(abs(Q(y)) > abs(Q_limit))
                        bus_type(y) = 0;
                        setQ = y;
                        cap_set = 1;
                        %Set the Q value of the bus
                        Q(setQ) = 1.75;
                    end
                end
            end
        end
    else
        %Run Power Flow Until Convergance
        while (met ~= true && counter < numIterations)
            [V, delta, counter, met, delta_y] = iteration(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, print);
        end
        if(met == true)
           find_current_system_info(Y_lines, V, N, line_bus_connections, numLines) 
           [real_loss, reactive_loss] = find_losses(delta_y, N)
        end
    end
   
    %Return to compare with cap bank
    if (cap_set == 1)
        testing_var_control = 0;
        disp("Running Newton Raspson with Capacitor Bank Incorporated")
        print = 1;
        newton_rasp_power_flow(Ybus, V, delta, P, Q, N, criteria, counter, bus_type, Y_lines, line_bus_connections, numLines, cap_set, numIterations, met, testing_var_control, print);
    end 
   
end

function currents = find_current_system_info(Y_lines, V, N, line_bus_connections, numLines)
    currents = find_current_directions(Y_lines, V, N, line_bus_connections, numLines);
end

function [real_loss, reactive_loss] = find_losses(delta_y, N)

    real_loss = 0;
    reactive_loss = 0;
    
    P_x = delta_y(1:N,1);
    Q_x = delta_y(N+1:N*2,1);
        
    for i = 1:N
       real_loss = real_loss + P_x(i); 
       reactive_loss = reactive_loss + Q_x(i);
    end

end

function currents = find_current_directions(Ylines, data, N, line_bus_connections, numLines)

    currents = zeros(1,N);

    for i = 1:numLines
        count = 0;
        while count < 1
            for k = 1:N
                if(line_bus_connections(k,i)==1)
                    bus1 = k;
                    count = count + 1;
                end
            end
        end
        
        count = 0;
        while count < 1
            for h = 1:N
                if(line_bus_connections(h,i)==1 && h ~= bus1)
                    bus2 = h;
                    count = count + 1;
                end
            end
        end
        
        currents(i) = find_current(data(bus1), data(bus2), Ylines(i)); %1-2
        %disp(currents(i))
    end
end 

function I = find_current(V1, V2, Y) 
    I = (V1-V2)*Y;
end 

%Milestone 4---------------------------------------------------------
function DC_power_flow()

end

function fast_decoupled_power_flow()

end

%Milestone 5---------------------------------------------------------
function [Zbus_pos, Zbus_neg, Zbus_zero] = find_Zbus_for_sequences(numLines, N, line_bus_connections, bus_connections, Y_lines, bus_types, TX_connections, numTX)
    
    numGen = 0;
    gen_buses = zeros(1, numGen);
    
    for i = 1:N
       if(bus_types(i) == 1 || bus_types(i) == 2)
          numGen = numGen + 1; 
          gen_buses(i) = i;
       end
    end
    numLines = numLines + numGen;
    
    bus_connections_pos = find_bus_connections_for_pos_seq(bus_connections, gen_buses);
    line_bus_connections_pos = find_line_bus_connections_pos_seq(line_bus_connections);
    Y_lines_pos = find_Y_lines_for_pos_seq(Y_lines, numLines);

    diagMatrix = find_diag_elements(line_bus_connections_pos, Y_lines_pos, N, numLines);
    off_diag = find_off_diag_elements(line_bus_connections_pos, Y_lines_pos, numLines, N);
    Ybus = find_Ybus_matrix(diagMatrix, off_diag, bus_connections_pos, N);
%     Ybus_table = array2table(Ybus);
%     disp('Positive Sequence Network Ybus')
%     disp(Ybus_table)
    Zbus_pos = inv(Ybus);
    Zbus_table = array2table(Zbus);
    disp('Positive Sequence Network Zbus')
    disp(Zbus_table)

    bus_connections_zero = find_bus_connections_for_zero_seq(bus_connections);
    line_bus_connections_zero = find_line_bus_connections_zero_seq(line_bus_connections);
    Y_lines_zero = find_Y_lines_for_zero_seq(Y_lines, numLines);

    diagMatrix = find_diag_elements(line_bus_connections_zero, Y_lines_zero, N, numLines);
    off_diag = find_off_diag_elements(line_bus_connections_zero, Y_lines_zero, numLines, N);
    Ybus = find_Ybus_matrix(diagMatrix, off_diag, bus_connections_zero, N);
%     Ybus_table = array2table(Ybus);
%     disp('Zero Sequence Network Ybus')
%     disp(Ybus_table)
    Zbus_zero = inv(Ybus);
    Zbus_table = array2table(Zbus);
    disp('Zero Sequence Network Zbus')
    disp(Zbus_table)

    bus_connections_neg = find_bus_connections_for_neg_seq(bus_connections);
    line_bus_connections_neg = find_line_bus_connections_neg_seq(line_bus_connections);
    Y_lines_neg = find_Y_lines_for_neg_seq(Y_lines, numLines);

    diagMatrix = find_diag_elements(line_bus_connections_neg, Y_lines_neg, N, numLines);
    off_diag = find_off_diag_elements(line_bus_connections_neg, Y_lines_neg, numLines, N);
    Ybus = find_Ybus_matrix(diagMatrix, off_diag, bus_connections_neg, N);
%     Ybus_table = array2table(Ybus);
%     disp('Negative Sequence Network Ybus')
%     disp(Ybus_table)
    Zbus_neg = inv(Ybus);
    Zbus_table = array2table(Zbus);
    disp('Negative Sequence Network Zbus')
    disp(Zbus_table)
end

function bus_connections_pos = find_bus_connections_for_pos_seq(bus_connections)
    bus_connections_pos = bus_connections;
end
function line_bus_connections_pos = find_line_bus_connections_pos_seq(line_bus_connections, numGen, gen_buses, N)
   
    m = zeros(numGen, N);
    count = 1;
    for i = 1:N
       if(gen_buses(i) ~= 0)
           m(count, gen_buses(i)) = 1;
           count = count + 1;
       end 
    end
    
    line_bus_connections_pos = [line_bus_connections; m];
end

function Y_lines_pos = find_Y_lines_for_pos_seq(Y_lines, numLines)
    Y_lines_pos = Y_lines;
end

function bus_connections_neg = find_bus_connections_for_neg_seq(bus_connections)
    bus_connections_neg = bus_connections;
end
function line_bus_connections_neg = find_line_bus_connections_neg_seq(line_bus_connections, numGen, gen_buses, N)
    m = zeros(numGen, N);
    count = 1;
    for i = 1:N
       if(gen_buses(i) ~= 0)
           m(count, gen_buses(i)) = 1;
           count = count + 1;
       end 
    end
    
    line_bus_connections_neg = [line_bus_connections; m];
end

function Y_lines_neg = find_Y_lines_for_neg_seq(Y_lines, numLines)
    Y_lines_neg = Y_lines;
end

function bus_connections_zero = find_bus_connections_for_zero_seq(bus_connections, TX_connections, numTX)
    
    for k = 1:numTX
        bus1 = TX_connections(1, k);
        bus2 = TX_connections(2, k);
        bus_connections(bus1, bus2) = 0;
        bus_connections(bus2, bus1) = 0;
    end
    
    bus_connections_zero = bus_connections;
end
function line_bus_connections_zero = find_line_bus_connections_zero_seq(line_bus_connections)

    m = zeros(numGen, N);
    count = 1;
    for i = 1:N
       if(gen_buses(i) ~= 0)
           m(count, gen_buses(i)) = 1;
           count = count + 1;
       end 
    end
    
    line_bus_connections_zero = [line_bus_connections; m];
end
function Y_lines_zero = find_Y_lines_for_zero_seq(Y_lines, numLines)
    Y_lines_zero = Y_lines; 
end

%Milestone 6----------------------------------------------------------
function [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = single_line_to_gnd(V_f, Z_f, Z0, Z1, Z2, bus)
    
    
    seq_fault_currents = zeros(1,3);
    
    seq_fault_currents(1) = V_f / (Z0 + Z1 + Z2 + 3*Z_f);
    seq_fault_currents(2) = seq_fault_currents(1);
    seq_fault_currents(3) = seq_fault_currents(2);
    
    seq_fault_voltages = [0; V_f; 0] - [Z0 0 0; 0 Z1 0; 0 0 Z2]*[seq_fault_currents(1); seq_fault_currents(2); seq_fault_currents(3)];
    seq_fault_voltages = seq_fault_voltages.';
%     seq_fault_voltages(1) = Z1 * seq_fault_currents(1);
%     seq_fault_voltages(2) = -1 * Z2 * seq_fault_currents(2);
%     seq_fault_voltages(3) = -1 * Z0 * seq_fault_currents(3);

    ph_fault_currents = zeros(1,3);
    ph_fault_voltages = zeros(1,3);
    
    ph_fault_currents(1) = 3 * seq_fault_currents(1);
    ph_fault_currents(2) = 0;
    ph_fault_currents(3) = 0;
    
    ph_fault_voltages(1) = Z_f * ph_fault_currents(1);
    ph_fault_voltages(2) = Z_f * ph_fault_currents(2);
    ph_fault_voltages(3) = Z_f * ph_fault_currents(3);
    
    disp(['Single Line to Ground Fault at Bus ' num2str(bus)])
    disp(' ')
    %[seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = single_line_to_gnd(V_f, Z_f, Z0, Z1, Z2);
    disp('Sequence Fault Current: ')
    disp(['I0: ' num2str(seq_fault_currents(1)) 'A'])
    disp(['I1: ' num2str(seq_fault_currents(2)) 'A'])
    disp(['I2: ' num2str(seq_fault_currents(3)) 'A'])
    disp(' ')
    disp('Sequence Fault Voltages:')
    disp(['V0: ' num2str(seq_fault_voltages(1)) 'V'])
    disp(['V1 ' num2str(seq_fault_voltages(2)) 'V'])
    disp(['V2 ' num2str(seq_fault_voltages(3)) 'V'])
    disp('Phase Fault Current: ')
    disp(['Ia: ' num2str(ph_fault_currents(1)) 'A'])
    disp(['Ib: ' num2str(ph_fault_currents(2)) 'A'])
    disp(['Ic: ' num2str(ph_fault_currents(3)) 'A'])
    disp(' ')
    disp('Phase Fault Voltages:')
    disp(['Va: ' num2str(ph_fault_voltages(1)) 'V'])
    disp(['Vb: ' num2str(ph_fault_voltages(2)) 'V'])
    disp(['Vc: ' num2str(ph_fault_voltages(3)) 'V'])
        
end

function [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = line_to_line_fault(V_f, Z_f, Z0, Z1, Z2)
    
    seq_fault_currents = zeros(1,3);
    seq_fault_voltages = zeros(1,3);
    
    seq_fault_currents(1) = V_f / (Z1 + Z2 + Z_f);
    seq_fault_currents(2) = -1 * seq_fault_currents(1);
    seq_fault_currents(3) = 0;
       
    seq_fault_voltages = Z_f * seq_fault_currents(1);
    
    ph_fault_currents = zeros(1,3);
    ph_fault_voltages = zeros(1,3);
    
    a = -0.5+i*0.8660254038;
    ph_fault_currents = zeros(1,3);
    ph_fault_currents(1) = 0;
    ph_fault_currents(2) = (a^2 - a) * seq_fault_currents(1);
    ph_fault_currents(3) = -1 * ph_fault_currents(2);
    
    ph_fault_voltages = Z_f * ph_fault_currents(2);
    
    disp(' ')
    disp(['Line to Line Fault at Bus ' num2str(bus)])
    %[seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = line_to_line_fault(V_f, Z_f, Z0, Z1, Z2);
    disp('Sequence Fault Current: ')
    disp(['I0: ' num2str(seq_fault_currents(1)) 'A'])
    disp(['I1: ' num2str(seq_fault_currents(2)) 'A'])
    disp(['I2: ' num2str(seq_fault_currents(3)) 'A'])
    disp(' ')
    disp('Sequence Fault Voltage (V1 - V2):')
    disp(['V: ' num2str(seq_fault_voltages) 'V'])
    disp('Phase Fault Current: ')
    disp(['Ia: ' num2str(ph_fault_currents(1)) 'A'])
    disp(['Ib: ' num2str(ph_fault_currents(2)) 'A'])
    disp(['Ic: ' num2str(ph_fault_currents(3)) 'A'])
    disp(' ')
    disp('Phase Fault Voltage (Vb - Vc):')
    disp(['V: ' num2str(ph_fault_voltages) 'V'])
    
end

function [sub_tran_fault_current, E] = symmetrical_fault(V_f, Zbus, bus, N)
    
    sub_tran_fault_current = V_f / Zbus(bus, bus);
    
    E = zeros(1,N); 
    
    for i = 1:N
        E(i) = (1-(Zbus(i,bus)/Zbus(bus,bus))) * V_f;
    end 
    
    disp(' ')
    disp(['Symmetrical Fault at Bus ' num2str(bus)])
    %[seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = line_to_line_fault(V_f, Z_f, Z0, Z1, Z2);
    disp('Sequence Fault Current: ')
    disp(['I0: ' num2str(seq_fault_currents(1)) 'A'])
    disp(['I1: ' num2str(seq_fault_currents(2)) 'A'])
    disp(['I2: ' num2str(seq_fault_currents(3)) 'A'])
    disp(' ')
    disp('Sequence Fault Voltage (V1 - V2):')
    disp(['V: ' num2str(seq_fault_voltages) 'V'])
    disp('Phase Fault Current: ')
    disp(['Ia: ' num2str(ph_fault_currents(1)) 'A'])
    disp(['Ib: ' num2str(ph_fault_currents(2)) 'A'])
    disp(['Ic: ' num2str(ph_fault_currents(3)) 'A'])
    disp(' ')
    disp('Phase Fault Voltage (Vb - Vc):')
    disp(['V: ' num2str(ph_fault_voltages) 'V'])
    
end

function [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = double_line_to_gnd(V_f, Z_f, Z0, Z1, Z2, bus)
    
    seq_fault_currents = zeros(1,3);
    seq_fault_voltages = zeros(1,3);
    
    seq_fault_currents(1) = V_f / (Z1 + Z2 + Z_f);
    seq_fault_currents(2) = -1 * seq_fault_currents(1);
    seq_fault_currents(3) = 0;
       
    seq_fault_voltages = Z_f * seq_fault_currents(1);
    
    ph_fault_currents = zeros(1,3);
    ph_fault_voltages = zeros(1,3);
    
    a = -0.5+i*0.8660254038;
    ph_fault_currents = zeros(1,3);
    ph_fault_currents(1) = 0;
    ph_fault_currents(2) = (a^2 - a) * seq_fault_currents(1);
    ph_fault_currents(3) = -1 * ph_fault_currents(2);
    
    ph_fault_voltages = Z_f * ph_fault_currents(2);
    
    disp(' ')
    disp(['Double Line to Ground Fault at Bus ' num2str(bus)])
    %[seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = line_to_line_fault(V_f, Z_f, Z0, Z1, Z2);
    disp('Sequence Fault Current: ')
    disp(['I0: ' num2str(seq_fault_currents(1)) 'A'])
    disp(['I1: ' num2str(seq_fault_currents(2)) 'A'])
    disp(['I2: ' num2str(seq_fault_currents(3)) 'A'])
    disp(' ')
    disp('Sequence Fault Voltage (V1 - V2):')
    disp(['V: ' num2str(seq_fault_voltages) 'V'])
    disp('Phase Fault Current: ')
    disp(['Ia: ' num2str(ph_fault_currents(1)) 'A'])
    disp(['Ib: ' num2str(ph_fault_currents(2)) 'A'])
    disp(['Ic: ' num2str(ph_fault_currents(3)) 'A'])
    disp(' ')
    disp('Phase Fault Voltage (Vb - Vc):')
    disp(['V: ' num2str(ph_fault_voltages) 'V'])
    
end

function perform_fault_analysis(V_f, Z_f, Zbus0, Zbus1, Zbus2)
    
    %Select the Bus to Fault At
    bus = 1;
    
    %Find the Equivalent Impendace for Each Sequence
    Z0 = Zbus0(bus, bus);
    Z1 = Zbus1(bus, bus);
    Z2 = Zbus2(bus, bus);
    
    [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = line_to_line_fault(V_f, Z_f, Z0, Z1, Z2, bus);
    [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = single_line_to_gnd(V_f, Z_f, Z0, Z1, Z2, bus);
    [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = symmetrical_fault(V_f, Z_f, Z0, Z1, Z2, bus);
    [seq_fault_currents, seq_fault_voltages, ph_fault_currents, ph_fault_voltages] = double_line_to_gnd(V_f, Z_f, Z0, Z1, Z2, bus);
   
end

%Helper Functions----------------------------------------------------
function num = convert_to_rect(mag, ang)
    [real, imag] = pol2cart(ang * (pi/180), mag);
    num = real + i*imag;
end

%Return Voltage Variable Functions
function Vp = get_Vp(Va, Vb, Vc)
    Vp = [Va; Vb; Vc];
end

function Vs = get_Vs(V0, V1, V2)
    Vs = [V0; V1; V2];
end

function [Va, Vb, Vc] = get_phase_voltage_components(Vp)
    Va = Vp(1);
    Vb = Vp(2);
    Vc = Vp(3);
end

function [V0, V1, V2] = get_sequence_voltage_components(Vs)
    V0 = Vs(1);
    V1 = Vs(2);
    V2 = Vs(3);
end

%Return Voltage Variable Functions
function Ip = get_Ip(Ia, Ib, Ic)
    Ip = [Ia; Ib; Ic];
end

function Is = get_Is(I0, I1, I2)
    Is = [I0; I1; I2];
end

function [Ia, Ib, Ic] = get_phase_current_components(Ip)
    Ia = Ip(1);
    Ib = Ip(2);
    Ic = Ip(3);
end

function [I0, I1, I2] = get_sequence_current_components(Is)
    I0 = Is(1);
    I1 = Is(2);
    I2 = Is(3);
end

%Transformation Functions
function Vp = find_phase_voltages_from_sequence(A, Vs)
    Vp = A * Vs;
end 

function Vs = find_sequence_voltages_from_phase(Vp,a)
    M = [1 1 1; 1 a a^2; 1 a^2 a];
    Vs = (1/3) * M * Vp;
end 

function Ip = find_phase_currents_from_sequence(A, Is)
    Ip = A * Is;
end 

function Is = find_sequence_currents_from_phase(Ip,a)
    M = [1 1 1; 1 a a^2; 1 a^2 a];
    Is = (1/3) * M * Ip;
end 
