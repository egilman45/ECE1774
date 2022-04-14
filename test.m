

clear 
clc

%------DEFINITIONS AND DATA COLLECTION--------%

w = 377;
j = 1i;

%Find the Number of Nodes and Lines in the System
[N, numLines] = find_lines_nodes();

%Find Bus Connections
bus_connections = find_bus_connections(N);

%Have User Define the Connections
line_bus_connections = find_line_connections(N, numLines);

%line_bus_connections = [1 1 0 0 0; 1 0 1 1 0; 0 1 1 0 1; 0 0 0 1 1];

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

%Bus Power Information
P = [0 0 0 0 1 1.1 1]; 
Q = [0 0 0 0 -.65 -.50 -.70]; 

%Flat Start Information
V_o = [1.0 0.93198 0.93648 1.0 0.91795 0.91373 0.90778];
delta = [0 (pi/180)*-4.5 (pi/180)*-4 (pi/180)*2.11 (pi/180)*-5.19 (pi/180)*-5.56 (pi/180)*-5.75];

N = 7;
flow_type = 1;

[J] = find_jacobian(V_o, Ybus, delta, N, bus_type, flow_type);

[newP, newQ, P_x, Q_x] = find_P_Q(Ybus, V_o, delta, N, P, Q);

delta_y = [P_x.'; Q_x.'];
disp(delta_y)

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

function [N, numLines] = find_lines_nodes()
    
%     prompt = {'How Many Nodes in the System? : ', 'How Many Lines in the System? : '};
%     dlgtitle = 'General System Data Collection';
%     dims = [1 35];
%     answer = inputdlg(prompt, dlgtitle,dims);
%     N = str2num(answer{1});
%     numLines = str2num(answer{2});
    N = 7;
    numLines = 7;
    
    
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

function conductor = find_conductor_info()
%     prompt = {'What is resistance (r)(ohm/mi)? : ', 'What is the Geometric Mean Radius (GMR)(ft)? : ', 'What is the diameter (d)(ft)? : ', 'What is distance between subconductors (d)(ft)? : ', 'How many subconductors in the bundle? : '};
%     dlgtitle = 'Conductor Type Data';
%     dims = [1 100];
%     definput = {'', ''};
%     answer = inputdlg(prompt, dlgtitle,dims);
%     conductor.r = str2num(answer{1});
%     conductor.GMR = str2num(answer{2});
%     conductor.diam = str2num(answer{3});
%     conductor.d = str2num(answer{4});
%     conductor.nb = str2num(answer{5});
    
    conductor.r = 0.385; %ohm/mi
    conductor.GMR = 0.0217; %ft
    conductor.diam = 0.642/12; %ft 
    conductor.d = 1.5; %ft
    conductor.nb = 2; 
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

function [J] = find_jacobian(V_o, Ybus, delta, N, bus_type, flow_type) 

    if(flow_type == 1)
        J1 = find_J1(V_o, Ybus, delta, N);
        J2 = find_J2(V_o, Ybus, delta, N);
        J3 = find_J3(V_o, Ybus, delta, N);
        J4 = find_J4(V_o, Ybus, delta, N);

        J = [J1 J2; J3 J4];
        
        J_t = array2table(J);
        disp(J_t)

        %Find slack bus 
        for i = 1:N
           if (bus_type(i) == 1)
              %Delete row and column for slack bus
               J(i, :) = [];

               J(:, i) = [];
               J((i+N-1), :) = [];

               J(:, (i+N-1)) = [];
           end
        end

        count = 1;

        %Find PV and delete the row and column for each
        for q = 1:N
           if(bus_type(q) == 2)
               J((q+N-2), :) = [];

               J(:, (q+N-2)) = [];
           end
           count = count + 1;
        end
    end
    

end

function [newP, newQ, P_x, Q_x] = find_P_Q(Ybus, V_o, delta, N, P, Q)
    
    newP = zeros(N,1);
    newQ = zeros(N,1);

       
    for k = 1:N
        result = 0; 
        for n = 1:N
            num = abs(Ybus(k,n)) * V_o(n) * cos(delta(k) - delta(n) - angle(Ybus(k,n)));
            result = result + num; 
        end
            newP(k,1) = V_o(k) * result;
    end
   
    for q = 1:N
        result = 0;
        for w = 1:N
            num = abs(Ybus(q,w)) * V_o(w) * sin(delta(q) - delta(w) - angle(Ybus(q,w)));
            result = result + num;
        end
            newQ(q,1) = V_o(q) * result;
    end

    
    P_x = newP.';
    Q_x = newQ.';
    newP = P - P_x;
    newQ = Q - Q_x;
    
end

function [real_loss, reactive_loss] = find_losses(P_x, Q_x, N)

    real_loss = 0;
    reactive_loss = 0;

    for i = 1:N
       real_loss = real_loss + P_x(i); 
       reactive_loss = reactive_loss + Q_x(i);
    end

end