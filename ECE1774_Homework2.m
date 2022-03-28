clear
clc
% ECE1774 Homework 2

%Problem 1

%Positive Sequenece Impedance

%Find Ybus and Invert to get Zbus
%Define Given Impedances
% impedance.generator = ;
% impedance.transformer = .1;
% impedance.line1 = 0.03;
% impedance.line2 = 0.03;
% impedance.line3 = 0.03;
% impedance.line4 = 0.03;
% impedance.load_bus3 = 0.05;
% impedance.load_bus4 = 0.05;
% impedance.load_bus5 = 0.05;

N = 5;
numLines = 9;

% line_bus_connections = [0 1 1 0 0;
%                         0 1 0 1 0;
%                         0 0 1 0 1;
%                         0 0 0 1 1;
%                         1 1 0 0 0;
%                         1 0 0 0 0
%                         0 0 1 0 0;
%                         0 0 0 0 1;
%                         0 0 0 1 0];
                    
line_bus_connections = [0 0 0 0 1 1 0 0 0;
                        1 1 0 0 1 0 0 0 0;
                        1 0 1 0 0 0 1 0 0;
                        0 1 0 1 0 0 0 0 1;
                        0 0 1 1 0 0 0 1 0];
                    
bus_connections = [1 1 0 0 0;
                   1 1 1 1 0;
                   0 1 1 0 1;
                   0 1 0 1 1;
                   0 0 1 1 1];
                    
Y_lines1 = [1/0.03 1/0.03 1/0.03 1/0.03 1/0.1 1/0.1 1/0.05 1/0.05 1/0.05];
Y_lines2 = [1/0.03 1/0.03 1/0.03 1/0.03 1/0.1 1/0.13 1/0.05 1/0.05 1/0.05];
Y_lines0 = [1/0.09 1/0.09 1/0.09 1/0.03 0 1/0.08 1/0.2 1/0.2 1/0.2];
                        
%Positive Sequence

%Find the Diagonal Elements of the Ybus
diagMatrix1 = find_diag_elements(line_bus_connections, Y_lines1, N, numLines);

%Find Outside Diaglonal Elements of Ybus
off_diag1 = find_off_diag_elements(bus_connections, Y_lines1, numLines, N);

%Find the Ybus
Ybus = find_Ybus_matrix(diagMatrix1, off_diag1, bus_connections, N);

Zbus1 = zeros(N,N);

for i = 1:N
   for j = 1:N
       Zbus1(i,j) = 1/(Ybus1(i,j));
   end
end

%Negative Sequence Impedance

%Zero Sequence Impedance


%%%%%%%%%%Functions%%%%%%%%%%%%%
function Zbus = find_Zbus(Ybus)
    Zbus = inv(Ybus);
end

function diagMatrix = find_diag_elements(connection, Y_lines, N, numLines)

    j = 1i;

    matrix = zeros(N, 1);
    sum = 0;
    for i = 1:N
        bus_info = connection(i,:);
        element = (bus_info .* Y_lines);
        for k = 1:numLines
            sum = sum + element(k);
        end
        matrix(i) = sum;
        sum = 0;
    end
    
    diagMatrix = matrix;

end

function off_diag = find_off_diag_elements(line_bus_connections, Y_lines, numLines, N)
    
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

