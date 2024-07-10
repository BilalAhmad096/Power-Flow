clear all;
clc;
                %% Power Flow analysis using Gauss Seidel%%
%Variables:
before = memory;
tic
% Your MATLAB code here


% %% Voltage values:
[V1] = randomvalues(1.00,1.04);
[V2]=randomvalues(0.99,1.03);
[V3]=randomvalues(0.98,1.02);
[V4]=randomvalues(1.00,1.04);
[V5]=randomvalues(0.96,1.00);
[V6]=randomvalues(1.01,1.05);
[V7]=randomvalues(0.99,1.03);
[V8]=randomvalues(0.97,1.01);
[V9]=randomvalues(0.98,1.02);

%% Active Power values:
[P2]=randomvalues(0.98,1.02);
[P3]=randomvalues(1.48,1.52);
[P4]=randomvalues(0.71,0.75);
[P5]=randomvalues(0.98,1.02);
[P6]=randomvalues(0.98,1.02);
[P7]=randomvalues(0.98,1.02);
[P8]=randomvalues(0.98,1.02);
[P9]=randomvalues(0.98,1.02);


%% Reactive power values:
% [Q3]=randomvalues(0.28,0.32);
% [Q8]=randomvalues(0.28,0.32);

%% results
Result=[];

%% no of iterations:

   % The program section to time. 

z=2;

for loop=1:z
    % Information about the bus matrix
      % no   V         Ang.    Pg         Qg         PL         QL    Type
      % (1) (2)        (3)     (4)        (5)        (6)        (7) 
    bus=[1  V1(loop)  0.000   0.00       0.00       0.00       0.00      1;
         2  V2(loop)  0.000   P2(loop)   0.00       0.00       0.00      2;
         3  V3(loop)  0.000   0.00       0.00       P3(loop)   Q3(loop)  3;
         4  V4(loop)  0.000   P4(loop)   0.00       0.00       0.00      2;
         5  V5(loop)  0.000   P5(loop)   0.00       0.00       0.00      2;            
         6  V6(loop)  0.000   P6(loop)   0.00       0.00       0.00      2;      
         7  V7(loop)  0.000   P7(loop)   0.00       0.00       0.00      2;
         8  V8(loop)  0.000   0.00       0.00       P8(loop)   Q8(loop)  3;
         9  V9(loop)  0.000   P9(loop)   0.00       0.00       0.00      2];
 

%        from   To       R          X
%        (1)   (2)      (3)        (4)
   line =[1     4	     0	      0.0576;
          4	    5	     0        0.0920;
          5     6	     0        0.1700;
          3     6	     0	      0.0586;
          6     7	     0        0.1008;
          7     8	     0        0.0720;
          8     2        0	      0.0625;
          8     9	     0	      0.1610; 
          9     4	     0	      0.0850];

% ++++++++++++++++++++++ Bus addmitance matrix Ybus +++++++++++++++++++++++
    [Y] = ybusmatrix(line);

% +++++++++++++++++++++++ Start iterative process +++++++++++++++++++++++++
    nbuses=length(bus(:,1)); % number of buses of the electric power system
    V=bus(:,2); Vprev=V; % Initial bus voltages
    Theta=bus(:,3); % Initial bus angles
% Net power (Generation - Load)
    P=bus(:,4)-bus(:,6);     
    Q=bus(:,5)-bus(:,7);
    tolerance=1; 
    iteration=0;
    while (tolerance > 1e-8)
        for k=2:nbuses
            PYV=0;
            for i=1:nbuses
                if k ~= i
                    PYV = PYV + Y(k,i)* V(i);  % Vk * Yik
                end
            end
            if bus(k,8)==2 % PV bus
                % Estimate Qi at each iteration for the PV buses
                Q(k)=-imag(conj(V(k))*(PYV + Y(k,k)*V(k)));
            end
            V(k) = (1/Y(k,k))*((P(k)-j*Q(k))/conj(V(k))-PYV); % Compute bus voltages
            if bus(k,8) == 2 % For PV buses, the voltage magnitude remains same, but the angle changes
                V(k)=abs(Vprev(k))*(cos(angle(V(k)))+j*sin(angle(V(k))));
            end
        end
        iteration=iteration+1;
        tolerance = max(abs(abs(V) - abs(Vprev)));
        Vprev=V;
    end
% +++++++++++++++++++++++++++++ Power flow ++++++++++++++++++++++++++++++++
% currents at each node
    I=Y*V;
% Power at each node
    S=V.*conj(I); % Complex power
    for k=1:nbuses
        if bus(k,8)==1
        % Real and reactive generation at the Slack bus
            Pgen(k)=real(S(k));
            Qgen(k)=imag(S(k));
        end
        if bus(k,8)==2
        % Real and reactive generation at the PV buses
            Pgen(k)=real(S(k))+bus(k,6);
            Qgen(k)=imag(S(k))+bus(k,7);
        end
        if bus(k,8)==3
            Pgen(k)=0;
            Qgen(k)=0;
        end
    end
% calculate the line flows and power losses
    FromNode=line(:,1);
    ToNode=line(:,2);
    nbranch = length(line(:,1)); % number of branches
% Define admmitance of lines
    r = line(:,3);
    rx = line(:,4);
    z = r + j*rx;
    y = ones(nbranch,1)./z;
% Define complex power flows
    Ss = V(FromNode).*conj((V(FromNode) - V(ToNode)).*y); % complex flow of the sending buses
    Sr = V(ToNode).*conj((V(ToNode) - V(FromNode)).*y); % complex low of the receiving buses

% Define active and reactive power flows
    Pij=real(Ss);
    Qij=imag(Ss);
    Pji=real(Sr);
    Qji=imag(Sr);

% Active power lossess
    P_loss=sum(Pij+Pji);

% Reactive power lossess
    Q_loss=sum(Qij+Qji);

% Extract values to CSV files:
    P1out=Pgen(1);
    P2out=Pgen(2);
    P3out=Pgen(3);
    P4out=Pgen(4);
    P5out=Pgen(5);
    P6out=Pgen(6);
    P7out=Pgen(7);
    P8out=Pgen(8);
    P9out=Pgen(9);


    Q1out=Qgen(1);
    Q2out=Qgen(2);
    Q3out=Qgen(3);
    Q4out=Qgen(4);
    Q5out=Qgen(5);
    Q6out=Qgen(6);
    Q7out=Qgen(7);
    Q8out=Qgen(8);
    Q9out=Qgen(9);


    Voltage=abs(V);
    V1out=Voltage(1);
    V2out=Voltage(2);
    V3out=Voltage(3);
    V4out=Voltage(4);
    V5out=Voltage(5);
    V6out=Voltage(6);
    V7out=Voltage(7);
    V8out=Voltage(8);
    V9out=Voltage(9);


    Delta2=(180/pi)*angle(V(2));
    Delta3=(180/pi)*angle(V(3));
    Delta4=(180/pi)*angle(V(4));
    Delta5=(180/pi)*angle(V(5));
    Delta6=(180/pi)*angle(V(6));
    Delta7=(180/pi)*angle(V(7));
    Delta8=(180/pi)*angle(V(8));
    Delta9=(180/pi)*angle(V(9));


    Result=[Result;P1out P2out -P3(loop) P4out P5out P6out P7out -P8(loop) P9out Q1out Q2out -Q3(loop) Q4out Q5out Q6out Q7out -Q8(loop) Q9out Delta2 Delta3 Delta4 Delta5 Delta6 Delta7 Delta8 Delta9  V1out   V2out   V3out V4out V5out V6out V7out V8out V9out];

    
end
timeElapsed = toc
after = memory;
memoryUsed = after.MemUsedMATLAB - before.MemUsedMATLAB;
disp(['Memory used by your code: ' num2str(memoryUsed) ' bytes']);

% writematrix(Result);
% type 'Result.txt';
%% 
% +++++++++++++++++++++++++++ Print results +++++++++++++++++++++++++++++++
disp('                Gauss Seidel Load-Flow Analysis')
disp('               Report of Power Flow Calculations ')
disp(' ')
fprintf('Number of iterations        : %g \n', iteration);
fprintf('Total real power losses     : %g.\n',P_loss);
fprintf('Total reactive power losses : %g.\n\n',Q_loss);
disp('                                      Generation             Load')
disp('     Bus      Volts     Angle      Real  Reactive      Real  Reactive ')
ywz=[bus(:,1)    abs(V)  (180/pi)*angle(V)  Pgen'  Qgen'  bus(:,6)  bus(:,7)];
disp(ywz)

disp('                      Line Flows                     ')
disp('    Line     From Bus   To Bus     Real    Reactive   ')
l=1:1:length(line(:,1));
xy=[l' FromNode ToNode Pij Qij];
disp(xy)
abs(V);
