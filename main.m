function main(timeStep, endTime)
    starttime = tic;
    DOF = 2;
    lengthAcross = 354; %m
    %Data in files follows format:
    %x,y
    x = 0:lengthAcross/512:lengthAcross; % increments in approx 25 cm
    y = x.*0;
%     y = x.*(x-lengthAcross)./4876.8;
    nodesData = [transpose(x) transpose(y)];
    
    %Retrieve number of nodes
    numNodes = size(nodesData,1);
    nodeToEvaluate = (numNodes-1)/2;
    %node1,node2
    sctrData = [transpose(1:(numNodes-1)) transpose(2:numNodes)];

    %Initialize mechanical property constants
    %propsData = csvread(strcat('props', num2str(questionNum), '.txt'));
    E = 200E9; %Pa  
    area = 1.487E-3; %m^2
    density = 7394; % kg/m^3
    dampening = 37.8; % 37.8 Nm/s
    momentOfInertia = 3.2E-7; % m^4
    gamma = 3/2; % for implicit
    beta = 8/5; % for implicit

    
    %Set End Nodes as fixed, all others to free
    nodeTypes = getNodeTypes(numNodes);

    %Retrieve number of elements and element information
    numElements = size(sctrData,1);

    %calculate and store all data necessary to build global matrices and theta for each element
    c = zeros(numElements,1);
    m = zeros(numElements,1);
    k = zeros(numElements,1);
    theta = zeros(numElements,1);
 
    % set unchanging element constants
    for i = 1:numElements
        node1 = sctrData(i,1);
        node2 = sctrData(i,2);
        x1 = nodesData(node1,1);
        y1 = nodesData(node1,2);
        x2 = nodesData(node2,1);
        y2 = nodesData(node2,2);
%         theta(i) = atan2d(y2-y1,x2-x1);
        length = sqrt(((x2-x1)^2)+((y2-y1)^2));
        k(i) = (area*E/length);
        c(i) = dampening;
        m(i) = length*area*density;
    end

    % Generated Gravitational forces to apply
    appliedForces = zeros((numNodes-2),3);
    for i = 1:numNodes-2
        appliedForces(i,1) = i+1;
        appliedForces(i,3) = (-9.81)*m(i); %N
    end
       
    %Create displacements and forces matrix; second column holds wether or not values
    %are known
    %UCurr = calculateU(numNodes, nodeTypes);
    UCurr = csvread('new.csv');
    F = calculateF(numNodes, nodeTypes, appliedForces);
    
    %UCurr = UCurr(:,1);
    VCurr = zeros(size(UCurr));
    aCurr = VCurr;
    F = F(:,1);
    
    % Create [M]
    M = calculateM(numNodes, numElements, sctrData, m);
    
            %size64 ensures the value passed to zeros() is an int
        size64 = int64(endTime/timeStep + 1);
        keepTrackU = zeros(size64, DOF + 1);
%         keepTrackVel = zeros(size64, DOF + 1);
%         keepTrackAcc = zeros(size64, DOF + 1);
    nodeToEvaluateX = nodeToEvaluate*2-1;
    nodeToEvaluateY = nodeToEvaluate*2;
    disp(nodeToEvaluateX);
    disp(nodeToEvaluateY);
    
    i = 1;
    for (j = 0: timeStep: endTime)
        disp(j);
        % set changing element constants
        for i2 = 1:numElements
            node1 = sctrData(i2,1);
            node2 = sctrData(i2,2);
            x1 = nodesData(node1,1) + UCurr(node1*2-1);
            y1 = nodesData(node1,2) + UCurr(node1*2);
            x2 = nodesData(node2,1) + UCurr(node2*2-1);
            y2 = nodesData(node2,2) + UCurr(node2*2);
            theta(i2) = atan2(y2-y1,x2-x1);
            %F = recalculateForce(F, i2, j, timeStep);
%             length = sqrt(((x2-x1)^2)+((y2-y1)^2));
%             k(i) = (area*E/length);
%             c(i) = dampening;
%             m(i) = length*area*density;
        end

        %Create [C],[Keff]
        C = calculateK(numNodes, numElements, sctrData, c, theta);
        K = calculateK(numNodes, numElements, sctrData, k, theta);

        [rows, cols] = size(K);

%         [Forganize, Uorganize, Morganize, Corganize, Korganize, knowns, order] = organizeMatrix(F, U, M, C, K); 

        Atilde = 2*M/(beta*(timeStep^2)) + K + 2*C*gamma/(beta*timeStep);
        Btilde = 2*M/(beta*(timeStep^2)) + 2*C*gamma/(beta*timeStep);
        Ctilde = 2*M/(beta*timeStep) + C*(2*gamma/beta - 1);
        Dtilde = M*(1 - beta)/beta + timeStep*C*((gamma - 1) + gamma*(1 - beta)/beta);

        %make matrices to record each iteration 
%         keepTrackU(1, 2) = UCurr(nodeToEvaluateX);
%         keepTrackU(1, 3) = UCurr(nodeToEvaluateY);
    %     keepTrackVel(1, 2) = keepTrackU(1, 2)/timeStep;
    %     keepTrackVel(1, 3) = keepTrackU(1, 3)/timeStep; 
    %     keepTrackAcc(1, 2) = keepTrackVel(1,2)/timeStep; 
    %     keepTrackAcc(1, 3) = keepTrackVel(1,3)/timeStep;

%         if questionNum == 2
%             Forganize(order(nodeToEvaluateX)) = j;
%         elseif questionNum == 3
%             %applying a sinusoid displacement to node 5 in the y-direction
%             omega = 0.1;
%             applyNode = 5;
%             Uorganize(nodeToDisplaceY) = 50 * sin(omega*j);
%         end

        FAtNextTimeStep = F + Btilde*UCurr + Ctilde*VCurr + Dtilde*aCurr; 

        UNext = solveMatrix(Atilde, FAtNextTimeStep, UCurr);
        aNext = (2/(beta*timeStep))*((UNext - UCurr)/timeStep)-((2/(beta*timeStep))*VCurr)- (((1-beta)/beta)*aCurr);
        VNext = VCurr + timeStep*((1-gamma)*aCurr + gamma*aNext);
        
        keepTrackU(i + 1, 1) = j;
        keepTrackU(i + 1, 2) = UNext(nodeToEvaluateX);
        keepTrackU(i + 1, 3) = UNext(nodeToEvaluateY);
%         keepTrackVel(i + 1, 1) = j;
%         keepTrackVel(i + 1, 2) = VorganizeNext(nodeToEvaluateX);
%         keepTrackVel(i + 1, 3) = VorganizeNext(nodeToEvaluateY);
%         keepTrackAcc(i + 1, 1) = j;
%         keepTrackAcc(i + 1, 2) = aorganizeNext(nodeToEvaluateX);
%         keepTrackAcc(i + 1, 3) = aorganizeNext(nodeToEvaluateY);
        
        UCurr = UNext; 
        VCurr = VNext; 
        aCurr = aNext; 
        
        i = i + 1;
    end
    
    %Commented disp cuts off y displacement
    %disp(keepTrackU(:,1:2));
    %disp(keepTrackVel(:,1:2));
    %disp(keepTrackAcc(:,1:2));
%     disp(keepTrackU);
    plot(keepTrackU(:,1), keepTrackU(:,3));
%     disp(keepTrackVel);
%     disp(keepTrackAcc);

%     % Reorganize to original Node Order
%     Usolve = zeros(size(Uorganize,1),1);
%     for i = 1:size(Uorganize,1)
%         Usolve(order(i)) = UCurr(i);
%     end

    %     disp(UCurr);
    disp(nodesData(:,1));
    disp(nodesData(:,2));
    sctrData = [transpose(sctrData(:,1)) ; transpose(sctrData(:,2))];
    disp(sctrData);
    postprocesser(nodesData(:,1),nodesData(:,2),sctrData,UCurr);

    csvwrite('old.csv', UCurr);
    disp('IMPLICIT COMPUTATION TIME (s):');
    disp(toc(starttime));
end

function nodeTypes = getNodeTypes(numNodes)
    % returns an array contaiNevermind ning data for each node type for the system
    % 1 - free, 2 - fixed, 3 - roller in x, 4 - roller in y
    nodeTypes = zeros(numNodes,1) + 1;
    nodeTypes(1) = 2;
    nodeTypes(numNodes) = 2;
end

function M = calculateM(numNodes, numElements, sctrData, m)
    DOF = 2;
    M = zeros(numNodes*DOF,numNodes*DOF);
    
    for i = 1:numElements
        node1 = sctrData(i,1);
        node2 = sctrData(i,2);
        M(node1*DOF-1,node1*DOF-1) = m(i)/2;
        M(node1*DOF,node1*DOF) = m(i)/2;
        M(node2*DOF-1,node2*DOF-1) = m(i)/2;
        M(node2*DOF,node2*DOF) = m(i)/2;
    end
end

function K = calculateK(numNodes, numElements, sctrData, k, theta)
    DOF = 2;
    K = zeros(numNodes*DOF,numNodes*DOF);
    
    for i = 1:numElements
        % need to add 1 as node numbers start at 0
        node1 = sctrData(i,1);
        node2 = sctrData(i,2);
        
        a = cos(theta(i));
        b = sin(theta(i));

        K(node1*2-1,node1*2-1) = K(node1*2-1,node1*2-1)+k(i)*(a^2);
        K(node1*2-1,node1*2) = K(node1*2-1,node1*2)+k(i)*(a*b);
        K(node1*2,node1*2-1) = K(node1*2,node1*2-1)+k(i)*(a*b);
        K(node1*2,node1*2) = K(node1*2,node1*2)+k(i)*(b^2);

        K(node1*2-1,node2*2-1) = K(node1*2-1,node2*2-1)-k(i)*(a^2);
        K(node1*2-1,node2*2) = K(node1*2-1,node2*2)-k(i)*(a*b);
        K(node1*2,node2*2-1) = K(node1*2,node2*2-1)-k(i)*(a*b);
        K(node1*2,node2*2) = K(node1*2,node2*2)-k(i)*(b^2);
  
        K(node2*2-1,node1*2-1) = K(node2*2-1,node1*2-1)-k(i)*(a^2);
        K(node2*2-1,node1*2) = K(node2*2-1,node1*2)-k(i)*(a*b);
        K(node2*2,node1*2-1) = K(node2*2,node1*2-1)-k(i)*(a*b);
        K(node2*2,node1*2) = K(node2*2,node1*2)-k(i)*(b^2);
    
        K(node2*2-1,node2*2-1) = K(node2*2-1,node2*2-1)+k(i)*(a^2);
        K(node2*2-1,node2*2) = K(node2*2-1,node2*2)+k(i)*(a*b);
        K(node2*2,node2*2-1) = K(node2*2,node2*2-1)+k(i)*(a*b);
        K(node2*2,node2*2) = K(node2*2,node2*2)+k(i)*(b^2);

    end
    
end

function U = calculateU(numNodes, nodeTypes)
    % For node Types:
    %1 - free, 2 - fixed, 3 - roller in x, 4 - roller in y
    DOF = 2;
    U = zeros(numNodes*DOF,2);
    
    for i = 1:numNodes
        if nodeTypes(i) == 2
            U(i*2-1,1) = 0;
            U(i*2,1) = 0;
            U(i*2-1,2) = 1;
            U(i*2,2) = 1;
        elseif nodeTypes(i) == 3
            U(i*2,1) = 0;
            U(i*2,2) = 1;
        elseif nodeTypes(i) == 4
            U(i*2-1,1) = 0;
            U(i*2-1,2) = 1;
        end
    end
end

function F = calculateF(numNodes, nodeTypes, appliedForces)
    % For node Types:
    %1 - free, 2 - fixed, 3 - roller in x, 4 - roller in y
    DOF = 2;
    [numAppliedForces,col] = size(appliedForces);
    
    F = zeros(numNodes*DOF,2);
        
    for i = 1:numNodes
        if nodeTypes(i) == 1
            F(i*2-1,1) = 0;
            F(i*2,1) = 0;
            F(i*2-1,2) = 1;
            F(i*2,2) = 1;
        elseif nodeTypes(i) == 2
            %none
        elseif nodeTypes(i) == 3
            F(i*2,1) = 0;
            F(i*2,2) = 1;
        elseif nodeTypes(i) == 4
            F(i*2-1,1) = 0;
            F(i*2-1,2) = 1;
        end
    end

    if numAppliedForces > 0
        for i = 1:numAppliedForces
            node = appliedForces(i,1);
            F(node*2-1,1) = appliedForces(i,2);
            F(node*2,1) = appliedForces(i,3);
        end
    end
end

%This function solves the matrix once it has been rearranged.
%
%This function takes 4 paramenters: -the "A", or global stiffness matrix   
%                                   -the "X", or global displacement matrix
%                                   -the "B", or force matrix 
%                                   
%It will return one matrix fully solved:  the "X" (displacement matrix)
%                                         AX = B
function X = solveMatrix(A,B,X)
    
% !!! NEED TO EXPLAIN IN POWERPOINT WHATS GOING ON HERE

    [rows,cols] = size(A);

    % if it is a mixed solution, then split up the "A" matrix into Ke, Kef,
    % Kfe, Kf matrices to solve. Then at the end, it will combine Xe and Xf
    % matrices as well as the Be and Bf matrices to return the full X and B
    % matrices. 
    Ke = A([3:rows-2],[3:cols-2]);
    Kef = [A([3:rows-2],[1:2]) A([3:rows-2],[cols-1:cols])];
%     Kfe = A([e+1:n], [1:e]);
%     Kf = A([e+1:n], [e+1:n]);

    Xe = [X([1:2],1); X([rows-1:rows],1)];
    Bf = B([3:rows-2], 1);

    % Use LUD to solve for Xf
    C = Bf - (Kef * Xe);

% %          % In the form of Kf*Xf = C
%          [L, U, Cprime] = LUD(Ke, C);
% %         % Now we have U*Xf = Cprime. We can solve for Xf using back
% %          % substitution
%          Xf = solveWithBackSubstitution(U, Cprime);

    Xf = gaussSeidel(Ke, C);
%     Be = (Ke * Xe) + (Kef * Xf);
%     B([1:e]) = Be;
    X([3:rows-2]) = Xf;
end

function [L, U, Cprime] = LUD(A, C)
    %Get size of A (m should = n)
    [Rows,Columns] = size(A);
    
    %Declare empty L and U matrices
    L = zeros(Rows, Columns);
    U = eye(Rows, Columns);
    Cprime = zeros(Rows, 1);
    
    %Outer loop iterates from i2j1 to the number of Rows
    %This is the number of times rows/column chunks in the L/U matrix must be filled in
    for i = 1:Rows
        %This loop handles L matrix
        for j = i:Rows
            %summationResult is the result of summing L(i,k)*U(k,j)
            summationResult = 0;
            k = 1;
            
            %This while loop calculates the values in the summation in the
            %formula for each value of L
            while (k <= i - 1)
                summationResult = summationResult + (L(j,k)*U(k,i));
                k = k + 1;
            end           
            
            %calculate the value of L(j,i)
            L(j,i) = A(j,i) - summationResult; 
        end
        
        %This loop handles U matrix
        for j = i:Columns
            %summationResult is summation of L(i,k)*U(k,j)
            summationResult = 0;
            k = 1;
            
            %calculate summationResult
            while (k <= i - 1)
                summationResult = summationResult + (L(i,k)*U(k,j));
                k = k + 1;
            end
            
            %calculate value of U(i,j)
            U(i,j) = (A(i,j) - summationResult)/L(i,i);
        end
    end
    
    Cprime(1) = C(1)/L(1,1);
    
    for i = 2:Rows
        summationResult = 0;
        k = 1;
        
        %using formula from a textbook (Applied Numerical Analysis,Gerald/Wheatly
        while (k <= i - 1)
            summationResult = summationResult + (L(i,k)*Cprime(k));
            k = k + 1;
        end
        
        Cprime(i) = (C(i) - summationResult)/L(i,i); 
    end
end

function Xf = solveWithBackSubstitution(U, Cprime)
    % create Xf which is the same size of 
    rows = size(Cprime,1);
    Xf = zeros(rows,1);
    
    for i = rows:-1:1
        sum = Cprime(i);
        for j = rows:-1:i
            sum = sum - (Xf(j)*U(i,j));
        end
        Xf(i) = sum;
    end

end

function X = gaussSeidel(K, F)
    %In the form of F=kx
    %Get size of array (length of K, F, and X)
    N = length(K);
    
    %Use gpuArray
    %K = gpuArray(Kcpu);
    %F = gpuArray(Fcpu);
    
    %'currentGuess' is the most up to date values for X - change it below
    %to modify the initial guess.
    currentGuess = zeros(N,1);
    
    %'relativeTolerance' is the tolerance value below which iterations for
    %the Gauss-Seidel method are halted. Change the tolerance below.
    relativeTolerance = 0.01;
    
    %'toleranceFlag' is used to signal that iterations should be stopped. 
    %When the flag is set, the program completes the current iteration and
    %then exits. So long as one value of x[i] is not below the tolerance at
    %an iteration, the program will continue.
    toleranceFlag = false;
    while(toleranceFlag ~= true)
        previousGuess = currentGuess;
        for i = 1:N            
            toleranceFlag = true;
            summation = 0;

            for j = 1:N
                if j ~= i
                    summation = summation + (K(i,j)*currentGuess(j)); 
                end
            end
            
            currentGuess(i) = (F(i)-summation)/K(i,i);
            error = abs(currentGuess(i)-previousGuess(i))/abs(previousGuess(i));
                    
            if error >= relativeTolerance
                toleranceFlag = false; 
            end 
        end 
    end
    X = currentGuess;
end