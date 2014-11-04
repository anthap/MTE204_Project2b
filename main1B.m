function main(questionNum, solveMethod, timeStep, endTime, nodeToEvaluate)

    %Begin by importing truss information from a CSV file
    %Data is stored in the arrays below and accessed throughout

    %Data in files follows format:
    %x,y
    nodesData = csvread(strcat('nodes', num2str(questionNum), '.txt'));
    %node1,node2
    sctrData = csvread(strcat('sctr', num2str(questionNum), '.txt'));
    %E,Area,Density,Dampening
    propsData = csvread(strcat('props', num2str(questionNum), '.txt'));
 
    %Retrieve number of nodes
    numNodes = size(nodesData,1);
    
    %Get data on node types based on the question
    nodeTypes = getNodeTypes(questionNum, numNodes);
    
    %Get data on applied forces based on the question
    appliedForces = getAppliedForces(questionNum);

    %Retrieve number of elements and element information
    numElements = size(propsData,1);
    
    %calculate and store all data necessary to build global matrices and theta for each element
    c = zeros(numElements,1);
    m = zeros(numElements,1);
    k = zeros(numElements,1);
    theta = zeros(numElements,1);
    
    for i = 1:numElements
        node1 = sctrData(i,1);
        node2 = sctrData(i,2);
        x1 = nodesData(node1,1);
        y1 = nodesData(node1,2);
        x2 = nodesData(node2,1);
        y2 = nodesData(node2,2);
        theta(i) = atan2d(y2-y1,x2-x1);
        length = sqrt(((x2-x1)^2)+((y2-y1)^2));
        E = propsData(i,1);
        area = propsData(i,2);
        k(i) = (area*E/length);
        c(i) = propsData(i,4);
        m(i) = length*area*propsData(i,3);
    end
       
    %Create displacements and forces matrix; second column holds wether or not values
    %are known
    U = calculateU(numNodes, nodeTypes);  
    F = calculateF(numNodes, nodeTypes, appliedForces);

    %Create [M],[C],[Keff]  
    M = calculateM(numNodes, numElements, sctrData, m);
    C = calculateK(numNodes, numElements, sctrData, c, theta);
    K = calculateK(numNodes, numElements, sctrData, k, theta);
    
    if solveMethod == 1
        solveExplicit(questionNum,M,C,K,F,U,timeStep,endTime,nodeToEvaluate);
    elseif solveMethod == 2
        solveImplicit(questionNum,M,C,K,F,U,timeStep,endTime,nodeToEvaluate);
    end
end

function nodeTypes = getNodeTypes(questionNum, numNodes)
    % returns an array containing data for each node type for a given
    % question
    % 1 - free, 2 - fixed, 3 - roller in x, 4 - roller in y
    nodeTypes = zeros(numNodes,1);
    
    if questionNum == 1
        nodeTypes(1) = 2;
        nodeTypes(2) = 1;
    elseif questionNum == 2
        nodeTypes(1) = 2;
        nodeTypes(2) = 1;
        nodeTypes(3) = 1;
        nodeTypes(4) = 1;
        nodeTypes(5) = 1;
        nodeTypes(6) = 1;
        nodeTypes(7) = 1;
        nodeTypes(8) = 1;
        nodeTypes(9) = 1;
        nodeTypes(10) = 1;
        nodeTypes(11) = 1;
    elseif questionNum == 3
        nodeTypes(1) = 4;
        nodeTypes(2) = 4;
        nodeTypes(3) = 1;
        nodeTypes(4) = 1;
        %while node 5 is in fact a y-roller, we designated it as a fixed
        %node as its displacement in y is known (it is applied to the system)
        nodeTypes(5) = 2;
    end
end

function appliedForces = getAppliedForces(questionNum)
    %Applied Forces matrix in format of:
    %Node#,Force in X, Force in Y
    if questionNum == 1
        appliedForces = zeros(1,3);
        appliedForces(1,1) = 2;
        appliedForces(1,2) = 10;
        appliedForces(1,3) = 0;
    elseif questionNum == 2
        %time variant, will be calculated for each time iteration
        appliedForces = zeros(1,3);
        appliedForces(1,1) = 11;
        appliedForces(1,2) = 1;
        appliedForces(1,3) = 0;
    elseif questionNum == 3
        appliedForces = int16.empty(0,0);
    end
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
        
        a = cosd(theta(i));
        b = sind(theta(i));

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

function solveImplicit(questionNum,M,C,K,F,U,timeStep,endTime,nodeToEvaluate);
    
    starttime = cputime;

    gamma = 3/2;
    beta = 8/5;
    DOF = 2; 
    
    [rows, cols] = size(K);
    %size64 ensures the value passed to zeros() is an int
    size64 = int64(endTime/timeStep + 1);
    keepTrackU = zeros(size64, DOF + 1);
    keepTrackVel = zeros(size64, DOF + 1);
    keepTrackAcc = zeros(size64, DOF + 1);
    
    [Forganize, Uorganize, Morganize, Corganize, Korganize, knowns, order] = organizeMatrix(F, U, M, C, K); 
    Uorganize = Uorganize(:,1);
    Forganize = Forganize(:,1);
    
    Atilde = 2*Morganize/(beta*(timeStep^2)) + Korganize + 2*Corganize*gamma/(beta*timeStep);
    Btilde = 2*Morganize/(beta*(timeStep^2)) + 2*Corganize*gamma/(beta*timeStep);
    Ctilde = 2*Morganize/(beta*timeStep) + Corganize*(2*gamma/beta - 1);
    Dtilde = Morganize*(1 - beta)/beta + timeStep*Corganize*((gamma - 1) + gamma*(1 - beta)/beta);

   %find X position of node to apply force
    for n = 1:rows
        if order(n) == (nodeToEvaluate*2 - 1)
            nodeToEvaluateX = order(n); 
        end
    end
    
    %find Y position of node to apply force
    for n = 1:rows
        if order(n) == nodeToEvaluate*2
            nodeToEvaluateY = order(n); 
        end
    end
    
    %set initial U,V,a

    UorganizeCurr = zeros(rows, 1); %THIS COULD be a matrix passed in if IC's aren't zero. 
    VorganizeCurr = zeros(rows, 1);
    aorganizeCurr = zeros(rows, 1); 
    
    %make matrices to record each iteration 
    keepTrackU(1, 2) = UorganizeCurr(nodeToEvaluateX);
    keepTrackU(1, 3) = UorganizeCurr(nodeToEvaluateY);
    keepTrackVel(1, 2) = keepTrackU(1, 2)/timeStep;
    keepTrackVel(1, 3) = keepTrackU(1, 3)/timeStep; 
    keepTrackAcc(1, 2) = keepTrackVel(1,2)/timeStep; 
    keepTrackAcc(1, 3) = keepTrackVel(1,3)/timeStep;
    
    %Find Y position of node to apply displacement (for Q3)
    for n = 1:rows
        if order(n) == 10
            nodeToDisplaceY = n; 
        end
    end

    i = 1;
    for (j = 0: timeStep: endTime)
        if questionNum == 2
            Forganize(order(nodeToEvaluateX)) = j;
        elseif questionNum == 3
            %applying a sinusoid displacement to node 5 in the y-direction
            omega = 0.1;
            applyNode = 5;
            Uorganize(nodeToDisplaceY) = 50 * sin(omega*j);
        end

        ForganizeAtNextTimeStep = Forganize + Btilde*UorganizeCurr + Ctilde*VorganizeCurr + Dtilde*aorganizeCurr; 
        [UorganizeNext, ForganizeAtNextTimeStep] = solveMatrix(Atilde, ForganizeAtNextTimeStep, Uorganize, knowns);  
        
        aorganizeNext = (2/(beta*timeStep))*((UorganizeNext - UorganizeCurr)/timeStep)-((2/(beta*timeStep))*VorganizeCurr)- (((1-beta)/beta)*aorganizeCurr);
        VorganizeNext = VorganizeCurr + timeStep*((1-gamma)*aorganizeCurr + gamma*aorganizeNext);
        
        keepTrackU(i + 1, 1) = j;
        keepTrackU(i + 1, 2) = UorganizeNext(nodeToEvaluateX);
        keepTrackU(i + 1, 3) = UorganizeNext(nodeToEvaluateY);
        keepTrackVel(i + 1, 1) = j;
        keepTrackVel(i + 1, 2) = VorganizeNext(nodeToEvaluateX);
        keepTrackVel(i + 1, 3) = VorganizeNext(nodeToEvaluateY);
        keepTrackAcc(i + 1, 1) = j;
        keepTrackAcc(i + 1, 2) = aorganizeNext(nodeToEvaluateX);
        keepTrackAcc(i + 1, 3) = aorganizeNext(nodeToEvaluateY);
        
        UorganizeCurr = UorganizeNext; 
        VorganizeCurr = VorganizeNext; 
        aorganizeCurr = aorganizeNext; 
        
        i = i + 1;
    end
    
    %Commented disp cuts off y displacement
    %disp(keepTrackU(:,1:2));
    %disp(keepTrackVel(:,1:2));
    %disp(keepTrackAcc(:,1:2));
    disp(keepTrackU);
    disp(keepTrackVel);
    disp(keepTrackAcc);
    disp('IMPLICIT COMPUTATION TIME (s):');
    disp(cputime - starttime);
end

function F = invert(M)
    determinant = det(M);
    [i, j] = size(M);
    trans = zeros(j,i);
    for a = 1:i
        for b = 1:j
            trans(a,b) = M(b,a);
        end
    end
    F = (1/determinant)*trans;
end

function [Forganize, Uorganize, Morganize, Corganize, Korganize, knowns, order] = organizeMatrix(Fmatrix, Umatrix, Mmatrix, Cmatrix, Kmatrix)
    knowns = 0; 
    [row, col] = size(Kmatrix); 
    for i = 1:row;
        if Umatrix(i,2) == 1;
            knowns = knowns + 1;
        end
    end        
    
    %assume the matrices are always the same height (A)
    Uorganize = Umatrix;
    Forganize = Fmatrix;
    Morganize = Mmatrix;
    Corganize = Cmatrix;
    Korganize = Kmatrix;
    order = transpose(1:row);
    for i = 1:row
        if Forganize(i,2) == 1
            j = i + 1;
            while j <= row
                if Forganize(j,2) == 0;
                    Forganize([i j],:) = Forganize([j i],:);
                    Morganize([i j],:) = Morganize([j i],:);
                    Corganize([i j],:) = Corganize([j i],:);
                    Korganize([i j],:) = Korganize([j i],:);
                    order([i j],:) = order([j i],:);
                    break
                end
                j = j + 1;
            end
            
        end
    end
    
    for i = 1:row
        if Uorganize(i,2) == 0;
            j = i+1;
            while j <= row 
                if Uorganize(j,2) == 1;
                    Uorganize([i j],:) = Uorganize([j i],:);
                    Morganize(:,[i j]) = Morganize(:,[j i]);
                    Corganize(:,[i j]) = Corganize(:,[j i]);
                    Korganize(:,[i j]) = Korganize(:,[j i]);
                    break
                end
                j = j + 1;
            end
        end            
    end
end

function [keepTrackU,keepTrackVel,keepTrackAcc] = solveExplicit(questionNum, M, C, K, F, U, timestep, endTime, nodeToEvaluate) 
    
    starttime = cputime;

    DOF = 2; 
    
    [rows, cols] = size(K);
    %size64 ensures the value passed to zeros() is an int
    size64 = int64(endTime/timestep + 1);
    keepTrackU = zeros(size64, DOF + 1);
    keepTrackVel = zeros(size64, DOF + 1);
    keepTrackAcc = zeros(size64, DOF + 1);

    %organize matrix function
    [Forganize, Uorganize, Morganize, Corganize, Korganize, knowns, order] = organizeMatrix(F, U, M, C, K); 
    Uorganize = Uorganize(:,1);
    Forganize = Forganize(:,1);
    
    A = (Morganize/(timestep^2)) + (Corganize/(2*timestep));
    B = Korganize - ((2*Morganize)/(timestep^2));
    D = (Morganize/(timestep^2)) - (Corganize/(2*timestep));
    
    %find X position
    for n = 1:rows
        if order(n) == (nodeToEvaluate*2 - 1)
            nodeToEvaluateX = order(n); 
        end
    end
    
    %find Y position 
    for n = 1:rows
        if order(n) == nodeToEvaluate*2
            nodeToEvaluateY = order(n); 
        end
    end
    
    %Set up variable for previous U starting with all zeros
    UorganizePrev = zeros(rows,1);
    UorganizeCurr = Uorganize;

    %keeptrackthings of special node to keeptrack of
    keepTrackU(1, 2) = UorganizeCurr(nodeToEvaluateX);
    keepTrackU(1, 3) = UorganizeCurr(nodeToEvaluateY);
    keepTrackVel(1, 2) = keepTrackU(1, 2)/timestep;
    keepTrackVel(1, 3) = keepTrackU(1, 3)/timestep; 
    keepTrackAcc(1, 2) = keepTrackVel(1,2)/timestep; 
    keepTrackAcc(1, 3) = keepTrackVel(1,3)/timestep;
    
    %store all the values in a matrix with # of col = num of timesteps
    i = 1;
  
    for j = 0:timestep:endTime
        %Apply question specific values
        if questionNum == 2
            Forganize(order(nodeToEvaluateX)) = j;
        end
       
        ForganizeAtTimestep = Forganize - B*UorganizeCurr - D*UorganizePrev;    
        
        [UorganizeNext, ForganizeAtTimestep] = solveMatrix(A, ForganizeAtTimestep, Uorganize, knowns);  
        
        keepTrackU(i + 1, 1) = j;
        keepTrackU(i + 1, 2) = UorganizeNext(nodeToEvaluateX);
        keepTrackU(i + 1, 3) = UorganizeNext(nodeToEvaluateY);
        keepTrackVel(i + 1, 1) = j;
        keepTrackVel(i + 1, 2) = (keepTrackU(i + 1, 2) - keepTrackU(i, 2))/timestep;
        keepTrackVel(i + 1, 3) = (keepTrackU(i + 1, 3) - keepTrackU(i, 3))/timestep; 
        keepTrackAcc(i + 1, 1) = j; 
        keepTrackAcc(i + 1, 2) = (keepTrackVel(i + 1, 2) - keepTrackVel(i, 2))/timestep; 
        keepTrackAcc(i + 1, 3) = (keepTrackVel(i + 1, 3) - keepTrackVel(i, 3))/timestep;

        UorganizePrev = UorganizeCurr;
        UorganizeCurr = UorganizeNext;
        i = i + 1;
    end
    
    %Commented disp cuts off y displacement
    %disp(keepTrackU(:,1:2));
    %disp(keepTrackVel(:,1:2));
    %disp(keepTrackAcc(:,1:2));
    disp(keepTrackU);
    disp(keepTrackVel);
    disp(keepTrackAcc);
    disp('EXPLICIT COMPUTATION TIME (s):');
    disp(cputime - starttime);
end

%This function solves the matrix once it has been rearranged.
%
%This function takes 4 paramenters: -the "A", or global stiffness matrix   
%                                   -the "X", or global displacement matrix
%                                   -the "B", or force matrix 
%                                   
%It will return one matrix fully solved:  the "X" (displacement matrix)
%                                         
function [X, B] = solveMatrix(A,B,X,e)

    [m,n] = size(A);
    [o,p] = size(B);
    
    %error handling in case the dimensions given aren't right:
    if m~=n
        X = false; 
    end
    i = 1;
   
    if p ~= 1 && o~=n
        B = false; 
    end

    %if it is not a mixed solution, Find the "X" matrix or "B" matrix
    %accordingly. Returns both of the matrices. 
    if e == n
        B = A*X;
        
    % if it is a mixed solution, then split up the "A" matrix into Ke, Kef,
    % Kfe, Kf matrices to solve. Then at the end, it will combine Xe and Xf
    % matrices as well as the Be and Bf matrices to return the full X and B
    % matrices. 

    else
        
        Ke = A([1:e],[1:e]);
        Kef = A([1:e],[e+1:n]);
        Kfe = A([e+1:n], [1:e]);
        Kf = A([e+1:n], [e+1:n]);

        Xe = X([1:e],1);
        Bf = B([e+1:n], 1);

        % Use LUD to solve for Xf
        C = Bf - (Kfe * Xe);
        
        % In the form of Kf*Xf = C
        [L, U, Cprime] = LUD(Kf, C);
        
        % Now we have U*Xf = Cprime. We can solve for Xf using back
        % substitution
        Xf = solveWithBackSubstitution(U, Cprime);
        
        Be = (Ke * Xe) + (Kef * Xf); 

        B([1:e]) = Be;
        X([e+1:n]) = Xf;
    end

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
