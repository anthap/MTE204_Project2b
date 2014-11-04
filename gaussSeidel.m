function X = gaussSeidel(K, F)

    %Get size of array (length of K, F, and X)
    N = length(K);
    
    %'currentGuess' is the most up to date values for X - change it below
    %to modify the initial guess
    currentGuess = zeros(N,1);
    
    %'relativeTolerance' is the tolerance value below which iterations for
    %the Gauss-Seidel method are halted. 
    relativeTolerance = 0.01;
    toleranceFlag = false;
    
    while(toleranceFlag ~= true)
        
        previousGuess = currentGuess;
        
        for i = 1:N
            
            summation = 0;
            
            for j = 1:N
                
                if j ~= i
                    
                    summation = summation + (K(i,j)*currentGuess(j));
                    
                end
                
            end
            
            currentGuess(i) = (F(i)-summation)/K(i,i);
            
            error = abs(currentGuess(i)-previousGuess(i))/abs(previousGuess(i));
                    
            if (error < relativeTolerance)
                
                toleranceFlag = true;
                
            end
            
        end
        
    end
    
    X = currentGuess;
 
end
