% take in the iteration number and time step and current applied force
% matrix; outputs the applied forces with Nik walking on it

function FWalk = NikWalks(Fapplied, j, timeStep)

    t = timeStep * j; 
    i = floor(t / 1.6);
    % Nik doesn't start walking on the first element and gets off before
    % last element
    if (i == 0 || i == 354)
        FWalk = Fapplied;
        
    else
        %using a sine function to model Nik's weight (88kg/863.28N) 
        %1.6s the time it takes for Nik to take a step 
        Fprev = (863.28/2)*sin(pi/1.6 * t) + (863.28/2);
        Fnow = (-863.28/2)*sin(pi/1.6 * t) + (863.28/2); 

        FWalk = Fapplied;
        Fwalk(i+1) = Fnow + Fapplied(i+1); 
        Fwalk(i) = Fprev + Fapplied(i);
    end
end