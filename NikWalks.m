function FWalk = NikWalks(Fapplied, j, timeStep)

    t = timeStep * j; 
    i = floor(t / 1.6);
    % Nik doesn't start walking on the first element and gets off before
    % last element
    if (i == 0 || i == 512)
        FWalk = Fapplied;
        
    else
        %using a sine function to model Nik's weight (88kg/863.28N) 
        %1.6s the time it takes for Nik to take a step 
        Fprev = (863.28/2)*sin(pi/1.6 * t) + (863.28/2);
        Fnow = (-863.28/2)*sin(pi/1.6 * t) + (863.28/2); 

        FWalk = Fapplied;
        Fwalk((i*2),1) = Fnow + Fapplied((i*2),1); 
        Fwalk(((i-1)*2),1) = Fprev + Fapplied(((i-1)*2),1);
    end
end