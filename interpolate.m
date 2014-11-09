function new = interpolate(old)
    
    new = zeros((length(old)-1)*2,1);

    for i = 1:length(new)/2

        disp(i);
        
        if (mod(i,2))
            
            new(i*2-1) = old(i);
            new(i*2) = old(i+1);
            
        else
            
            new(i*2-1) = (old((i+1))+old(i-1))/2;
            new(i*2) = (old(i)+old(i+2))/2;
            
        end
        
    end
    
    csvwrite('new.csv', new);
            

end