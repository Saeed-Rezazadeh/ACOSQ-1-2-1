function [SDR , Distortion , T , codebook] = COSQ_4(f  , y_1 , y_2 , y_3 , Pr_z , T , codebook , delta )
y_1_2 = (y_1 - 1) * 2 + y_2 ;
FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    T_u = zeros(length(T) , 1) ;
    parfor u_index = 1 : length(T)
        d_4 = zeros(2 , 1) ;
        summation = 0 ;
        
        u = T(u_index , 1) ;
        x_prime = T(u_index , 2 + y_1) ;
        if (mod(x_prime - 1, 2) == 0)
            x_3 = 1 ;
        else
            x_3 = 2 ;
        end
        for x_4 = 1 : 2
            for y_4 = 1 : 2
                summation = summation ...
                    + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) ...
                    * (u - codebook(y_4)) ^ 2 ;
            end
            d_4(x_4) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min(d_4) ;
        T_u(u_index , 1) = partition_index ;
    end
    T(: , 5) = T_u ;
    
    %% Optimal Centroids
    for y_4 = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for x_4 = 1 : 2
            u_index = find (T(: , 5) == x_4) ;
            
            for u_i = 1 : length(u_index)
                x_prime = T(u_index(u_i) , 2 + y_1) ;
                if (mod(x_prime - 1, 2) == 0)
                    x_3 = 1 ;
                else
                    x_3 = 2 ;
                end
                numerator = numerator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1)...
                    * T(u_index(u_i) , 1) * f(u_index(u_i)) ;
                denominator = denominator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1)...
                    * f(u_index(u_i)) ;
                
            end
        end
        codebook(y_4) = numerator / denominator ; 
    end
    %% Distortion
    D(2) = distortion_3(f , y_1 , y_2 , y_3 , codebook , delta , Pr_z , T) ;
    fprintf (FileID , 'Overall D = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR = %f\n' , SDR) ;
fclose (FileID) ;
end