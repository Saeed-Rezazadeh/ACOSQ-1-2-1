function [SDR , Distortion , T , codebook] = COSQ_2(f , y_1 , Pr_z , T , codebook , delta )

FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    T_u = zeros(length(T) , 1) ;
    parfor u_index = 1 : length(T)
        d_2 = zeros(4 , 1) ;
        summation = 0 ;
        
        u = T(u_index , 1) ;
        x_1 = T(u_index , 2) ;
        for x_2 = 1 : 2
            for x_3 = 1 : 2
                x_prime = (x_2 - 1) * 2 + x_3 ;
                for y_2 = 1 : 2
                    for y_3 = 1 : 2
                        y_prime = (y_2 - 1) * 2 + y_3 ;
                        
                        summation = summation + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 )...
                            * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) * (u - codebook(y_prime)) ^ 2 ;
                    end
                end
                d_2(x_prime) = summation ;
                summation = 0 ;
            end
        end
        [~ , partition_index] = min(d_2) ;
        T_u(u_index , 1) = partition_index ;
    end
    T(: , 3) = T_u ;
    %% Optimal Centroids
    for y_2 = 1 : 2
        for y_3 = 1 : 2
            y_prime = (y_2 - 1) * 2 + y_3 ;
            
            numerator = 0 ;
            denominator = 0 ;
            for x_2 = 1 : 2
                for x_3 = 1 : 2
                    x_prime = (x_2 - 1) * 2 + x_3 ;
                    u_index = find (T(: , 3) == x_prime) ;
                    for u_i = 1 : length(u_index)
                        x_1 = T(u_index(u_i) , 2) ;
                        
                        numerator = numerator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                            * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                            * T(u_index(u_i) , 1) * f(u_index(u_i)) ;
                        
                        denominator = denominator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                            * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                            * f(u_index(u_i)) ;
                        
                    end
                end
            end
            codebook(y_prime) = numerator / denominator ; 
        end
    end
    %% Distortion
    D(2) = distortion_2(f , y_1 , codebook , delta , Pr_z , T) ;
    fprintf (FileID , 'Overall D = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR = %f\n' , SDR) ;
fclose (FileID) ;
end