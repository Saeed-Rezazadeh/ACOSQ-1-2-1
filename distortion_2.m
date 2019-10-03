function Distortion = distortion_2(f ,  y_1 , codebook , delta , Pr_z , T)
summation = 0 ;
for x_2 = 1 : 2
    for x_3 = 1 : 2
        x_prime = (x_2 - 1) * 2 + x_3 ;
        for y_2 = 1 : 2
            for y_3 = 1 : 2
                y_prime = (y_2 - 1) * 2 +y_3 ;
                
                u_index = find(T(: , 3) == x_prime) ;
                for u_i = 1 : length(u_index)
                    x_1 = T(u_index(u_i) , 2) ;
                    summation = summation + delta * Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                        * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                        * (T(u_index(u_i) , 1) - codebook(y_prime)) ^ 2 * f(u_index(u_i)) ;
                    
                end
            end
        end
    end
end
Distortion = summation ;
end