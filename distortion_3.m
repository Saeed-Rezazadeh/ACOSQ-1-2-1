function Distortion = distortion_3(f ,  y_1 , y_2 , y_3 , codebook , delta , Pr_z , T)
summation = 0 ; 
for x_4 = 1 : 2 
    for y_4 = 1 : 2 
        u_index = find (T(: , 5) == x_4) ; 
        for u_i = 1 : length(u_index) 
            x_prime = T(u_index(u_i) , 2 + y_1) ; 
            if (mod(x_prime - 1, 2 ) == 0) 
                x_3 = 1 ;
            else
                x_3 = 2 ; 
            end 
            summation = summation + delta ...
                * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1)... 
                * f(u_index(u_i)) * (T(u_index(u_i) , 1) - codebook(y_4)) ^ 2 ;
        end 
    end
end 
Distortion = summation ; 
end