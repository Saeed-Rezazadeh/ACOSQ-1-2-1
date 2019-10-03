function [ f_u_given_y_1_y_2] = generate_pdf_step_3(Pr_z , f_u_given_y_1 , T , y_1 , y_2 ,y_3 , delta)
%% numerator
numerator = zeros (length(T) , 1) ;
for u_index = 1 : length(T)
    x_1 = T(u_index , 2 ) ;
    
    x_prime = T(u_index , 2 + y_1) ;
    binary_x_prime = de2bi(x_prime - 1 , log2(4) , 'left-msb') ;
    
    x_2 = binary_x_prime(1)  + 1 ;
    x_3 = binary_x_prime(2)  + 1 ;
    
    numerator(u_index) = Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
        * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
        * f_u_given_y_1(u_index) ;
end

%% denominator
denominator = 0 ;
for x_2 = 1 : 2
    for x_3 = 1 : 2
        x_prime = (x_2 - 1) * 2 + x_3 ;
        
        u_index_x_prime = find(T(: , 2 + y_1) == x_prime) ;
        for u_i = 1 : length(u_index_x_prime)
            x_1 = T(u_index_x_prime(u_i) , 2) ;
            
            denominator = denominator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 )...
                * delta * f_u_given_y_1(u_index_x_prime(u_i)) ;
        end
    end
end
f_u_given_y_1_y_2 = numerator ./ denominator ;
f_u_given_y_1_y_2 = f_u_given_y_1_y_2 ./ (sum(f_u_given_y_1_y_2) * delta ) ;
end