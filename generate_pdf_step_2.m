function [f_u_given_y_1] = generate_pdf_step_2(y_1 , T , Pr , f , delta )
%% numerator
numerator = zeros(length(T) , 1) ;
for u_index = 1 : length(T)
    x_1 = T(u_index , 2) ;
    numerator (u_index) = Pr(x_1 , y_1) * f(u_index) ;
end

%% denominator
summation = 0 ; 
for x_1 = 1 : 2
    u_index_x = find(T(: , 2) == x_1) ;
    summation = summation + Pr(x_1 , y_1) * sum(f(u_index_x)) * delta ;
end
denominator = summation ;
f_u_given_y_1 = numerator ./ denominator ;
f_u_given_y_1 = f_u_given_y_1 ./ (sum(f_u_given_y_1) * delta ) ;
end