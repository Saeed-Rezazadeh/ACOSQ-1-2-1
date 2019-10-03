function Probability_y_1 = Pr_y_1 (y_1 , f , T , delta , Pr)

summation = 0 ;
for x_1 = 1 : 2
    u_index = find (T(: , 2) == x_1) ;
    summation = summation + Pr(x_1 , y_1) * delta * sum(f(u_index)) ;
end
Probability_y_1 = summation ;
end