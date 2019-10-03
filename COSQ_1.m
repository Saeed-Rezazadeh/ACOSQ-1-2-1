function [SDR , Distortion , T , codebook] = COSQ_1 (Pr , f , T ,  codebook , delta)

FileID = fopen ('Results.txt' , 'a') ;
D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal partitions
    for u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros (1 , 2) ;
        u = T(u_index , 1) ;
        for x_1 = 1 : 2
            for y_1 = 1 : 2
                summation = summation + Pr(x_1 , y_1) * (u - codebook (y_1)) ^ 2 ;
            end
            temp (x_1) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min (temp) ;
        T_u (u_index) = partition_index ;
    end
    T (: , 2) = T_u ;
    
    %% Optimal Centroids
    parfor y_1 = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for x_1 = 1 : 2
            u_index = find (T (: , 2) == x_1) ;
            u = T(u_index , 1) ;
            numerator = numerator + Pr (x_1 , y_1) * sum (u .* f(u_index)) ;
            denominator = denominator + Pr (x_1 , y_1) * sum (f(u_index)) ;
        end
        codebook(y_1) = numerator / denominator ;
    end
    %% Distortion 
    D(2) = distortion_1 (f , codebook , delta , Pr , T ) ;
    fprintf (FileID , 'Overall D = %f\n' ,D(2)) ; 
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR = %f\n' , SDR) ; 
fclose (FileID) ; 
end