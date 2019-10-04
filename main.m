%% The script corresponds to the Algorithm 4 with r = ( 1 2 1)
clc
clear all
close all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;

%% Source distributions
sigma = 1 ; % the standard devision of the source
u = 0 ; % the source's mean value
alpha = 500000  ;  % the size of the training set to find the initial codebook using the so called splitting algorithm


%% Channel's cross-over probability epsilon
epsilon = unique ([10^-5 : 2 * 10^-5 : 10^ -4 , 10^-4 : 10^ -4 : 10^ -3 , 10^-3 ,0.005 0.01  0.05  0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.
SIZE = length(epsilon) ;
noise =  [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE , SIZE : -1 : 1] ;


% The variable resolution determines the accuracy of the Riemann summation
resolution = 2 ^ 11 ;


%% Initialize parameters
SDR_1 = zeros (length(noise) , 1) ;
Final_SDR_2 = zeros(length(noise) , 1) ;

SDR_4 = zeros (length(noise)  , 1) ;
Final_SDR_4 = zeros(length(noise) , 1) ;


Probability_y_1 = zeros (2 , 1) ;
Probability_y_1_y_2_y_3 = zeros (4 , 1) ;
[Training_set , T , delta] = initialization (sigma , u , alpha , resolution) ;
%% Noise Correlation
for noise_Delta = [0 5 10]
    D_2 = zeros (2 , 1) ;
    D_4 = zeros (4 , 1) ;
    
    for k = 1 : length(noise)
        i = noise(k) ;
        % Compute the source pdf. We herein consider a zero-mean
        % unit-variance Gaussian source distribution.
        u = T(: , 1) ;
        f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
        f = f./ (sum(f) * delta ) ;
        
        Pr_z = [(1 - epsilon(i) + noise_Delta) / (1 + noise_Delta)  , epsilon(i) / (1 + noise_Delta) ;
            (1 - epsilon(i)) / (1 + noise_Delta)  , (epsilon(i) + noise_Delta) / (1 + noise_Delta)] ;
        
        Pr_1 = [1 - epsilon(i) , epsilon(i) ;
            epsilon(i) , 1 - epsilon(i)] ;
        
        %% COSQ for bit 1 (step 1)
        
        if (k == 1)
            % Using the splitting algorithm to find the initial codebook.
            [~  , codebook] = kmeans (Training_set , 2 , 'MaxIter',1000 , 'OnlinePhase','on') ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon.
            load ('codebook_rate_1.mat' , 'codebook')
        end
        % The first step of the ACOSQ described in Section 4.1 and
        % Algorithm 4. In this step a 1-bit COSQ is designed for the source pdf f as computed in line 51.
        [SDR_1(k) ,  D_1 , T , codebook] = COSQ_1(Pr_1 , f , T(: , 1) ,  codebook , delta) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon. In other words, the codebook obtained at every step is used to initialize
        % the quantizer in the SAME step for the new epsilon.
        save ('codebook_rate_1' , 'codebook' ) ;
        fprintf (FileID , '\nrate = 1\n') ;
        %% COSQ for bit 2 (step 2)
        for y_1 = 1 : 2
            % Based on (4.8) in the Thesis compute the source conditional pdf given the received
            % sequence y_1 where y_1 is the channel output corresponding to
            % the transmitted sequence in the first step.
            [f_u_given_y_1] = generate_pdf_step_2(y_1 , T , Pr_1 , f , delta ) ;
            
            for inner_noise = 1 : k
                i = noise(inner_noise) ;
                Pr_z = [(1 - epsilon(i) + noise_Delta) / (1 + noise_Delta)  , epsilon(i) / (1 + noise_Delta) ;
                    (1 - epsilon(i)) / (1 + noise_Delta)  , (epsilon(i) + noise_Delta) / (1 + noise_Delta)] ;
                
              
                
                if (inner_noise == 1)
                    % Find the initial codebook for the conditional source
                    % pdf using splitting algorithm.
                    codebook = init_codebook(4 , f_u_given_y_1 , delta , T , alpha) ;
                else
                    % We slightly increase the channel's cross-over probability,
                    % setting the codebook from the system with small epsilon as
                    % the initial state of the system with new epsilon.
                    Data = ['codebook_y_1_' num2str(y_1)] ;
                    load (Data) ;
                end
                % The second step of the ACOSQ described in Section 4.1 and
                % Algorithm 4. In this step a 2-bit COSQ is designed for the conditional source pdf f_u_given_y_1
                % as computed in line 85.
                [SDR_2 , D_2(y_1) , hold_T , codebook] = COSQ_2(f_u_given_y_1 , y_1 , Pr_z , T(: , 1 : 2) , codebook , delta ) ;
            end
            T (: , 2 + y_1) = hold_T(: , 3);
            % Compute the probability of P(Y_1 = y_1) for y_1 = 0 , 1.
            Probability_y_1(y_1) = Pr_y_1 (y_1 , f , T(: , [1 2]) , delta , Pr_1) ;
            Data = ['codebook_y_1_' num2str(y_1)] ;
            save (Data , 'codebook' ) ;
        end
        Probability_y_1 = Probability_y_1 ./ sum(Probability_y_1) ;
        % As mentioned in the Thesis, the ultimate distortion at every step is the
        % weighted sum of the conditional distortion functions given the
        % received sequence y_1.
        Final_D_2 = sum(D_2 .* Probability_y_1) ;
        % Compute the SDR value in the following way as the source is zero
        % mean and unit variance.
        Final_SDR_2(k) = 10 * log10(1 / Final_D_2) ;
        fprintf (FileID , 'Final SDR_2 = %4.2f\n' ,  Final_SDR_2(k)) ;
        fprintf (FileID , '\nrate = 3\n') ;
        
        %% COSQ for rate 4 (step 3)
        codebook_4 = [] ;
        for y_1 = 1 : 2
            for y_2 = 1 : 2
                for y_3 = 1 : 2
                    y_1_2_3 = (y_1 - 1) * 4 + (y_2 - 1) * 2 + y_3  ;
                    
                    [f_u_given_y_1] = generate_pdf_step_2(y_1 , T , Pr_1 , f , delta ) ;
                    % According to (4.8) compute the conditional pdf given the received sequence
                    % y_1y_2y_3 where y_2y_3 is the corresponding channel output for
                    % the sequence transmitted in the second step.
                    [f_u_given_y_1_y_2_y_3] = generate_pdf_step_3(Pr_z , f_u_given_y_1 , T , y_1 , y_2 , y_3 , delta) ;
                    for inner_noise = 1 : k
                        i = noise(inner_noise) ;
                        Pr_z = [(1 - epsilon(i) + noise_Delta) / (1 + noise_Delta)  , epsilon(i) / (1 + noise_Delta) ;
                            (1 - epsilon(i)) / (1 + noise_Delta)  , (epsilon(i) + noise_Delta) / (1 + noise_Delta)] ;
                        
                        
                        if (inner_noise == 1)
                            % Find the initial codebook for the conditional source
                            % pdf using splitting algorithm.
                            codebook = init_codebook(2 , f_u_given_y_1_y_2_y_3 , delta , T , alpha) ;
                        else
                            Data = ['codebook_y_1_2_3_' num2str(y_1_2_3)] ;
                            load (Data) ;
                        end
                        % The last step of the ACOSQ described in Section 4.1 and
                        % Algorithm 4. In this step a 1-bit COSQ is designed for the conditional source pdf f_u_given_y_1_y_2_y_3
                        % as computed in line 138.
                        [SDR_4(k) , D_4(y_1_2_3) , hold_T , codebook] = COSQ_4(f_u_given_y_1_y_2_y_3 , y_1 , y_2 , y_3 , Pr_z , T(: , 1 : 4) , codebook , delta ) ;
                    end
                    codebook_4 = [codebook_4 ; codebook ] ;
                    T(: , 4 + y_1_2_3) = hold_T(: , 5) ;
                    
                    Probability_y_1_y_2_y_3(y_1_2_3) = Pr_y_1_y_2_y_3( y_1, y_2 , y_3 , Pr_1 , Pr_z , f  , T , delta) ;
                    
                    Data = ['codebook_y_1_2_3_' num2str(y_1_2_3)] ;
                    save (Data , 'codebook' ) ;
                    
                    % Store the ultimate partition indexes and codebook to compute
                    % experimental results.
                    Data = ['T\T_k_' num2str(k) '_delta_' num2str(noise_Delta)] ;
                    save (Data , 'T' ,'codebook_4') ;
                end
            end
        end
        Probability_y_1_y_2_y_3 = Probability_y_1_y_2_y_3 ./ sum(Probability_y_1_y_2_y_3) ;
        Final_D_4 = sum(D_4 .* Probability_y_1_y_2_y_3) ;
        Final_SDR_4(k) = 10 * log10(1 / Final_D_4) ;
        fprintf (FileID , 'Final SDR = %4.2f\n' ,  Final_SDR_4(k)) ;
        fprintf (FileID , '\nrate = 4\n') ;
        fprintf (FileID , 'noise = %d\n' ,  i) ;
    end
    fprintf (FileID , '\nSDR for rate 2\n') ;
    clear D_2
    D_2 = zeros (SIZE , 1) ;
    final_SDR_2 = zeros (SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var = Final_SDR_2(index) ;
        hold_var = hold_var(:) ;
        final_SDR_2(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_2(i) = 10 ^(- final_SDR_2(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_2(i)) ;
        fprintf (FileID , '\nSDR_2 = %7.4f' , final_SDR_2(i)) ;
    end
    fprintf (FileID , '\nSDR for rate 4\n') ;
    clear D_4 ;
    D_4 = zeros(SIZE , 1) ;
    final_SDR_4 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var = Final_SDR_4(index) ;
        hold_var = hold_var(:) ;
        final_SDR_4(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_4(i) = 10 ^(- final_SDR_4(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_4(i)) ;
        fprintf (FileID , '\nSDR_4 = %7.4f\n' , final_SDR_4(i)) ;
    end
    % Store the ultimate SDR values for the given channel parameters 
    Data = ['ACOSQ_1_2_1_delta_' num2str(noise_Delta)] ;
    save (Data , 'final_SDR_4' , 'Final_SDR_4' , 'final_SDR_2' , 'Final_SDR_2' , 'epsilon') ;
    clear D_4 D_2;
end
