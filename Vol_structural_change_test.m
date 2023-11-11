
  
function [CUSUM1,QS1,hat_omega2,tau_x_set] = Vol_structural_change_test(y,band)
% PURPOSE:
%     Obtain CUMSUM and QS tests in our paper
%   
% USAGE:
%     [CUSUM1,QS1] = Vol_structural_change_test(y,band)
% 
% INPUTS:
%     y: input time series, row vector
%     band: parameter controlling the values for bandwidth.
%         band == -1: using plug-in method to choose bandwidth
%         band ~= -1: bandwidth is cv*T^(-1/5)       
%
% OUTPUTS:
%     CUSUM1: CUMSUM  test statistic in the paper
%     QS1 : QS test statistic in the paper  
%     hat_omega2 : estimated long run variance \omega^2
%     tau_x_set : nonparametric volatility estimators 
%
% edited by Ruike Wu @ Xiamen University 

    T = length(y);
    m = T^(1/3);
    y2 = abs(y);
    
    if band == -1
        h0 = plugin_tau_basic(y);
    else
        h0 = band*T^(-0.2);
    end


    %    %local constant
    h = h0;
    s_set = 1:T;
    for  x = 1 :T
        ker_set = kernel(T,x,s_set,h0);
        hat_tau_x_set = ker_set.*y2;
        hat_tau_x = sum(hat_tau_x_set)/sum(ker_set);
        tau_x_set(x)= hat_tau_x;
    end
    sigma2_hat = tau_x_set;
    
    
    hat_et = y2 - sigma2_hat;
    hat_gamma_L = (sum(hat_et.^2)/T)*autocorr(hat_et,T-1);
    s_set2 = 1:T-1;
    ker_set_out = Bartlett_kernel((1/m).*s_set2);
    hat_omega2_part2 = sum(ker_set_out.*hat_gamma_L(2:T));
    hat_omega2 = hat_gamma_L(1) + 2*hat_omega2_part2;
    
        
    
    tilde_v = y2 - mean(y2);
    SET_r_set = [];
    for M = 1 :T
        tilde_v_r = tilde_v(1:M);
        ST_r = (T^(-0.5))*abs(sum(tilde_v_r));
        SET_r  = ST_r/sqrt(hat_omega2);
        SET_r_set(M) = SET_r;
    end
    CUSUM1 = max(SET_r_set);
    QS1 = mean(SET_r_set.^2);       
end


