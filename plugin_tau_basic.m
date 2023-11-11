function h_plug= plugin_tau_basic(y)
    T = length(y);
    y2 = abs(y);
    % under H0
    sigma0_hat = mean(y2);
    hat_et = y2 - sigma0_hat;
    hat_v  =  hat_et;
    
    c = 2;
    hat_v_lag = hat_v(1:T-1);
    hat_v_no_lag = hat_v(2:T);
    cons_col =  ones(1,length(hat_v_lag));
    X = [cons_col ;hat_v_lag];
    hat_a = ((X*X')^(-1))*X*hat_v_no_lag';
    hat_a(1)  = [];
    omega2_ini = (1/(1-hat_a))^2;
    
    h_plug = c*((omega2_ini)^(1/5))*T^(-0.2);
end