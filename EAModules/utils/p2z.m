% input array of p value
% output array of z score z_i = invert normal CDF of (1-p_i)
function basic_z_score = p2z(array_p_value);
for i = 1: length(array_p_value)
    basic_z_score(i) = norminv(1-array_p_value(i), 0 ,1); %mu = 0, sigma = 1
    if (basic_z_score(i) == inf)
        fprintf('pval of the %i is too small so the z-score = inf, set new p val 9.999e-16\n',i)
        array_p_value(i) = 9.9999999999999999999999999e-16;
        basic_z_score(i) = norminv(1-array_p_value(i), 0 ,1); %mu = 0, sigma = 1
    end
end


