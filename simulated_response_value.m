function input = simulated_response_value(theta)
% generate a set of simulated RSP data with SIG(offset) = 5
input.RP = 100;
input.RN = 40;
input.OP = 90;
input.ON = input.OP + 180;
input.SIG = 70;
input.SIG_offset = 5;
input.offset = input.RP/5;

if (theta >0 && theta <=360 && rem(theta,1)==0)
    input.data = input.RP*exp(-1/2.*(theta-input.OP).^2/input.SIG.^2)+input.RN*exp(-1/2.*(theta-input.ON).^2/input.SIG.^2)+input.offset+input.SIG_offset*randn(1);
else
    error = 'input error. This is not a valid value of theta.'
end
