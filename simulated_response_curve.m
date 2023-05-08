function input = simulated_response_curve(interval)

% generate a set of simulated RSP data with SIG(offset) = 5
input.RP = 100;
input.RN = 40;
input.OP = 90;
input.ON = input.OP + 180;
input.SIG = 70;
input.SIG_offset = 5;
input.offset = input.RP/5;
theta = 0:interval:359;
if (interval >=0 && interval <=90 && rem(interval,1)==0)
    input.data = input.RP*exp(-0.5.*(theta-input.OP).^2/input.SIG.^2)+input.RN*exp(-0.5.*(theta-input.ON).^2/input.SIG.^2)+input.offset+input.SIG_offset*randn(1,360/interval);
else
    error = 'input error. This is not a valid value of theta.'
end