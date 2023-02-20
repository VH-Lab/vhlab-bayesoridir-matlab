function input = response_value_data(theta)
% generate a set of simulated RSP data with SIG(offset) = 5
RP = 100;
RN = 40;
OP = 90;
ON = OP + 180;
SIG = 70;
offset = RP/5;
SIG_offset = 5;

if (theta >=0 && theta <360 && rem(theta,1)==0)
    input(theta).RSP = RP*exp(-1/2.*(theta-OP).^2/SIG.^2)+RN*exp(-1/2.*(theta-ON).^2/SIG.^2)+offset+SIG_offset*randn(1);

else
    error = 'input error. This is not a valid value of theta'
end
