% generate simulate data
parameters = data.generate_fixed_parameters(7,6,270,40,3);
ang = 0:45:359;
data1 = data.generate_simulate_data(ang,parameters);
[data2,noise] = data.generate_noise(data1,5);

%plotting
figure(1),hold on
plot(ang,data2,'r*') 
plot(ang,data1,'k.')
ylim([0,15])
xlabel('theta')
ylabel('response magnitude')
title('simulate tuning curve (+50% noise)')