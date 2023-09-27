% generate simulate data
parameters = data.generate_fixed_parameters(7,6,270,40,3);
ang = 0:45:359;
data1 = data.generate_simulate_data(ang,parameters);
[data2,noise] = data.generate_noise(data1,5);
figure(1),hold on
plot(ang,data2(:,1),'r')
plot(ang,data2(:,2),'g')
plot(ang,data2(:,3),'b')
plot(ang,data2(:,4),'y')
plot(ang,data2(:,5),'c')
parameters = simulation.generate_fixed_parameters(7,6,270,40,3);
ang = 0:45:359;
data1 = simulation.generate_simulate_data(ang,parameters);
[data2,noise] = simulation.generate_noise(data1,1);
figure(1),hold on
plot(ang,data2,'r')
plot(ang,data1,'k')