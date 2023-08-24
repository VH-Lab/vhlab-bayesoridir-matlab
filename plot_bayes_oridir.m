function plot_bayes_oridir(output_struct)

figure(),plot(output_struct.marginal_likelihood.Rp.values,output_struct.marginal_likelihood.Rp.likelihoods),
xlabel("Rpref(Hz)"),
ylabel("probability of Rpref");
1
figure(),plot(output_struct.marginal_likelihood.Rn.values,output_struct.marginal_likelihood.Rn.likelihoods),
xlabel("alpha"),
ylabel("probability of alpha");
2
figure(),plot(output_struct.marginal_likelihood.theta_pref.values,output_struct.marginal_likelihood.theta_pref.likelihoods),
xlabel("θpref"),
ylabel("probability of θpref");
3
figure(),plot(output_struct.marginal_likelihood.sigma.values,output_struct.marginal_likelihood.sigma.likelihoods),
xlabel("σ"),
ylabel("probability of σ");
4
figure(),plot(output_struct.marginal_likelihood.Rsp.values,output_struct.marginal_likelihood.Rsp.likelihoods),
xlabel("Rsp"),
ylabel("probability of Rsp");
5
figure(),plot(output_struct.descriptors.oi.values,output_struct.descriptors.oi.likelihoods),
xlabel("Orientation Index"),
ylabel("probability of Orientation Index");
end