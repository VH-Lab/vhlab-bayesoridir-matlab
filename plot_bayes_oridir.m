function plot_bayes_oridir(bayes_structure);

plot(bayes_structure.vrp,bayes_structure.p_rp,'b');

plot(vrp-20,p_rp,'b'),
xlabel("Rpref(Hz)"),
ylabel("probability of Rpref"),

plot(vrn-20,p_rn,'b'),
xlabel("Rnull(Hz)"),
ylabel("probability of Rnull"),

plot(vang,p_ang,'b'),
xlabel("θpref"),
ylabel("probability of θpref"),

plot(vsig,p_sig,'b'),
xlabel("σ"),
ylabel("probability of Rpref"),