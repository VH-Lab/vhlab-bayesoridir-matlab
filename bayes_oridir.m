function [out] = bayes_oridir(angles, responses, varargin)
% BAYES_ORIDIR - estimates orientation and direction fit parameters in a Bayesian manner
%
% OUT = BAYES_ORIDIR(ANGLES, RESPONSES, ...)
%
% Given a set of ANGLES (such as [0:30:360-30]) and a set of RESPONSES to each
% angle, where RESPONSES(i,j) is the ith response to a grating presented at
% ANGLES(j), return a Bayesian estimation of fit parameters to the function:
%
% R(THETA) = RSP+RP*exp(-vlt.math.angdiff(OP-ANGLES).^2/(2*SIG^2))+RN*exp(-vlt.math.angdiff(180+OP-ANGLES).^2/(2*SIG^2));
%
% OUT is an output structure with fields:
%  P_RP - posterior probability of RP
%  P_RN - posterior probability of RN
% % (fill in)
%  V_OP - values of OP that are examined (normally linspace(0,360,360) )
%  V_RP - values of RP that are examined (normally )
%  LIKELIHOOD - Joint posterior likelihood (RSP, RP, RN, OP, SIG)



  
  

  
p_rp = sum(sum(sum(lik,4),3),2);
p_rp = p_rp/sum(p_rp);
p_rn = sum(sum(sum(lik,4),3),1);
p_rn = p_rn/sum(p_rn);
p_ang = squeeze(sum(sum(sum(lik,4),2),1));
p_ang = p_ang/sum(p_ang);
p_sig = squeeze(sum(sum(sum(lik,3),2),1));
p_sig = p_sig/sum(p_sig);

out.p_rp = p_rp;
out.r_rn = p_rn;
 % . . 
   
