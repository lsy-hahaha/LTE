function [ pathloss_enb_to_ue ] = LTE_pathloss_enb_to_ue( ditance_enb_to_ue )
%LTE_PATHLOSS_ENB_TO_UE Summary of this function goes here
%   Detailed explanation goes here
distance_enb_to_ueinkm= ditance_enb_to_ue/1000;
Prob=min(0.018./distance_enb_to_ueinkm,1).*(1-exp(-distance_enb_to_ueinkm/0.063))+exp(-distance_enb_to_ueinkm/0.063);
PLlos=103.4+24.2*log10(distance_enb_to_ueinkm);
PLnlos=131.1+42.8*log10(distance_enb_to_ueinkm);
pathloss_enb_to_ue=Prob.*PLlos+(1-Prob).*PLnlos;
end

