function [ prb_snr ] = calculate_PRB_SNR( sc_snr,num_prb)

nUE=length(sc_snr);
for ue_i=1:nUE
    for prb_i=1:num_prb
        prb_snr(prb_i,ue_i)=mean(sc_snr{1,ue_i}(1+(prb_i-1)*12:prb_i*12));
    end
end


end

