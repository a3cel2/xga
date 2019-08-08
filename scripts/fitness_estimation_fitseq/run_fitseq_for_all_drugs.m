drugs = ["CTRL","itraconazole","methotrexate","tamoxifen","beauvericin","benomyl","bisantrene","camptothecin","cisplatin","colchicine","cycloheximide","fluconazole","imatinib","ketoconazole","miconazole","mitoxantrone","valinomycin"];
mating_types = ["alpha","A"];

for i = 1:length(drugs)
    for j = 1:length(mating_types)
        drug = drugs(i);
        mating_type = mating_types(j);
        t_seq_vec = transpose(csvread(strcat("input_data/avail_timepoints_",drug,"_",mating_type,".csv")));
        BC_num_mat_original = csvread(strcat("input_data/fitseq_input_",drug,"_",mating_type,".csv")); 
        cell_depth = [];
        file_name = strcat('fitseq_output/',char(drug),'_',char(mating_type));
        [x_estimate_result, r_estimate_result, x_mean_estimate_result] =  FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name);
    end
end
