function  [orig_corr,pval,sig] = permutation_correlate_nonmatch(data1,data2,iter,corr_type)

%if only comparing with one dataset, just create empty vector data2

    size1 = size(data1,1); % data measure
    size2 = size(data2,1); % individual difference measure
     
    null_corr = zeros(size(data1,2),iter);
    
    orig_corr = corr(data1,data2,'type',corr_type,'rows','complete');
    rand_vec = randi(size2,size2,iter);
    for x = 1:iter
        %rand_vec = rand(size2,1);
        %[~,sort_rand] = sort(rand_vec);
        %data2_rand = data2(sort_rand); 
        %data2_rand = data2(rand_vec(:,x));
        null_corr(:,x) = corr(data1,data2(rand_vec(:,x)),'type',corr_type,'rows','complete');
    end
    
    pval_right = ((iter - sum(orig_corr(:)'>null_corr'))/iter)';
    pval_left = ((sum(orig_corr(:)'<null_corr'))/iter)';
    pval_combo = (horzcat(pval_left,pval_right));
    pval = min(pval_combo')';
    sig = double(pval<0.05 & pval~=0);
end