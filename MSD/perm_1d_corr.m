function pval = perm_1d_corr(data1,data2,iter)
  size1 = size(data1,1);
  size2 = size(data2,1);
  if size1 ~= size2
    sprintf('%s','size mismatch')
  end
  null_corr = zeros(iter,1);
  orig_corr = corr(data1,data2,'rows','pairwise');
  for x = 1:iter
    rand_vec = rand(size1,1);
    [~,sort_rand] = sort(rand_vec);
    data1_rand = data1(sort_rand); 
    null_corr(x,1) = corr(data1_rand,data2,'rows','pairwise');  
  end
%   pos_thr = prctile(null_corr,95);
%   neg_thr = prctile(null_corr,5);
%   sig = double(orig_corr > pos_thr) + double(orig_corr < neg_thr);
  pval1 = (iter-sum(orig_corr>null_corr))/iter;
  pval2 = sum(orig_corr>null_corr)/iter;
  temp = horzcat(pval1,pval2);
  pval = min(temp);
end