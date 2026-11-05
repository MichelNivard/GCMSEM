function check_sld4m_storage(upstream_dir, work_dir)
% Check every native returned estimate and deletion against ordinary arrays.
  addpath(upstream_dir); rng(8451);
  n = 2000; p = 6; B = 20;
  annot = [ones(n,1), double(rand(n,p-1) > .6)];
  ell2 = 1 + rand(n,p)*5; ell4 = .5 + rand(n,p)*3;
  chisq = (randn(n,1) .* sqrt(1 + ell2*(.1:.1:.6)')).^2;
  idx = (1:n)'; idx(7:37:end) = []; chisq = chisq(idx);
  w2 = 1./ell2(idx,1); w4 = 1./ell4(idx,1);
  disk = cell(1,3); arrays = {ell2,ell4,annot};
  for k=1:3
    path = fullfile(work_dir, sprintf('storage-test-%d.f64',k));
    fid = fopen(path,'wb','ieee-le'); fwrite(fid,arrays{k}','double'); fclose(fid);
    % scipy.io.savemat stores dimension metadata as int64. The adapter must
    % still return double sums: integer arithmetic would round native ratios.
    disk{k} = SLD4MDiskMatrix(path,int64(n),int64(p),sum(arrays{k}));
    assert(isa(sum(disk{k}),'double'));
  end
  expected = cell(1,14); observed = cell(1,14);
  [expected{:}] = SLD4M(chisq,ell2,ell4,idx,annot,w2,w4,B);
  [observed{:}] = SLD4M(chisq,disk{1},disk{2},idx,disk{3},w2,w4,B);
  max_difference = 0;
  for k=1:14
    assert(isequal(size(expected{k}),size(observed{k})));
    assert(isequal(isnan(expected{k}),isnan(observed{k})));
    difference = max(abs(expected{k}(:)-observed{k}(:)));
    max_difference = max(max_difference,difference);
    assert(difference <= 1e-11*max(1,max(abs(expected{k}(:)))));
  end
  fprintf('All 14 native outputs agree; maximum absolute difference %.17g\n',max_difference);
end
