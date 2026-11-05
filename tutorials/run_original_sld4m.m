function run_original_sld4m(reference_dir, input_file, upstream_dir, output_file)
% Follow the released example.m using the author's complete reference inputs.
% Input SNP indices are matched to the original RSIDs in original row order.
% SLD4M.m is unmodified; SLD4MDiskMatrix changes storage, not the arithmetic.
  addpath(upstream_dir);
  r = load(fullfile(reference_dir, 'reference_metadata.mat'));
  x = load(input_file);
  ell2 = SLD4MDiskMatrix(fullfile(reference_dir, 'ell2.f64'), r.nrows, r.ncols, r.ell2_sums);
  ell4 = SLD4MDiskMatrix(fullfile(reference_dir, 'ell4.f64'), r.nrows, r.ncols, r.ell4_sums);
  annot = SLD4MDiskMatrix(fullfile(reference_dir, 'annot.f64'), r.nrows, r.ncols, r.annot_sums);
  assert(isa(sum(annot),'double'), 'Reference totals must use double precision');
  no_blocks = 100;
  fprintf('Original S-LD4M: %d matched SNPs; %d annotations; %d blocks\n', ...
      numel(x.chisq), r.ncols, no_blocks);
  tic;
  [ma,ma_se,h,h_se,ma_enrich,ma_enrich_se,h_enrich,h_enrich_se,ma_jk,h_jk] = ...
      SLD4M(x.chisq,ell2,ell4,x.idx_ref,annot,r.EW2(x.idx_ref),r.EW4(x.idx_ref),no_blocks);
  elapsed_seconds = toc;
  save(output_file,'ma','ma_se','h','h_se','ma_enrich','ma_enrich_se', ...
      'h_enrich','h_enrich_se','ma_jk','h_jk','no_blocks','elapsed_seconds','-v7');
  fprintf('Native base-category Ma = %.10g (SE %.10g); %.1f seconds\n',ma(1),ma_se(1),elapsed_seconds);
end
