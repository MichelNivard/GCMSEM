function run_sld4m(input_file, upstream_dir, output_file, input_mode)
% Execute the original SLD4M.m without changing its regressions or jackknife.
% input_file is produced by prepare_polygenicity.py from the released scores.
if nargin < 4, input_mode = 'beta'; end
if ~strcmp(input_mode,'beta') && ~strcmp(input_mode,'zscore'), error('input_mode must be beta or zscore'); end
addpath(upstream_dir);
x = load(input_file);
% Dense annotation input also works with Octave broadcasting.
x.annot = full(x.annot);
p = size(x.ell2,2); nt = size(x.beta,2);
ma = zeros(p,nt); ma_se = ma; genetic_variance = ma; genetic_variance_se = ma;
ma_jk = zeros(p,x.no_blocks,nt); variance_jk = ma_jk;
for t=1:nt
  fprintf('SLD4M trait %d/%d; tags=%d; annotations=%d\n',t,nt,size(x.beta,1),p);
  % Per-genotype-standardized squared effects are an upstream-supported input.
  % Both regressions use the same fixed reference weights as GCMSEM.
  if strcmp(input_mode,'zscore')
    chisq = (x.beta(:,t)./x.se(:,t)).^2;
  else
    chisq = x.beta(:,t).^2;
  end
  [a,ae,h,he,~,~,~,~,aj,hj] = SLD4M(chisq,x.ell2,x.ell4,...
    x.idxldsc,x.annot,x.weights,x.weights,x.no_blocks,1,eye(p));
  ma(:,t)=a;ma_se(:,t)=ae;genetic_variance(:,t)=h;genetic_variance_se(:,t)=he;
  ma_jk(:,:,t)=aj;variance_jk(:,:,t)=hj;
end
save(output_file,'ma','ma_se','genetic_variance','genetic_variance_se',...
  'ma_jk','variance_jk','input_mode','-v7');
fprintf('Saved %s\n',output_file);
end
