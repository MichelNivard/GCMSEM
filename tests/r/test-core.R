library(GCMSEM)
checks<-0L
check<-function(x,msg){checks<<-checks+1L;if(!isTRUE(x))stop('FAILED: ',msg);cat('ok',checks,msg,'\n')}
near<-function(x,y,tol=1e-7)max(abs(x-y))<tol
fails<-function(expr)inherits(try(force(expr),silent=TRUE),'try-error')
m<-gcmsem_model('y ~ b*x',c('x','y'),gaussian_observed=FALSE)
p<-m$start;p[]<-c(.8,.3,1.2,.5,.1,.7,.4)
v<-implied_moments(m,p)
check(near(v[c('x|x','x|y','y|y')],c(.8,.32,.628)),'analytic one-way covariance')
check(near(v['x|x|y'],.12),'analytic xxy cumulant')
check(near(v['x|y|y|y'],1.2*.4^3),'analytic xyyy cumulant')
J<-implied_moments(m,p,jacobian=TRUE)$jacobian
fd<-sapply(seq_along(p),function(i){u<-w<-p;u[i]<-u[i]+1e-5;w[i]<-w[i]-1e-5;(implied_moments(m,u)-implied_moments(m,w))/(2e-5)})
check(near(J,fd,1e-7),'all exact derivatives match finite differences')
check(source_diagnostics(m,p)$necessary_conditions_pass,'known source moments feasible')
mm<-gcmsem_model('f =~ a*x + a*y + z\ny ~ b*x\nz ~ c*y\nx ~~ d*z\nk3(x,x,y) = free("e")\nk4(y,z,z,z) = free("h")',c('x','y','z'),gaussian_observed=FALSE,gaussian_latent=FALSE)
check(length(moment_keys(mm$observed))==31,'3-trait full cumulant vector')
check(sum(names(mm$start)=='a')==1,'equality loading one parameter')
jj<-implied_moments(mm,jacobian=TRUE)$jacobian
ff<-sapply(seq_along(mm$start),function(i){u<-w<-mm$start;u[i]<-u[i]+1e-5;w[i]<-w[i]-1e-5;(implied_moments(mm,u)-implied_moments(mm,w))/(2e-5)})
check(near(jj,ff,1e-7),'latent/mixed-source exact derivatives')
# Canonical tuple syntax makes repeated positions explicit, independently of labels.
a<-gcmsem_model('k3(x,x,y)=free("y"); k3(x,y,y)=free("x"); k4(y,x,y,x)=-0.2',c('x','y'))
ap<-a$start;ap['y']<-3;ap['x']<-2
av<-implied_moments(a,ap)
check(near(av[c('x|x|y','x|y|y','x|x|y|y')],c(3,2,-.2)),
      'explicit repeated indices remain unambiguous when labels equal variable names')
anon<-gcmsem_model('k3(x,x,y); k3(x,y,y)=free(); k4(x,x,y,y)=free()',c('x','y'))
check(all(c('k3(x,x,y)','k3(x,y,y)','k4(x,x,y,y)')%in%names(anon$start)),
      'standalone and free() create anonymous source cumulants')
eq<-gcmsem_model('k3(x,x,y)=free("s"); k3(x,y,y)=free("s")',c('x','y'))
eqp<-eq$start;eqp['s']<-.75
eqv<-implied_moments(eq,eqp)
check(sum(names(eq$start)=='s')==1 && near(eqv[c('x|x|y','x|y|y')],rep(.75,2)),
      'quoted repeated labels impose exact equality')
single_quote<-gcmsem_model("k3(x,x,y)=free('s')",c('x','y'))
check('s'%in%names(single_quote$start),'single-quoted labels supported without evaluation')
fixed_cumulants<-gcmsem_model('k3(x,x,y)=0; k3(x,y,y)=-1.25e-2; k4(x,x,y,y)=+2e-1',c('x','y'))
check(near(implied_moments(fixed_cumulants)[c('x|x|y','x|y|y','x|x|y|y')],c(0,-.0125,.2)),
      'finite signed scientific numeric literals fix exact source cumulants')
for(bad_syntax in c('k3(x,y)', 'k3(x,x,y,y)', 'k4(x,x,y)', 'k3(x,x,y,)',
                     'k3(x,,y)', 'k3(x,x+y,y)', 'k3(x,x,y)=s',
                     'k3(x,x,y)=free(s)', 'k3(x,x,y)=1/2',
                     'k3(x,x,y)=Inf', 'k3(x,x,y)=1e999',
                     'k3(x,x,y)=free("s","t")','k3(x,x,y)=free("bad label")',
                     'k3(x,x,y)=free("s")oops','k3(x,x,y)=',
                     'k3(x,x,unknown)')) {
  check(fails(gcmsem_model(bad_syntax,c('x','y'))),paste('reject malformed cumulant:',bad_syntax))
}
check(fails(gcmsem_model('k3(x,x,y); k3(y,x,x)',c('x','y'))),
      'duplicate symmetric explicit cumulants rejected')
.gcmsem_parser_executed<-FALSE
injection<-paste0('k3(x,x,y)=free(assign(".gcmsem_parser_executed",TRUE,envir=.GlobalEnv))')
check(fails(gcmsem_model(injection,c('x','y'))) && !.gcmsem_parser_executed,
      'model text never executes expressions')
rm(.gcmsem_parser_executed)
for(old in c('x ~-- s*y','x ~~- s*y')) {
  err<-try(gcmsem_model(old,c('x','y')),silent=TRUE)
  expected<-if(grepl('~--',old,fixed=TRUE))'k3(x, x, y) = free("s")' else 'k3(x, y, y) = free("s")'
  check(inherits(err,'try-error') && grepl(expected,as.character(err),fixed=TRUE),
        paste('retired shorthand gives exact migration:',old))
}
legacy_messages<-character()
legacy<-withCallingHandlers(gcmsem_model('x ~~~ s*y*y; x ~~~~ t*x*y*y',c('x','y')),
  warning=function(w){legacy_messages<<-c(legacy_messages,conditionMessage(w));invokeRestart('muffleWarning')})
canonical<-gcmsem_model('k3(x,y,y)=free("s"); k4(x,x,y,y)=free("t")',c('x','y'))
lp<-legacy$start;lp[c('s','t')]<-c(.8,1.2)
check(near(implied_moments(legacy,lp),implied_moments(canonical,lp)),
      'deprecated repeated-tilde models preserve numerical meaning')
check(length(legacy_messages)==2 &&
      any(grepl('k3(x, y, y) = free("s")',legacy_messages,fixed=TRUE)) &&
      any(grepl('k4(x, x, y, y) = free("t")',legacy_messages,fixed=TRUE)),
      'legacy warnings spell out exact canonical expansions')
check(fails(gcmsem_model('x ~~ z',c('x','y'))),'unknown variable rejected')
check(fails(gcmsem_model('x ~~ y\ny ~~ x',c('x','y'))),'duplicate symmetric terms rejected')
check(fails(gcmsem_model('x ~ x',c('x','y'))),'self regression rejected')
check(fails(implied_moments(m,keys='y|x')),'noncanonical key rejected')
V<-diag(.001,length(v));dimnames(V)<-list(names(v),names(v))
f<-gcmsem(m,v,V,n_starts=2,verbose=FALSE)
check(near(coef(f),p,2e-4),'fit known moments/path')
check(f$diagnostics$inference_valid,'interior identified inference valid')
check(all(is.finite(parameter_estimates(f)$se)),'sandwich SEs finite')
qs<-lengths(strsplit(names(v),'|',fixed=TRUE));unit<-1e-3
vv<-v*unit^qs;VV<-V*outer(unit^qs,unit^qs)
fs<-gcmsem(m,vv,VV,n_starts=2,verbose=FALSE)
check(near(coef(fs)['b'],coef(f)['b'],2e-5),'tiny-unit path invariance')
check(near(vcov(fs)['b','b'],vcov(f)['b','b'],1e-6),'tiny-unit SE invariance')
check(fs$diagnostics$inference_valid,'tiny-unit inference valid')
fg<-gcmsem_model('y ~ b*x',c('x','y'))
f2<-gcmsem(fg,v,V,orders=2,n_starts=1,verbose=FALSE)
check(length(fitted(f2))==3 && near(coef(f2)['b'],.4,1e-5),'observed order subsetting of full V')
W<-V;W[1,2]<-W[2,1]<-.0003
fw<-gcmsem(m,v,W,n_starts=1,verbose=FALSE)
check(near(coef(fw),coef(f),1e-5),'DWLS ignores offdiagonal V for estimates')
check(max(abs(vcov(fw)-vcov(f)))>1e-7,'sandwich uses offdiagonal V')
Vbad<-V;Vbad[1,2]<-Vbad[2,1]<-.1
check(fails(gcmsem(m,v,Vbad,verbose=FALSE)),'indefinite V rejected')
u<-gcmsem_model('y ~ b*x\nx ~ c*y',c('x','y'))
fu<-gcmsem(u,v,V,orders=2,n_starts=1,verbose=FALSE)
check(!fu$diagnostics$inference_valid && all(is.na(parameter_estimates(fu)$se)),'underidentified SEs withheld')
bad<-gcmsem_model('k4(x,x,y,y)=20',c('x','y'))
check(!source_diagnostics(bad)$necessary_conditions_pass,'impossible cross-kurtosis rejected')
fixed<-gcmsem_model('x ~~ 1*x','x')
fv<-implied_moments(fixed);vf<-diag(.01,3);dimnames(vf)<-list(names(fv),names(fv))
fixedfit<-gcmsem(fixed,fv,vf,n_starts=1,verbose=FALSE)
check(fixedfit$objective==0 && nrow(parameter_estimates(fixedfit))==0,'fixed model supported')
cat('PASS:',checks,'checks\n')
# Compiled tensor contractions use the independently implemented R kernel as oracle.
ri<-implied_moments(mm,jacobian=TRUE,backend='R')
ni<-implied_moments(mm,jacobian=TRUE,backend='native')
check(near(ri$estimate,ni$estimate,1e-12)&&near(ri$jacobian,ni$jacobian,1e-12),'compiled kernel agrees with R oracle')
# Deterministic fixed-zero sources are harmless when every incident cumulant is zero.
zero<-gcmsem_model('f =~ 0*x; f ~~ 0*f', 'x')
check(source_diagnostics(zero)$necessary_conditions_pass && source_diagnostics(zero)$covariance_min_eigen>0,'fixed-zero source removed from boundary diagnostics')
# Independent finite differences for the exact source-PSD penalty gradient.
pm<-gcmsem_model('x ~~ c*y; k3(x,x,y)=free("s"); k3(x,y,y)=free("t"); k4(x,x,y,y)=free("h")',c('x','y'),gaussian_observed=FALSE)
pp<-pm$start;pp[c('h','c','s','t')]<-c(6,.2,.3,.4)
ps<-GCMSEM:::.constraint_setup(pm,pm$start)
pa<-GCMSEM:::.source_penalty(pm,pp,ps)
pfd<-sapply(seq_along(pp),function(i){a<-b<-pp;a[i]<-a[i]+1e-5;b[i]<-b[i]-1e-5;(GCMSEM:::.source_penalty(pm,a,ps)$value-GCMSEM:::.source_penalty(pm,b,ps)$value)/(2e-5)})
check(max(abs(pa$gradient-pfd)/(1+abs(pfd)))<1e-7,'source PSD penalty exact derivative')
repeatm<-gcmsem_model('y ~ 0*x',c('x','y'),gaussian_observed=FALSE)
repeatp<-repeatm$start;repeatp[c('k4(x,x,x,x)','k4(y,y,y,y)')]<--4
repeats<-GCMSEM:::.constraint_setup(repeatm,repeatm$start)
repeata<-GCMSEM:::.source_penalty(repeatm,repeatp,repeats)
repeatfd<-sapply(seq_along(repeatp),function(i){a<-b<-repeatp;a[i]<-a[i]+1e-5;b[i]<-b[i]-1e-5;(GCMSEM:::.source_penalty(repeatm,a,repeats)$value-GCMSEM:::.source_penalty(repeatm,b,repeats)$value)/(2e-5)})
check(max(abs(repeata$gradient-repeatfd)/(1+abs(repeatfd)))<1e-7,'PSD derivative valid at repeated eigenvalues')
# Inconsistent empirical target is projected to a necessary-feasible boundary.
cm<-gcmsem_model('x ~~ v*x','x',gaussian_observed=FALSE)
ct<-c('x|x'=1,'x|x|x'=0,'x|x|x|x'=-4)
cv<-diag(.1,3);dimnames(cv)<-list(names(ct),names(ct))
cf<-gcmsem(cm,ct,cv,n_starts=1,verbose=FALSE,feasibility='penalty')
check(cf$diagnostics$source$necessary_conditions_pass && !cf$diagnostics$inference_valid,'optional penalized boundary fit flagged for inference')
check(cf$objective>0 && cf$constraint_penalty<1e-12,'raw DWLS Q separated from constraint penalty')
cu<-1e-3^(2:4);cfu<-gcmsem(cm,ct*cu,cv*outer(cu,cu),n_starts=1,verbose=FALSE,feasibility='penalty')
check(abs(cf$objective-cfu$objective)<1e-4,'PSD constrained fit unit invariance')
checkfit<-suppressWarnings(gcmsem(cm,ct,cv,n_starts=1,verbose=FALSE))
check(!checkfit$diagnostics$source$necessary_conditions_pass && checkfit$objective<1e-8,'default unconstrained higher-order Heywood estimate retained')
cat('PASS total:',checks,'checks\n')
# Python bundle interoperability and mandatory metadata checks.
tmp<-tempfile('gcmsem-bundle-')
write.table(data.frame(key=names(v),order=qs,estimate=unname(v),se=sqrt(diag(V))),paste0(tmp,'.moments.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(key=names(v),V,check.names=FALSE),paste0(tmp,'.sampling_cov.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('{"schema_version":"1.0","moment_type":"cumulant","traits":["x","y"]}',file=paste0(tmp,'.metadata.json'))
bundle<-read_moments(tmp)
check(near(bundle$moments$estimate,v)&&near(bundle$V,V),'versioned Python bundle read')
cat('{"moment_type":"cumulant"}',file=paste0(tmp,'.metadata.json'))
check(fails(read_moments(tmp)),'missing schema version rejected')
cat('{"schema_version":"1.0"}',file=paste0(tmp,'.metadata.json'))
check(fails(read_moments(tmp)),'missing cumulant convention rejected')
unlink(paste0(tmp,c('.moments.tsv','.sampling_cov.tsv','.metadata.json')))
cat('PASS all:',checks,'checks\n')

# User-requested unconstrained Heywood handling: formal estimates are not censored.
check(checkfit$diagnostics$formal_se_available && all(is.finite(parameter_estimates(checkfit)$se)),
      'higher-order Heywood retains computable formal sandwich SE')
check(all(parameter_estimates(checkfit)$se_status=='formal_only_flagged'),
      'formal-only uncertainty explicitly labelled')
nt<-ct;nt['x|x']<- -1
messages<-character()
nf<-withCallingHandlers(gcmsem(cm,nt,cv,n_starts=3,verbose=FALSE),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
check(coef(nf)['v']<0 && abs(coef(nf)['v']+1)<1e-6,'negative source variance retained without clipping')
check(nf$diagnostics$formal_se_available && !nf$diagnostics$source_admissible,
      'classical Heywood retains formal SE and distinct admissibility flag')
check(nf$diagnostics$source$covariance_min_eigen<0,
      'negative variance remains in covariance eigenvalue diagnostic')
check(any(grepl('CLASSICAL HEYWOOD CASE',messages,fixed=TRUE)),'classical Heywood emits prominent warning')
messages<-character()
hf<-withCallingHandlers(gcmsem(cm,ct,cv,n_starts=3,verbose=FALSE),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
check(any(grepl('HIGHER-ORDER HEYWOOD CASE',messages,fixed=TRUE)),'higher-order Heywood emits prominent warning')
check(hf$objective==min(hf$starts$objective[hf$starts$convergence==0]),
      'unconstrained selection uses minimum converged raw Q')
check(hf$objective<cf$objective && !hf$diagnostics$source_admissible,
      'inadmissible lower-Q estimate preferred over admissible penalty fit by default')
nuf<-suppressWarnings(gcmsem(cm,nt*cu,cv*outer(cu,cu),n_starts=1,verbose=FALSE))
check(abs(coef(nuf)['v']/1e-6+1)<1e-5,'negative variance tiny-unit invariance')
cat('PASS final:',checks,'checks\n')
