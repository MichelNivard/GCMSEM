# Independent implementation of a sparse linear cumulant SEM.
.name_ok <- function(x) grepl("^[A-Za-z][A-Za-z0-9_.]*$", x)
.sorted_tuple <- function(x) paste(sort(x), collapse = ",")

.numeric_literal <- function(x) grepl("^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$", x)

# This parser deliberately never uses parse(), eval(), or expression substitution.
# Cumulant arguments are variable identifiers, and labels are quoted identifiers.
.parse_cumulant_statement <- function(line) {
  m <- regexec("^k([34])\\s*\\(([^()]*)\\)\\s*(?:=\\s*(.+))?$", line, perl=TRUE)
  bits <- regmatches(line,m)[[1L]]
  if(!length(bits))stop("Malformed cumulant statement. Use k3(x, x, y) = free(\"s\") or k4(x, x, y, y) = 0: ",line)
  q <- as.integer(bits[2L]); inside <- bits[3L]
  name_pattern <- "[A-Za-z][A-Za-z0-9_.]*"
  argument_pattern <- paste0("^\\s*",name_pattern,"\\s*(,\\s*",name_pattern,"\\s*){",q-1L,"}$")
  if(!grepl(argument_pattern,inside,perl=TRUE))
    stop("k",q," requires exactly ",q," plain variable identifiers separated by commas: ",line)
  tuple <- trimws(strsplit(inside,",",fixed=TRUE)[[1L]])
  rhs <- if(length(bits)>=4L)trimws(bits[4L]) else ""
  modifier <- NULL
  if(nzchar(rhs) && !grepl("^free\\s*\\(\\s*\\)$",rhs,perl=TRUE)) {
    if(.numeric_literal(rhs)) {
      if(!is.finite(as.numeric(rhs)))stop("Fixed cumulant must be a finite numeric literal: ",rhs)
      modifier <- rhs
    } else {
      lm <- regexec("^free\\s*\\(\\s*([\"'])([A-Za-z][A-Za-z0-9_.]*)\\1\\s*\\)$",rhs,perl=TRUE)
      lb <- regmatches(rhs,lm)[[1L]]
      if(!length(lb))stop("Cumulant RHS must be free(), free(\"label\"), or a finite numeric literal; expressions and unquoted labels are not accepted: ",rhs)
      modifier <- lb[3L]
    }
  }
  list(lhs=tuple[1L],op=paste0("k",q),rhs=rhs,tuple=tuple,modifier=modifier,line=line)
}

.canonical_cumulant_text <- function(tuple,modifier=NULL) {
  rhs <- if(is.null(modifier))"free()" else if(.numeric_literal(modifier))modifier else paste0('free("',modifier,'")')
  paste0("k",length(tuple),"(",paste(tuple,collapse=", "),") = ",rhs)
}

#' Define a genomic comoment structural equation model
#' @export
gcmsem_model <- function(syntax, observed, latent = character(),
                         gaussian_observed = TRUE, gaussian_latent = TRUE) {
  if (!is.character(observed) || !length(observed) || anyDuplicated(observed) ||
      any(!.name_ok(observed))) stop("observed must contain unique simple variable names")
  if (!is.character(syntax) || length(syntax) != 1L) stop("syntax must be one character string")
  lines <- trimws(unlist(strsplit(gsub("#[^\n]*", "", syntax), "[;\n]")))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) stop("Empty model syntax")
  parsed <- lapply(lines, function(line) {
    if(grepl("^k[34]\\s*\\(",line,perl=TRUE))return(.parse_cumulant_statement(line))
    m <- regexpr("=~|~~~~|~~~|~--|~~-|~~|~", line, perl = TRUE)
    if (m[1] < 1L) stop("Unsupported model statement: ", line)
    op <- regmatches(line, m)
    lhs <- trimws(substr(line, 1L, m[1] - 1L))
    rhs <- trimws(substr(line, m[1] + attr(m, "match.length"), nchar(line)))
    if (!.name_ok(lhs) || !nzchar(rhs)) stop("Malformed model statement: ", line)
    list(lhs=lhs, op=op, rhs=rhs, line=line)
  })
  latent <- unique(c(latent, vapply(Filter(function(z) z$op == "=~", parsed), `[[`, "", "lhs")))
  if (any(!.name_ok(latent)) || any(latent %in% observed)) stop("Observed and latent names must be valid and disjoint")
  nodes <- c(observed, latent); n <- length(nodes)
  rows <- list()
  put <- function(kind, tuple, modifier = NULL, default_start = 0, automatic = FALSE) {
    if (any(!tuple %in% nodes)) stop("Undeclared variable(s): ", paste(setdiff(tuple,nodes), collapse=", "), ". Declare observed or latent variables.")
    ii <- match(tuple, nodes)
    if (kind == "cumulant") ii <- sort(ii)
    id <- paste(kind, .sorted_tuple(if(kind=="path") seq_along(ii) else ii), collapse=":")
    if(kind=="path") id <- paste("path",paste(ii,collapse=","),sep=":")
    fixed <- !is.null(modifier) && .numeric_literal(modifier)
    if (!is.null(modifier) && !fixed && !.name_ok(modifier)) stop("Invalid parameter modifier: ",modifier)
    label <- if (is.null(modifier) || fixed) "" else modifier
    value <- if(fixed) as.numeric(modifier) else default_start
    if(!is.finite(value)) stop("Non-finite parameter")
    old <- which(vapply(rows,function(r)r$id==id,logical(1)))
    if(length(old) && !rows[[old]]$automatic) stop("Duplicate model term: ",paste(tuple,collapse=","))
    row <- list(id=id,kind=kind,tuple=ii,order=if(kind=="path")0L else length(ii),
                fixed=fixed,value=value,label=label,automatic=automatic)
    if(length(old)) rows[[old]] <<- row else rows[[length(rows)+1L]] <<- row
  }
  gauss_nodes <- function(arg, pool, nm) {
    if(is.logical(arg) && length(arg)==1L && !is.na(arg)) return(if(arg)pool else character())
    if(is.character(arg) && all(arg %in% pool)) return(arg)
    stop(nm," must be TRUE, FALSE, or variable names within its node class")
  }
  gaussian <- c(gauss_nodes(gaussian_observed,observed,"gaussian_observed"),
                gauss_nodes(gaussian_latent,latent,"gaussian_latent"))
  for(v in nodes) {
    put("cumulant",c(v,v),if(v %in% latent)"1" else NULL,1,TRUE)
    for(q in 3:4) put("cumulant",rep(v,q),if(v %in% gaussian)"0" else NULL,
                      if(q==3)0 else 0.5,TRUE)
  }
  for(z in parsed) {
    if(z$op %in% c("k3","k4")) {
      start<-if(length(z$tuple)==4L && length(unique(z$tuple))==1L)0.5 else 0
      put("cumulant",z$tuple,z$modifier,start)
      next
    }
    terms <- trimws(unlist(strsplit(z$rhs,"(?<![eE])\\+",perl=TRUE)))
    if(any(!nzchar(terms))) stop("Empty RHS term in ",z$line)
    for(term in terms) {
      bits <- trimws(strsplit(term,"\\*")[[1]])
      if(z$op %in% c("~","=~")) {
        if(!length(bits)%in%c(1L,2L)) stop("Regression/loading requires [modifier*]variable: ",term)
        mod <- if(length(bits)==2L)bits[1] else NULL
        var <- tail(bits,1L)
        if(var=="1")stop("Mean/intercept models are not fitted by GCMSEM 0.1; choose orders 2:4")
        if(z$lhs==var)stop("Self regression/loading is not allowed")
        put("path",if(z$op=="~")c(z$lhs,var) else c(var,z$lhs),mod,0.2)
      } else {
        q <- switch(z$op,"~~"=2L,"~~~"=3L,"~~~~"=4L,"~--"=3L,"~~-"=3L)
        need <- if(z$op%in%c("~--","~~-"))1L else q-1L
        if(!length(bits)%in%c(need,need+1L))stop("Operator ",z$op," requires ",need," RHS variables and an optional modifier: ",term)
        mod <- if(length(bits)==need+1L)bits[1] else NULL
        vars <- tail(bits,need)
        tuple <- switch(z$op,"~--"=c(z$lhs,z$lhs,vars),"~~-"=c(z$lhs,vars,vars),c(z$lhs,vars))
        if(z$op%in%c("~--","~~-"))
          stop("Operator ",z$op," has been retired because repeated-index direction is ambiguous. Replace this term with ",.canonical_cumulant_text(tuple,mod),".")
        if(z$op%in%c("~~~","~~~~"))
          warning("Higher-cumulant operator ",z$op," is deprecated. This term means ",.canonical_cumulant_text(tuple,mod),"; use that explicit syntax instead.",call.=FALSE)
        start <- if(q==2L && length(unique(tuple))==1L)1 else if(q==4L && length(unique(tuple))==1L)0.5 else 0
        put("cumulant",tuple,mod,start)
      }
    }
  }
  # Shared labels represent one parameter, including across different orders.
  free_names <- character(); starts <- numeric(); lower <- numeric()
  for(j in seq_along(rows)) {
    r <- rows[[j]]
    if(r$fixed) {rows[[j]]$parameter <- 0L; next}
    nm <- if(nzchar(r$label))r$label else if(r$kind=="path")paste0(nodes[r$tuple[1]],"~",nodes[r$tuple[2]]) else paste0("k",r$order,"(",paste(nodes[r$tuple],collapse=","),")")
    id <- match(nm,free_names)
    lo <- if(r$kind=="cumulant" && r$order==2L && length(unique(r$tuple))==1L)0 else -Inf
    if(is.na(id)) {free_names<-c(free_names,nm);starts<-c(starts,r$value);lower<-c(lower,lo);id<-length(free_names)} else lower[id]<-max(lower[id],lo)
    rows[[j]]$parameter <- id
  }
  names(starts)<-names(lower)<-free_names
  model <- list(syntax=syntax, observed=observed,latent=latent,nodes=nodes,rows=rows,
                start=starts,lower=lower,upper=setNames(rep(Inf,length(starts)),free_names),
                gaussian_defaults=gaussian,version="0.1.0")
  class(model)<-"gcmsem_model"; model
}

#' @export
print.gcmsem_model <- function(x,...) {
  cat("GCMSEM model:",length(x$observed),"observed,",length(x$latent),"latent variables;",length(x$start),"free parameters\n")
  cat("Observed:",paste(x$observed,collapse=", "),"\n")
  cat("Higher source cumulants default to zero for:",paste(x$gaussian_defaults,collapse=", "),"\n")
  invisible(x)
}

.resolve_par <- function(model,par) {
  if(is.null(par))return(model$start)
  if(!is.numeric(par) || length(par)!=length(model$start) || any(!is.finite(par)))stop("par must provide all finite free parameters")
  if(!is.null(names(par))) {
    if(anyDuplicated(names(par)) || !setequal(names(par),names(model$start)))stop("Named par must match model parameter names")
    par<-par[names(model$start)]
  }
  setNames(as.numeric(par),names(model$start))
}
