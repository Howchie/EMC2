## Script to test limits of bivariate normals
devtools::load_all(); library(microbenchmark)
n=50
ks = seq(-3,3,length.out=n)
hs = seq(-3,3,length.out=n)
ps = seq(-1,1,length.out=n)

# Create function for test
pbvn = function(hs,ks,ps,algorithm){
  if(!all(sapply(list(hs,ks,ps), function(x) length(x) == length(hs))))
  {stop("length of inputs must be equal")}
  out = numeric(length(hs)); count=0
  for (k in ks) {
    for (h in hs) {
      for (p in ps) {
        count=count+1
        if (algorithm=="Tsay"){out[count]=pbvn_tsay(h,k,p)}
        else if (algorithm=="Drezner"){out[count]=pbvn_drezner(h,k,p)}
        else if (algorithm=="Genz"){out[count]=pbvn_tvpack(h,k,p)}
        else{stop("Invalid algorithm specified")}
      }
    }
  }
  return(out)
}

bench_pbvn = microbenchmark(tsay = pbvn(hs,ks,ps,"Tsay"),
                            drezner = pbvn(hs,ks,ps,"Drezner"),
                            genz = pbvn(hs,ks,ps,"Genz"),times=20)
