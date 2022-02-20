# Code to reproduce the analyses shown in Figures 1B and 1C of Rajkumar et al. 2022
# Author: Mathieu Lajoie (mathieu.lajoie2@mcgill.ca)
# February 18, 2022

# Moving average used to estimate P(mutation | TMB) i.e. P(y == 1 | x) 
# y must be an integer binary vector giving the mutation status in each sample.
# x must be a numeric vector giving the TMB of each sample.
mavg = function(y, x, k = 201, pc = 1, pad = TRUE){
  
  checkmate::assertInteger(y, any.missing = FALSE, lower = 0,upper = 1)
  checkmate::assertNumeric(x, finite = TRUE, any.missing = FALSE, len = length(x))
  
  # window size must be odd
  checkmate::assertTRUE(k %% 2 == 1)
  
  n = length(x)
  o = order(x)
  res = stats::filter(y[o], filter = rep(1/k,k), method = "conv", sides = 2) + pc/k
  
  if(pad){
    pad.size = floor(k/2)
    res[1:pad.size] = res[pad.size+1]
    res[(n-pad.size+1):n] = res[n-pad.size]
  }
  
  res = res[order(o)]
  res = res * sum(y,na.rm = TRUE)/sum(res,na.rm = TRUE) # re-scale so that expected nb mutations equals observed nb.
  
  if(!all.equal(sum(res),sum(y))){
    print(paste("Nb mut",sum(res)))
    print(paste("Sum prob",sum(y)))
    warning("expected nb events don't match density in mavg()")
  }
  res
}

# Simulate co-occurring mutations while taking TMB into account.
# Returns a data.frame with expected number of co-occ. mutation, a two-tailed p-value and 95% confidence interval. 
# An extra row with observed number is added for subsequent plotting purpose.
# m1 and m2 must be integer binary vectors giving the mutation status in each sample for gene 1 and gene 2.
# p1 and p2 must be numeric vectors given the mutation probability in each sample (according to their TMB) for gene 1 and gene 2.
mc_sim_cooc = function(m1, m2, p1, p2, nb_sim = 1e5, return.all.data = FALSE, show.plot = TRUE){
  
  checkmate::assertInteger(m1, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assertInteger(m2, lower = 0, upper = 1, any.missing = FALSE, len = length(m1))
  checkmate::assertNumeric(p1, lower = 0, upper = 1, any.missing = FALSE, len = length(m1))
  checkmate::assertNumeric(p2, lower = 0, upper = 1, any.missing = FALSE, len = length(m1))
  
  nb1 = sum(m1);nb1
  nb2 = sum(m2);nb2
  nb.samp = length(m1)
  
  checkmate::assertTRUE(all.equal(sum(p1), nb1)) # sanity check
  checkmate::assertTRUE(all.equal(sum(p2), nb2))
  
  nb.inter = sum(m1 & m2); nb.inter
  idx = 1:length(m1)
  
  # we must use *rates* instead of probs for sampling
  r1 = -log(1-p1)
  r2 = -log(1-p2) 
  
  res = replicate(nb_sim, expr = {
    s1 = sample(idx, size = nb1, replace = FALSE, prob = r1)
    s2 = sample(idx, size = nb2, replace = FALSE, prob = r2)
    sum(s1 %in% s2)
  })
  
  if(show.plot){
    hist(res, breaks = 1e3, xlim=c(0,max(c(res,nb.inter+1))), main = "", xlab = "nb cooc")
    abline(v = nb.inter,col="red")
  }
  
  # p-value formula from:
  # Davison, Anthony Christopher, and David Victor Hinkley. Bootstrap methods and their application. No. 1. Cambridge university press, 1997.
  pval.low = (sum(res <= nb.inter)+1)/(nb_sim+1)  
  pval.high = (sum(res >= nb.inter)+1)/(nb_sim+1)
  sres = sort(res)
  tprob = table(factor(res,levels = 0:max(c(nb.inter,res))))/nb_sim
  psim = tprob[as.character(nb.inter)] 
  
  sim.cooc = mean(res)
  sim.ci95.low = sres[nb_sim*0.025]
  sim.ci95.high = sres[nb_sim*0.975]
  
  pval.low = c(NA, pval.low)
  pval.high = c(NA, pval.high)
  pval.two.tailed = pmin(1, pmin(pval.low, pval.high)*2 )
  
  dfres = data.frame(test = c("obs","sim"),
                     nb.cooc = c(nb.inter, sim.cooc),
                     ci95.low = c(NA, sim.ci95.low),
                     ci95.high = c(NA, sim.ci95.high),
                     pval.low,
                     pval.high,
                     pval.two.tailed)
  
  if(return.all.data){
    return(list(df = dfres, sim = res, density = tprob))
  }else{
    return(dfres)
  }
  
}
