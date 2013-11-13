function metacommunity_input(in_param)
# stochastic competition model - n species with logistic growth in p patches

#turn off the broadcasting warning
warning("off", "Octave:broadcast");

npatches = 5;

nspecies = 5;
#all parameters are the same for all patches
r_mean = 1.5;
r_range = 0.8;

r = sort(unifrnd(r_mean - r_range/2, r_mean + r_range/2,1,nspecies)); #intrinsic growth rates; environmental stochasticity could be modelled by adding random noise to r values

intra = 0.01; #the strength of the intraspecific comeptition (the same for all species)
inter_intra = 1.0; #the average ratio of interspecific and intraspecific competition coefficient (this multiplies nondiagonal elements of the competition matrix)
a_var = 5.0; #interspecific variation in a(i,j) values
#increasing the a_var leads to increasing the mean a(i,j) value! This can confound the effects of increasing the a(i,j) variation.

growth_competition_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off

a = zeros(nspecies,nspecies);
for i = 1:nspecies
  a(i,i) = intra; #set intraspecific competition coefficients;remember that K(i)=r(i)/a(i,i)
    for j = 1:nspecies
      if i != j
	a(i,j) = a(i,i)*(r(i)/r(j))^a_var;
      endif
  endfor
endfor


if growth_competition_trade_off == 0
  b = reshape(a(randperm(numel(a))),nspecies,nspecies); #permute the a(i,j)_values randomly; need to restore the diagonal after this operation
  xx = [1:nspecies^2](b==a(1,1)); #indices of the elements that are equal to a(i,i)
  for ii = 1:length(xx)
    for jj = ([0:nspecies-1]*(nspecies+1)+1)
      if xx(ii) == jj
	xx(ii) = 999;
      endif
    endfor
  endfor
  xx(xx == 999) = [];
  for jj = ([0:nspecies-1]*(nspecies+1)+1)
    if b(jj) != a(1,1)
      b(xx(end)) = b(jj);
      xx(end) = [];
      b(jj) = a(1,1);
    endif
  endfor
  a = b;
  clear b;
endif

#now multiply non-diagonal elements by inter_intra
aa = a.*inter_intra;
for i = 1:nspecies
  aa(i,i) = intra; #reset intraspecific competition coefficients
endfor

K = r./(diag(a)');

m0_mean = 0.1; #mean species basal dispersal rate

#competition-colonization trade-off
m0_var = 5.0; #interspecific variation in m0 values; 1.0 = direct proportionality with r/r_mean when the trade-off is on
#I use r as an indirect measure of competition strength (species with high r are week competitors and hence strong dispersers)

competition_colonization_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off
m0 = m0_mean*((r/r_mean).^m0_var);

if competition_colonization_trade_off == 0
  m0 = m0(randperm(length(m0))); #permute the m0_values randomly
endif

#save all parameters
save(in_param, "npatches", "nspecies", "r_mean", "r_range", "r", "intra", "inter_intra", "a_var", "growth_competition_trade_off", "a", "K", "m0_mean", "m0_var", "competition_colonization_trade_off", "m0");

endfunction
