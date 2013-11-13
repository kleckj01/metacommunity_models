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

a_var = 5.0; #interspecific variation in a(i,j) values
#increasing the a_var leads to increasing the mean a(i,j) value! This can confound the effects of increasing the a(i,j) variation.

growth_competition_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off

a = zeros(nspecies,nspecies);
for i = 1:nspecies
  a(i,i) = 0.01; #set intraspecific competition coefficients;remember that K(i)=r(i)/a(i,i)
    for j = 1:nspecies
      if i != j
	a(i,j) = a(i,i)*(r(i)/r(j))^a_var;
      endif
  endfor
endfor


if competition_colonization_trade_off = 0
  b = reshape(a(randperm(numel(a))),nspecies,nspecies); #permute the a(i,j)_values randomly; need to restore the diagonal after this operation
  xx = [1:nspecies^2](b==a(1,1)); #indices of the elements that are equal to a(i,i)
  for ii = 1:length(xx)
    for jj in [0:nspecies-1]*(nspecies+1)+1
      if ii == jj
	xx(ii) = [];
      endif
    endfor
  endfor
  a = b;
  rm(b);
endif

#I should use mean(a(i,j)) value for the analysis of the results

K = r./(diag(a)');

m0_mean = 0.1; #mean species basal dispersal rate

#competition-colonization trade-off
m0_var = 5.0; #strength of the competition-colonization trade-off; 0.0 = no trade-off, 1.0 = direct proportionality with competition

competition_colonization_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off
m0 = m0_mean*((r/r_mean).^m0_var);

if competition_colonization_trade_off = 0
  m0 = m0(randperm(length(m0))); #permute the m0_values randomly
endif

#I use r as an indirect measure of competition strength (species with high r are week competitors and hence strong dispersers)
#mean(m0) differs from m0_mean! I should use mean(m0) for the analysis of the results

##GEOMETRIC MEAN AND MEDIAN ARE UNAFFECTED AND THESE ARE MORE SUITABLE BECAUSE M0 AND A(I,J) ARE RATES SO THEY HAVE A MULTIPLICATIVE EFFECT! ARITHMETIC MEAN IS INAPPROPRIATE IN SUCH CASES, I CAN ALSO LOG-TRANSFORM THE VALUES AND CALCULATE ARITHMETIC MEAN FROM THE LOG-TRANSFORMED VALUES.


#save all parameters
save(in_param, "npatches", "nspecies", "r_mean", "r_range", "r", "a_var", "growth_competition_trade_off", "inter_intra", "a", "K", "m0_mean", "m0_var", "competition_colonization_trade_off", "m0");

endfunction
