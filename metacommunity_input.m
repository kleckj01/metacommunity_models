function metacommunity_input(in_param)
# stochastic competition model - n species with logistic growth in p patches
#all parameters are the same for all patches

#turn off the broadcasting warning
warning("off", "Octave:broadcast");

npatches = 10;

nspecies = 10;

r_mean = 1.5;
r_range = 2.0;

r = sort(unifrnd(r_mean - r_range/2, r_mean + r_range/2,1,nspecies)); #intrinsic growth rates

intra = 0.01; #the strength of intraspecific comeptition (the same for all species); remember that K(i)=r(i)/a(i,i)
inter_intra = 1.0; #the average ratio of interspecific and intraspecific competition coefficient (this multiplies nondiagonal elements of the competition matrix)
a_var = 2.0; #interspecific variation of a(i,j) values

growth_competition_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off
#The shape of the trade-off: linear on log-log scale

a = zeros(nspecies,nspecies);
for i = 1:nspecies
  a(i,i) = intra; #set intraspecific competition coefficients
    for j = 1:nspecies
      if i != j
	a(i,j) = a(i,i)*(r(i)/r(j))^a_var;
      endif
  endfor
endfor

if growth_competition_trade_off == 0
  r = r(randperm(length(r))); #permute the r_values randomly to destroy the correlation between r and a(i,j), this keeps the original strength of the competition hierarchy, so that some species are stronger but can have large r at the same time
endif

#now multiply non-diagonal elements by inter_intra
aa = a.*inter_intra;
for i = 1:nspecies
  aa(i,i) = intra; #reset intraspecific competition coefficients
endfor

#calculate species competitive strength as an arithmeitc mean of a(i,j) elemenents that describe its effect on all other species
for i = 1:nspecies
  comp = aa(:,i);
  comp(i) = []; #remove the intraspecific term
  comp_strength(i) = mean(comp);
endfor

#The shape of the growth-competition trade-off: linear on log-log scale

K = r./(diag(a)');

m0_mean = 0.1; #mean species basal dispersal rate

m0_var = 2.0; #interspecific variation in m0 values

competition_colonization_trade_off = 1; #1 means that the trade-off is on, 0 means that it is off
#The shape of the competition_colonization trade-off: linear on log-log scale

m0 = m0_mean*(((comp_strength./mean(comp_strength)).^-1).^m0_var);

if competition_colonization_trade_off == 0
  m0 = m0(randperm(length(m0))); #permute the m0_values randomly
endif

#disturbance (its effect is the same for all species)
dist_freq = 0.2; #the rate of disturbance events

dist_extent = 0.1; #the proportion of patches affected by disturbance
#I have to make sure that the number of patches times the proportion to be disturbed is an integer!
if (npatches*dist_extent-round(npatches*dist_extent)) != 0
  "WARNING: the number of patches to be disturbed is not an integer!"
endif
  
dist_intensity = 1.0; #the proportion of individuals of all species killed by the disturbance in the disturbed patch/-es

#comparing different combinations of dist_freq, dist_magnitude and dist_extent allows to disentangle the role of overall increase in mortality
#(~very frequent disturbances of very low intensity; i.e. a small proportion of individuals dies in all patches every few time steps) and true disturbance
#(i.e., event that occurs only now and then but has larger intensity so it is not just background mortality)
#the product of the three parameters could be a good measure of the overall increase in mortality rate;

#save all parameters
save(in_param, "npatches", "nspecies", "r_mean", "r_range", "r", "intra", "inter_intra", "a_var", "growth_competition_trade_off", "a", "K", "comp_strength", "m0_mean", "m0_var", "competition_colonization_trade_off", "m0", "dist_freq", "dist_extent", "dist_intensity");

endfunction
