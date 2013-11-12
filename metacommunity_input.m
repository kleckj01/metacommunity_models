function metacommunity_input(in_param)
# stochastic competition model - n species with logistic growth in p patches

#turn off the broadcasting warning
warning("off", "Octave:broadcast");

npatches = 5;

nspecies = 5;
#all parameters are the same for all patches
r_mean = 1.5;
r_range = 0.8;

r = linspace(r_mean - r_range/2, r_mean + r_range/2, nspecies); #intrinsic growth rates; environmental stochasticity could be modelled by adding random noise to r values

asym = 1.5; #a(2,1) = a(2,2)*asym; a(2,3) = a(2,2)/asym, etc.; 1.0 is the symetric competition case
#increasing the asym leads to increasing the mean a(i,j) value; this can confound the effects of increasing the overall competition strength and incresing competition asymmetry
#solution could be to multiply non-diagonal elements of the matrix a by a certain number to make sure that mean a(i,j) is constant (equal to a(i,i))
#BEWARE: it changes the shape of the asymetry (but the result seems reasonable)

a = zeros(nspecies,nspecies);
for i = 1:nspecies
  a(i,i) = 0.01; #set intraspecific competition coefficients;remember that K(i)=r(i)/a(i,i)
  j = i;
  while (j < nspecies)
    a(i,j+1) = a(i,j)/asym;
    j = j+1;
  end
  j = i;
  while (j > 1)
    a(i,j-1) = a(i,j)*asym;
    j = j-1;
  end
endfor

#recalculate the a(i,j) values to keep mean(a(i,j)) = a(i,i)
div = (sum(sum(a-diag(diag(a),nspecies,nspecies)))/(nspecies^2-nspecies))/a(1,1);
a = a/div;
for i = 1:nspecies
  a(i,i) = 0.01;
endfor

K = r./(diag(a)');

m0_mean = 0.1; #mean species basal dispersal rate

#competition-colonization trade-off - I NEED TO MAKE AN ASSUMPTION ON THE FUNCTIONAL FORM OF THE TRADE-OFF.
#I use r as an indirect measure of competition strength (species with high r are week competitors and hence strong dispersers)
cc_strength = 5.0; #strength of the competition-colonization trade-off; 0.0 = no trade-off, 1.0 = direct proportionality with competition
m0 = m0_mean*((r/r_mean).^cc_strength);

#save all parameters
save(in_param, "npatches", "nspecies", "r_mean", "r_range", "r", "asym", "inter_intra", "a", "K", "m0_mean", "cc_strength", "m0");

endfunction
