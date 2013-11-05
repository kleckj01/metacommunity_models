function metacommunity(logg)

#start recording output using diary
#diary diary1.txt
diary (logg);

#turn off the broadcasting warning
warning ("off", "Octave:broadcast");

# stochastic competition model - 2 species with logistic growth in a single patch
tic() #set timer
clear
set(0,'DefaultAxesFontSize',18);

npatches = 5;

patches = 1:npatches;

nspecies = 5;
#all parameters are the same for all patches
r_mean = 1.5;
r_range = 0.4;

r = linspace(r_mean - r_range, r_mean + r_range, nspecies); #intrinsic growth rates; environmental stochasticity could be modelled by adding random noise to r values
K = 100;
asym = 1.0; #a(2,1) = a(2,2)*asym; a(2,3) = a(2,2)/asym, etc.

a = zeros(nspecies,nspecies);
for i = 1:nspecies
  a(i,i) = r(i)/K; #set intraspecific competition coefficients so that K(i)=r(i)/a(i,i)
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

m0_mean = 0.2; #mean species basal dispersal rate

#competition-colonization trade-off - I NEED TO MAKE AN ASSUMPTION ON THE FUNCTIONAL FORM OF THE TRADE-OFF.
cc_strength = 0.0; #strength of the competition-colonization trade-off; 0.0 = no trade-off, 1.0 = direct proportionality with competition
cc = diag(a)';
m0 = m0_mean*((cc/mean(cc)).^cc_strength);

t(1)=0;
t_max=10;

x = zeros(npatches,nspecies);

 for i = 1:npatches
   for j = 1:nspecies
     x(i,j) = round(K/nspecies); #initial abundances at K with equal abundance of all species
    endfor
 endfor
 
i=1;

while(t(i)<t_max)
  #rates for patch
  #patch = x, species = y
  b(:,:)=[x(:,:,i).*r]; #birth rates
  for p = 1:npatches
    d(p,:)=[x(p,:,i).*(sum((a(:,:).*x(p,:,i))'))]; #death rates
  endfor
  m(:,:)=[x(:,:,i).*m0]; #emigration rates
  
  for p = 1:npatches
    rates(p)=sum(b(p,:))+sum(d(p,:))+sum(m(p,:)); #sum of all rates for individual patches 
  endfor
  
  u1=rand;
  u2=rand;
  
  t(i+1)=-log(u1)/sum(rates)+t(i); #interevent times follow the exponential distribution
  x(:,:,i+1)=x(:,:,i);
    
#select a patch with a probability proportional to the total rate of events in the patch
  randVal = rand(1);
  select_patch = min(find(cumsum(rates)/sum(rates) >= randVal));
  rates = 0;
  
#select a species within the selected patch
  s_rates = b(select_patch,:)+d(select_patch,:)+m(select_patch,:); #vector of rates for individual species in patch p
  randVal = rand(1);
  select_species = min(find(cumsum(s_rates)/sum(s_rates) >= randVal));
  s_rates = 0;

#select a process for this species
  p_rates = [b(select_patch,select_species),d(select_patch,select_species),m(select_patch,select_species)]; #vector of rates for the selected species in patch p
  randVal = rand(1);
  process = min(find(cumsum(p_rates)/sum(p_rates) >= randVal));

  if(process==1)
      x(select_patch,select_species,i+1) = x(select_patch,select_species,i) + 1; #birth
    elseif(process==2)
      x(select_patch,select_species,i+1) = x(select_patch,select_species,i) - 1; #death
    elseif(process==3)
      x(select_patch,select_species,i+1) = x(select_patch,select_species,i) - 1; #migration from one patch to the other
      p_select = patches;#select patch for immigration; in this case selection is random (all patches have the same probability of receiving the immigrant)
      p_select(select_patch) = []; #exclude the patch where the migrant originates
      if(length(p_select) > 1)
	immigration_patch = p_select(round(unifrnd(1,length(p_select)))); #select randomly one patch
	else
	immigration_patch = p_select;
      end
      x(immigration_patch,select_species,i+1) = x(immigration_patch,select_species,i) + 1; #immigration
  end # the end of the if loop
  
i=i+1; #end of the while loop
end

#final abundances of all species in all patches
x(:,:,length(x(1,1,:)))


toc() #reports the number of seconds elapsed


#plot in time
# colours = [1:5,1:5];
# lwd = [1,1,1,1,1,2,2,2,2,2];
# 
# for pa = 1:npatches
#   figure()
#   dat = reshape(x(pa,:,:),nspecies,length(x(pa,:,:)));
#   for sp = 1:nspecies
#     plot(t,dat(sp,:),int2str(colours(sp)),'LineWidth',lwd(sp));
#     hold on;
#   endfor
#   axis([0,t_max,0,max(max(max(x)))*1.1]);
#   title('Patch')
#   xlabel('t');
#   ylabel('x1, x2');
# endfor

#end recording
diary off

endfunction
