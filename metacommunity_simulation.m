function metacommunity_simulation(in_param,output)
# stochastic competition model - n species with logistic growth in n patches

#turn off the broadcasting warning
warning("off", "Octave:broadcast");

tic() #set timer

#load input parameters
load(in_param)

#set duration
t(1)=0;
t_max=20;

x = zeros(npatches,nspecies);

#the following starts with all patches occupied by all species
# for i = 1:npatches
#   x(i,:) = round(K./nspecies); #initial abundances at K with equal abundance of all species
# endfor

#the following starts wiht one patch filled, others empty
 for i = 1
   x(i,:) = round(K./nspecies); #initial abundances at K with equal abundance of all species
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
      p_select = 1:npatches;#select patch for immigration; in this case selection is random (all patches have the same probability of receiving the immigrant)
      p_select(select_patch) = []; #exclude the patch where the migrant originates
      if(length(p_select) > 1)
	immigration_patch = p_select(randperm(length(p_select))(1)); #select randomly one patch
	else
	immigration_patch = p_select;
      end
      x(immigration_patch,select_species,i+1) = x(immigration_patch,select_species,i) + 1; #immigration
  end # the end of the if loop
  
i=i+1; #end of the while loop
end

#final abundances of all species in all patches
save(output, "t", "x");

toc() #reports the number of seconds elapsed

endfunction
