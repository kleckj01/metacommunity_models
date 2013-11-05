function plot_dynamics(t,data)

#plot in time
 colours = [1:5,1:5];
 lwd = [1,1,1,1,1,2,2,2,2,2];
 npatches = length(data(:,1,1));
 nspecies = length(data(1,:,1));

 for pa = 1:npatches
   figure()
   dat = reshape(data(pa,:,:),npatches,length(data(pa,:,:)));
   for sp = 1:nspecies
     plot(t,dat(sp,:),int2str(colours(sp)),'LineWidth',lwd(sp));
     hold on;
   endfor
   axis([0,t(end),0,100]);
   title(cstrcat('Patch ',int2str(pa)));
   xlabel('Time');
   ylabel('Abundance');
 endfor

endfunction
