function plot_dynamics(t,x)

#plot in time
 colours = [1:5,1:5];
 lwd = [1,1,1,1,1,2,2,2,2,2];
 npatches = length(x(:,1,1));
 nspecies = length(x(1,:,1));
 tmax = t(end);

 for pa = 1:npatches
   figure()
   dat = reshape(x(pa,:,:),nspecies,length(x(pa,:,:)));
   for sp = 1:nspecies
     plot(t,dat(sp,:),int2str(colours(sp)),'LineWidth',lwd(sp));
     hold on;
   endfor
   axis([0,tmax,0,200]);
   title('Patch')
   xlabel('t');
   ylabel('x1, x2');
 endfor

endfunction



