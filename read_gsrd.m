clear all; close all; clc
clear all
r=1;

path='/Users/affe_prime/Desktop/opts_reac_1.txt';

fid=fopen(path);
a=fscanf(fid,'%f');

y_s=a(1); x_s=a(2);
t_end=a(3); t_int=a(4);
F_in_a=a(5); F_in_b=a(6);
beta=a(7); n_cs=a(8);
tau_v=a(9); tau_c=a(10);
gr_pr=a(11);

for t=1:t_end/t_int
    for j=1:x_s
        eps(:,j,t)=a(12+(t-1)*(3+n_cs)*y_s*x_s+(j-1)*y_s:11+(t-1)*(3+n_cs)*y_s*x_s+j*y_s);
        u_x(:,j,t)=a(12+((t-1)*(3+n_cs)+1)*y_s*x_s+(j-1)*y_s:11+((t-1)*(3+n_cs)+1)*y_s*x_s+j*y_s);
        u_y(:,j,t)=a(12+((t-1)*(3+n_cs)+2)*y_s*x_s+(j-1)*y_s:11+((t-1)*(3+n_cs)+2)*y_s*x_s+j*y_s);
        for k=1:n_cs
            psi(:,j,k,t)=a(12+((t-1)*(3+n_cs)+2+k)*y_s*x_s+(j-1)*y_s:11+((t-1)*(3+n_cs)+2+k)*y_s*x_s+j*y_s);
        end
    end
end

clear a

%%
close all
figure('units','normalized','outerposition',[0.3 0.05 0.3 0.8])
set(gcf,'color','w')
colormap jet

vid=0;
if vid
    writerObj=VideoWriter(['/Users/affe_prime/Desktop/vid.avi']);
    writerObj.FrameRate=15;
    writerObj.Quality=100;
    open(writerObj);
end

l_p=70;
u_p=-30;

for t=1:size(eps,3)-1
    for i=1:n_cs/4
        subplot(4,n_cs/4,i)
        imagesc(psi(:,:,i+n_cs/2,t))
        daspect([1 1 1])
        
        subplot(4,n_cs/4,i+n_cs/4)
        imagesc(psi(:,:,i,t)-psi(:,:,i+n_cs/4,t))
        daspect([1 1 1])
        
        subplot(4,n_cs/4,i+n_cs/2)
        imagesc(psi(:,:,i+3*n_cs/4,t))
        daspect([1 1 1])
    end
    
    subplot(4,1,4)
    imagesc(eps(:,:,t))
    colorbar
    daspect([1 1 1])
    text(l_p,u_p,['t=' num2str((t+1)*t_int)],'color','k','fontsize',20,'interpreter','tex')
    if vid
        frame=getframe(gcf);
        writeVideo(writerObj,frame);
    end
    pause(0.0000001)
end
if vid
    close(writerObj);
end