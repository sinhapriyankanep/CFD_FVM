%%2D part
clc
close all
clear all
%% given parameters
L_x=2;
L_y=1;
N_x=20;
N_y=10;
rho=1;
D=0.1;
u=1;
v=1;
dx=L_x/N_x;
dy=L_y/N_y;
Pe_x= (rho*u*dx)/D; 
Pe_y= (rho*v*dy)/D;
%% discretization
x=linspace(dx/2,L_x-(dx/2),N_x);
y=linspace(dy/2,L_y-(dy/2),N_y);

c=zeros(N_y,N_x);
uc=zeros(N_y,N_x);
c_left= zeros(1,N_y);
c_left= zeros(1,N_y);
%% boundary condition
% Gaussian initial case
sigma=L_y/5;
co=1;
yc=L_y/2;
p1=2*(sigma^2);

   for  row=1:N_y

           c_left(row)=co * exp(-((y(row) - yc)^2) / p1);
           
         
           uc_left(row)=exp( - ((y(row) - yc)^2) / p1 ); 
   end
uc(:,1)=uc_left;
c(:,1)=c_left;
    c_new=c;
    uc_new=uc;
epsilon=1e-06;
residual=1;
iter_count=0;
u_iter_count=0;
u_residual=1;


          
            
%% looping (using CDS)

while residual > epsilon
   % iter_count = iter_count + 1;
        
    for i=2:N_y-1 %interior points
        for j=2:N_x-1
         %CDS
          aE=(D/dx)-(u/2);
             aW=(D/dx) +(u/2) ;
             aN=(D/dy) - (v/2);
             aS=(D/dy) +(v/2);
             aP= ((2*D)/dx) + ((2*D)/dy);
             c_new(i,j)=  (aE*c(i,j+1)+aW*c(i,j-1)+aN*c(i-1,j)+aS*c(i+1,j))/aP;

             
            
           

        end
    end

    
    %for lefyt

     for i=2:N_y-1
        for j=1


              % CDS
              aE=D/dx;
              aN=(D/dy)- (v/2);
              aS=(D/dy) +(v/2);
              bE_L= ((2*D)/dy) +(u/2) ;
              aP= (u/2) + ((4*D)/dy) + (D/dx);
              c_new(i,j)=(aE*c_new(i,j+1)+aW*0+aN*c_new(i-1,j)+aS*c_new(i+1,j)+bE_L*c_left(i))/aP;
          
        end
    end

     %for right boundary

     for j=20
        for i=2: N_y-1
            % upwind
              aW=(D/dx) + (u/2);
             aE=0;
             aS=(D/dy) + (v/2);
              aN=(D/dy) - (v/2);
              aP= (u/2) + (D/dx) + ((2*D/dy));
              %bE_rb=(2*D-F)*phi_right;
              c_new(i,j)=(aE*0+aW*c(i,j-1)+aN*c(i-1,j)+aS*c(i+1,j))/aP;

              

        end
    end

     %for top boundary

     for i=1
        for j=2:N_x-1
             
% CDS
             aW=(D/dx) +(u/2);
             aE=(D/dx) - (u/2);
             aN=0;
              aS=(D/dy) +(v/2);
              aP=((2*D)/dx) + (D/dy) + (v/2);
              %bE_tb=(2*D-F)*phi_top;
              c_new(i,j)=(aE*c(i,j+1)+aW*c(i,j-1)+aN*0+aS*c(i+1,j))/aP;

             
        end
    end


  %for bottom boundary

     for i=N_y
        for j=2:N_x-1
              
 
%cds
             aW=(u/2);
             aE=(D/dx) - (u/2);
             aN=(D/dy) - (v/2);
              aS=0;
              aP=((2*D)/dx) ;
              %bE_bb=(2*D+ F)*phi_bot;
              c_new(i,j)=(aE*c(i,j+1)+aW*c(i,j-1)+aS*0+aN*c(i-1,j))/aP;



        end
    end


% corners
% top left corner
for i=1
    for j=1
        
% cds
        aP= (u/2) + (v/2) + ((3*D)/dx) + (D/dy);
        aN=0;
        aW=0;
        aS=(D/dy) +(v/2);
        aE=(D/dx) - (u/2);
        be_tlc=((2*D)/dx) + u;
        c_new(i,j)=(aE*c(i,j+1)+aS*c(i+1,j)+be_tlc.*c_left(i))/aP;



%
    end
end


% topright corner
for i=1
    for j=N_x
%cds
            aE=0;
        aN=0;
        aP= (u/2) +(v/2) + (D/dx) + (D/dy);
        aS=(D/dy) +(v/2);
        aW=(D/dx) + (u/2);
        %be_trc=(2*D-F)*phi_right+(2*D-F)*phi_top;
        c_new(i,j)=(aW*c(i,j-1)+aS*c(i+1,j))/aP;

        
       
    end
end

%bottom left corner
for i=N_y
    for j=1


          % cds
        aE=(D/dx) - (u/2) ;
        aN=(D/dy) - (v/2);
        aW=0;
        aS=0;
       aP=(u/2) - (v/2) +(3*D) + (D/dy);
        blc=((2*D)/dx) + u;
        c_new(i,j)=(aE*c(i,j+1)+aN*c(i-1,j)+blc.*c_left(i))/aP;
       
        



    end
end

%bottom right corner
for i=N_y
    for j=N_x
     % cds
        aE=0;
        aN=(D/dy) - (v/2) ;
        aW=(D/dx) + (u/2);
        aS=0;
        aP= (u/2) - (v/2) + (D/dx) + (D/dy);
        %be_brc=(2*D-F)*phi_right+(2*D + F)*phi_bot;
        c_new(i,j)=(aN*c(i-1,j)+aW*c(i,j-1))/aP;




    end
end

    iter_count=iter_count+1;
 residual=sum(abs(c_new(:)-c(:)));
   %u_residual=sum(abs(uc_new(:)-uc(:)));
   % semilogy(iter_count,residual,'o')
  %  hold on
    c=c_new;
    %uc=uc_new;
end

%% Looping using UPWIND
while u_residual> epsilon

    for i=2:N_y-1 %interior points
        for j=2:N_x-1
         

             
             %  upwind
             u_aE=(D/dx) ;
             u_aW=(D/dx) +u ;
             u_aN=(D/dy);
             u_aS=(D/dy) +v;
             u_aP= u+v+ ((2*D)/dx) + ((2*D)/dy);
             uc_new(i,j)=  (u_aE*uc(i,j+1)+u_aW*uc(i,j-1)+u_aN*uc(i-1,j)+u_aS*uc(i+1,j))/u_aP;

           

        end
    end

    
    %for lefyt

     for i=2:N_y-1
        for j=1


              
              % upwind
              u_aE=D/dx;
              u_aN=(D/dy);
              u_aS=(D/dy) +v;
              ubE_L= ((2*D)/dx) + u ;
              u_aP=u+ v+((3*D)/dx)+ ((2*D)/dy);
              uc_new(i,j)=(u_aE*uc_new(i,j+1)+u_aW*0+u_aN*uc_new(i-1,j)+u_aS*uc_new(i+1,j)+ubE_L*uc_left(i))/u_aP;
        end
    end

     %for right boundary

     for j=20
        for i=2: N_y-1
            
              % upwind
              u_aW=(D/dx) + u;
             u_aE=0;
             u_aS=(D/dy) +v;
              u_aN=(D/dy);
              u_aP= u+v+ (D/dx) +((2*D)/dy);
              %bE_rb=(2*D-F)*phi_right;
              uc_new(i,j)=(u_aE*0+u_aW*uc(i,j-1)+u_aN*uc(i-1,j)+u_aS*uc(i+1,j))/u_aP;

        end
    end

     %for top boundary

     for i=1
        for j=2:N_x-1
             


              % upwind
              u_aW=(D/dx) +u;
             u_aE=(D/dx);
             u_aN=0;
              u_aS=(D/dy) +v;
              u_aP=u+v+ ((2*D)/dx)+ (D/dy); % changed
              %bE_tb=(2*D-F)*phi_top;
              uc_new(i,j)=(u_aE*uc(i,j+1)+u_aW*uc(i,j-1)+u_aN*0+u_aS*uc(i+1,j))/u_aP;
        end
    end


  %for bottom boundary

     for i=N_y
        for j=2:N_x-1
              
 

%upwind
             u_aW=(D/dx)+u;
             u_aE=(D/dx);
             u_aN=(D/dy);
              u_aS=0;
              u_aP=u+((2*D)/dx) + (D/dy); %changed
              %bE_bb=(2*D+ F)*phi_bot;
              uc_new(i,j)=(u_aE*uc(i,j+1)+u_aW*uc(i,j-1)+u_aS*0+u_aN*uc(i-1,j))/u_aP;

        end
    end


% corners
% top left corner
for i=1
    for j=1
        

% upwind
        u_aP= u+v + ((3*D)/dx) + (D/dy);
        u_aN=0;
        u_aW=0;
        u_aS=(D/dy) +v;
        u_aE=(D/dx);
        u_be_tlc=((2*D)/dx) + u;
        uc_new(i,j)=(u_aE*uc(i,j+1)+u_aS*uc(i+1,j)+u_be_tlc.*uc_left(i))/u_aP;


    end
end


% topright corner
for i=1
    for j=N_x


        
        %upwind
            u_aE=0;
        u_aN=0;
        u_aP= u+v + (D/dx) + (D/dy);
        u_aS=(D/dy) +v;
        u_aW=(D/dx) + u;
        %be_trc=(2*D-F)*phi_right+(2*D-F)*phi_top;
        uc_new(i,j)=(u_aW*uc(i,j-1)+u_aS*uc(i+1,j))/u_aP;
    end
end

%bottom left corner
for i=N_y
    for j=1


      
       
        % upwind
        u_aE=(D/dx);
        u_aN=(D/dy);
        u_aW=0;
        u_aS=0;
       u_aP=u +((3*D)/dx) + (D/dy);
        u_blc=((2*D)/dx) +u;
        uc_new(i,j)=(u_aE*uc(i,j+1)+u_aN*uc(i-1,j)+u_blc.*uc_left(i))/u_aP;



    end
end

%bottom right corner
for i=N_y
    for j=N_x
     


% upwind
        u_aE=0;
        u_aN=(D/dy) ;
        u_aW=(D/dx) + u ;
        u_aS=0;
        u_aP= u+ (D/dx) + (D/dy);
        %be_brc=(2*D-F)*phi_right+(2*D + F)*phi_bot;
        uc_new(i,j)=(u_aN*uc(i-1,j)+u_aW*uc(i,j-1))/u_aP;



    end
end

    u_iter_count=u_iter_count+1;
 %residual=sum(abs(c_new(:)-c(:)));
   u_residual=sum(abs(uc_new(:)-uc(:)));
   % semilogy(iter_count,residual,'o')
  %  hold on
    %c=c_new;
    uc=uc_new;


end




figure
contourf(x,L_y-y,c,256,'EdgeColor','none')
%colorbar
cb = colorbar;                     % Get colorbar handle
ylabel(cb, 'Concentration')  % Label the
ylabel('y')
xlabel('x')
title('CDS- Concentration contour')
figure
contourf(x,L_y-y,uc,128,'EdgeColor','none')
cb1 = colorbar;                     % Get colorbar handle
ylabel(cb1, 'Concentration')
ylabel('y')
xlabel('x')
title('Upwind- Concentration contour')



%% plotting
c_left_boundary = c(:,1);
c_right_boundary = c(:,end);

uc_left_boundary = uc(:,1);
uc_right_boundary = uc(:,end);
%single
figure;
plt1=plot( y, c_left_boundary,'m--', 'LineWidth', 2); hold on;
plt2=plot( y, c_right_boundary, 'r--', 'LineWidth', 2);
ylabel('concentration')
xlabel('y')
title('CDS-concentration ')
legend([plt1 plt2],{'Left boundary', 'Right boundary'}) ;
grid on
% Plot comparison
figure;
pl1=plot( y, c_left_boundary,'m--', 'LineWidth', 2); hold on;
pl2=plot( y, c_right_boundary, 'r--', 'LineWidth', 2);
pl3=plot( y, uc_left_boundary,'c-', 'LineWidth', 0.5); 
pl4=plot( y, uc_right_boundary, 'g-', 'LineWidth', 0.5);
ylabel('concentration')
xlabel('y')
title('Comparison of Left and Right Boundary Concentration')
legend([pl1 pl2 pl3 pl4],{'CDS-Left boundary', 'CDS-Right boundary','Upwind-Left boundary','Upwind-Right boundary'}) ;
grid on

%figure;

% ylabel('Concentration')
% xlabel('y')
% title('Upwind Comparison of Left and Right Boundary Concentration')
% legend('Left Boundary (Inlet)', 'Right Boundary (Outlet)', 'Location', 'best')
% grid on