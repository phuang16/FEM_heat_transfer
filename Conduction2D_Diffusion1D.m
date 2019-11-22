close all; clear all; clc

% Tissue Parameters
rho = 1000*1e-9;% kg/mm^3
c = 1000; % J/(kg C)
u = 0; % velocity(mm/s)
k = 0.00055;% diffusion term (W/mm degC)
%h = 500*1e6; % heat tranfer coefficient (W/mm^2 k)
%Tn = 37; % degC
%q = 500*10^6; % heat flow

SAR_all = [0, 100]; % Specific absorption rate (W/kg)
timepts_all = [30, 60, 120]; % time points (sec)
metabolic = 0; % metabolic rate

for iSAR = 1:length(SAR_all)
    SAR = SAR_all(iSAR);
    
    % in y direction
    L1 = 20;% length (mm)
    nel = 10;% number of elements
    nmz = nel + 1;%number of nodes
    
    % in x direction
    L2 = 2 ;
    nel2 = 10;% number of elements
    nmz2 = nel2 + 1; %number of nodes
    
    ele_tot = nmz*nmz2;
    %% coordinate for nodes and triangluar elements
    coor_y = [0 : L1/nel : L1 ];  %  coordinates in y
    coor_x = [0 : L2/nel2 : L2 ];  %  coordinates in x
    
    % === coordinate for all nodes ===
    for iel = 1:nmz2
        node_loc_all((iel-1)*nmz+1:iel*nmz, :) = [coor_x(iel)*ones(nmz,1), coor_y'] ;
    end
    
    % == order of the nodes in each element ===
    % node_order_all = [ 1, 2, 5;   % corresponding to local nodes i, j, k of an triangular element
    %                    2, 3, 6;
    %                    3, 4, 7;
    %                    5, 6, 9;
    %                    6, 7, 10;
    %                    7, 8, 11;
    %                    2, 6, 5;
    %                    3, 7, 6;
    %                    4, 8, 7;
    %                    6, 10, 9;
    %                    7, 11, 10;
    %                    8, 12, 11];
    node_order_all = zeros(nmz*nmz2, 3);
    for i = 1:nel2
        node_order_all((i-1)*nel+1:i*nel, 1) =  (i-1)*nmz+1:i*nmz-1 ; %i
        node_order_all((i-1)*nel+1:i*nel, 2) =  (i-1)*nmz+2:i*nmz ; %j
        node_order_all((i-1)*nel+1:i*nel, 3) =  (i)*nmz+1:(i+1)*nmz-1 ; %k
        
        node_order_all((nel2+i-1)*nel+1: (nel2+i)*nel, 1) =  (i-1)*nmz+2:i*nmz ; %i
        node_order_all((nel2+i-1)*nel+1: (nel2+i)*nel, 2) =  (i)*nmz+2:(i+1)*nmz ; %j
        node_order_all((nel2+i-1)*nel+1: (nel2+i)*nel, 3) =  (i)*nmz+1:(i+1)*nmz-1 ; %k
    end
    
    %% find the elements with boundary edge
    clear  boundary_upper boundary_lower
    for  i = 1:nel2
        ele_boundary_upper(i) = (nel2+i)*nel ;
        ele_boundary_lower(i) = nel*(i-1)+1 ;
    end
    
    %% === Compute matrix at each element ===
    % If all element size is the same:
    Area_total  = node_loc_all(end,:) - node_loc_all(1,:);
    Area = Area_total(1)*Area_total(2)/(ele_tot);
    
    fe = zeros(3,1);
    fe_all = zeros(3, 1, ele_tot);
    
    for iele = 1:ele_tot
        coorde = node_loc_all(node_order_all(iele, :), :);
        
        % Me
        Me = rho*c*Area/12*[2, 1, 1;, 1, 2, 1; 1, 1, 2];
        Me_all(:,:,iele) = Me;
        
        % Diffusion 
        bi = coorde(2,2)-coorde(3,2);
        bj = coorde(3,2)-coorde(1,2);
        bk = coorde(1,2)-coorde(2,2);
        ci = coorde(3,1)-coorde(2,1);
        cj = coorde(1,1)-coorde(3,1);
        ck = coorde(2,1)-coorde(1,1);
        
        B = [bi^2, bi*bj, bi*bk; bj*bi, bj^2, bj*bk; bk*bi, bk*bj, bk^2];
        C = [ci^2, ci*cj, ci*ck; cj*ci, cj^2, cj*ck; ck*ci, ck*cj, ck^2];
        
        Ke = (k/(4*Area))*(B+C);
        Ke_all(:,:,iele) = Ke;
        
        % Boundary
        if ( any(iele== ele_boundary_lower) == 1) % lower bound, boundary edge i-k
            Tn = 37;
            L_boundary = coorde(3,1)- coorde(1,1);
            fe = (L_boundary*k/(4*Area))*([bi*Tn + bk*Tn; 0; bi*Tn + bk*Tn] + [ci*Tn + ck*Tn; 0; ci*Tn + ck*Tn]) ;
            fe_all(:,:,iele) = fe;
        end
        if ( any(iele== ele_boundary_upper) == 1) % upper bound, boundary edge i-j
            Tn = 35;
            L_boundary = coorde(2,1)- coorde(1,1);
            fe = (L_boundary*k/(4*Area))*([bi*Tn + bj*Tn; bi*Tn + bj*Tn; 0] + [ci*Tn + cj*Tn; ci*Tn + cj*Tn; 0]);
            fe_all(:,:,iele) = -fe;
        end
        
    end
    
    for iele = 1:ele_tot
        fe_all(:,:,iele) = fe_all(:,:,iele) + (rho*SAR + metabolic )*(Area/3)*ones(3, 1);
    end
    
    %% ==== Assembly =====
    global_M_allele = zeros(ele_tot, ele_tot);
    global_K_allele = zeros(ele_tot, ele_tot);
    global_F_allele = zeros(ele_tot, 1);
    
    for iele = 1:ele_tot
        global_M_iele = zeros(ele_tot, ele_tot);
        global_K_iele = zeros(ele_tot, ele_tot);
        global_F_iele = zeros(ele_tot, 1);
        
        global_M_iele(node_order_all(iele, :),node_order_all(iele, :)) =  Me_all(:,:,iele);
        global_M_allele(:,:,iele)  = global_M_iele;
        global_K_iele(node_order_all(iele, :),node_order_all(iele, :)) = Ke_all(:,:,iele);
        global_K_allele(:,:,iele)  = global_K_iele;
        
        global_F_iele(node_order_all(iele, :)) =  fe_all(:,:,iele);
        global_F_allele(:,iele) = global_F_iele;
    end
    
    global_M = sum(global_M_allele, 3);
    global_K = sum(global_K_allele, 3);
    global_F = sum(global_F_allele, 2);
    
    %% load in matrix from 1D diffusion, mapped to global matrix
    Me_1D_diffus = dlmread('diffusion1D_20mm_10ele_Me.dat');
    A_1D_diffus = dlmread('diffusion1D_20mm_10ele_A.dat');
    
    node_vessel = nmz*nel2+1:nmz*(nel2+1);
    
    global_M_1D_diffus= zeros(ele_tot, ele_tot);
    global_M_1D_diffus(node_vessel,node_vessel) =  Me_1D_diffus;
    
    global_A_1D_diffus = zeros(ele_tot, ele_tot);
    global_A_1D_diffus(node_vessel,node_vessel) =  -A_1D_diffus;
    
    global_M_TissueVessel = global_M + global_M_1D_diffus;
    global_A_TissueVessel = global_K + global_A_1D_diffus;
    
    %%  find the nodes at the upper and lower boundary
    
    node_loc_upperele = node_order_all( ele_boundary_upper, :);
    node_boundary_upper = node_loc_upperele(:,1:2);
    node_boundary_upper = unique(node_boundary_upper(:));
    
    node_loc_lowerele = node_order_all( ele_boundary_lower, :);
    node_boundary_lower = node_loc_lowerele(:,[1 3]);
    node_boundary_lower = unique(node_boundary_lower(:));
    
    %% solve for T
    IS = zeros(nmz,nmz2); IS([1 nmz],:) = 1; Boundary = IS(:); % nodes with Boundary condition temperature
    T = zeros(ele_tot, 1); T(node_boundary_lower) = 37; T(node_boundary_upper) = 35;
    
    %%
    deltat = 1e-3; %integration step
    for it = 1:length(timepts_all)
        time = timepts_all(it);
        %Minv=inv(global_M(find(~Boundary),find(~Boundary)));
        Minv=inv(global_M_TissueVessel(find(~Boundary),find(~Boundary)));
        
        for i=0:deltat:time
            % T(find(~Boundary))=T(find(~Boundary))+deltat* Minv*(-global_K(find(~Boundary),:)*T + global_F(~Boundary)) ;
            T(find(~Boundary))=T(find(~Boundary))+deltat* Minv*(-global_A_TissueVessel(find(~Boundary),:)*T + global_F(~Boundary)) ;
            
        end
        
        T_mat = zeros(nmz, nmz2);
        for i = 1:nmz2
            T_mat(:,i) = T( (i-1)*nmz +1 : i*nmz);
        end
        
        %%
        % === plotting===
        [X, Y] = meshgrid(coor_x, coor_y);

        figure(iSAR); subplot(1, length(timepts_all), it)
        surf(X, Y, T_mat);
        caxis([34, 39]);
        h = colorbar;
        ylabel(h, 'Temperature (degC)',  'fontsize', 12)
        view(0,90)
        ylabel('y (mm)', 'fontsize', 14); xlabel('x (mm)', 'fontsize', 14)
        title(strcat('SAR = ', num2str(SAR), ' W/kg, ', num2str(timepts_all(it)), ' sec'), 'fontsize', 14)
        
        T_all(iSAR, :) = T_mat(:);
        
    end
    
end

