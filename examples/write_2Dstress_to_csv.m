function convert_2Dstress_to_csv(faultprop,stresses,depidx,outdir,apx)
    % get initial stress, pressure, and stress changes with respect to
    % initial time step
    [sne0, tau0] = stresses.get_initial_stress();
    [dsne, dtau] = stresses.get_stress_changes();
    [p_init] = stresses.get_initial_pressure();


    % define matrix dimensions
    nr_of_pillars=size(stresses.get_initial_stress(),1);
    pillar_lengths=size(sne0{1,1},1);
    
    % initialize empty matrices to write out later
    data_ens_out=zeros(pillar_lengths,nr_of_pillars);
    data_tau_out=zeros(pillar_lengths,nr_of_pillars);
    data_pressure_out=zeros(pillar_lengths,nr_of_pillars);
    data_dtau_out=zeros(pillar_lengths,nr_of_pillars);
    data_dens_out=zeros(pillar_lengths,nr_of_pillars);
    data_depths_out=zeros(pillar_lengths,nr_of_pillars);
    data_dists_out=zeros(pillar_lengths,nr_of_pillars);
    
    % initialize along strike pillar distance array
    pillar_dists=[0];
    for p=2:1:nr_of_pillars
        coord_x0=faultprop.x_coor(p-1);
        coord_x1=faultprop.x_coor(p);
        coord_y0=faultprop.y_coor(p-1);
        coord_y1=faultprop.y_coor(p);
    
        distance=sqrt((coord_x1-coord_x0)^2+(coord_y1-coord_y0)^2);
        pillar_dists(p)=pillar_dists(p-1)+distance;
    
    end

    for p=1:1:nr_of_pillars
        % initial effective normal stress (ens) and shear stress (tau) profiles

        data_ens_out(:,p)=sne0{p,1};
        data_tau_out(:,p)=tau0{p,1};
        data_pressure_out(:,p)=p_init{p,1};
        
        data_dtau_out(:,p)=dtau{p,1}(:,depidx);
        data_dens_out(:,p)=dsne{p,1}(:,depidx);

        data_depths_out(:,p)=stresses.y+stresses.ensemble{p,1};
        data_dists_out(:,p)=pillar_dists(p);
    end

    writematrix(data_ens_out,strcat(outdir,'ens',apx,'.csv'));
    writematrix(data_tau_out,strcat(outdir,'tau',apx,'.csv'));
    writematrix(data_pressure_out,strcat(outdir,'pressure',apx,'.csv'));
    writematrix(data_dens_out,strcat(outdir,'dens',apx,'.csv'));
    writematrix(data_dtau_out,strcat(outdir,'dtau',apx,'.csv'));
    writematrix(data_depths_out,strcat(outdir,'depths',apx,'.csv'));
    writematrix(data_dists_out,strcat(outdir,'pillar_dists',apx,'.csv'));
end