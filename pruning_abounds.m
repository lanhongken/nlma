function simulations = pruning_abounds( M_, options_, shock_sequence, simul_length, pruning_order, pruning_type, use_cached_nlma_values, initial_state, call_back, call_back_arg )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pruning_abounds.m
%
% This file implements all of the second and third order pruning algorithms
% documented in "Pruning in Perturbation DSGE Models" by Hong Lan and 
% Alexander Meyer-Gohde.
%
% The pruning algorithms are called by setting the option pruning_type:
%
% For second order approximations (options_.order=2 in Dynare)
% pruning_type = 
%'kim_et_al' :  The second order algorithm of 
%               KIM, J., S. KIM, E. SCHAUMBURG, AND C. A. SIMS (2008): 
%               "Calculating and Using Second-Order Accurate Solutions of 
%               Discrete Time Dynamic Equilibrium Models",Journal of 
%               Economic Dynamics and Control, 32(11), 3397-414. This algorithm
%				is identical to the series expansion at the determinisitic 
%				steady state.
%
%
%'lan_meyer-gohde' : The second order algorithm of LAN, H., AND A. 
%                    MEYER-GOHDE(2013): "Solving DSGE Models with a 
%                    Nonlinear Moving Average", Journal of Economic 
%                    Dynamics and Control, 37(12), 2643-2667.
%					 This algorithm is identical to a series expansion at the 
%					 matched stochastic steady state and the second order 
%					 algorithm of DEN HAAN, W. J., AND J.
%                     DE WIND (2012): "Nonlinear and Stable Perturbation 
%                     Based Approximations", Journal of Economic Dynamics 
%                     and Control, 36(10), 1477-497.
%
%'lan_meyer-gohde_det' : The second order nonlinear moving average with
%					 a deterministic past.
%
% For third order approximations (options_.order=3 in Dynare)
% pruning_type = 
%'andreasen' : The third order algorithm of ANDREASEN, M. M. (2012): "On 
%              the Effects of Rare Disasters and Uncertainty Shocks for 
%              Risk Premia in Non-Linear DSGE Models", Review of Economic 
%              Dynamics, 15(3), 295-316. This algorithm
%				is identical to the series expansion at the determinisitic 
%				steady state.
%
%'fernandez-villaverde_et_al' : The third order algorithm of 
%               FERNANDEZ-VILLAVERDE, J., P. A. GUERRO N-QUINTANA, J. 
%               RUBIO-RAMI REZ, AND M. URIBE (2011): "Risk Matters: The 
%               Real Effects of Volatility Shocks", American Economic
%               Review, 101(6), 2530-61.
%
%'den_haan_de_wind' : The third order algorithm of DEN HAAN, W. J., AND J.
%                     DE WIND (2012): "Nonlinear and Stable Perturbation 
%                     Based Approximations", Journal of Economic Dynamics 
%                     and Control, 36(10), 1477-497.
%
%'den_haan_de_wind_alt' : The alterante third order algorithm of DEN HAAN, W. J., AND J.
%                     DE WIND (2012): "Nonlinear and Stable Perturbation 
%                     Based Approximations", Journal of Economic Dynamics 
%                     and Control, 36(10), 1477-497. This uses the highest
%					  available stationary series in the recursion for the order
%					  in question
%
%'lan_meyer-gohde' : The second order algorithm of LAN, H., AND A. 
%                    MEYER-GOHDE(2013): "Solving DSGE Models with a 
%                    Nonlinear Moving Average", Journal of Economic 
%                    Dynamics and Control, 37(12), 2643 - 2667. This algorithm 
%				     is identical to a series expansion at the 
%					 matched stochastic steady state
%
%'lan_meyer-gohde_det' : The third order nonlinear moving average with
%					 a deterministic past.
%
%THIS VERSION: 1.1.0 March 13, 2014
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
global oo_
%--------------------------------------------------------------------------
% 0.1. Dynare version check
%--------------------------------------------------------------------------
  % Starting in Dynare 4.4.0, the following fields are no longer in oo_.dr, 
  % they can be found in M_
    [numeric_version] = return_dynare_version(dynare_version);
    if numeric_version >= 4.4 
        nstatic = M_.nstatic;
        nspred = M_.nspred; % note M_.nspred = M_.npred+M_.nboth;
        % nboth = M_.nboth;
        % nfwrd = M_.nfwrd;
    else
        nstatic = oo_.dr.nstatic;
        nspred = oo_.dr.npred;
        % nboth = oo_.dr.nboth;
        % nfwrd = oo_.dr.nfwrd;
    end

  % Build state variable selector
    select_state = nstatic+1:nstatic+nspred;

%--------------------------------------------------------------------------
% 0.2. Set up initial value of the state and pre-allocate memory
%--------------------------------------------------------------------------
    
    simul_length_p1 = simul_length + 1;
    
    simulation_first = zeros(M_.endo_nbr,simul_length_p1);
    if pruning_order >= 2
        simulation_second = zeros(M_.endo_nbr,simul_length_p1);
        if pruning_order >= 3
            simulation_third = zeros(M_.endo_nbr,simul_length_p1);
            simulation_first_sigma_2 = zeros(M_.endo_nbr,simul_length_p1);
        end
    end
    
    if nargin < 7
        use_cached_nlma_values = 0;
    end
    if nargin >= 8
        simulation_first( :, 1 ) = initial_state.first( oo_.dr.order_var );
        if pruning_order >= 2
            simulation_second( :, 1 ) = initial_state.second( oo_.dr.order_var );
            if pruning_order >= 3
                simulation_third( :, 1 ) = initial_state.third( oo_.dr.order_var );
                simulation_first_sigma_2( :, 1 ) = initial_state.first_sigma_2( oo_.dr.order_var );
            end
        end
    end
    
%--------------------------------------------------------------------------
% 1. Simulate first order solution for all the algorithms
%--------------------------------------------------------------------------
  if pruning_order == 1
    E = shock_sequence;
    for t=2:simul_length_p1
      simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t-1);
       if nargin >= 10
          call_back(call_back_arg);
       end
    end
    simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
    simulations.constant=oo_.dr.ys;
    simulations.total=simulations.first+repmat(simulations.constant,[1 simul_length]);    
  end

%--------------------------------------------------------------------------
% 2. Simulate second order pruned solutions
%--------------------------------------------------------------------------
  if pruning_order == 2
    
     % Kim et al's second order pruned solution
       if strcmp(pruning_type,'kim_et_al')
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014, in nlma_simul.m, shock_sequence is oo_.exo_simul' by default.
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*(oo_.dr.ghuu*exe+oo_.dr.ghs2);
          for t=2:simul_length
               simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
               exe=alt_kron(E(:,t),E(:,t));
               sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
               sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
               simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                       +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe+oo_.dr.ghs2);
               if nargin >= 9
                  call_back(call_back_arg);
               end
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.second=simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.constant=oo_.dr.ys;
          simulations.total=simulations.second+simulations.first+repmat(simulations.constant,[1 simul_length]);
       end
    
     % Den haan and de Wind's second order pruned solution
       if strcmp(pruning_type,'den_haan_de_wind')
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe = alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*oo_.dr.ghuu*exe;
          for t=2:simul_length
               simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
               exe=alt_kron(E(:,t),E(:,t));
               sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
               sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
               simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                       +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe);
               if nargin >= 9
                  call_back(call_back_arg);
               end
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.second=simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.constant=oo_.dr.ys+(1/2)*oo_.dr.ghs2(oo_.dr.inv_order_var);
          simulations.total=simulations.second+simulations.first+repmat(simulations.constant,[1 simul_length]);
       end
      
     % Lan and Meyer-gohde's second order solution
       if strcmp(pruning_type, 'lan_meyer-gohde')
            if use_cached_nlma_values
                ghs2_nlma = oo_.dr.ghs2_nlma;
            else
              % Compute nlma's y_sigma^2
                ghs2_state_nlma = (eye(nspred)-oo_.dr.ghx(nstatic+1:nstatic+nspred,:))\(oo_.dr.ghs2(nstatic+1:nstatic+nspred,:));
                ghs2_nlma = [ oo_.dr.ghx(1:nstatic,:)*ghs2_state_nlma+oo_.dr.ghs2(1:nstatic,:)
                              ghs2_state_nlma
                              oo_.dr.ghx(nstatic+nspred+1:M_.endo_nbr,:)*ghs2_state_nlma+oo_.dr.ghs2(nstatic+nspred+1:M_.endo_nbr,:)]; 
              % Save results
                oo_.dr.ghs2_state_nlma = ghs2_state_nlma;
                oo_.dr.ghs2_nlma = ghs2_nlma;
            end
          % Simulation
            E = shock_sequence;
            for t = 2:simul_length_p1
                simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t-1);
                exe = alt_kron(E(:,t-1),E(:,t-1));
                sxe = alt_kron(simulation_first(select_state,t-1),E(:,t-1));
                sxs = alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
                simulation_second(:,t) = oo_.dr.ghx*simulation_second(select_state,t-1)...
                                         +(1/2)*( oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe );
               if nargin >= 9
                  call_back(call_back_arg);
               end
            end
            simulations.first = simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
            simulations.second = simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
            simulations.constant = oo_.dr.ys + 0.5*ghs2_nlma(oo_.dr.inv_order_var,:);
            simulations.total = simulations.second + simulations.first...
                                 +repmat( simulations.constant,[1 simul_length] );
       end
     % Lan and Meyer-gohde's second order solution with a deterministic past
       if strcmp(pruning_type, 'lan_meyer-gohde_det')
          % Simulation
            E = shock_sequence;
            simulation_first(:,1) = oo_.dr.ghu*E(:,1);
            exe = alt_kron(E(:,1),E(:,1));
            simulation_second(:,1) = (1/2)*oo_.dr.ghuu*exe;
            for t = 2:simul_length
                simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
                exe = alt_kron(E(:,t),E(:,t));
                sxe = alt_kron(simulation_first(select_state,t-1),E(:,t));
                sxs = alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
                simulation_second(:,t) = oo_.dr.ghx*simulation_second(select_state,t-1)...
                                         +(1/2)*( oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe );
            end
            simulations.first = simulation_first(oo_.dr.inv_order_var,:);
            simulations.second = simulation_second(oo_.dr.inv_order_var,:);
            simulations.total = simulations.second + simulations.first...
                                 +repmat( oo_.dr.ys + 0.5*oo_.dr.ghs2(oo_.dr.inv_order_var,:),[1 simul_length] );
       end
  end

%--------------------------------------------------------------------------
% 3. Simulate third order pruned solutions
%--------------------------------------------------------------------------
  if pruning_order==3
     if options_.pruning == 0
         % This seperates the risk correction of the first order coefficients 
         % from the first order coefficients themselves and seperates out 
         % the blocks of the third order policy function ghxxx etc.
            oo_ = full_block_dr_new(oo_,M_,options_);
     end
     
     % Andreasen's third order pruned solution
       if strcmp(pruning_type,'andreasen')
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*(oo_.dr.ghuu*exe+oo_.dr.ghs2);
          simulation_first_sigma_2(:,1)=(1/2)*(oo_.dr.ghuss*E(:,1));
          simulation_third(:,1)=(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe);
          for t=2:simul_length
              simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe+oo_.dr.ghs2);
              simulation_first_sigma_2(:,t)=oo_.dr.ghx*simulation_first_sigma_2(select_state,t-1)+(1/2)*(oo_.dr.ghxss*simulation_first(select_state,t-1)+oo_.dr.ghuss*E(:,t-1));
              simulation_third(:,t)=oo_.dr.ghx*simulation_third(select_state,t-1)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t-1),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t-1))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe))...
                                     +oo_.dr.ghxx*alt_kron(simulation_second(select_state,t-1),simulation_first(select_state,t-1))+oo_.dr.ghxu*alt_kron(simulation_second(select_state,t-1),E(:,t-1));
               if nargin >= 9
                  call_back(call_back_arg);
               end
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.second=simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.first_sigma_2=simulation_first_sigma_2(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.third=simulation_third(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.constant=oo_.dr.ys;
          simulations.total=simulations.third+simulations.first_sigma_2...
                            +simulations.second+simulations.first+repmat(simulations.constant,[1 simul_length]);
       end
	  % Third order series expansion with accumulated constant risk correction
       if strcmp(pruning_type,'series_expansion_det')
	       % Compute nlma's y_sigma^2
           ghs2_state_nlma = (eye(oo_.dr.npred)-oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:))\(oo_.dr.ghs2(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:));
           ghs2_nlma = [ oo_.dr.ghx(1:oo_.dr.nstatic,:)*ghs2_state_nlma+oo_.dr.ghs2(1:oo_.dr.nstatic,:)
                         ghs2_state_nlma
                         oo_.dr.ghx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*ghs2_state_nlma+oo_.dr.ghs2(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)]; 
         % Save results
           oo_.dr.ghs2_state_nlma = ghs2_state_nlma;
           oo_.dr.ghs2_nlma = ghs2_nlma;
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*(oo_.dr.ghuu*exe);
          simulation_first_sigma_2(:,1)=(1/2)*(oo_.dr.ghuss*E(:,1));
          simulation_third(:,1)=(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe);
          for t=2:simul_length
              simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe);
              simulation_first_sigma_2(:,t)=oo_.dr.ghx*simulation_first_sigma_2(select_state,t-1)+(1/2)*(oo_.dr.ghxss*simulation_first(select_state,t-1)+oo_.dr.ghuss*E(:,t));
              simulation_third(:,t)=oo_.dr.ghx*simulation_third(select_state,t-1)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe))...
                                     +oo_.dr.ghxx*alt_kron(simulation_second(select_state,t-1),simulation_first(select_state,t-1))+oo_.dr.ghxu*alt_kron(simulation_second(select_state,t-1),E(:,t));
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,:);
          simulations.second=simulation_second(oo_.dr.inv_order_var,:);
          simulations.first_sigma_2=simulation_first_sigma_2(oo_.dr.inv_order_var,:);
          simulations.third=simulation_third(oo_.dr.inv_order_var,:);
          simulations.total=simulations.third+simulations.first_sigma_2...
                            +simulations.second+simulations.first+repmat( oo_.dr.ys + 0.5*ghs2_nlma(oo_.dr.inv_order_var,:),[1 simul_length] );
       end      
     % Den haan and de Wind's third order pruned solution  %% -> AMG 2.2.16
       if strcmp(pruning_type,'den_haan_de_wind')
                    % Compute nlma's y_sigma^2
           ghs2_state_nlma = (eye(oo_.dr.npred)-oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:))\(oo_.dr.ghs2(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:));
           ghs2_nlma = [ oo_.dr.ghx(1:oo_.dr.nstatic,:)*ghs2_state_nlma+oo_.dr.ghs2(1:oo_.dr.nstatic,:)
                         ghs2_state_nlma
                         oo_.dr.ghx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*ghs2_state_nlma+oo_.dr.ghs2(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)]; 
         % Compute nlma's y_sigma^2e and y_sigma^2y^state
           % y_sigma^2e           
             ghuss_nlma = oo_.dr.ghuss + oo_.dr.ghxu*alt_kron(ghs2_state_nlma,eye(M_.exo_nbr));
           % y_sigma^2y^state
             ghxss_nlma = oo_.dr.ghxss + oo_.dr.ghxx*alt_kron(ghs2_state_nlma,eye(oo_.dr.npred));
         % Save results
           oo_.dr.ghs2_state_nlma = ghs2_state_nlma;
           oo_.dr.ghs2_nlma = ghs2_nlma;
           oo_.dr.ghuss_nlma = ghuss_nlma;
           oo_.dr.ghxss_nlma = ghxss_nlma;
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*oo_.dr.ghuu*exe+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          simulation_third(:,1)=(1/2)*oo_.dr.ghuu*exe+(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          for t=2:simul_length
              simulation_first(:,t)=(oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_first(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= (oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_second(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe);
               simulation_third(:,t)=(oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_third(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t-1)...
                                     +(1/2)*(oo_.dr.ghxx*alt_kron(simulation_second(select_state,t-1),simulation_second(select_state,t-1))+2*oo_.dr.ghxu*alt_kron(simulation_second(select_state,t-1),E(:,t-1))+oo_.dr.ghuu*exe)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t-1),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t-1))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe));
               if nargin >= 9
                  call_back(call_back_arg);
               end
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.second=simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.third=simulation_third(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.constant=oo_.dr.ys + 0.5*ghs2_nlma(oo_.dr.inv_order_var);
          simulations.total=simulations.third+repmat(simulations.constant,[1 simul_length]);
       end
	   
     % Den haan and de Wind's alternative third order pruned solution  %% -> AMG 2.2.16
       if strcmp(pruning_type,'den_haan_de_wind_alt')         
           % Compute nlma's y_sigma^2
           ghs2_state_nlma = (eye(oo_.dr.npred)-oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:))\(oo_.dr.ghs2(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:));
           ghs2_nlma = [ oo_.dr.ghx(1:oo_.dr.nstatic,:)*ghs2_state_nlma+oo_.dr.ghs2(1:oo_.dr.nstatic,:)
                         ghs2_state_nlma
                         oo_.dr.ghx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*ghs2_state_nlma+oo_.dr.ghs2(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)]; 
         % Compute nlma's y_sigma^2e and y_sigma^2y^state
           % y_sigma^2e           
             ghuss_nlma = oo_.dr.ghuss + oo_.dr.ghxu*alt_kron(ghs2_state_nlma,eye(M_.exo_nbr));
           % y_sigma^2y^state
             ghxss_nlma = oo_.dr.ghxss + oo_.dr.ghxx*alt_kron(ghs2_state_nlma,eye(oo_.dr.npred));
         % Save results
           oo_.dr.ghs2_state_nlma = ghs2_state_nlma;
           oo_.dr.ghs2_nlma = ghs2_nlma;
           oo_.dr.ghuss_nlma = ghuss_nlma;
           oo_.dr.ghxss_nlma = ghxss_nlma;
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*oo_.dr.ghuu*exe+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          simulation_third(:,1)=(1/2)*oo_.dr.ghuu*exe+(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,1);
          for t=2:simul_length
              simulation_first(:,t)=(oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_first(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= (oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_second(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe);
				sxe=alt_kron(simulation_second(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_second(select_state,t-1),simulation_second(select_state,t-1));
               simulation_third(:,t)=(oo_.dr.ghx+(1/2)*oo_.dr.ghxss)*simulation_third(select_state,t-1)+(oo_.dr.ghu+(1/2)*oo_.dr.ghuss)*E(:,t)...
                                     +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_second(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t))+oo_.dr.ghxuu*alt_kron(simulation_second(select_state,t-1),exe));
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,:);
          simulations.second=simulation_second(oo_.dr.inv_order_var,:);
          simulations.third=simulation_third(oo_.dr.inv_order_var,:);
          simulations.total=simulations.third+repmat( oo_.dr.ys + 0.5*ghs2_nlma(oo_.dr.inv_order_var,:),[1 simul_length] );
       end
       
     % Fernandez-villaverde et al's third order pruned solution
       if strcmp(pruning_type,'fernandez-villaverde_et_al')
          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*(oo_.dr.ghuu*exe+oo_.dr.ghs2);
          simulation_first_sigma_2(:,1)=(1/2)*(oo_.dr.ghuss*E(:,1));
          simulation_third(:,1)=(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe);
          for t=2:simul_length
              simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe+oo_.dr.ghs2);
              simulation_first_sigma_2(:,t)=oo_.dr.ghx*simulation_first_sigma_2(select_state,t-1)+(1/2)*(oo_.dr.ghxss*simulation_first(select_state,t-1)+oo_.dr.ghuss*E(:,t-1));
              simulation_third(:,t)= oo_.dr.ghx*simulation_third(select_state,t-1)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t-1),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t-1))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe));
               if nargin >= 9
                  call_back(call_back_arg);
               end
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.second=simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.first_sigma_2=simulation_first_sigma_2(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.third=simulation_third(oo_.dr.inv_order_var,2:simul_length_p1);
          simulations.constant=oo_.dr.ys;
          simulations.total=simulations.third+simulations.first_sigma_2+simulations.second+simulations.first+repmat(simulations.constant,[1 simul_length]);
       end
      
    % Lan and Meyer-Gohde's third order solution   
      if strcmp(pruning_type, 'lan_meyer-gohde')
           if use_cached_nlma_values
               ghs2_nlma = oo_.dr.ghs2_nlma;
               ghuss_nlma = oo_.dr.ghuss_nlma;
               ghxss_nlma = oo_.dr.ghxss_nlma;
           else
             % Compute nlma's y_sigma^2
               ghs2_state_nlma = (eye(nspred)-oo_.dr.ghx(nstatic+1:nstatic+nspred,:))\(oo_.dr.ghs2(nstatic+1:nstatic+nspred,:));
               ghs2_nlma = [ oo_.dr.ghx(1:nstatic,:)*ghs2_state_nlma+oo_.dr.ghs2(1:nstatic,:)
                             ghs2_state_nlma
                             oo_.dr.ghx(nstatic+nspred+1:M_.endo_nbr,:)*ghs2_state_nlma+oo_.dr.ghs2(nstatic+nspred+1:M_.endo_nbr,:)]; 
             % Compute nlma's y_sigma^2e and y_sigma^2y^state
               % y_sigma^2e           
                 ghuss_nlma = oo_.dr.ghuss + oo_.dr.ghxu*alt_kron(ghs2_state_nlma,eye(M_.exo_nbr));
               % y_sigma^2y^state
                 ghxss_nlma = oo_.dr.ghxss + oo_.dr.ghxx*alt_kron(ghs2_state_nlma,eye(nspred));
             % Save results
               oo_.dr.ghs2_state_nlma = ghs2_state_nlma;
               oo_.dr.ghs2_nlma = ghs2_nlma;
               oo_.dr.ghuss_nlma = ghuss_nlma;
               oo_.dr.ghxss_nlma = ghxss_nlma;
           end
         % Simulation
           E = shock_sequence;
           for t = 2:simul_length_p1
               simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t-1);
               exe = alt_kron(E(:,t-1),E(:,t-1));
               sxe = alt_kron(simulation_first(select_state,t-1),E(:,t-1));
               sxs = alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
               simulation_second(:,t) = oo_.dr.ghx*simulation_second(select_state,t-1)...
                                        +(1/2)*( oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe );
               simulation_first_sigma_2(:,t) = oo_.dr.ghx*simulation_first_sigma_2(select_state,t-1)...
                                              +(1/2)*(ghuss_nlma*E(:,t-1)+ghxss_nlma*simulation_first(select_state,t-1));
               simulation_third(:,t) = oo_.dr.ghx*simulation_third(select_state,t-1)...
                                       +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t-1),exe))...
                                       +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t-1))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe))...
                                       +oo_.dr.ghxx*alt_kron(simulation_second(select_state,t-1),simulation_first(select_state,t-1))...
                                       +oo_.dr.ghxu*alt_kron(simulation_second(select_state,t-1),E(:,t-1));
               if nargin >= 9
                  call_back(call_back_arg);
               end
           end
           simulations.first = simulation_first(oo_.dr.inv_order_var,2:simul_length_p1);
           simulations.second = simulation_second(oo_.dr.inv_order_var,2:simul_length_p1);      
           simulations.first_sigma_2 = simulation_first_sigma_2(oo_.dr.inv_order_var,2:simul_length_p1);
           simulations.third = simulation_third(oo_.dr.inv_order_var,2:simul_length_p1);
           simulations.constant = oo_.dr.ys + 0.5*ghs2_nlma(oo_.dr.inv_order_var,:);
           simulations.total = simulations.third +simulations.first_sigma_2 + simulations.second + simulations.first...
                               +repmat( simulations.constant,[1 simul_length] );
      end
     % Lan and Meyer-gohde's third order solution with a deterministic past
       if strcmp(pruning_type, 'lan_meyer-gohde_det')

          %E=oo_.exo_simul';
          E = shock_sequence; % -> HL, Sept 01, 2014
          simulation_first(:,1)=oo_.dr.ghu*E(:,1);
          exe=alt_kron(E(:,1),E(:,1));
          simulation_second(:,1)=(1/2)*(oo_.dr.ghuu*exe);
          simulation_first_sigma_2(:,1)=(1/2)*(oo_.dr.ghuss*E(:,1));
          simulation_third(:,1)=(1/6)*oo_.dr.ghuuu*alt_kron(E(:,1),exe);
          for t=2:simul_length
              simulation_first(:,t)=oo_.dr.ghx*simulation_first(select_state,t-1)+oo_.dr.ghu*E(:,t);
              exe=alt_kron(E(:,t),E(:,t));
              sxe=alt_kron(simulation_first(select_state,t-1),E(:,t));
              sxs=alt_kron(simulation_first(select_state,t-1),simulation_first(select_state,t-1));
              simulation_second(:,t)= oo_.dr.ghx*simulation_second(select_state,t-1)...
                                      +(1/2)*(oo_.dr.ghxx*sxs+2*oo_.dr.ghxu*sxe+oo_.dr.ghuu*exe);
              simulation_first_sigma_2(:,t)=oo_.dr.ghx*simulation_first_sigma_2(select_state,t-1)+(1/2)*(oo_.dr.ghxss*simulation_first(select_state,t-1)+oo_.dr.ghuss*E(:,t));
              simulation_third(:,t)=oo_.dr.ghx*simulation_third(select_state,t-1)...
                                     +(1/6)*(oo_.dr.ghxxx*alt_kron(simulation_first(select_state,t-1),sxs)+oo_.dr.ghuuu*alt_kron(E(:,t),exe))...
                                     +(1/2)*(oo_.dr.ghxxu*alt_kron(sxs,E(:,t))+oo_.dr.ghxuu*alt_kron(simulation_first(select_state,t-1),exe))...
                                     +oo_.dr.ghxx*alt_kron(simulation_second(select_state,t-1),simulation_first(select_state,t-1))+oo_.dr.ghxu*alt_kron(simulation_second(select_state,t-1),E(:,t));
          end
          simulations.first=simulation_first(oo_.dr.inv_order_var,:);
          simulations.second=simulation_second(oo_.dr.inv_order_var,:);
          simulations.first_sigma_2=simulation_first_sigma_2(oo_.dr.inv_order_var,:);
          simulations.third=simulation_third(oo_.dr.inv_order_var,:);
          simulations.total=simulations.third+simulations.first_sigma_2...
                            +simulations.second+simulations.first+repmat( oo_.dr.ys + 0.5*oo_.dr.ghs2(oo_.dr.inv_order_var,:),[1 simul_length] );
       end      
      
  end
