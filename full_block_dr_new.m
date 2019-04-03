function oo_=full_block_dr_new(oo_,M_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% full_block_dr_new.m
%
% This file rearranges the third-order coeffiecient oo_.dr.g_3 from Dynare 
% into blocks oo_.dr.ghxxx, oo_.dr.ghxxu, oo_.dr.ghxuu, oo_.dr.ghuuu,
% oo_.dr.ghxss, and oo_.dr.ghuss analogously to the presentation of second
% order coefficients as oo_.dr.g_2 as oo_.dr.ghxx, oo_.dr.ghxu, oo_.dr.ghuu,
% and oo_.dr.ghs2
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


if options_.order>=3 && options_.pruning == 0 %make sure the solution is a third-order unpruned solution (otherwise, everything is already provided by Dynare)
%% Recover the certainty equivalent decision rule   
    options_.order = 1;
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;%not always set by dynare for higher order perturbations
    end
    exo_simul = oo_.exo_simul;
    oo_.exo_simul = zeros(M_.maximum_lag+2,M_.exo_nbr);% these two lines are needed to prevent dynare from using the simulation to evaluate the derivatives
    exo_det_simul = oo_.exo_det_simul;
    oo_.exo_det_simul = [];
    [numeric_version] = return_dynare_version(dynare_version);
    if numeric_version<4.3 %Starting in Dynare 4.3.0, dr1 is removed. 
        [dr,~,~,~,~] = dr1(oo_.dr,0,M_,options_,oo_);
    elseif numeric_version>=4.3 %replace with stochastic_solvers in 4.3
        [dr,~] = stochastic_solvers(oo_.dr,0,M_,options_,oo_);   
    else
        disp('Error, no certainty equivalent dr calculated');
    end
    assert( numeric_version >= 4.4 );

    % nstatic = M_.nstatic;
    nspred = M_.nspred; % note M_.nspred = M_.npred+M_.nboth;
    % nboth = M_.nboth;
    % nfwrd = M_.nfwrd;

    oo_.dr.ghx=dr.ghx;
    oo_.dr.ghu=dr.ghu;
    oo_.exo_simul = exo_simul;
    oo_.exo_det_simul = exo_det_simul;
    %[oo_.dr,~,~,~,~] = dr1(oo_.dr,0,M_,options_,oo_);
    %[oo_.dr,~] =stochastic_solvers(oo_.dr,0,M_,options_,oo_);
    options_.order=3;
%% First unravel the folded decision rule
    b=zeros(nspred+M_.exo_nbr,(nspred+M_.exo_nbr)^2);
    a=tril(ones(nspred+M_.exo_nbr));
    i=find(a);
    a(i)=1:length(i);
    b(:,1:nspred+M_.exo_nbr)=a+a'- diag(diag(a));
    for j=2:nspred+M_.exo_nbr
        a = tril(ones(nspred+M_.exo_nbr-j+1));
        i = find(a);
        a(i) = (1:length(i))+max(max(b));
         b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))=[b(:,(1:(j-1))*(nspred+M_.exo_nbr)-(nspred+M_.exo_nbr-j)),[zeros(j-1,nspred+M_.exo_nbr-j+1);a]];
        b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))=tril(b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr)));
        b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))=  b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))+ b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))'- diag(diag(b(:,(j-1)*(nspred+M_.exo_nbr)+1:j*(nspred+M_.exo_nbr))));
    end
    G_3_select=b(:);
    G3=zeros(M_.endo_nbr,(nspred+M_.exo_nbr)^3);
    G3(:,G_3_select>0)=oo_.dr.g_3(:,G_3_select(G_3_select>0));
%% Then convert the three-fold Kronecker product to a three-fold block Kronecker product with blocks x_{t-1} and u_t [see Koning, Neudecker, and Wansbeek (1991)]
    temp=0:(nspred+M_.exo_nbr)^2:(nspred+M_.exo_nbr)^2*(nspred+M_.exo_nbr-1);
    G3=G3(:,repmat([vec(reshape(1:(nspred*(nspred+M_.exo_nbr)), (nspred+M_.exo_nbr), nspred)');vec(reshape((nspred*(nspred+M_.exo_nbr))+1:(nspred*(nspred+M_.exo_nbr))+(M_.exo_nbr*(nspred+M_.exo_nbr)),(nspred+M_.exo_nbr), M_.exo_nbr)')],[(nspred+M_.exo_nbr),1])+vec(temp(ones(1,(nspred+M_.exo_nbr)^2),:)));
    %G3=G3(:,repmat([vec(reshape(1:(nspred*(nspred+M_.exo_nbr)), (nspred+M_.exo_nbr), nspred)');vec(reshape((nspred*(nspred+M_.exo_nbr))+1:(nspred*(nspred+M_.exo_nbr))+(M_.exo_nbr*(nspred+M_.exo_nbr)),(nspred+M_.exo_nbr), M_.exo_nbr)')],[(nspred+M_.exo_nbr),1])+vec(repmat(0:(nspred+M_.exo_nbr)^2:(nspred+M_.exo_nbr)^2*(nspred+M_.exo_nbr-1),[(nspred+M_.exo_nbr)^2 1])));
    G3=G3(:,[vec(reshape(1:(nspred*(nspred+M_.exo_nbr)^2), (nspred+M_.exo_nbr)^2, nspred)');vec(reshape((nspred*(nspred+M_.exo_nbr)^2)+1:(nspred*(nspred+M_.exo_nbr)^2)+(M_.exo_nbr*(nspred+M_.exo_nbr)^2),(nspred+M_.exo_nbr)^2, M_.exo_nbr)')]);
%% Then pick off the desired blocks for third order terms and calculate the uncertainty correction of the first-order terms
    oo_.dr.ghxxx=6*G3(:,1:nspred^3);
    oo_.dr.ghxxu=6*G3(:,1+nspred^3+2*M_.exo_nbr*nspred^2+M_.exo_nbr^2*nspred:nspred^3+3*M_.exo_nbr*nspred^2+M_.exo_nbr^2*nspred);
    oo_.dr.ghxuu=6*G3(:,1+nspred^3+3*M_.exo_nbr*nspred^2+2*M_.exo_nbr^2*nspred:nspred^3+3*M_.exo_nbr*nspred^2+3*M_.exo_nbr^2*nspred);
    oo_.dr.ghuuu=6*G3(:,1+nspred^3+3*M_.exo_nbr*nspred^2+3*M_.exo_nbr^2*nspred:end);
    %oo_.dr.ghs2x=2*(oo_.dr.g_1(:,1:nspred)-oo_.dr.ghx); 
    %oo_.dr.ghs2u=2*(oo_.dr.g_1(:,1+nspred:end)-oo_.dr.ghu);
    oo_.dr.ghxss = 2*(oo_.dr.g_1(:,1:nspred)-oo_.dr.ghx); %<---- Notation change: use ghxss to align with Dynare HL, March 2014
    oo_.dr.ghuss=2*(oo_.dr.g_1(:,1+nspred:end)-oo_.dr.ghu); %<---- Notation change: use ghuss to align with Dynare HL, March 2014
end
