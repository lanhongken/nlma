function [ deflect, M, oo ] = compute_deflected_linear_approximation( M, options, oo, ModeRSS1OrMean2 )

% Copyright: Alexander Meyer-Gohde 2016, Tom Holden 2018 (only slight modifications from TH)

% You are free to use/modify/redistribute this program so long as original authorship credit is given and you in no way impinge on its free distribution.
% This software is provided as is with no guarantees of any kind.

if ModeRSS1OrMean2 > 2
    disp('Incompatible mode')
    deflect=[];
    return
end

[numeric_version] = return_dynare_version(dynare_version);
assert( numeric_version >= 4.4 );

nstatic = M.nstatic;
nspred = M.nspred; % note M_.nspred = M_.npred+M_.nboth;
% nboth = M_.nboth;
nfwrd = M.nfwrd;

if options.order>=3
    if options.pruning == 0
        oo = full_block_dr_new(oo,M,options);
    end
end
if options.order>=2
    if isempty(options.qz_criterium)==1
        options.qz_criterium=1+1e-6;
    end
    state_var =  lyapunov_symm(oo.dr.ghx(nstatic+1:nstatic+nspred,:),oo.dr.ghu(nstatic+1:nstatic+nspred,:)*M.Sigma_e*oo.dr.ghu(nstatic+1:nstatic+nspred,:)',2-options.qz_criterium,1e-12,options.lyapunov_complex_threshold);
end

if options.order>=1
    deflect.y=oo.dr.ys(oo.dr.order_var);
    deflect.y_x=oo.dr.ghx;
    deflect.y_u=oo.dr.ghu;
end
if options.order>=2
    if ModeRSS1OrMean2==2
        x_2=(eye(nspred)-oo.dr.ghx(nstatic+1:nstatic+nspred,:))\(oo.dr.ghuu(nstatic+1:nstatic+nspred,:)*vec(M.Sigma_e)+oo.dr.ghxx(nstatic+1:nstatic+nspred,:)*vec(state_var)+oo.dr.ghs2(nstatic+1:nstatic+nspred,:));
        y_2=[oo.dr.ghx(1:nstatic,:)*x_2+oo.dr.ghuu(1:nstatic,:)*vec(M.Sigma_e)+oo.dr.ghxx(1:nstatic,:)*vec(state_var)+oo.dr.ghs2(1:nstatic,:);x_2;oo.dr.ghx(nstatic+nspred+1:M.endo_nbr,:)*x_2+oo.dr.ghuu(nstatic+nspred+1:M.endo_nbr,:)*vec(M.Sigma_e)+oo.dr.ghxx(nstatic+nspred+1:M.endo_nbr,:)*vec(state_var)+oo.dr.ghs2(nstatic+nspred+1:M.endo_nbr,:)];
        deflect.y=deflect.y+0.5*y_2;
    elseif ModeRSS1OrMean2==1
        x_2=(eye(nspred)-oo.dr.ghx(nstatic+1:nstatic+nspred,:))\(oo.dr.ghs2(nstatic+1:nstatic+nspred,:));
        y_2=[oo.dr.ghx(1:nstatic,:)*x_2+oo.dr.ghs2(1:nstatic,:);x_2;oo.dr.ghx(nstatic+nspred+1:M.endo_nbr,:)*x_2+oo.dr.ghs2(nstatic+nspred+1:M.endo_nbr,:)];
        deflect.x_2=x_2;
        deflect.y=deflect.y+0.5*y_2;
    end
end
if options.order==3 && ModeRSS1OrMean2 > 0
    deflect.y_x=deflect.y_x+0.5*(oo.dr.ghxx*kron(x_2,eye(nspred))+2*(oo.dr.g_1(:,1:nspred)-oo.dr.ghx));
    deflect.y_u=deflect.y_u+0.5*(oo.dr.ghxu*kron(x_2,eye(M.exo_nbr))+2*(oo.dr.g_1(:,1+nspred:end)-oo.dr.ghu));
end
deflect.y=deflect.y(oo.dr.inv_order_var);
deflect.y_y=[zeros(M.endo_nbr,nstatic),deflect.y_x,zeros(M.endo_nbr,nfwrd)];
deflect.y_y=deflect.y_y(oo.dr.inv_order_var,oo.dr.inv_order_var);
deflect.y_e=deflect.y_u(oo.dr.inv_order_var,:);
