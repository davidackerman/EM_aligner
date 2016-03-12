function [x2, R] = solve_AxB(K,Lm,options,d)
%% SOLVE
tic
R = []; 
%     in case of iterative methods you can use semilogy(R(:,1),R(:,2),'-o');
%     xlabel('Iteration Number');
%     ylabel('Relative Residual');
%     to view iterative method progress;
if strcmp(options.solver,'backslash--noreg')
    if options.verbose,disp('------------ Performing backslash -- no reg');end
    x2 = K\Lm;
end

%%%%%%%%%%%%%%%%% preconditioning
if ~strcmp(options.solver,'backslash') ...
        && ~strcmp(options.solver,'backslash--noreg') ...
        && ~strcmp(options.solver,'lsqlin')
    if options.use_ilu && isempty(options.L2)
        L2 = [];U2 = [];
        if options.verbose,disp('Calculating preconditioner');end
        %q = colamd(K);
        [L2,U2] = ilu(K,struct('type',options.ilu_type,'droptol',options.ilu_droptol, 'udiag', options.ilu_udiag));
    else
        L2 = options.L2;
        U2 = options.U2;
    end
    iL2 = L2;
    iU2 = U2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = -999;
if strcmp(options.solver, 'lsqlin')
    if options.verbose,disp('------------ Performing backslash x2 = K\Lm solution ');end
    x2 = lsqlin(A, b, B, d);
    
    
elseif strcmp(options.solver, 'backslash')
    if options.verbose,disp('------------ Performing backslash x2 = K\Lm solution ');end
    x2 = K\Lm;
    
elseif strcmp(options.solver, 'gmres')
    if options.verbose,disp('------------ Performing gmres');end
    %x1 = K\Lm; % sosi
    %[x2,flag,relres,iter,resvec] = gmres(K,Lm,options.restart,options.tol, options.maxit, L2, U2, x1);
    [x2,flag,relres,iter,resvec] = gmres(K,Lm,options.restart,options.tol, options.maxit,L2, U2, d);
    disp(flag);
    % flag = 0
    % gmres converged to the desired tolerance tol within maxit outer iterations.
    % flag = 1
    % gmres iterated maxit times but did not converge.
    % flag = 2
    % Preconditioner M was ill-conditioned.
    % flag = 3
    % gmres stagnated. (Two consecutive iterates were the same.)
elseif strcmp(options.solver, 'tfqmr')
    if options.verbose,disp('------------ Performing tfqmr');end
    [x2,flag,relres,iter,resvec] = tfqmr(K,Lm,options.tol, options.maxit,L2,U2,d);
    
elseif strcmp(options.solver, 'bicg')
    if options.verbose,disp('------------ Performing biconjugate gradients');end
    [x2,flag,relres,iter,resvec] = bicg(K,Lm,options.tol, options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'bicgstab')
    if options.verbose,disp('------------ Performing biconjugate gradients stabelized');end
    [x2,flag,relres,iter,resvec] = bicgstab(K,Lm,options.tol, options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'lsqr')
    if options.verbose,disp('------------ Performing lsqr');end
    [x2,flag,relres,iter,resvec] = lsqr(A,b,options.tol,options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'symmlq')
    if options.verbose,disp('------------ Performing symmlq');end
    [x2,flag,relres,iter,resvec] = symmlq(K,Lm,options.tol, options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'cgs')
    if options.verbose,disp('------------ Performing cgs');end
    [x2,flag,relres,iter,resvec] = cgs(K,Lm,options.tol, options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'pcg')
    if options.verbose,disp('------------ Performing pcg with incomplete cholesky preconditioning');end
    [x2,flag,relres,iter,resvec] = pcg(K,Lm,options.tol, options.maxit, L2, U2, d);
    
elseif strcmp(options.solver, 'minres')
    if options.verbose,disp('------------ Performing minres');end
    [x2,flag,relres,iter,resvec] = minres(K,Lm,options.tol, options.maxit, L2, U2, d);
end
if options.verbose
disp('-----------    Time ----------');
toc
disp('------------------------------');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sosi
% x1 = K\Lm;
% disp(num2str([d(100:130) x1(100:130) x2(100:130)]));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.verbose,
    disp('Flag (default is -999 for backslash)');disp(flag);
    disp('Done!');
end
if ~strcmp(options.solver,'backslash')...
        && ~strcmp(options.solver,'backslash--noreg')
    c1 = [0:numel(resvec)-1];
    c2 = resvec/norm(Lm);
    R  = [c1(:) c2(:)];
else 
    R = 0;
end

if options.debug && ~strcmp(options.solver,'backslash')...
        && ~strcmp(options.solver,'backslash--noreg')
    
    semilogy(R(:,1),R(:,2),'-o');
    xlabel('Iteration Number');
    ylabel('Relative Residual');
    warning on;
    drawnow;
    pause(3);
end

% % the svd solution
% [U, S, V] = svds(K);
% invS = 1./S;
% invS(isinf(invS)) = 0;
% x2 = U*invS*V'*Lm;


x2 = real(x2);


%% Error  ------- sosi -------- error progress as function of lambda needs to be done on A and b and NOT on K and Lm
% % % err = norm(K*d-Lm);
% % err = norm(A*d-b);
% % if options.verbose,disp('Error (pre-optimization  -- Ax=b):');
% %     disp(num2str(err));
% % end
% % % % err = norm(K*x2-Lm);
% % err = norm(A*x2-b);
% % if options.verbose,disp('Error (post-optimization -- Ax=b):');
% %     disp(num2str(err));
% % end


% err1 = norm(K*x1-Lm);
% if options.verbose,
%     disp('Error (optimal with backslash):');
%     disp(num2str(err1));
%     disp('Ratio error to optimal error');
%     disp(num2str(abs((err)/err1)));
% end
