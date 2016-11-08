delete('amg_matlab.out');
diary('amg_matlab.out');
diary on

tol = input('target convergence factor [0.5]: ');
if length(tol) == 0; tol = 0.5; end
fprintf(1,'tol = %g\n',tol);

ctol = input('coarsening aggressiveness [0.7]: ');
if length(ctol) == 0; ctol = 0.7; end
fprintf(1,'ctol = %g\n',ctol);

stol = input('sparsification tolerance [1e-4]: ');
if length(stol) == 0; stol = 1e-4; end
fprintf(1,'stol = %g\n',stol);

compute_cf = input('compute contraction factor (not necessary), yes - no [default]: ', 's');
if (max(strcmp(compute_cf, {'yes' 'y' 'YES' 'Y'})) ~= 1)
    compute_cf = 'no'; 
else
    compute_cf = 'yes';
end
fprintf(1,'%s\n',compute_cf);

global data A nullspace

tstart = tic;
[id A] = amg_import();
timport = toc(tstart);
[data, tcrs, tsmooth, tintp] = amg_setup(A,'new',ctol, tol, stol);
Abottom = full(data.A{length(data.n)});
fprintf(1,'A_%d = %g (',length(data.n),Abottom);
if Abottom<1e-9
  nullspace = 1;
else
	nullspace = 0;
	fprintf(1,'no ');
end
fprintf(1,'nullspace)\n');
tic;
amg_export(data,id,nullspace);
texport = toc;

% Compute contraction factor only if asked explicitely
if (strcmp(compute_cf, 'yes') == 1)
    fprintf(1,'Done!\n\nComputing largest eigenvalues of iteration matrix:\n');
    lambda = eigs(@Efun,length(A(:,1)))
    if nullspace>0
	    rho = abs(lambda(2));
    else
        rho = abs(lambda(1));
    end
	fprintf(1,'Error contraction factor: %g\n', rho);
end
ttot = toc(tstart);

fprintf(1,'Timers:\n');
fprintf(1,' Total time = %g seconds\n', ttot);
fprintf(1,' Import time = %g seconds\n', timport);
fprintf(1,' Coarsening time = %g seconds\n', tcrs);
fprintf(1,' Smoother time = %g seconds\n', tsmooth);
fprintf(1,' Interpolation time = %g seconds\n', tintp);
fprintf(1,' Export time = %g seconds\n', texport);
trest = ttot-timport-tcrs-tsmooth-tintp-texport;
fprintf(1,' Rest = %g seconds\n', trest);

diary off