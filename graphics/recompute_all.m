%load('../data/matrix_laplace2d_grid256.mat', 'A', 'y0');
%load('../data/matrix_advection_diffusion_grid256.mat', 'A', 'y0');
%load('../data/matrix_diffusion_grid256.mat', 'A', 'y0');

%name_list = ;

%grid_list = ;

for name = {'diffusion-2'}%{'diffusion', 'advection_diffusion', 'laplace2d'}
    for grid = {64, 128, 256}
        load(strcat('../data/matrix_', name{1}, '_grid', int2str(grid{1}), '.mat'), 'A', 'y0');
        ch = 0;
        m = 49;
        tol = 1e-16;
        tolY = 0.0;
        B = y0;
        E = speye(size(B, 1));
        LE = chol(E, 'lower');
        opts.tol=1e-2;
        s1=eigs(-A,E,1,'lm',opts);
        s2=eigs(-A,E,1,'sm',opts);
        [Z, resnorm, Zall] = rksm(A,E,LE,B,m,tol,s1,s2,ch,tolY);

        resnorm = reshape(resnorm, 1, numel(resnorm));
        save(strcat('rksm_', name{1}, '_grid',  int2str(grid{1}), '.mat'), 'Z', 'resnorm');

        [Z,resnorm]=kpik(A,E,LE,B,m,tol,tolY);
        resnorm = reshape(resnorm, 1, numel(resnorm));
        save(strcat('kpik_', name{1}, '_grid',  int2str(grid{1}), '.mat'), 'Z', 'resnorm');
        disp(strcat(name{1}, '_grid',  int2str(grid{1}), ' is handled.'))
    end
end
disp('Finish!')