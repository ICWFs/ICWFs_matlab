function [x0,y0]=mcsampling(PDaux,N_traj_aux,threshold)

%define the possible values for the (x,y) pair
PDaux = PDaux/sum(PDaux(:));

row_vals = [1:size(PDaux,1)]'*ones(1,size(PDaux,2));  %all x values
col_vals = ones(size(PDaux,1),1)*[1:size(PDaux,2)];  %all y values

%convert your 2D problem into a 1D problem
PDaux = PDaux(:);
PDaux(PDaux<threshold*max(PDaux(:))) = 1E-12*max(PDaux(:));
row_vals = row_vals(:);
col_vals = col_vals(:);

%calculate your fake 1D CDF, assumes sum(A(:))==1
CDF = cumsum(PDaux);%cumsum(PDaux); %remember, first term out of of cumsum is not zero

%because of the operation we're doing below (interp1 followed by ceil)
%we need the CDF to start at zero
CDF = [0; CDF(:)];

%generate random values
rand_vals = rand(N_traj_aux,1);  %spans zero to one

%look into CDF to see which index the rand val corresponds to
out_val = interp1(CDF,[0:1/(length(CDF)-1):1],rand_vals); %spans zero to one
ind = ceil(out_val*length(PDaux));

%using the inds, you can lookup each pair of values
ind(isnan(ind)) = [];

x0 = row_vals(ind);
y0 = col_vals(ind);