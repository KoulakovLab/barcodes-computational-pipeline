%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeGiN 'registercl_v012'                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
function [ REGISTEREDcl ] = registercl_v012( bigCL,smallCL, plot )
%%
% Input Variables
  bCluster = bigCL;        % Big Cluster
  sCluster = smallCL;      % Small Cluster needs some massaging

  bClsz = size( bCluster );
  sClsz = size( sCluster );
%  sCluster = [ sCluster; zeros( (bClsz(1)-sClsz(1)),sClsz(2) ) ];

  bCluster = bigCL(:);     % Big Cluster
  sCluster = sCluster(:);  % Small Cluster

  plt = plot;              % whether or not to plot the results
%
% Code Variables
  values = ones( bClsz(1)*bClsz(2),1 );   % Any elements of 'values' that have duplicate values of 'bCluster' and 'sCluster' are added together.

  nlls = find( sCluster==0 );
  values(nlls) = 0;   % Any elements of 'values' that are zero are ignored, along with the corresponding values of 'bCluster' and 'sCluster'.
  sCluster(nlls) = 1; % See null 'values' above (line 23).

  nmbrBCL = max( bCluster );   % number of clusters in bCluster
  nmbrSCL = max( sCluster );   % number of clusters in sCluster
%
% Cluster Registration
  regcl = sparse( bCluster,sCluster,values, nmbrBCL,nmbrSCL );   % ~1s

  REGISTEREDcl = regcl;   % Returning the values of the registered clusters (output)
%
% Plotting the results...
  if ( plt == 1 )
    figure; hold on;
      spy( REGISTEREDcl );
      % axis image;
      ylabel('Big Clusters'); xlabel('Small Clusters');
      ttl3 = sprintf('Registered: %d by %d clusters', nmbrBCL,nmbrSCL);
      title( ttl3 );
    hold off;
  end
%
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eNd 'registercl_v012'                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





































