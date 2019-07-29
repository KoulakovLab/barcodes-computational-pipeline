%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeGiN 'pipeline_v001'                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
function [ STATISTICS ] = pipeline_v001( ConnMat,CORA,CORD,Ci,ClstrC,uC1clstrmp , OUTPUT )
%%
% Input Variables
  C1 = ConnMat;
  CORA1 = CORA;
  CORD1 = CORD;
  C1i = Ci;
  ClusterC1 = ClstrC;
  uC1clm = uC1clstrmp;
  STATISTICS = OUTPUT;   % sparse matrix, size( ClusterC1 )
% Pipeline
  % random donor selection
  ndnrs = size( C1,1 );      % number of donors
  temp = randperm( ndnrs )';
  randind_01 = temp(1:2:end);
  randind_02 = temp(2:2:end);
  rD01 = C1(randind_01,:);   % random DONORS 01
  rD02 = C1(randind_02,:);   % random DONORS 02
  % Clustering
  [ClusterRD01,rD01shuff, rD01i,rD01j,rD01clm,rD01clmmod] = lumpclusters_v006( rD01,CORA1,CORD1(randind_01, randind_01), 0.05,0 );
  [ClusterRD02,rD02shuff, rD02i,rD02j,rD02clm,rD02clmmod] = lumpclusters_v006( rD02,CORA1,CORD1(randind_02, randind_02), 0.05,0 );
  % Unique clusters
  urD01clm = uniqueclustermap_v001( rD01clm );
  urD02clm = uniqueclustermap_v001( rD02clm );
  % Registration
  invRD01i = zeros( size(rD01i) );              % maps C1shuff into rD01shuff
  for cntr=1:length(rD01i)
    temp = randind_01(rD01i);                   % length(find( C1( randind_01(rD01i), rD01j) - C1(C1i(invRD01i), C1j) ))  ==  0  ==  length(find( rD01shuff - C1shuff(invRD01i, :) ))
    invRD01i(cntr) = find( C1i==temp(cntr) );   % C1i(invRD01i) = randind_01(rD01i)
  end
  invRD02i = zeros( size(rD02i) );              % maps C1shuff into rD02shuff
  for cntr=1:length(rD02i)
    temp = randind_02(rD02i);                   % length(find( C1( randind_02(rD02i), rD02j) - C1(C1i(invRD02i), C1j) ))  ==  0  ==  length(find( rD02shuff - C1shuff(invRD02i, :) ))
    invRD02i(cntr) = find( C1i==temp(cntr) );   % C1i(invRD02i) = randind_02(rD02i)
  end
  burD01clm = zeros( size(uC1clm) );
  burD02clm = zeros( size(uC1clm) );
  burD01clm( invRD01i, :) = urD01clm;   % Compare Clusters: C1shuff & rD01shuff :: C1shuff(invRD01i,:) = rD01shuff
  burD02clm( invRD02i, :) = urD02clm;   % Compare Clusters: C1shuff & rD02shuff :: C1shuff(invRD02i,:) = rD02shuff
  ccD01reg = registercl_v012( uC1clm,burD01clm, 0 );
  ccD02reg = registercl_v012( uC1clm,burD02clm, 0 );
  % Alignment
  nBIGcl = size( ccD01reg,1 );
  nS1cl  = size( ccD01reg,2 );
  nS2cl  = size( ccD02reg,2 );
  [valtmp1, indtmp1] = max( ccD01reg , [] , 1 );
  [valtmp2, indtmp2] = max( ccD02reg , [] , 1 );
  mtchS1 = [];
  mtchS2 = [];
  jjj    = [];
  parfor cntr=1:nBIGcl
    s1ind = find( indtmp1==cntr );
    s2ind = find( indtmp2==cntr );
    if ( length(s1ind) > length(s2ind) )
	  jjj = [ jjj; cntr * ones( length(s2ind),1 ) ];   % mapping into BIG clusters
	  iii = [1:length(s2ind)];
      mtchS1 = [mtchS1 ; s1ind(iii)' ];
      mtchS2 = [mtchS2 ; s2ind' ];
	else
	  jjj = [ jjj; cntr * ones( length(s1ind),1 ) ];   % mapping into BIG clusters
	  iii = [1:length(s1ind)];
      mtchS1 = [mtchS1 ; s1ind' ];
      mtchS2 = [mtchS2 ; s2ind(iii)' ];
    end
  end
  % Pipeline Output
  cocoa = ClusterRD01(mtchS1);
  cacao = ClusterRD02(mtchS2);
  ntrsct = intersect( find(cocoa), find(cacao) );
  STATISTICS( jjj(ntrsct) ) = STATISTICS( jjj(ntrsct) ) + 1;
%
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eNd 'pipeline_v001'                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





































