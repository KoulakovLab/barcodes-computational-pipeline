%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeGiN 'pipeline_v002'                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%
function [ STATISTICS ] = pipeline_v002( ConnMat,CORA,CORD,Cj,ClstrC,uC1clstrmp , OUTPUT )
%%
% Input Variables
  C1 = ConnMat;
  CORA1 = CORA;
  CORD1 = CORD;
  C1j = Cj;
  ClusterC1 = ClstrC;
  uC1clmA = uC1clstrmp;
  STATISTICS = OUTPUT;   % sparse matrix, size( ClusterC1 )
% Pipeline
  % random acceptor selection
  nccptrs = size( C1,2 );
  temp = randperm( nccptrs );
  randina_11 = temp(1:2:end);
  randina_12 = temp(2:2:end);
  rA11 = C1(:,randina_11);   % random ACCEPTORS 11
  rA12 = C1(:,randina_12);   % random ACCEPTORS 12
  % Clustering
  [ClusterRA11,rA11shuff, rA11i,rA11j,rA11clm,rA11clmmod] = lumpclusters_v006( rA11,CORA1(randina_11, randina_11),CORD1, 0.05,0 );   % ~5s
  [ClusterRA12,rA12shuff, rA12i,rA12j,rA12clm,rA12clmmod] = lumpclusters_v006( rA12,CORA1(randina_12, randina_12),CORD1, 0.05,0 );   % ~4s
  % Unique clusters
%  ttuC1clmA  = uniqueclustermap_v001(   C1clm' );
%  uC1clmA  = ttuC1clmA';
  tturA11clm = uniqueclustermap_v001( rA11clm' );
  tturA12clm = uniqueclustermap_v001( rA12clm' );
  urA11clm  = tturA11clm';
  urA12clm  = tturA12clm';
  % Registration
  invRA11j = zeros( size(rA11j) );              % maps C1shuff into rA11shuff
  for cntr=1:length(rA11j)
    temp = randina_11(rA11j);                   % length(find( C1( rA11i, randina_11(rA11j)) - C1(C1i, C1j(invRA11j)) ))  ==  0  ==  length(find( rA11shuff - C1shuff(:, invRA11j) ))
    invRA11j(cntr) = find( C1j==temp(cntr) );   % C1j(invRA11j)' = randina_11(rA11j)
  end
  invRA12j = zeros( size(rA12j) );              % maps C1shuff into rA12shuff
  for cntr=1:length(rA12j)
    temp = randina_12(rA12j);                   % length(find( C1( rA12i, randina_12(rA12j)) - C1(C1i, C1j(invRA12j)) ))  ==  0  ==  length(find( rA12shuff - C1shuff(:, invRA12j) ))
    invRA12j(cntr) = find( C1j==temp(cntr) );   % C1j(invRA12j)' = randina_12(rA12j)
  end
  burA11clm = zeros( size(uC1clmA) );
  burA12clm = zeros( size(uC1clmA) );
  burA11clm(:, invRA11j) = urA11clm;   % Compare Clusters: C1shuff & rA11shuff :: C1shuff(:,invRA11j) = rA11shuff
  burA12clm(:, invRA12j) = urA12clm;   % Compare Clusters: C1shuff & rA12shuff :: C1shuff(:,invRA12j) = rA12shuff
  ccA11reg = registercl_v012( uC1clmA,burA11clm, 0 );
  ccA12reg = registercl_v012( uC1clmA,burA12clm, 0 );
  % Alignment
  nBIGcl = size( ccA11reg,1 );
  nS1cl  = size( ccA11reg,2 );
  nS2cl  = size( ccA12reg,2 );
  [valtmp1, indtmp1] = max( ccA11reg , [] , 1 );
  [valtmp2, indtmp2] = max( ccA12reg , [] , 1 );
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
  cocoa = ClusterRA11';
  cacao = ClusterRA12';
  cocoa = cocoa(mtchS1);
  cacao = cacao(mtchS2);
  ntrsct = intersect( find(cocoa), find(cacao) );
  STATISTICS( jjj(ntrsct) ) = STATISTICS( jjj(ntrsct) ) + 1;
%
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eNd 'pipeline_v002'                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





































