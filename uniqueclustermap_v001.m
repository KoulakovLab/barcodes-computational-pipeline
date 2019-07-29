%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeGiN                                                                        %
% uniqueclustermap_v001.m                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ UniqueCLusterMap ] = uniqueclustermap_v001( CLusterMap )
%%
clm = CLusterMap;
szclm = size( clm );
uniqueclm = zeros( szclm(1) , szclm(2) );
%
uniqueclm(:,1) = clm(:,1) - 1;   % the first column receives the first list of clusters

for cntr=2:szclm(2)
  tempBFR = clm(:,cntr-1);
  tempCRRNT = clm(:,cntr);  

  if ( all((tempCRRNT - tempBFR) == 0) )
    uniqueclm(:,cntr) = uniqueclm(:,cntr-1);
  else
    tempSHIFT = uniqueclm(:,cntr-1);
    uniqueclm(:,cntr) = tempSHIFT(end) + uniqueclm(:,1);
  end
end
%
UniqueCLusterMap = uniqueclm;
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eNd                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


