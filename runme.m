%%%
% 
%%%

clear
close all
load data

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Calibri')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultAxesFontWeight', 'bold')
set(0,'DefaultAxesLineWidth', 2)
% Change default text fonts.
set(0,'DefaultTextFontname', 'Calibri')
set(0,'DefaultTextFontSize', 12)
set(0,'DefaultTextFontWeight', 'bold')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline: try *3*: DONORS
% 17:55 Tuesday, November 11, 2014
tic;
  [ClusterC1,C1shuff, C1i,C1j,C1clm,C1clmmod] = lumpclusters_v006( C1,CORA1,CORD1, 0.05,0 );
toc;
%%%

szClstrC1 = size( ClusterC1 );
pvalueDONORS = sparse( szClstrC1(1),szClstrC1(2) );

tic;
  for cntr=1:100
    temp = pvalueDONORS;
    [ pvalueDONORS ] = pipeline_v001( C1,CORA1,CORD1,C1i,ClusterC1,uC1clm , temp );
  end
toc;   % ~20.0114 hours :: each iteration = ~3.7096 min  (using PARFOR = 1.9571 min/iteration)
clear ans cntr temp szClstrC1

save('data.mat', 'pvalueDONORS','-append');

figure, hold on;
  spy( pvalueDONORS );
  xlabel('Acceptors'); ylabel('Donors');
  title(sprintf('Donors pValue, 100 iterations, %d nonzero elements', length(find(pvalueDONORS))));
hold off;


tic;
  parfor cntr=1:100
    temp = pipeline_v001( C1,CORA1,CORD1,C1i,ClusterC1,uC1clm , sparse( 456,786 ) );
    pvalueDONORSv02 = pvalueDONORSv02 + temp;
  end
toc;   % ~3.5787 hours :: each iteration = 1.9571 min
clear ans cntr temp szClstrC1
save('data.mat', 'pvalueDONORSv02','-append');


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
pvalueDONORSv03 = sparse( 456,786 );
pvalueACCEPTORSv03 = sparse( 456,786 );
tempD = sparse( 456,786 );
tempA = sparse( 456,786 );

tic;
  parfor cntr=1:300
    tempD = pipeline_v001( C1,CORA1,CORD1,C1i,ClusterC1,uC1clm , sparse( 456,786 ) );
    tempA = pipeline_v002( C1,CORA1,CORD1,C1j,ClusterC1,uC1clmA , sparse( 456,786 ) );
    pvalueDONORSv03 = pvalueDONORSv03 + tempD;
    pvalueACCEPTORSv03 = pvalueACCEPTORSv03 + tempA;
  end
toc;   % ~3.5787 hours
clear ans cntr tempD tempA szClstrC1
save('data.mat', 'pvalueDONORSv03','pvalueACCEPTORSv03','-append');

%%%

pval = [0:0.001:1]';
abnDNRS   = zeros( size(pval) );
abnCCPTRS = zeros( size(pval) );
abnTTL    = zeros( size(pval) );
tic;
for cntr=1:length(pval)
  abnDNRS(cntr)   = length(find( pvalDNRSttl   >= pval(cntr) ));
  abnCCPTRS(cntr) = length(find( pvalCCPTRSttl >= pval(cntr) ));
  abnTTL(cntr)    = length(find( pvalTTL       >= pval(cntr) ));
end
toc;   % 1.831074 seconds
clear ans cntr

figure, hold on;
  plot( pval(2:end),abnDNRS(2:end)   , 'r' );
  plot( pval(2:end),abnCCPTRS(2:end) , 'b' );
  plot( pval(2:end),abnTTL(2:end)    , 'g' );
  xlabel('p-value'); ylabel('Abundance (number of clusters >= pvalue)');
hold off;

save('data.mat', 'pvalDNRSttl','pvalCCPTRSttl','pvalTTL','-append');
save('data.mat', 'pval','abnDNRS','abnCCPTRS','abnTTL','-append');

figure, hold on;
  plot( pval(2:end),abnTTL(2:end) , 'b' , 'LineWidth',2 );
  plot( [0.16],[658] , 'r*','MarkerSize',5 );
  xlabel('p-value'); ylabel('Abundance (number of clusters >= pvalue)');
  title({'Number of Nodes','Donors + Acceptors; 500 iterations'});
hold off;




cnnstrDNRS   = zeros( size(pval) );
cnnstrCCPTRS = zeros( size(pval) );
cnnstrTTL    = zeros( size(pval) );
tic;
for cntr=1:length(pval)
  cnnstrDNRS(cntr)   = sum(ClusterC1(find( pvalDNRSttl   >= pval(cntr) )));
  cnnstrCCPTRS(cntr) = sum(ClusterC1(find( pvalCCPTRSttl >= pval(cntr) )));
  cnnstrTTL(cntr)    = sum(ClusterC1(find( pvalTTL       >= pval(cntr) )));
end
toc;   % 2.492616 seconds
clear ans cntr

figure, hold on;
  plot( pval(2:end),cnnstrDNRS(2:end)   , 'r' );
  plot( pval(2:end),cnnstrCCPTRS(2:end) , 'b' );
  plot( pval(2:end),cnnstrTTL(2:end)    , 'g' );
  xlabel('p-value'); ylabel('Connection Strength');
hold off;

save('data.mat', 'cnnstrDNRS','cnnstrCCPTRS','cnnstrTTL','-append');

figure, hold on;
  plot( pval(2:end),cnnstrTTL(2:end) , 'b' , 'LineWidth',2 );
  plot( [0.142 ; 0.16],[7.619e5 ; 2.395e5] , 'r*','MarkerSize',5 );
  xlabel('p-value'); ylabel('Sum Connection Strength');
  title({'Sum Conn Strength for Clusters > p-value','Donors + Acceptors; 500 iterations'});
hold off;


figure, hold on;
  plot( abnTTL(2:end),cnnstrTTL(2:end) , 'b' , 'LineWidth',2 );
  plot( [abnTTL(2); abnTTL(end); 658],[cnnstrTTL(2); cnnstrTTL(end); 2.395e5] , 'r*','MarkerSize',5 );
  xlabel('Node Count'); ylabel('Sum Connection Strength');
  title({'Sum Conn Strength vs Node Count','Donors + Acceptors; 500 iterations'});
hold off;

figure, hold on;
  plot3( abnTTL(2:end),cnnstrTTL(2:end),pval(2:end) , 'b' , 'LineWidth',2 );
  xlabel('Node Count'); ylabel('Sum Connection Strength'); zlabel('p-Value');
  title({'Sum Conn Strength vs Node Count','Donors + Acceptors; 500 iterations'});
  grid on; view( [35 , 64] );
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline: try *3*: ACCEPTORS
% 10:22 Wednesday, November 12, 2014

szClstrC1 = size( ClusterC1 );
pvalueACCEPTORS = sparse( szClstrC1(1),szClstrC1(2) );

tic;
  for cntr=1:100
    temp = pvalueACCEPTORS;
    [ pvalueACCEPTORS ] = pipeline_v002( C1,CORA1,CORD1,C1j,ClusterC1,uC1clmA , temp );
  end
toc;   % ~21.4240 hours :: each iteration = ~12.8076 min
clear ans cntr temp szClstrC1

save('data.mat', 'pvalueACCEPTORS','-append');

figure, hold on;
  spy( pvalueACCEPTORS );
  xlabel('Acceptors'); ylabel('Donors');
  title(sprintf('Acceptors pValue, 100 iterations, %d nonzero elements', length(find(pvalueACCEPTORS))));
hold off;


tic;
  parfor cntr=1:100
    temp = pipeline_v002( C1,CORA1,CORD1,C1j,ClusterC1,uC1clmA , sparse( 456,786 ) );
    pvalueACCEPTORSv02 = pvalueACCEPTORSv02 + temp;
  end
toc;   % ~2.2374 hours :: each iteration = 1.9571 min
clear ans cntr temp szClstrC1
save('data.mat', 'pvalueACCEPTORSv02','-append');


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Venn diagram
length(find( (cocoa > 0) .* (cacao > 0) ))   % same as 'ntrsct' (INTERSECT) on line 211 above (jjj).
length(find( (cocoa > 0) .* (cacao == 0) ))
length(find( (cocoa == 0) .* (cacao > 0) ))
length(find( (cocoa == 0) .* (cacao == 0) ))

length(find( cocoa > 0 ))
length(find( cacao > 0 ))

%%%

pvalueDONORSv06    = sparse( 456,786 );
pvalueACCEPTORSv06 = sparse( 456,786 );
tempD = sparse( 456,786 );
tempA = sparse( 456,786 );

tic;
  parfor cntr=1:500
    tempD = pipeline_v001( C1,CORA1,CORD1,C1i,ClusterC1,uC1clm  , sparse( 456,786 ) );
    tempA = pipeline_v002( C1,CORA1,CORD1,C1j,ClusterC1,uC1clmA , sparse( 456,786 ) );
    pvalueDONORSv06    = pvalueDONORSv06    + tempD;
    pvalueACCEPTORSv06 = pvalueACCEPTORSv06 + tempA;
  end
toc;   % 3.1599 days
clear ans cntr tempD tempA szClstrC1
save('data.mat', 'pvalueDONORSv06','pvalueACCEPTORSv06','-append');

pvalDNRSttl   = pvalueDONORSv06 / max(max( pvalueDONORSv06 ));
pvalCCPTRSttl = pvalueACCEPTORSv06 / max(max( pvalueACCEPTORSv06 ));
pvalTTL = ( pvalDNRSttl + pvalCCPTRSttl ) / 2;

pval = [0:0.001:1]';
abnDNRS   = zeros( size(pval) );
abnCCPTRS = zeros( size(pval) );
abnTTL    = zeros( size(pval) );

tic;
parfor cntr=1:length(pval)
  abnDNRS(cntr)   = length(find( pvalDNRSttl   >= pval(cntr) ));
  abnCCPTRS(cntr) = length(find( pvalCCPTRSttl >= pval(cntr) ));
  abnTTL(cntr)    = length(find( pvalTTL       >= pval(cntr) ));
end
toc;   % 1.831074 seconds
clear ans cntr

cnnstrDNRS   = zeros( size(pval) );
cnnstrCCPTRS = zeros( size(pval) );
cnnstrTTL    = zeros( size(pval) );
tic;
parfor cntr=1:length(pval)
  cnnstrDNRS(cntr)   = sum(ClusterC1(find( pvalDNRSttl   >= pval(cntr) )));
  cnnstrCCPTRS(cntr) = sum(ClusterC1(find( pvalCCPTRSttl >= pval(cntr) )));
  cnnstrTTL(cntr)    = sum(ClusterC1(find( pvalTTL       >= pval(cntr) )));
end
toc;   % 2.492616 seconds
clear ans cntr

figure, hold on;
  plot( abnTTL(2:end),cnnstrTTL(2:end) , 'b' , 'LineWidth',2 );
%  plot( [abnTTL(2); abnTTL(end); 658],[cnnstrTTL(2); cnnstrTTL(end); 2.395e5] , 'r*','MarkerSize',5 );
  xlabel('Node Count'); ylabel('Sum Connection Strength');
  title({'Sum Conn Strength vs Node Count','Donors + Acceptors; 500 iterations'});
hold off;

figure, hold on;
  plot3( abnTTL(2:end),cnnstrTTL(2:end),pval(2:end) , 'b' , 'LineWidth',2 );
  xlabel('Node Count'); ylabel('Sum Connection Strength'); zlabel('p-Value');
  title({'Sum Conn Strength vs Node Count','Donors + Acceptors; 500 iterations'});
  grid on; view( [-5 , 34] );
hold off;


%%%

% Connection Matrices Plots...
cmap5 = [
        1 , 1      , 1 ;
        1 , 0.8431 , 0 ;
	    0 , 0.3922 , 0 ;
		0 , 0      , 1 ;
		1 , 0      , 0
];

cmap5p = [
        0.8275 , 0.8275 , 0.8275 ;
        1 , 0.8431 , 0 ;
	    0 , 0.3922 , 0 ;
		0 , 0      , 1 ;
		1 , 0      , 0
];

cmap10 = [
        1 , 1      , 1 ;
        1 , 1      , 1 ;
        1 , 0.8431 , 0 ;
        1 , 0.8431 , 0 ;
	    0 , 0.3922 , 0 ;
	    0 , 0.3922 , 0 ;
		0 , 0      , 1 ;
		0 , 0      , 1 ;
		1 , 0      , 0
		1 , 0      , 0
];

cmap1000 = [
         repmat( [1 , 1      , 1] , 200,1 ) ;
         repmat( [1 , 0.8431 , 1] , 200,1 ) ;
         repmat( [1 , 0.3922 , 1] , 200,1 ) ;
         repmat( [1 , 0      , 1] , 200,1 ) ;
         repmat( [1 , 0      , 0] , 200,1 ) ;
];



cocoa = (pvalueACCEPTORSv07 + pvalueDONORSv06) .* logical(ClusterC1);
szcc = size( cocoa );
figure, hold on;
  for cntrPVL=1:1000
    [ tempx , tempy ] = ind2sub( szcc , find( cocoa == cntrPVL ));
    if ( cntrPVL < 200 )
	  clr = [ 0.8275 , 0.8275 , 0.8275 ];   % grey
	elseif ( cntrPVL < 400 )
	  clr = [ 1 , 0.8431 , 0 ];   % yellow
	elseif ( cntrPVL < 600 )
	  clr = [ 0 , 0.3922 , 0 ];   % green
	elseif ( cntrPVL < 800 )
	  clr = [ 0 , 0      , 1 ];   % blue
	else
	  clr = [ 1 , 0      , 0 ];   % red
	end
    plot( tempy,tempx, 'LineStyle','none', 'Marker','.','MarkerFaceColor',clr,'MarkerEdgeColor',clr,'MarkerSize',3 );
%    plot( tempy,tempx, 'LineStyle','none', 'Marker','o','MarkerFaceColor',clr,'MarkerEdgeColor',clr,'MarkerSize',2 );
  end
  caxis( [0,1] ), colormap( cmap5p ); colorbar,
  axis ij, axis image,
  title({'Random Sampling ConnMat Overlap:', '1000 samples, 1201 nonzero elements'});
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
hold off;
clear ans szcc tempx tempy clr cntrPVL cntrCCPTR cntrDNR


pvalDNRSttl   = pvalueDONORSv06 / 500;
pvalCCPTRSttl = pvalueACCEPTORSv07 / 500;
pvalTTL       = 0.001 * (pvalueACCEPTORSv07 + pvalueDONORSv06) .* logical(ClusterC1);
pval = [0:0.001:1]';
abnDNRS   = zeros( size(pval) );
abnCCPTRS = zeros( size(pval) );
abnTTL    = zeros( size(pval) );
tic;
for cntr=1:length(pval)
  abnDNRS(cntr)   = length(find( pvalDNRSttl   >= pval(cntr) ));
  abnCCPTRS(cntr) = length(find( pvalCCPTRSttl >= pval(cntr) ));
  abnTTL(cntr)    = length(find( pvalTTL       >= pval(cntr) ));
end
toc;   % 1.831074 seconds
clear ans cntr

figure, hold on;
  plot( pval(2:end),abnDNRS(2:end)   , 'r' , 'LineWidth',2 );
  plot( pval(2:end),abnCCPTRS(2:end) , 'g' , 'LineWidth',2 );
  plot( pval(2:end),abnTTL(2:end)    , 'b' , 'LineWidth',2 );
  xlabel('p-value'); ylabel('Abundance (number of clusters >= pvalue)');
  title({'Number of Nodes','Donors + Acceptors; 1000 iterations'});
hold off;

figure, % hold on;
  semilogy( pval(2:end),abnTTL(2:end) , 'b' , 'LineWidth',2 );
  xlabel('p-value'); ylabel('Abundance (number of clusters >= pvalue)');
  title({'Number of Nodes','Donors + Acceptors; 1000 iterations'});
hold off;



figure, hold on;
  spy(ClusterC1);
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
  title('Original Cluster ConnMat: 1909 nonzero elements');
hold off;

figure,   % hold on;
%  spy(pvalueDONORSv06);
  imagesc( pvalueDONORSv06 ./ 500 );
  axis image,
  colormap( cmap5 );
  colorbar,
  title({'Random Donor Sampling ConnMat:', '500 samples, 754 nonzero elements'});
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
hold off;

figure,   % hold on;
%  spy(pvalueACCEPTORSv06);
  imagesc( pvalueACCEPTORSv07 ./ 500 );
  axis image,
  colormap( cmap5 );
  colorbar,
  title({'Random Acceptor Sampling ConnMat:', '500 samples, 1082 nonzero elements'});
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
hold off;


figure,   % hold on;
%  spy(pvalueACCEPTORSv07);
  imagesc( pvalueACCEPTORSv07 ./ 100 );
  axis image,
  colormap( cmap5 );
  colorbar,
  title({'Random Acceptor Sampling ConnMat:', '100 samples, 978 nonzero elements'});
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
hold off;



cocoa = (pvalueDONORSv06 ./ 500) .* logical(ClusterC1) ;
cacao = (pvalueACCEPTORSv07 ./ 500) .* logical(ClusterC1) ;

figure, hold on;
  imagesc( 0.001*(pvalueACCEPTORSv07 + pvalueDONORSv06) .* logical(ClusterC1) );
  colormap( cmap10 ); colorbar,
  spy( 0.001*(pvalueACCEPTORSv07 + pvalueDONORSv06) .* logical(ClusterC1) , 'k' );
  axis image,
  title({'Random Sampling ConnMat Overlap:', '500 samples, 1201 nonzero elements'});
  xlabel('Acceptor Clusters'); ylabel('Donor Clusters');
hold off;

cocoa = logical(pvalueACCEPTORSv07 + pvalueDONORSv06) .* ClusterC1;
sum(sum(cocoa)) / sum(sum(ClusterC1))   % 0.9058 :: ~91% of the total connection strength is recovered from 1201 nonzero elements (0.6291 percent of the total number of nonzero elements, 1909)



%%% See line 200 %%%

%pvalueDONORSv07    = sparse( 456,786 );
%tempD = sparse( 456,786 );
pvalueACCEPTORSv07 = sparse( 456,786 );
tempA = sparse( 456,786 );

% matlabpool open   % matlabpool('size')

tic;
  parfor cntr=1:500
%  for cntr=1
%    tempD = pipeline_v001( C1,CORA1,CORD1,C1i,ClusterC1,uC1clm  , sparse( 456,786 ) );
%    pvalueDONORSv07    = pvalueDONORSv07    + tempD;
    tempA = pipeline_v002p( C1,CORA1,CORD1,C1j,ClusterC1,uC1clm , sparse( 456,786 ) );
    pvalueACCEPTORSv07 = pvalueACCEPTORSv07 + tempA;
  end
toc;   % ~4.7h
figure, spy( pvalueACCEPTORSv07 , 5 );
clear ans cntr tempD tempA szClstrC1

save('data.mat', 'pvalueACCEPTORSv07','-append');

% matlabpool close




















































