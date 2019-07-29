%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load 'edit.mat' before you run this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%

function [ ClusterConnMat, ConnMatShuffled, rowsshuff,colsshuff,clmask,clmaskmod ] = lumpclusters_v006( ConnMat,CORRACCEPT,COORDONOR, wsp,pp )
%%
% Variables
C1 = ConnMat;
CORA1 = CORRACCEPT;   % corr( ConnMat ) = corr_matr( C1 )
CORD1 = COORDONOR;    % corr( ConnMat' ) = corr_matr( C1' )
wwsspp = wsp;
print = pp;
%
% Acceptors
[ cla1 ] = watershedc2( (CORA1>0.05).*CORA1, wwsspp );
[ vvv, aaa ] = values2( cla1 );
CCCA = zeros( size(C1) );

iii2 = [];
ddd = [];

[saaa, ind] = sort( aaa, 'descend' );
vvv = vvv( ind );

C1A = C1; % matrix to be compressed;

ic = 1;

%%% CLUSTER COUNTING PROBLEM %%%
for j=1:length(vvv)
  ii = find( cla1==vvv(j) );
  iii2 = [iii2; ii];
%  ddd = [ ddd; j*ones(size(ii)) ];

  if j==1
    ibcla1 = ii;
  end

  if vvv(j)>0
    ddd = [ ddd; ic*ones(size(ii)) ];
    C1A(:,ic)=sum(C1(:,ii),2);
    ic = ic + 1;
  else
    ddd = [ ddd; [ic:(ic+length(ii)-1)]' ];
    C1A(:,ic:(ic+length(ii)-1)) = C1(:,ii);
    ic = ic + length(ii);
  end
end
%%% CLUSTER COUNTING PROBLEM %%%

C1A = C1A(:,1:(ic-1));
CBA1 = CORA1(ibcla1,ibcla1);

CCCA = (CORA1(iii2, iii2));
ddda = ddd;
ddd = ddd*ones(size(ddd'));
ddd = ddd + ddd';
ddd = mod(ddd,3);

Cshuff = C1(:,iii2);   % acceptors index shuffling

if ( print==1 )
  figure(1)
    imagesc(CBA1);
    axis image; colorbar;

  figure(2)
    imagesc(CCCA+ddd*0.01);
    axis image; colorbar

  figure(3)
    spy(Cshuff); axis image;
    title('Within cluster connectivity shuffled');
    xlabel('acceptors');
    ylabel('donors');

  figure(4)
    spy(C1A);
end

%
% Donors
[ cld1 ] = watershedc2( (CORD1>0.05).*CORD1, wwsspp );
[ vvv, aaa ] = values2(cld1);
CCCD = zeros( size(C1A) );

iii3 = [];
ddd = [];

[saaa, ind] = sort(aaa, 'descend');
vvv = vvv(ind);

C1AD = C1A; % matrix to be compressed;
ic = 1;

%%% CLUSTER COUNTING PROBLEM %%%
for j=1:length(vvv)
  ii = find(cld1==vvv(j));
  iii3 = [iii3; ii];
%  ddd = [ddd; j*ones(size(ii))];

  if vvv(j)>0
    ddd = [ddd; ic*ones(size(ii))];
    C1AD(ic,:) = sum(C1A(ii,:),1);
    ic = ic + 1;
  else
    ddd = [ ddd; [ic:(ic+length(ii)-1)]' ];
    C1AD(ic:(ic+length(ii)-1),:) = C1A(ii,:);
    ic = ic + length(ii);
  end
end
%%% CLUSTER COUNTING PROBLEM %%%

C1AD = C1AD(1:(ic-1),:);
CCCD = (CORD1(iii3, iii3));

dddd = ddd;
ddd = ddd*ones(size(ddd'));
ddd = ddd+ddd';
ddd = mod(ddd,3);

Cshuff = Cshuff(iii3,:);   % donors index shuffling

if ( print==1 )
  figure(5)
    imagesc(CCCD+ddd*0.01);
    axis image;
    colorbar;

  figure(6)
    spy(Cshuff);
    axis image;
    title('Within cluster connectivity shuffled');
    xlabel('acceptors');
    ylabel('donors');
end

ddd = dddd*ones(size(ddda'))+ones(size(dddd))*ddda';
%
dddmod = mod(ddd,3);
dddmod(1,1)=10;

if ( print==1 )
  figure(7); hold on;
    imagesc(dddmod*.1);
    m = gray;
    colormap(m(end:(-1):1,:));
    axis image;
    spy(Cshuff);
  hold off;

  figure(8); hold on;
    spy(C1AD);
  hold off;
end
%
% Returning the values of the function
ClusterConnMat = C1AD;
ConnMatShuffled = Cshuff;
rowsshuff = iii3;
colsshuff = iii2;
clmask = ddd;
clmaskmod = dddmod;
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
