function M = MeshImproveObtuses( M , L , maxIts )

  if nargin < 3, maxIts = 1e4; end
  it = 0;
  while 1, it = it+1; if it > maxIts, break; end
    A = meshQuality( M ,'angles');
    [t,s] = ind2sub( size(A) , argmax( A ) ); %disp([t,s]);
    if A(t,s) < L, break; end
    disp(A(t,s));

    switch s
      case 1, abc = M.tri(t,[2,3,1]);
      case 2, abc = M.tri(t,[3,1,2]);
      case 3, abc = M.tri(t,[1,2,3]);
    end
    a = abc(1); b = abc(2); c = abc(3);
    aa = M.xyz(a,:);
    bb = M.xyz(b,:);
    cc = M.xyz(c,:);

    x = bb-aa;
    y = cc-aa;
    dd = aa + dot( x , y ) / dot( x , x ) * x;
    d = size( M.xyz ,1)+1; M.xyz(d,:) = dd;

    ts = find( sum( ismember( M.tri ,[a,b] ) ,2) > 1 );

    for tt = ts(:).'
      T0 = M.tri( tt ,:); T1 = T0;
      T0( T0 == b ) = d;
      T1( T1 == a ) = d;
      
      M.tri( tt    ,:) = T0;
      M.tri( end+1 ,:) = T1;
      for f = fieldnames( M ,'^tri.+').', f = f{1};
        M.(f)( end+1 ,:,:,:,:) = M.(f)( tt ,:,:,:,:);
      end
    end

  end

  M = MeshTidy( M ,0,true);

end
