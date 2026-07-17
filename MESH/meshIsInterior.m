function v = meshIsInterior( M , G , ALGs , FULL )
%{

m = ndmat( [0 1],[0 1],[0 1],[0 1],[0 1],0 );
m = unique( m ,'rows');
m( all( ~m(:,1:3) ,2) & any( m(:,4:5) ,2) ,:) = [];
% m(end+1,end) = 0;
T = [];
for c = 1:size(m,1)
  m(c,:)
  meshIsInterior( RV , G , m(c,:) );
  tic
  meshIsInterior( RV , G , m(c,:) );
  t = toc;
  T(c,1) = t;
end
[~,ord] = sort( T ); T = T(ord); m = m(ord,:); [m,T]


%}



  VERBOSE = false;

  if 0
    err = false;
    
    if ~err
      B = MeshBoundary( M );
      if isfield( B , 'tri' ) && ~isempty( B.tri )
        err = true;
      end
    end
    if ~err
      if CheckSelfIntersections( M )
        err = true;
      end
    end
    if err
      warning('Mesh shouldn''t look watertight');
    end
  end

  if nargin < 3, ALGs = [1,2,2.1,4,5.1,5.1,7,8]; end

  if isnumeric( ALGs ), ALGs = num2cell( ALGs(:).' ); end
  if ischar( ALGs ) && isrow( ALGs ),    ALGs = { ALGs }; end
  for a = ALGs(:).', a = a{1};
    if ischar(a), a = lower(a); end
    switch a
      case { 1 ,'bb','boundingbox' }
      case { 2 ,'mb','mini','miniball' }
      case { 2.1 ,'mb2d','miniball2d' }
      case { 3 ,'el','ellipse' }
      case { 3.01 ,'aael','axisalignedellipse' }
      case { 3.1 ,'el2d','ellipse2d' }
      case {4,'ico','icosahedrom'}
      case {4.1,'ico1','icosahedrom1'}
      case {5,'seeds'}
      case {5.1,'mseeds','ms','movilseeds'}
      case {7,'mesh'}
      case {7.01,'pseudonormal'}
      case {7.1,'meshfull'}
      case {8,'tetgen'}
        if isempty( which('tetgen') ),error('Unavailable function ''tetgen''.'); end
      case {9,'winding','wn'}
      case {10,'fastwinding','fwn','fw'}
      otherwise
        error('Unknown algorithm.');
    end
  end

  if nargin < 4, FULL = false; end

  M = struct('xyz',double(M.xyz),'tri',double(M.tri));
  M = MeshTidy( M ,0,true);

  CH = [];  %the convex hull
  T  = [];  %the tetgen
  LOC = [];
  MBC = []; MBR = [];  %center and radius of the hull miniball (cached by getMB)
  FWN = [];            %fast-winding face clusters: dipoles+radii (cached by getFWN)
  FADJ = [];           %nT x 3 edge-neighbor faces, per corner (cached by getFADJ)
  function addNormalsToM_tri()
    if ~isfield( M , 'triNORMALS' )
      M = MeshFixCellOrientation( M );
      M.triNORMALS = meshNormals( M );
    end
  end
  function addNormalsToM_xyz()
    if ~isfield( M , 'xyzNORMALS' )
      addNormalsToM_tri();
      try,    M.xyzNORMALS = meshNormals( M , 'best' );
      catch,  M.xyzNORMALS = meshNormals( M , 'angle' );   %the vertex weighting WITH the sign-test guarantee (Baerentzen-Aanaes); 'area' has counterexamples
      end
    end
  end
  function addNormalsToCH_tri()
    getCH();
    if ~isfield( CH , 'triNORMALS' )
      CH.triNORMALS = meshNormals( CH );
    end
  end
  function addNormalsToCH_xyz()
    if ~isfield( CH , 'xyzNORMALS' )
      addNormalsToCH_tri();
      CH.xyzNORMALS = meshNormals( CH , 'best' );
    end
  end
  function getCH()
    if isempty( CH )
      CH.xyz = double( M.xyz );
      CH.tri = double( convhulln( CH.xyz ) );
      CH.tri = CH.tri( : ,[1,3,2]);
      CH = MeshTidy( CH ,0,true);
    end
  end
  function getMB()
    if isempty( MBR )
      getCH();
      [ MBC , MBR ] = miniball( CH.xyz );  MBC = MBC(:).';
    end
  end
  function getFWN()
  %cluster the faces on a uniform grid over their centroids and store, per
  %cluster: the area-weighted center, the DIPOLE (sum of face vector areas), a
  %bounding radius, and the grouped face list -- the far-field data of the fast
  %winding number (Barill et al. 2018, first order).
    if ~isempty( FWN ), return; end
    addNormalsToM_tri();          %ensures a CONSISTENT face orientation first
    tP1 = M.xyz( M.tri(:,1) ,:);  tP2 = M.xyz( M.tri(:,2) ,:);  tP3 = M.xyz( M.tri(:,3) ,:);
    AN  = cross( tP2 - tP1 , tP3 - tP1 , 2 ) / 2;    %oriented face vector areas
    fa  = sqrt( sum( AN.^2 ,2) );                    %face areas
    CT  = ( tP1 + tP2 + tP3 ) / 3;                   %face centroids
    gN  = min( 48 , max( 16 , round( 16 * ( size(M.tri,1) / 2e4 )^(1/3) ) ) );
                  %clusters-per-axis GROWS with the mesh (16 up to ~20k faces) so the
                  %near-exact shell stays thin: balances the far-field cost
                  %O(nQueries x nClusters) against the exact near work per query.
    bbw = [ min( CT ,[],1) ; max( CT ,[],1) ];
    cw  = max( diff( bbw ,1,1) , eps ) / gN;
    [ ~ , ~ , cid ] = unique( max( min( floor( ( CT - bbw(1,:) ) ./ cw ) , gN-1 ) , 0 ) ,'rows');
    ncl = max( cid );
    sa  = max( accumarray( cid , fa , [ncl,1] ) , eps );
    fC  = [ accumarray( cid , fa.*CT(:,1) , [ncl,1] ) , ...
            accumarray( cid , fa.*CT(:,2) , [ncl,1] ) , ...
            accumarray( cid , fa.*CT(:,3) , [ncl,1] ) ] ./ sa;
    fN  = [ accumarray( cid , AN(:,1) , [ncl,1] ) , ...
            accumarray( cid , AN(:,2) , [ncl,1] ) , ...
            accumarray( cid , AN(:,3) , [ncl,1] ) ];
    dv  =           sqrt( sum( ( tP1 - fC(cid,:) ).^2 ,2) );
    dv  = max( dv , sqrt( sum( ( tP2 - fC(cid,:) ).^2 ,2) ) );
    dv  = max( dv , sqrt( sum( ( tP3 - fC(cid,:) ).^2 ,2) ) );
    fR  = accumarray( cid , dv , [ncl,1] , @max );   %cluster bounding radius
    [ ~ , ord_ ] = sort( cid );
    ptr_ = cumsum( [ 1 ; accumarray( cid , 1 , [ncl,1] ) ] );
    FWN  = struct( 'C',fC , 'N',fN , 'R',fR , 'ord',ord_ , 'ptr',ptr_ , ...
                   'P1',tP1 , 'P2',tP2 , 'P3',tP3 );
  end
  function getT()
    if isempty( T )
      T = tetgen( M );
      T = struct( 'xyz' , double( T.xyz ) , 'tri' , double( T.tri ) );
    end
  end
  function generateLocalizer()
    if isempty( LOC )
      vtkClosestElement( M );
      LOC = onCleanup( @()vtkClosestElement( [] , [] ) );
    end
  end
  function getFADJ()
  %nT x 3 face adjacency: FADJ(f,j) is the face sharing the edge OPPOSITE corner j
  %of face f (0 if none, e.g. a boundary edge). The edge opposite corner j joins
  %the OTHER two corners -- matching the barycentric convention bc(:,j) ~ 0 <=> the
  %closest point lies on that edge. Built from meshEsuE (edge-based element
  %adjacency) once and cached.
    if ~isempty( FADJ ), return; end
    Tl  = M.tri;  nTl = size( Tl ,1);
    A   = meshEsuE( M , true , 'e' );        %A{f} = faces edge-adjacent to f
    eop = [ 2 3 ; 3 1 ; 1 2 ];               %corners of the edge opposite corner j
    FADJ = zeros( nTl , 3 );
    for f = 1:nTl
      nbs = A{f};  if isempty( nbs ), continue; end
      for j = 1:3
        u = Tl( f , eop(j,1) );  v = Tl( f , eop(j,2) );
        for g = nbs(:).'
          if any( Tl(g,:)==u ) && any( Tl(g,:)==v ), FADJ(f,j) = g; break; end
        end
      end
    end
  end

  % center mesh and queries by a COMMON shift: the inside/outside test is
  % translation-invariant, and keeping the coordinates near the origin avoids the
  % catastrophic cancellation of distances / dot-products when the data lives far
  % from the origin (world or image coordinates with a large offset). The mean is
  % used as the shift (the miniball center would do as well, but it would force
  % the convex-hull computation even when no algorithm needs it).
  cXYZ  = mean( M.xyz , 1 );
  M.xyz = M.xyz - cXYZ;
  G = double(G) - cXYZ;
  v = zeros( size(G,1) , 1 ,'int8'); G(:,4) = ( 1:size(G,1) ).';

  w = ismember( G(:,1:3) , M.xyz ,'rows');
  G( w ,:) = [];

  for a = ALGs(:).', a = a{1};
    if isempty( G ), break; end
    if ischar(a), a = lower(a); end

    if VERBOSE
      nG = size(G,1);
      fprintf('Going with ALG: %s\n', uneval( a ) );
      fprintf('    %d in G  ', nG );
      START = tic;
    end
    switch a
      case { 1 ,'bb','boundingbox' }
        getCH();
        BB = meshBB( CH );
        for d = 1:3
          w = ( G(:,d) < BB(1,d) ) | ( G(:,d) > BB(2,d) );
          v( G(w,4) ) = -1; G( w ,:) = [];
        end

      case { 2 ,'mb','mini','miniball' }
        getMB();
        w = sum( ( G(:,1:3) - MBC ).^2 ,2) > MBR^2;
        v( G(w,4) ) = -1; G( w ,:) = [];

      case { 2.1 ,'mb2d','miniball2d' }
        getCH();
        [C,R] = miniball( CH.xyz(:,[1,2]) ); C = C(:).';
        w = sum( ( G(:,[1,2]) - C ).^2 ,2) > R^2;
        v( G(w,4) ) = -1; G( w ,:) = [];

        [C,R] = miniball( CH.xyz(:,[1,3]) ); C = C(:).';
        w = sum( ( G(:,[1,3]) - C ).^2 ,2) > R^2;
        v( G(w,4) ) = -1; G( w ,:) = [];

        [C,R] = miniball( CH.xyz(:,[2,3]) ); C = C(:).';
        w = sum( ( G(:,[2,3]) - C ).^2 ,2) > R^2;
        v( G(w,4) ) = -1; G( w ,:) = [];

      case { 3 ,'el','ellipse' }
        %MINIMUM-VOLUME enclosing ellipsoid (Loewner-John MVEE) of the hull
        %vertices, via minEllipsoid (Khachiyan; enclosure GUARANTEED, rotation
        %included): the tightest ellipsoidal prune. Anything outside it is
        %outside the mesh. Conservative >1+1e-12 per minEllipsoid's help.
        getCH();
        [ EA , EC ] = minEllipsoid( CH.xyz );
        w = sum( ( ( G(:,1:3) - EC ) * EA ) .* ( G(:,1:3) - EC ) ,2) > 1 + 1e-12;
        v( G(w,4) ) = -1; G( w ,:) = [];

      case { 3.01 ,'aael','axisalignedellipse' }
        %AXIS-ALIGNED enclosing ellipsoid (fminsearch over per-axis scales +
        %miniball; the OLD ALG 3): kept for reference -- usually LOOSER (no
        %rotation) and slower than the MVEE above.
        getCH();
        [C,R] = miniball( CH.xyz , 'ellipse' );
        w = sum( ( G(:,1:3) - C ).^2 .* ( 1./(R.^2) ) ,2) > 1;
        v( G(w,4) ) = -1; G( w ,:) = [];

      case { 3.1 ,'el2d','ellipse2d' }
        getCH();
        [C,R] = miniball( CH.xyz .* [1,1,0] , 'ellipse' );
        w = sum( bsxfun( @times , bsxfun( @minus , G(:,[1,2]) , C(:,[1,2]) ).^2 , 1./(R(:,[1,2]).^2) ) ,2) > 1;
        v( G(w,4) ) = -1; G( w ,:) = [];

        [C,R] = miniball( CH.xyz .* [1,0,1] , 'ellipse' );
        w = sum( bsxfun( @times , bsxfun( @minus , G(:,[1,3]) , C(:,[1,3]) ).^2 , 1./(R(:,[1,3]).^2) ) ,2) > 1;
        v( G(w,4) ) = -1; G( w ,:) = [];
        
        [C,R] = miniball( CH.xyz .* [0,1,1] , 'ellipse' );
        w = sum( bsxfun( @times , bsxfun( @minus , G(:,[2,3]) , C(:,[2,3]) ).^2 , 1./(R(:,[2,3]).^2) ) ,2) > 1;
        v( G(w,4) ) = -1; G( w ,:) = [];

      case {4,'ico','icosahedrom'}
        getCH();
        S = meshFacesCenter( sphereMesh(0) );
        S = [ S ; S .* [1,1,-1] ];
        S( S(:,3) > 0 ,:) = [];
        S = unique( S ,'rows');
        for d = 1:size(S,1)
          CHD = CH.xyz * S(d,:).';
          GD  = G(:,1:3) * S(d,:).';

          m = max( CHD ); w = GD > m; GD( w ,:) = [];
          v( G(w,4) ) = -1; G( w ,:) = []; 

          m = min( CHD ); w = GD < m;
          v( G(w,4) ) = -1; G( w ,:) = []; 
        end
      
      case {4.1,'ico1','icosahedrom1'}
        getCH();
        S = meshFacesCenter( sphereMesh(1) );
        S = [ S ; S .* [1,1,-1] ];
        S( S(:,3) > 0 ,:) = [];
        S = unique( S ,'rows');        
        for d = 1:size(S,1)
          CHD = CH.xyz * S(d,:).';
          GD  = G(:,1:3) * S(d,:).';

          m = max( CHD ); w = GD > m; GD( w ,:) = [];
          v( G(w,4) ) = -1; G( w ,:) = []; 

          m = min( CHD ); w = GD < m;
          v( G(w,4) ) = -1; G( w ,:) = []; 
        end

      case {5,'seeds'}
        bb = [ min( G(:,1:3) ,[],1) ; max( G(:,1:3) ,[],1) ];
        S = [ ndmat( bb(:,1) , bb(:,2) , bb(:,3 ) ) ; G( unique( round( linspace( 1 , size(G,1) , 500 ) ) ) ,1:3) ];

        [ io , r ] = isInterior( S );
        for s = 1:size(S,1)
          if isnan( r(s) ), continue; end
          w = sum( ( G(:,1:3) - S(s,:) ).^2 ,2) < 0.99*r(s)^2;
          v( G(w,4) ) = io(s); G( w ,:) = [];
        end

      case {5.1,'mseeds','ms','movilseeds'}
        bb = meshBB( M ); bb = 10*( bb - mean( bb ,1) ) + mean( bb ,1);


        addNormalsToM_tri();
        if ~exist( 'fc' , 'var' )
          getMB();
          off = 1e-3 * MBR;   %straddle offset RELATIVE to the object scale (not absolute)
          fc = meshFacesCenter( M );
          [~,S] = FarthestPointSampling( fc ,[],0,50);
          S = [ fc(S,:) + off * M.triNORMALS(S,:) ; fc(S,:) - off * M.triNORMALS(S,:) ];
        else
          S = G( unique( round( linspace( 1 , size(G,1) , 100 ) ) ) ,1:3);
        end

        [ io , r , cp ] = isInterior( S );
        for s = 1:size(S,1)
          if isnan( r(s) ), continue; end
%           plotMESH( M ,'[0.2]','nice'); hplotMESH( transform( sphereMesh(4) ,'s',r(s),'t',S(s,:) ) ,'b[0.5]','nice'); ze
          while 1
            if io(s) < 0, nS = S(s,:) + ( S(s,:) - cp(s,:) )*0.5;
            else,         nS = S(s,:) + ( S(s,:) - cp(s,:) )*0.2;
            end
            [ nio , nr ] = isInterior( nS );
            if isnan(nr) || nr < r(s), break; end
             S(s,:) = nS;
            io(s,:) = nio;
             r(s,:) = nr;
%             cp(s,:) = ncp;
            if nio < 0 && ( nS(1) < bb(1,1) || nS(1) > bb(2,1) || ...
                            nS(2) < bb(1,2) || nS(2) > bb(2,2) || ...
                            nS(3) < bb(1,3) || nS(3) > bb(2,3) )
              break
            end
          end
%           hplotMESH( transform( sphereMesh(4) ,'s',r(s),'t',S(s,:) ) ,'r[0.5]','nice'); ze
        end

%         plotMESH( M ,'[0.2]','nice');
%         for s = 1:size(S,1)
%           if io(s)<0, continue; end
%           hplotMESH( transform( sphereMesh(4) ,'s',r(s),'t',S(s,:) ) ,'r[0.5]','nice'); ze
%         end

        for s = find( io(:).' > 0 )
          w = sum( ( G(:,1:3) - S(s,:) ).^2 ,2) < 0.99*r(s)^2;
          v( G(w,4) ) = 1; G( w ,:) = [];
        end
        for s = find( io(:).' < 0 )
          w = sum( ( G(:,1:3) - S(s,:) ).^2 ,2) < 0.99*r(s)^2;
          v( G(w,4) ) = -1; G( w ,:) = [];
        end
        

      case {7,'mesh'}
        [ io , r ] = isInterior( G(:,1:3) );
        
        w = ~isnan(r);
        v( G(w,4) ) = io(w); G( w ,:) = [];

      case {7.01,'pseudonormal'}
        [ io , r ] = isInterior_pseudonormal( G(:,1:3) );

        w = ~isnan(r);
        v( G(w,4) ) = io(w); G( w ,:) = [];

      case {7.1,'meshfull'}
        [ io , r ] = isInterior_meshfull( G(:,1:3) );
        
        w = ~isnan(r);
        v( G(w,4) ) = io(w); G( w ,:) = [];

      case {8,'tetgen'}
        getT();
        w = tsearchn( T.xyz , T.tri , G(:,1:3) );
        w = ~isnan(w);
        v( G( w,4) ) =  1;
        v( G(~w,4) ) = -1;
        G( [ w | ~w ] ,:) = [];

      case {9,'winding','wn'}
        %generalized winding number (Jacobson et al. 2013): |w|~1 inside, ~0
        %outside (the abs makes it independent of the GLOBAL orientation; it only
        %needs a CONSISTENT one). Decides every remaining point, like tetgen.
        wn = windingNumber( G(:,1:3) );
        w  = abs( wn ) > 0.5;
        v( G( w,4) ) =  1;
        v( G(~w,4) ) = -1;
        G( [ w | ~w ] ,:) = [];

      case {10,'fastwinding','fwn','fw'}
        %fast winding number (Barill et al. 2018, first order): far face clusters
        %enter through their area-normal DIPOLE, near clusters are summed exactly.
        wn = fastWindingNumber( G(:,1:3) );
        w  = abs( wn ) > 0.5;
        v( G( w,4) ) =  1;
        v( G(~w,4) ) = -1;
        G( [ w | ~w ] ,:) = [];

%       case {6,'convexhull','ch'}
%         addNormalsToCH_tri();
% 
%         [ io , r ] = isInterior( CH , G(:,1:3) );
%         
%         w = io < 0 & ~isnan(r);
%         v( G(w,4) ) = -1; G( w ,:) = [];

%       case {6.1,'seedsconvexhull'}
%         addNormalsToCH_tri();
% 
%         ids = unique( round( linspace( 1 , size(G,1) , 100 ) ) );
%         CS = G( ids ,1:3);
% 
%         [ io , r ] = isInterior( CH , CS );
%         for c = 1:size(CS,1)
%           if isnan( r(c) ), continue; end
%           if io(c) >= 0, continue; end
% 
%           w = sum( bsxfun( @minus , G(:,1:3) , CS(c,:) ).^2 ,2) < 0.99*r(c)^2;
%           v( G(w,4) ) = io(c); G( w ,:) = [];
%         end

    end
    if VERBOSE
      fprintf('    took %f seconds for removing %d points.\n' , toc( START ) , nG - size(G,1) );
    end

  end

  if 0 && ALGs(6), try
    T = getT();
    Tv = meshQuality( T , 'volume' );
    for it = 1:20
      if ~any( Tv ), break; end
      [~,id] = max( Tv ); Tv( id ) = 0;

      C = mean( T.xyz( T.tri( id ,:) ,:) ,1);
      [~,~,R] = vtkClosestElement( M , C ); R = R*0.9999; R = R*R;
      w = sum( bsxfun( @minus , G(:,1:3) , C(:).' ).^2 ,2) < R;
      v( G(w,4) ) = true;
      G( w ,:) = [];
%       w = ~isnan( tsearchn( T.xyz( T.tri( id ,:) ,:) , [1 2 3 4 ] , G(:,1:3) ) );
%       v( G(w,4) ) = true;
%       G( w ,:) = [];
    end
  end; end
  if 0 && ALGs(7), try
    CH = getCH();
    CH = MeshFixCellOrientation( CH );
    CH.triNORMALS = meshNormals( CH );
    CH.xyzNORMALS = meshNormals( CH , 'uniform' );

    w = ~meshIsInterior_helper( CH , G(:,1:3) );
    v( G(w,4) ) = false;
    G( w ,:) = [];
  end; end

  if FULL && ~isempty( G )
    addNormalsToM_xyz();
    w = meshIsInterior_helper( M , G(:,1:3) );
  
    v( G( w,4) ) =  1;
    v( G(~w,4) ) = -1;
    v = v > 0;
  end

  function [ io , r , cp ] = isInterior( x )
    addNormalsToM_tri();
    generateLocalizer();
    [ e , cp , r , bc ] = vtkClosestElement( double( x ) );
    r( any( bc < 1e-6 ,2) ) = NaN;
    io = sign( dot( cp - x , M.triNORMALS( e ,:) ,2) );
  end

  function [ io , r , cp ] = isInterior_pseudonormal( x )
  % 'pseudonormal' (ALG 7.01) -- inside(+1)/outside(-1) from the sign of (cp-x).n,
  % where cp is the closest surface point and n is the ANGLE-WEIGHTED PSEUDONORMAL
  % at cp, obtained by BARYCENTRIC INTERPOLATION of the three vertex pseudonormals
  % of the closest face (M.xyzNORMALS, computed 'best' with an 'angle' fallback).
  %
  % Why the vertex normals MUST be angle-weighted (not uniform / area): the sign
  % test is valid at a vertex only if n lies inside the cone of the incident face
  % normals, and the incident-ANGLE weighting is the one that carries that
  % guarantee.
  %   [TW98] G. Thurmer & C. A. Wuthrich, "Computing vertex normals from polygonal
  %          facets", Journal of Graphics Tools 3(1):43-46, 1998.   (angle weights)
  %   [BA05] J. A. Baerentzen & H. Aanaes, "Signed Distance Computation Using the
  %          Angle Weighted Pseudonormal", IEEE Transactions on Visualization and
  %          Computer Graphics 11(3):243-253, 2005.  (proves the sign correctness
  %          of the pseudonormal at a face, an edge and a vertex).
  % Uniform / area weights can fall OUTSIDE that cone on skewed fans and then flip
  % the sign (measured on a spiky mesh: ~36/40 and ~25/40 vertex probes wrong, vs 0
  % for angle / 'best') -- this is why addNormalsToM_xyz uses 'best'.
  %
  % NEEDED CORRECTION (see ALG 7.1 'meshfull'): barycentrically interpolating the
  % two ENDPOINT vertex pseudonormals along an EDGE is NOT the exact edge
  % pseudonormal. [BA05] shows the exact normal at a point ON an edge is n_f1 + n_f2
  % -- the sum of the TWO incident FACE normals (the bisector / max-margin axis of
  % the dihedral wedge), which depends ONLY on those two faces. The interpolation
  % instead mixes in every OTHER face incident to the two endpoints, so on sharp
  % asymmetric creases it can tilt out of the edge's normal wedge and give the wrong
  % sign (measured: ~1/894 edge probes wrong even with 'best'; ALG 7.1 fixes it to
  % 0). Use 7.01 as the fast single-field variant, 7.1 when edges must be exact.
    addNormalsToM_xyz();
    generateLocalizer();
    [ e , cp , r , bc ] = vtkClosestElement( double( x ) );

    r( r < 1e-8 ) = NaN;
    io = sign( dot( cp - x , bc(:,1) .* M.xyzNORMALS( M.tri(e,1) ,:) + bc(:,2) .* M.xyzNORMALS( M.tri(e,2) ,:) + bc(:,3) .* M.xyzNORMALS( M.tri(e,3) ,:) ,2) );
  end

  function [ io , r , cp ] = isInterior_meshfull( x )
  % 'meshfull' (ALG 7.1) -- like 7.01 but EXACT on edges and vertices. It picks the
  % feature-based Baerentzen-Aanaes pseudonormal [BA05] from the barycentric
  % coordinates of the closest point cp on face e:
  %   * all bc >  tol -> cp INSIDE the face      -> n = face normal
  %   * one  bc ~= 0  -> cp on the OPPOSITE edge -> n = n_f(e) + n_f(neighbor), the
  %                      sum of the two incident FACE normals (the exact edge
  %                      pseudonormal; neighbor via meshEsuE / getFADJ)
  %   * two  bc ~= 0  -> cp at the CORNER (bc~1) -> n = angle-weighted vertex
  %                      pseudonormal (M.xyzNORMALS, 'best')            [TW98,BA05]
  % Only the SIGN of (cp-x).n is used, so the summed edge normal need not be unit.
  % This removes the residual ~0.1% edge misclassifications that 7.01 can make on
  % sharp creases, while leaving the face and vertex cases identical to the theory.
    addNormalsToM_xyz();          %'best' vertex pseudonormals (+ oriented face normals)
    getFADJ();
    generateLocalizer();
    [ e , cp , r , bc ] = vtkClosestElement( double( x ) );
    r( r < 1e-8 ) = NaN;

    tol = 1e-6;
    z   = bc < tol;                          %barycentric coords ~ 0
    nz  = sum( z , 2 );                       %0 -> face, 1 -> edge, 2 -> vertex

    N = M.triNORMALS( e ,:);                  %default (face interior): face normal

    we = find( nz == 1 );                     %EDGE hits
    if ~isempty( we )
      [ ~ , j ] = max( z(we,:) , [] , 2 );    %the single zero corner (edge opposite it)
      g  = FADJ( sub2ind( size(FADJ) , e(we) , j ) );   %neighbor face across that edge
      ok = g > 0;                             %boundary edge (g==0): keep face normal
      Ne = M.triNORMALS( e(we) ,:);
      Ne( ok ,:) = Ne( ok ,:) + M.triNORMALS( g(ok) ,:);
      N( we ,:) = Ne;
    end

    wv = find( nz == 2 );                     %VERTEX hits
    if ~isempty( wv )
      [ ~ , jv ] = max( bc(wv,:) , [] , 2 );  %the corner with bc ~ 1
      N( wv ,:) = M.xyzNORMALS( M.tri( sub2ind( size(M.tri) , e(wv) , jv ) ) ,:);
    end

    io = sign( dot( cp - x , N , 2 ) );
  end

  function w = windingNumber( x )
  %generalized winding number at the rows of x: the sum of the signed solid
  %angles of EVERY face, normalized by 4*pi. BRUTE FORCE (all faces per query)
  %but fully vectorized and row-blocked.
    addNormalsToM_tri();          %ensures a CONSISTENT face orientation first
    w = solidAngleW( x , M.xyz( M.tri(:,1) ,:) , M.xyz( M.tri(:,2) ,:) , M.xyz( M.tri(:,3) ,:) );
  end

  function w = fastWindingNumber( x )
  %fast winding number: single-level Barnes-Hut. A cluster whose bounding sphere
  %is farther than BETAFW times its radius enters through its DIPOLE
  %  w_far = ( N_c . (C_c - x) ) / ( 4*pi*|C_c - x|^3 ) ,   N_c = sum of vector areas
  %(the first-order far field of the exact solid angle); near clusters are summed
  %EXACTLY. Larger BETAFW = more accurate and slower; ~2.5 is plenty for the
  %0.5-threshold here (the approximation error is ~1e-3 in w).
    getFWN();
    BETAFW = 2.5;
    nQ_  = size( x , 1 );
    w    = zeros( nQ_ , 1 );
    nCl_ = size( FWN.C , 1 );
    blk_ = max( 1 , floor( 4e6 / nCl_ ) );
    for kk = 1:blk_:nQ_
      bq = kk : min( kk+blk_-1 , nQ_ );
      dX = FWN.C(:,1).' - x(bq,1);                   %block x nClusters
      dY = FWN.C(:,2).' - x(bq,2);
      dZ = FWN.C(:,3).' - x(bq,3);
      DD = dX.*dX + dY.*dY + dZ.*dZ;
      isfar = DD > ( BETAFW * FWN.R.' ).^2;
      iD3 = zeros( size(DD) );  iD3(isfar) = DD(isfar).^(-1.5);   %0 where near (no 0/0)
      w(bq) = ( ( dX .* iD3 ) * FWN.N(:,1) ...       %dipole far field, summed over clusters
              + ( dY .* iD3 ) * FWN.N(:,2) ...
              + ( dZ .* iD3 ) * FWN.N(:,3) ) / ( 4*pi );
      for cc = find( any( ~isfar , 1 ) )             %near clusters: exact solid angles
        qn = bq( ~isfar(:,cc) );
        tt = FWN.ord( FWN.ptr(cc) : FWN.ptr(cc+1)-1 );
        w(qn) = w(qn) + solidAngleW( x(qn,:) , FWN.P1(tt,:) , FWN.P2(tt,:) , FWN.P3(tt,:) );
      end
    end
  end

end


function w = solidAngleW( x , Q1 , Q2 , Q3 )
%sum over the triangles (Q1,Q2,Q3) of their signed solid angle seen from each
%row of x, normalized by 4*pi (van Oosterom & Strackee 1983):
%  tan(O/2) = det([a b c]) / ( |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b| )
%fully vectorized (block x triangles via implicit expansion), row-blocked so the
%peak memory stays bounded (~15 block-x-nT temporaries).
  nT = size( Q1 , 1 );
  nQ = size( x  , 1 );
  w  = zeros( nQ , 1 );
  blk = max( 1 , floor( 1e6 / max( nT , 1 ) ) );
  for k = 1:blk:nQ
    b  = k : min( k+blk-1 , nQ );
    ax = Q1(:,1).' - x(b,1);  ay = Q1(:,2).' - x(b,2);  az = Q1(:,3).' - x(b,3);
    bx = Q2(:,1).' - x(b,1);  by = Q2(:,2).' - x(b,2);  bz = Q2(:,3).' - x(b,3);
    cx = Q3(:,1).' - x(b,1);  cy = Q3(:,2).' - x(b,2);  cz = Q3(:,3).' - x(b,3);
    na = sqrt( ax.*ax + ay.*ay + az.*az );
    nb = sqrt( bx.*bx + by.*by + bz.*bz );
    nc = sqrt( cx.*cx + cy.*cy + cz.*cz );
    de = ax.*( by.*cz - bz.*cy ) + ay.*( bz.*cx - bx.*cz ) + az.*( bx.*cy - by.*cx );
    dn = na.*nb.*nc + ( ax.*bx + ay.*by + az.*bz ).*nc ...
                    + ( bx.*cx + by.*cy + bz.*cz ).*na ...
                    + ( cx.*ax + cy.*ay + cz.*az ).*nb;
    w(b) = sum( atan2( de , dn ) , 2 ) / ( 2*pi );   %Omega = 2*atan2 ; w = sum(Omega)/4pi
  end
end
