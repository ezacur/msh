function [v,c,IT] = meshVolume( M , varargin )
%MESHVOLUME  Volume, centroid and 2nd-moment matrix of a closed surface mesh.
%
%   V              = meshVolume( M )
%   [V,C]          = meshVolume( M )
%   [V,C,IT]       = meshVolume( M , ... )
%
%   For a CLOSED, consistently-oriented triangular surface mesh M (fields .xyz
%   and .tri) computes, via divergence-theorem polynomial-moment integrals over
%   the faces:
%     V  : the enclosed volume (scalar).
%     C  : the volume centroid,  C = [ integral(x) integral(y) integral(z) ]/V.
%     IT : the 3x3 matrix of SECOND MOMENTS,  IT(i,j) = integral( x_i*x_j dV ).
%
%   IMPORTANT about IT:
%     * IT is the SECOND-MOMENT (a.k.a. covariance-of-mass) matrix, NOT the
%       physical inertia tensor. They share eigenVECTORS (principal axes), so
%       IT is fine for aligning to principal axes; for the physical inertia
%       tensor use  trace(IT)*eye(3) - IT.
%     * IT is taken about the ORIGIN, unless 'center' is given (then about C).
%
%   OPTIONS (name flags, any order):
%     'center'            translate M to its centroid before the 2nd moments,
%                         so IT is about the centroid (C is still in the
%                         original frame).
%     'noorient'          do NOT call MeshFixCellOrientation first. By DEFAULT
%                         the mesh is re-oriented, so V is always returned
%                         POSITIVE; with 'noorient' V keeps its SIGN (negative
%                         if the faces wind inward).
%     'noclean'           do NOT call MeshTidy first (and skip the open-mesh
%                         check). By DEFAULT the mesh is cleaned and a warning
%                         is issued if it looks open.
%                         ( default clean+orient costs ~6x vs 'noclean','noorient';
%                           pass both when M is already clean and oriented. )
%     'byslicing' , N  |  'slice' , N
%                         alternative estimate by Cavalieri z-slices instead of
%                         the divergence theorem: N>0 = number of slices,
%                         N<0 = slice thickness. Cross-section area per slice is
%                         summed times the spacing. (Requires the 2D polygon
%                         library -- polygon/area/union.)
%
%   EXAMPLE (unit-sphere sanity check, should be ~1):
%     M.xyz = randn(1000,3);
%     M.xyz = bsxfun(@rdivide,M.xyz,sqrt(sum(M.xyz.^2,2)));
%     M.tri = convhulln( M.xyz );
%     meshVolume( M ) / (4/3*pi)
%
% See also MeshFixCellOrientation, MeshTidy, MeshBoundary, meshSlice,
%          MeshAlignByInertiaTensor.

  M = Mesh( M ,0);
  
  [varargin,~,sliceIntegration] = parseargs(varargin,'byslicing','slice' ,'$DEFS$',0);
  
  if ~~sliceIntegration
    N = sliceIntegration;
    
    cl = @(x)x([1:end,1],:);
    bb = meshBB( M );
    if N > 0, Zs = dualVector( linspace( bb(1,3) , bb(2,3) , N ) );
    else,     Zs = ( bb(1,3) + N/2 ):-N:( bb(2,3) - N/2 );
    end
    v = 0;
    c = [];
    for z = Zs(:).'
      T = getPlane([0 0 z;0 0 1]);
      try, C = meshSlice( M , T ); if isempty( C ), continue; end; catch, continue; end
      C = transform( C , minv(T) ) * eye(3,2);
      C = polyline( C ); %figure; pplot(C); title(uneval(z));
      
      
      A = [];
      for a = 1:C.np
        D = double( C(a) );
        if size(D,1) < 3, continue; end
        if ~isequal( D(1,:) , D(end,:) ), D(end+1,:) = D(1,:); end
        if isempty( A ), A = polygon(D);
        else,            A = union( A , polygon(D) );
        end
        if nargout > 1
          D(:,3) = 0;
          D = transform( D ,T);
          D = Mesh( D ,'ClosedContour');
          D = MeshAddField( D , 'xyzZ' , z );
          c = MeshAppend( c , D );
        end
      end
      if ~isempty( A ), v = v + area( A ); end
    end
    v = v * mean( diff(Zs) );
    
    return;
  end
  

  [varargin,CLEAN ] = parseargs(varargin,'noclean' ,'$FORCE$',{false,true});
  [varargin,ORIENT] = parseargs(varargin,'noorient','$FORCE$',{false,true});

  if CLEAN
    M = MeshTidy( M ,0,true);

    bounds = MeshBoundary( M );
    if isfield( bounds , 'tri' ) && size( bounds.tri ,1) > 0, warning('MESH look open. try with vtkFillHolesFilter'); end
  end
  

  if ORIENT
    M = MeshFixCellOrientation( M );
  end
  
  
  N = cross( ( M.xyz(M.tri(:,2),:) - M.xyz(M.tri(:,1),:) ) , ( M.xyz(M.tri(:,3),:) - M.xyz(M.tri(:,2),:) ) , 2 );
  A = sqrt( sum( N.^2 , 2 ) );
  N = bsxfun( @rdivide , N , A );
  A = A/2;
  
  %v = vtkMassProperties( M , 'GetVolume' );
  
  v =  sum(  ( M.xyz( M.tri(:,1),: ) + M.xyz( M.tri(:,2),: ) + M.xyz( M.tri(:,3),: ) ) .* N , 2 );
  
  w = ~isfinite( A ) | ~isfinite( v );
  v(w) = [];
  A(w) = [];
  
  
  v = 2*sum( v .* A )/( 3 * factorial(3) );

%   ISC = {};
%   disp([5 4 3]); 
%   disp( moment( [ 5 4 3] ) );
%   disp([3 4 5]); moment( [ 3 4 5] );

%   uneval([ moment([2 0 0]) moment([1 1 0]) moment([1 0 1]) ;
%     0               moment([0 2 0]) moment([0 1 1]) ;
%     0               0               moment([0 0 2]) ])

  if nargout > 1
    c = [ moment([1 0 0]) ,  moment([0 1 0]) , moment([0 0 1]) ];
    c = c/v;
  end
  
  
  [varargin,CENTER ] = parseargs(varargin,'center' ,'$FORCE$',{true,false});
  if CENTER
    M.xyz = bsxfun( @minus , M.xyz , c(:).' );
  end
  if nargout > 2
    IT = [ moment([2 0 0])/2 ,  moment([1 1 0])   , moment([1 0 1])   ;...
           0                 ,  moment([0 2 0])/2 , moment([0 1 1])   ;...
           0                 ,  0                 , moment([0 0 2])/2 ];
    IT = IT + IT.';
  end

  
  function mo = moment( pqr )
    p1 = pqr(1);  p2 = pqr(2);  p3 = pqr(3);
    P = p1 + p2 + p3;

    N_P = bsxfun( @rdivide , N , [p1 p2 p3]+1 );
    
    mo = zeros( size(M.tri,1) , 1 );
    
    for k1 = 0:p1,  for k2 = 0:p2,  for k3 = 0:p3

          K  = k1 + k2 + k3;
          
          Is = zeros( size(M.tri,1) , 1 );
          for j1 = 0:k1,   for j2 = 0:k2,   for j3 = 0:k3
                J  = j1 + j2 + j3;
                
%                 if any( cellfun( @(s) isequal(s,[j1 j2 j3 k1 k2 k3 P]) , ISC ) )
%                   disp( [j1 j2 j3 k1 k2 k3 P] );
%                 else
%                   ISC{end+1} = [j1 j2 j3 k1 k2 k3 P];
%                 end
                
                Is = Is + comb( [k1 k2 k3] , [j1 j2 j3] ) * factorial( J ) * factorial( K - J ) * ...
                              prod( bsxfun( @power , M.xyz( M.tri(:,1) , : ) , [    j1      j2      j3 ] ) , 2 ) .* ...
                              prod( bsxfun( @power , M.xyz( M.tri(:,2) , : ) , [ k1-j1 , k2-j2 , k3-j3 ] ) , 2 ) .* ...
                              sum( N_P .*  (   M.xyz( M.tri(:,1),: ) * ( 1 + J     ) + ...
                                               M.xyz( M.tri(:,2),: ) * ( 1 - J + K ) + ...
                                               M.xyz( M.tri(:,3),: ) * ( 1 + P - K ) ) , 2 );

          end,  end, end

          mo = mo + comb( [p1 p2 p3] , [k1 k2 k3] ) * factorial( P - K ) * Is .* ...
                     prod( bsxfun( @power , M.xyz( M.tri(:,3) , : ) , [ p1-k1 , p2-k2 , p3-k3 ] ) , 2 );
                   
    end,  end,  end

    mo = 2*sum( mo .* A )/( 3 * factorial( P + 3 ) );
    
    function C = comb( P , K )
      C = nchoosek( P(1) , K(1) ) * nchoosek( P(2) , K(2) ) * nchoosek( P(3) , K(3) );
    end
  end

end
