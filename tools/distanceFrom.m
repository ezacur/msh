function [d,cp,eid,bc] = distanceFrom( xyz , M , getBOUNDARY , QUIET )
%DISTANCEFROM  Distance from points to a mesh / wireframe / point cloud.
%
%   [d,cp,eid,bc] = distanceFrom( xyz , M )
%   [d,cp,...]    = distanceFrom( xyz , M , getBOUNDARY , QUIET )
%
%   Dispatches on M:
%     * triangle mesh (celltype 5) -> vtkClosestElement (batched with progress
%       for > 5e5 queries; QUIET silences it). eid = closest face, bc = its
%       barycentric coordinates (N x 3).
%     * wireframe (celltype 3)     -> d2Wireframe. eid = closest segment,
%       bc = [w1 w2] (N x 2).
%     * numeric point cloud (nx3)  -> closest point (vtk or knnsearch).
%     * polyline object            -> ClosestElement.
%
%   getBOUNDARY (triangle meshes and wireframes): if true, d is NEGATED at the
%   queries whose closest point lies ON the OPEN boundary of M -- for a triangle
%   surface (celltype 5) an open-boundary edge or vertex; for a wireframe
%   (celltype 3) a FREE END (degree-1 node) -- i.e. the query projects OUTSIDE
%   the surface patch / BEYOND the polyline's open end (used to reject boundary
%   matches in registration). |d| is unchanged; closed meshes and closed loops
%   are unaffected. The boundary data is memoized (persistent) per connectivity.
%
%   Non-finite query rows yield NaN outputs. Empty xyz / empty M return NaNs.
%
% See also d2Wireframe, vtkClosestElement, meshCelltype, MeshBoundary.

  persistent storedBOUNDARIES
  if isempty( storedBOUNDARIES ), storedBOUNDARIES = cell(0,3); end

  if isempty( xyz )
    n   = 0;
    d   = NaN(n,1);
    cp  = NaN(n,3);
    eid = NaN(n,1);
    bc  = NaN(n,3);
    return;
  end
  if isstruct( M ) && ( isempty( M.xyz ) || isempty( M.tri ) )
    n   = size( xyz ,1);
    d   = NaN(n,1);
    cp  = NaN(n,3);
    eid = NaN(n,1);
    bc  = NaN(n,3);
    return;
  end


  if nargin < 3 || isempty( getBOUNDARY ), getBOUNDARY = false; end
  if nargin < 4, QUIET = false; end

  CT = 0;
  if isstruct( M ), CT = meshCelltype( M ); end
  if ~( isscalar(CT) && ( CT == 3 || CT == 5 ) ), getBOUNDARY = false; end   %boundary rejection only for wireframes (3) & triangle surfaces (5)

  try, M.tri = double( M.tri ); end
  if getBOUNDARY
    b = find( cellfun( 'prodofsize' , storedBOUNDARIES(:,1) ) == numel( M.tri ) );
    b = b( arrayfun( @(c)isequal( storedBOUNDARIES{c,1} , M.tri ) , b ) );

    if ~isempty( b )
      b = b(1);
      BOUNDARYelements = storedBOUNDARIES{b,2};
      BOUNDARYedges    = storedBOUNDARIES{b,3};
    elseif CT == 3
      %wireframe: the boundary is the set of FREE ENDS (degree-1 nodes).
      %MeshBoundary returns exactly those node ids -> slot 3 holds NODE IDS here
      %(not edges); slot 2 flags segments touching one (only for the any()
      %early-out below).
      BOUNDARYedges    = unique( MeshBoundary( M.tri ) );
      BOUNDARYelements = any( ismember( M.tri , BOUNDARYedges ) ,2);

      b = size( storedBOUNDARIES ,1) + 1;
      storedBOUNDARIES{b,1} = M.tri;
      storedBOUNDARIES{b,2} = BOUNDARYelements;
      storedBOUNDARIES{b,3} = BOUNDARYedges;
    else
      %triangle surface: slot 3 holds the open-boundary edges plus [0 u] vertex
      %sentinels; slot 2 flags faces touching the boundary.
      BOUNDARYedges = MeshBoundary( M.tri );
      BOUNDARYedges = sort( BOUNDARYedges ,2);
      u = unique( BOUNDARYedges );
      BOUNDARYedges( end + (1:numel(u)) , 2 ) = u;
      BOUNDARYedges = sortrows( BOUNDARYedges );

      BOUNDARYelements = any( ismember(  M.tri , u ) ,2);

      b = size( storedBOUNDARIES ,1) + 1;
      storedBOUNDARIES{b,1} = M.tri;
      storedBOUNDARIES{b,2} = BOUNDARYelements;
      storedBOUNDARIES{b,3} = BOUNDARYedges;
    end
  end

  
  if getBOUNDARY && ~any( BOUNDARYelements ), getBOUNDARY = false; end
  getBC = getBOUNDARY || nargout > 3;
  
  if     0
  elseif isstruct( M )  &&  CT == 3

    [d, eid, cp, bc] = d2Wireframe( xyz , struct('xyz',double(M.xyz),'tri',double(M.tri) ) );

    if getBOUNDARY
      %negate d where the closest point falls on a FREE END (degree-1 node) of
      %the wireframe -- i.e. the query projects BEYOND the polyline's open end.
      %cp lands on a node exactly when one barycentric weight is ~0 (d2Wireframe
      %clamps t to the endpoint); that node is M.tri(eid,argmax(bc)). Flip d if
      %it is a free end. BOUNDARYedges holds the free-end node ids here.
      w        = isfinite( eid ) & all( isfinite( bc ) ,2);
      [~,ep]   = max( bc ,[],2 );                        %1 -> M.tri(eid,1), 2 -> M.tri(eid,2)
      onNode   = false( numel(eid) ,1);
      onNode(w) = min( bc(w,:) ,[],2 ) <= 1e-5;          %the other weight ~0 => cp sits on a vertex
      node     = zeros( numel(eid) ,1);
      node(w)  = M.tri( sub2ind( size(M.tri) , eid(w) , ep(w) ) );
      hit      = onNode & ismember( node , BOUNDARYedges );
      d( hit ) = -d( hit );
    end

  elseif isstruct( M )  &&  CT == 5

    nxyz = size( xyz , 1 );
    
    eid = NaN( nxyz , 1 );
    cp  = NaN( nxyz , 3 );
    d   = NaN( nxyz , 1 );
    if getBC
    bc  = NaN( nxyz , 3 );
    end
    gcw = []; try, gcw = getCurrentWorker(); end
    if ~isempty( gcw ) || nxyz < 5e5
      w = all( isfinite( xyz ) ,2);
      if any(w)


        if getBC, [ eid(w,1) , cp(w,:) , d(w,1) , bc(w,:) ] = vtkClosestElement( M , xyz(w,:) );
        else,     [ eid(w,1) , cp(w,:) , d(w,1)           ] = vtkClosestElement( M , xyz(w,:) );
        end
      end
    else
      bunchSize = 1e5;

      vtkClosestElement( [] , [] );
      vtkClosestElement( M ); CLEAN = onCleanup( @()vtkClosestElement( [] , [] ) );
      for e = 1:bunchSize:nxyz
        w = e + ( 0:bunchSize-1 ); w = w( w <= nxyz ); w = w( all( isfinite( xyz(w,:) ) ,2) );
        if any(w)
          if getBC, [ eid(w,1) , cp(w,:) , d(w,1) , bc(w,:) ] = vtkClosestElement( xyz(w,:) );
          else,     [ eid(w,1) , cp(w,:) , d(w,1)           ] = vtkClosestElement( xyz(w,:) );
          end
          if ~QUIET
            fprintf('(%9d - %9d   of   %9d)  %g %% done\n' , w(1) , w(end) , nxyz , w(end)/nxyz * 100 );
          end
        end
      end
      clearvars( 'CLEAN' );
    end
      
    
    if getBOUNDARY
      %negate d where the closest point lies ON the open boundary (edge or
      %vertex). Everything stays in FULL-row space: compacting through w and
      %then indexing d with subset positions used to flip the WRONG rows when
      %xyz carried non-finite points.
      w = isfinite( eid );
      e = false( numel(eid) ,1);
      e(w) = BOUNDARYelements( eid(w) );      %closest face touches the boundary
      if ~any( e ), return; end

      b = bc > 1e-5;                          %non-negligible barycentric coords
      b( ~e ,:) = false;                      %only boundary faces
      b( all(b,2) ,:) = false;                %strictly interior cp -> not on an edge

      B = zeros( numel(eid) , 3 );
      B(w,:) = M.tri( eid(w) ,:);
      B = B .* b;                             %node ids where the cp lies (0 elsewhere)
      B = sort( B ,2);
      B = B( : ,2:3);                         %[a b] edge, [0 u] vertex, [0 0] none

      k = find( B(:,2) );
      k = k( ismember( B(k,:) , BOUNDARYedges ,'rows' ) );
      d( k ) = -d( k );
    end
    
  elseif isnumeric( M )  &&  size( M ,2) == 3
    
    try
      [~,cp,d] = vtkClosestPoint( struct('xyz',double(M)) , double( xyz ) );
    catch
      [n,d] = knnsearch( double(M) , double(xyz) ,'K',1);
      cp = M(n,:);
    end
    
  elseif isa( M , 'polyline' )
    
    %[ ~ , cp , d ] = closestElement( M , xyz );
    [ ~ , cp , d ] = ClosestElement( double( M ) , double( xyz ) );
    %[ ~ , cp , d ] = closestElement( M , double(xyz) );
    
%     error('not implemented yet');
    
  end



end
