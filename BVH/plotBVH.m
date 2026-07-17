function h = plotBVH( B , M , depth0 )
%plotBVH  Interactive viewer of a BVH blob: the mesh + the node volumes.
%
%   plotBVH( B )           draws the blob's own geometry (B.X/B.Tri) and the
%                             node volumes (AABB wireframes or sphere circles)
%                             of ONE tree depth at a time, in WORLD space (the
%                             global frame is applied, so you see exactly what
%                             the queries see).
%   plotBVH( B , M )       draws mesh M instead (must match the blob).
%   plotBVH( B , M , d0 )  starts at depth d0 (default 1 = the root).
%   h = plotBVH( ... )     returns the figure handle.
%
%   KEYS:
%     up / down     go one level deeper / shallower
%     a             toggle ALL depths at once (color-coded by depth)
%     l             toggle leaf-only view (only leaf slots, at any depth)
%     f             toggle drawing the frame axes (build frame at the centroid;
%                   PCA-aligned when the cloud is anisotropic)
%     r             reset view
%
%   Leaf slots are drawn in RED, internal slots in BLUE (single-depth mode);
%   in all-depths mode each depth gets a color from the colormap. The title
%   shows depth / node count / slot count / volume type / leaf sizes.
%
% See also BVH, bvhClosestElement, bvhIntersectRay.

  if nargin < 2, M = []; end
  if nargin < 3 || isempty( depth0 ), depth0 = 1; end
  if ~isstruct( B ) || ~isfield( B , 'bounds4' )
    error('plotBVH:B','B must be a BVH blob.');
  end

  %node depths (root = 1), children from child4 > 0
  ch = double( B.child4 );
  nN = size( ch ,2);
  dep = zeros( 1 , nN );  dep(1) = 1;
  stack = 1;
  while ~isempty( stack )
    n = stack(end);  stack(end) = [];
    for k = 1:4
      c = ch(k,n);
      if c > 0, dep(c) = dep(n) + 1;  stack(end+1) = c; end %#ok<AGROW>
    end
  end
  maxDep = max( dep );

  %world-space transform of the frame
  Af = B.frame(1:3,1:3);  tf = B.frame(1:3,4).';

  S = struct( 'B',B, 'ch',ch, 'dep',dep, 'maxDep',maxDep, 'Af',Af, 'tf',tf, ...
              'depth',min(max(depth0,1),maxDep), 'all',false, 'leaves',false, ...
              'axes',true, 'hVol',[], 'hFrm',[] );

  hf = figure( 'Name' , sprintf('BVH blob: %d elems, %d nodes, %s' , ...
               B.nE , nN , upper(B.volume) ) , 'Color' , 'w' );
  ax = axes( 'Parent' , hf );  hold( ax , 'on' );

  %---- the mesh itself (world space) ----------------------------------------
  if isempty( M )
    Xw = B.X * Af.' + tf;
    T  = double( B.Tri );
  else
    Xw = double( M.xyz );  Xw(:,end+1:3) = 0;
    T  = double( M.tri );
  end
  k = sum( T > 0 ,2);
  w = k >= 3;                                     %triangles / tets (their faces)
  if any( w )
    F3 = T( w & k==3 , 1:3 );
    if any( k == 4 )
      T4 = T( k==4 , 1:4 );
      F3 = [ F3 ; T4(:,[1 2 3]) ; T4(:,[1 2 4]) ; T4(:,[1 3 4]) ; T4(:,[2 3 4]) ];
    end
    patch( 'Parent',ax, 'Vertices',Xw, 'Faces',F3, ...
           'FaceColor',[0.85 0.85 0.85], 'EdgeColor',[0.6 0.6 0.6], ...
           'FaceAlpha',0.35, 'EdgeAlpha',0.3 );
  end
  w = k == 2;                                     %segments
  if any( w )
    E = T( w , 1:2 );
    xs = [ Xw(E(:,1),1) , Xw(E(:,2),1) , NaN(size(E,1),1) ].';
    ys = [ Xw(E(:,1),2) , Xw(E(:,2),2) , NaN(size(E,1),1) ].';
    zs = [ Xw(E(:,1),3) , Xw(E(:,2),3) , NaN(size(E,1),1) ].';
    line( xs(:) , ys(:) , zs(:) , 'Parent',ax, 'Color',[0.3 0.3 0.3] );
  end
  w = k == 1;                                     %points
  if any( w )
    plot3( ax , Xw(T(w,1),1) , Xw(T(w,1),2) , Xw(T(w,1),3) , '.' , ...
           'Color',[0.3 0.3 0.3] , 'MarkerSize',6 );
  end

  axis( ax , 'equal' );  grid( ax , 'on' );  view( ax , 3 );
  rotate3d( hf , 'on' );

  guidata( hf , S );
  set( hf , 'KeyPressFcn' , @(src,ev) onKey( src , ev , ax ) );
  redraw( hf , ax );
  if nargout > 0, h = hf; end
end

%% --------------------------------------------------------------------------
function onKey( hf , ev , ax )
  S = guidata( hf );
  switch ev.Key
    case 'uparrow',   S.depth = min( S.depth+1 , S.maxDep );  S.all = false;
    case 'downarrow', S.depth = max( S.depth-1 , 1 );         S.all = false;
    case 'a',         S.all    = ~S.all;
    case 'l',         S.leaves = ~S.leaves;
    case 'f',         S.axes   = ~S.axes;
    case 'r',         view( ax , 3 );  axis( ax , 'equal' );
    otherwise, return;
  end
  guidata( hf , S );
  redraw( hf , ax );
end

function redraw( hf , ax )
  S = guidata( hf );
  delete( S.hVol( ishandle( S.hVol ) ) );  S.hVol = [];
  delete( S.hFrm( ishandle( S.hFrm ) ) );  S.hFrm = [];

  if S.all, deps = 1:S.maxDep; else, deps = S.depth; end
  cmap = lines( max( S.maxDep ,7) );
  nSlots = 0;  nNodes = 0;

  for d = deps
    nodes = find( S.dep == d );
    nNodes = nNodes + numel( nodes );
    slots = [];  isleaf = [];
    for n = nodes
      for k = 1:4
        c = S.ch(k,n);
        if c == 0, continue; end
        if S.leaves && c ~= -1, continue; end
        slots(end+1,:)  = [ n , k ];          %#ok<AGROW>
        isleaf(end+1,1) = ( c == -1 );        %#ok<AGROW>
      end
    end
    if isempty( slots ), continue; end
    nSlots = nSlots + size( slots ,1);

    if S.all
      colL = cmap( d ,:);  colI = cmap( d ,:);
    else
      colL = [0.85 0.1 0.1];  colI = [0.1 0.3 0.85];   %leaf red, internal blue
    end
    S.hVol = [ S.hVol , drawSlots( ax , S , slots( isleaf==1 ,:) , colL ) ];
    S.hVol = [ S.hVol , drawSlots( ax , S , slots( isleaf==0 ,:) , colI ) ];
  end

  if S.axes                                    %build-frame axes at the centroid
    L = 0.25 * max( max( S.B.X ,[],1) - min( S.B.X ,[],1) );
    o = S.tf;
    for a = 1:3
      v = S.Af(:,a).' * L;
      S.hFrm(end+1) = line( o(1)+[0 v(1)] , o(2)+[0 v(2)] , o(3)+[0 v(3)] , ...
        'Parent',ax, 'Color', [a==1, a==2, a==3]*0.8 , 'LineWidth',2 );
    end
  end

  if S.all, dtxt = sprintf('ALL (1..%d)', S.maxDep);
  else,     dtxt = sprintf('%d / %d', S.depth , S.maxDep);
  end
  if S.leaves, dtxt = [ dtxt , '  [leaves only]' ]; end
  title( ax , sprintf('depth %s   |   %d nodes, %d slots   |   %s, leaves [%g %g]   (keys: up/down, a, l, f, r)' , ...
         dtxt , nNodes , nSlots , upper(S.B.volume) , S.B.leafSize(1) , S.B.leafSize(2) ) , ...
         'FontName','monospaced' , 'FontSize',9 );
  guidata( hf , S );
end

function h = drawSlots( ax , S , slots , col )
%draw the volumes of the given (node,slot) pairs, transformed to WORLD space.
%kdop slots draw their AABB part (the hull of the 8 extra diagonal planes is
%not rendered); obb slots draw the rotated box from its axes.
  h = [];
  if isempty( slots ), return; end
  b4 = double( S.B.bounds4 );
  if S.B.vol == 6
    %capsules: draw the core segments (the radius is implied)
    n = size( slots ,1);
    XYZ = NaN( 3*n , 3 );
    for i = 1:n
      c = slots(i,1);  k = slots(i,2);
      P0 = [ b4( 4*0+k ,c) , b4( 4*1+k ,c) , b4( 4*2+k ,c) ];
      P1 = [ b4( 4*3+k ,c) , b4( 4*4+k ,c) , b4( 4*5+k ,c) ];
      Q  = [ P0 ; P1 ] * S.Af.' + S.tf;
      XYZ( 3*i-2:3*i-1 ,:) = Q;
    end
    h = line( XYZ(:,1) , XYZ(:,2) , XYZ(:,3) , 'Parent',ax, 'Color',col , 'LineWidth',1.5 );
  elseif S.B.vol == 5
    %RSS: draw the core rectangles (the swept radius is implied)
    n = size( slots ,1);
    XYZ = NaN( 6*n , 3 );
    for i = 1:n
      c = slots(i,1);  k = slots(i,2);
      A3 = [ b4( 4*0+k ,c) b4( 4*1+k ,c) b4( 4*2+k ,c) ;
             b4( 4*3+k ,c) b4( 4*4+k ,c) b4( 4*5+k ,c) ;
             b4( 4*6+k ,c) b4( 4*7+k ,c) b4( 4*8+k ,c) ];
      ul = b4( 4*9+k  ,c);  uh = b4( 4*10+k ,c);
      vl = b4( 4*11+k ,c);  vh = b4( 4*12+k ,c);
      w0 = b4( 4*13+k ,c);
      crn = [ ul vl w0 ; uh vl w0 ; uh vh w0 ; ul vh w0 ; ul vl w0 ] * A3;
      crn = crn * S.Af.' + S.tf;
      XYZ( 6*i-5:6*i-1 ,:) = crn;
    end
    h = line( XYZ(:,1) , XYZ(:,2) , XYZ(:,3) , 'Parent',ax, 'Color',col , 'LineWidth',1 );
  elseif S.B.vol ~= 1
    %box wireframes (aabb, kdop's aabb part, obb): 8 corners + corner walk
    n = size( slots ,1);
    path = [ 1 2 4 3 1 5 6 8 7 5 NaN 2 6 NaN 4 8 NaN 3 7 ];  %corner walk
    XYZ = NaN( numel(path)*n , 3 );
    cube = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 1 1 0 ; 0 0 1 ; 1 0 1 ; 0 1 1 ; 1 1 1 ];
    for i = 1:n
      c = slots(i,1);  k = slots(i,2);
      if S.B.vol == 3
        A3 = [ b4( 4*0+k ,c) b4( 4*1+k ,c) b4( 4*2+k ,c) ;
               b4( 4*3+k ,c) b4( 4*4+k ,c) b4( 4*5+k ,c) ;
               b4( 4*6+k ,c) b4( 4*7+k ,c) b4( 4*8+k ,c) ];   %rows = axes
        lo = [ b4( 4*9+k  ,c) , b4( 4*10+k ,c) , b4( 4*11+k ,c) ];
        hi = [ b4( 4*12+k ,c) , b4( 4*13+k ,c) , b4( 4*14+k ,c) ];
        crn = ( lo + cube .* ( hi - lo ) ) * A3;   %axis coords -> build space
      else
        lo = [ b4(   k ,c) , b4( 8+k ,c) , b4(16+k ,c) ];
        hi = [ b4( 4+k ,c) , b4(12+k ,c) , b4(20+k ,c) ];
        crn = lo + cube .* ( hi - lo );
      end
      crn = crn * S.Af.' + S.tf;               %build -> world
      seg = NaN( numel(path) ,3);
      w = ~isnan( path );
      seg(w,:) = crn( path(w) ,:);
      XYZ( (i-1)*numel(path)+(1:numel(path)) ,:) = seg;
    end
    XYZ = XYZ( 1:find( any(~isnan(XYZ),2) ,1,'last') ,:);
    h = line( XYZ(:,1) , XYZ(:,2) , XYZ(:,3) , 'Parent',ax, 'Color',col , 'LineWidth',1 );
  else
    %spheres: 3 great circles each (cheap wire representation)
    tt = linspace( 0 , 2*pi , 25 ).';
    ring = [ cos(tt) , sin(tt) , 0*tt ];
    n = size( slots ,1);
    XYZ = NaN( n*3*(numel(tt)+1) , 3 );
    row = 1;
    for i = 1:n
      c = slots(i,1);  k = slots(i,2);
      ce = [ b4( k ,c) , b4( 4+k ,c) , b4( 8+k ,c) ];
      r  = b4( 12+k ,c);
      if r < 0, continue; end
      for a = 1:3
        R = ring( : , circshift( 1:3 , a ) ) * r + ce;
        R = R * S.Af.' + S.tf;
        XYZ( row:row+numel(tt)-1 ,:) = R;
        row = row + numel(tt) + 1;
      end
    end
    XYZ = XYZ( 1:max(row-1,1) ,:);
    h = line( XYZ(:,1) , XYZ(:,2) , XYZ(:,3) , 'Parent',ax, 'Color',col );
  end
end
