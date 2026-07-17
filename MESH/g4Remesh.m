function [R, Q] = g4Remesh( M, el, Nits, varargin )
% G4REMESH High-quality isotropic remeshing using geometry3Sharp.
%
%   [R, Q] = G4REMESH(M, EL, NITS, OPTIONS...) performs iterative remeshing
%   on the input mesh M to achieve a target edge length EL over NITS iterations.
%
%   This function wraps the powerful C# library 'geometry3Sharp' (g4) to provide
%   robust remeshing with advanced constraints, projection, and annealing capabilities.
%
%   INPUTS:
%   -------
%   M      : Input mesh structure. Must contain:
%            .xyz - [N x 3] matrix of vertex coordinates.
%            .tri - [F x 3] matrix of triangle indices.
%
%   el     : Target Edge Length strategy. Controls the density of the result.
%            This argument is flexible to support "Annealing" (gradually reducing
%            edge length to converge stably):
%
%            1. Scalar (> 0) : Constant target length for all iterations.
%            2. Scalar (< 0) : Heuristic density. Calculates target length based
%                              on surface area divided by abs(el).
%            3. Vector [Nits]: Explicit schedule. You can pass a vector defining
%                              the target length for every iteration.
%                              Example: geospace(start_len, end_len, Nits)
%            4. Cell {A, B}  : Auto-Annealing. Generates a geometric progression
%                              from length A to length B over 75% of iterations,
%                              then holds B constant.
%            5. Cell {A,B,N} : Same as above, but interpolates over N steps.
%
%   Nits   : Integer. Number of iterations (passes). Typical range: 10-50.
%
%   OPTIONS (VARARGIN):
%   -------------------
%
%   --- PROJECTION (Crucial for shape preservation) ---
%   'auto'              : Automatically sets the input mesh M as the projection target.
%                         Prevents the mesh from shrinking/losing volume. (Recommended).
%   M_Struct (struct)   : Uses a different mesh structure as the projection target.
%   'AfterRefinement'   : Projects vertices only at the end of a pass (faster) rather
%                         than after every topological op (safer, default 'Inline').
%
%   --- BOUNDARY & CONSTRAINTS (Topology Control) ---
%   The following options control how the mesh boundaries (holes/edges) are handled.
%   
%   'MeshBoundary' (or 'SoftBoundary', 'SlidingBoundary') [RECOMMENDED FOR QUALITY]
%       "The Train Track". Converts boundary edges into smooth 3D mathematical curves.
%       Vertices are constrained to slide along these curves.
%       - Geometry: Preserved (projected on curve).
%       - Topology: Dynamic. Edges can split/collapse to match target density.
%       - Result: High-quality equilateral triangles at the boundary.
%
%   'FixBoundaryEdges' (or 'HardBoundary', 'FixBoundary') [RECOMMENDED FOR MATCHING]
%       "The Brick Wall". Locks all boundary vertices and edges in place.
%       - Geometry: Exact original positions.
%       - Topology: Static. No splitting or collapsing allowed on the boundary.
%       - Result: Exact boundary match, but potential for poor "needle" triangles
%         if the target density differs from the original boundary resolution.
%
%   'FixBoundaryNodes'
%       "The Fence Posts". Locks the specific vertices of the boundary in place,
%       but does NOT lock the edges between them (unless FixBoundaryEdges is also used).
%       - Geometry: Original nodes are fixed anchors.
%       - Topology: Edges between anchors CAN split.
%       - Use case: Combine with 'MeshBoundary' to keep corners fixed while
%         smoothing the segments between them.
%
%   'PreserveBoundaryLoops'
%       Alternative implementation of 'FixBoundaryEdges'. Tries to identify
%       topological loops. Less robust on non-manifold meshes.
%
% 
%   Nx1 fixed indexes
% 
%   Nx3 Matrix (Coordinates)
%       Pass a matrix of coordinates. The function finds the nearest vertices
%       and locks them (Fixed Position). Useful for internal landmarks.
%
%   'tri...' (Region Constraints)
%       If a string starts with 'tri' (e.g., 'triLeftAtrium'), the function looks
%       for M.triLeftAtrium. It extracts the boundary of that sub-region and
%       applies a "Sliding Curve" constraint (like MeshBoundary) to it.
%       Essential for preserving internal boundaries between tissue types.
%
%   --- ALGORITHM SETTINGS ---
%   'SmoothSpeedT', val : (0.0 to 1.0). Controls smoothing aggressiveness.
%                         - 0.05: Conservative (DEFAULT; preserves features).
%                         - 0.5 : Balanced.
%                         - 1.0 : Fluid (Best for 'MeshBoundary' + Projection).
%   'cotan'             : Use Cotangent-weighted smoothing (preserves curvature
%                         better without projection) instead of Uniform smoothing.
%
%   --- PERFORMANCE & DETERMINISM ---
%   'deterministic'     : (Default). Disables parallel processing. Ensures that
%   (or 'noEnable...')    running the code twice yields bit-exact same results.
%   'EnableParallelSmooth': Enables multi-threaded smoothing. Faster on large meshes
%                           but results may vary slightly between runs due to race conditions.
%
%   --- OUTPUT ---
%   'quiet'             : Suppresses console output.
%   'verbose'           : Force verbose output (Default).
%
%   ADVANCED: CONSTRAINT ORDER & LOGIC
%   ----------------------------------
%   Constraints are applied sequentially. If you combine constraints, the order
%   matters significantly (Last one wins/overwrites).
%
%   Scenario: You want high-quality boundaries ('MeshBoundary') but you need
%             specific corners to remain exactly at their original positions.
%
%   CORRECT ORDER:  ..., 'MeshBoundary', 'FixBoundaryNodes', ...
%       1. 'MeshBoundary' turns edges into sliding curves.
%       2. 'FixBoundaryNodes' hammers "nails" into the original vertices.
%       Result: Vertices are anchored, but the edges between them can curve and split.
%
%   INCORRECT ORDER: ..., 'FixBoundaryNodes', 'MeshBoundary', ...
%       1. Vertices are fixed.
%       2. 'MeshBoundary' overwrites the "Fixed" status with "Sliding".
%       Result: Anchors are lost; vertices will slide along the curve.
%
%   OUTPUTS:
%   --------
%   R : Resulting mesh structure.
%   Q : (Optional) Cell array containing the mesh state at every iteration.
%       Warning: Requesting Q significantly slows down execution.
%
if 0
  M = loadv('c:\Dropbox\mTools\Corify_tools\RESOURCES\MARKERS3D.mat','MARKERS3D');
  %M = MeshSubdivide( M ,'loop' );

  [M,ids] = MeshForceNode( M , M.lmk );
 %extent( M.xyz( abs(ids) ,:) - M.lmk )
  
  M = MeshReOrderNodes( M , unique( [ abs(ids(:).') , 1:size( M.xyz ,1) ] ,'stable' ) );
  extent( M.xyz( 1:128 ,:) - M.lmk )

  while 1, try, M = MeshFillSmallestHole( M ); catch, break; end; end

  el = 10;
  R = g4Remesh( M , el , 1000 , M , 1:128 );   % fix original nodes 1..128

  plotMESH( R )

  %%
end

  if isnumeric( el ) && isscalar(el) && el < 0
    if mod( -el , 1 ), error('g4Remesh:el','A negative EL (target triangle count) must be an integer.'); end
    el = sqrt( meshSurface( M )/abs( el ) * 2 / sind(60) );
  elseif iscell( el ) && numel( el ) == 2
    Nits = max( Nits , 30 );                                    % annealing needs enough iterations
    el = geospace( el{1} , el{2} , floor( Nits * 0.75 / 10 ) );
    el = repmat( el , 10 , 1 );
    el = el(:); el( end:Nits ,:) = el(end);
  elseif iscell( el ) && numel( el ) == 3
    plateauLength = el{3};
    Nits = max( Nits , ceil( 2 * plateauLength / 0.75 ) );      % ensure >=2 annealing levels
    el = geospace( el{1} , el{2} , floor( Nits * 0.75 / plateauLength ) );
    el = repmat( el , plateauLength , 1 );
    el = el(:); el( end:Nits ,:) = el(end);
  end
  el = double( el );

  Moriginal = MeshGenerateIDs( M , 'xyz' );          % stable per-vertex IDs (1..N) BEFORE tidy
  M  = MeshTidy( Mesh( Moriginal ) ,NaN,true);       % xyzID is carried through the renumbering
  M0 = Mesh( M ,0);
  G  = G4( M );

  %map an ORIGINAL (user-facing) vertex index -> tidied index (0 if MeshTidy dropped it)
  orig2tidy = zeros( size(Moriginal.xyz,1) , 1 );
  orig2tidy( M.xyzID ) = 1:numel( M.xyzID );
  R  = g4.Remesher(G);
  
  R.ProjectionMode                        = g4.('Remesher+TargetProjectionMode').Inline;
  R.SmoothType                            = g4.('Remesher+SmoothTypes').Uniform;
  R.AllowCollapseFixedVertsWithSameSetID  = true;
  R.EnableSmoothInPlace                   = false;
  R.EnableFlips                           = true;
  R.EnableSplits                          = true;
  R.EnableCollapses                       = true;
  R.EnableSmoothing                       = true;
  R.PreventNormalFlips                    = true;

  FVcoords = []; % FixedVerticesCoordinates (VertexConstraint added)
  constraintID         = 1;     % constraints
  SmoothSpeedT         = 0.05;  % antes era 0.25
  EnableParallelSmooth = [];
  VERBOSE              = true;
  while ~isempty( varargin )
    V = varargin{1}; varargin(1) = [];
    if 0
    elseif ischar( V ) && ( strcmpi( V , 'SmoothSpeedT' ) || strcmpi( V , 'Speed') )
      SmoothSpeedT = varargin{1}; varargin(1) = [];

    elseif ischar( V ) && strcmpi( V , 'quiet' )
      VERBOSE = false;

    elseif ischar( V ) && strcmpi( V , 'verbose' )
      VERBOSE = true;

    elseif ischar( V ) && ( strcmpi( V , 'noEnableParallelSmooth' ) || strcmpi( V , 'deterministic' ) )
      EnableParallelSmooth = false;   %to ensure determinism..

    elseif ischar( V ) && strcmpi( V , 'EnableParallelSmooth' )
      EnableParallelSmooth = true;    %by default it is true.

    elseif isstruct( V )
      P = G4( MeshTidy( Mesh( V ,0) ,NaN,true) );
      R.SetProjectionTarget( g4.MeshProjectionTarget(P) );

    elseif ischar( V ) && strcmpi( V ,'auto' )
      G0 = g4.DMesh3(G);
      R.SetProjectionTarget( g4.MeshProjectionTarget.Auto( G0 ) );

%     elseif isnumeric( V ) && isvector( V ) && all( V >= 1) && all( ~mod(V,1) ) %fix the points in a single constraint (they can be eventually be collapsed if AllowCollapseFixedVertsWithSameSetID is true)
%       error('esto hay que revisarlo, ahora utilizo R.Constraints.SetOrUpdate.');
%       if isempty(C), C = g4.MeshConstraints(); end
%       V = unique( V );
%       constraintID = constraintID + 1;
%       for v = 1:numel( V )
%         VC(end+1,:) = M.xyz( V(v) ,:);
%         C.SetOrUpdateVertexConstraint( V(v)-1 , g4.VertexConstraint( true , constraintID ) );
%       end

%     elseif iscell( V ) && numel(V) == 1 && isvector( V{1} )    %fix the points in different constraints (it is safer, even if AllowCollapseFixedVertsWithSameSetID is true)
%       error('esto hay que revisarlo, ahora utilizo R.Constraints.SetOrUpdate.');
%       if isempty(C), C = g4.MeshConstraints(); end
%       V = V{1};
%       for v = 1:numel( V )
%         constraintID = constraintID + 1;
%         VC(end+1,:) = M.xyz( V(v) ,:);
%         C.SetOrUpdateVertexConstraint( V(v)-1 , g4.VertexConstraint( true , constraintID ) );
%       end

    elseif isnumeric( V ) && isvector( V ) && all( V > 0 ) && ~any( mod( V ,1) )
      if isempty( R.Constraints ), R.SetExternalConstraints( g4.MeshConstraints() ); end
      V = V(:).';
      if any( V > numel( orig2tidy ) ), error('g4Remesh:constraint','Constrained vertex index out of range (original mesh has %d vertices).', numel(orig2tidy) ); end
      tv = orig2tidy( V );                                     % ORIGINAL index -> tidied index
      if any( tv == 0 ), error('g4Remesh:constraint','%d constrained vertices were removed by MeshTidy.', sum(tv==0) ); end
      FVcoords = [ FVcoords ; Moriginal.xyz( V ,:) ];
      for v = tv(:).'
        constraintID = constraintID + 1;
        R.Constraints.SetOrUpdateVertexConstraint( int32( v - 1 ) , g4.VertexConstraint( true , constraintID ) );
      end

    elseif isnumeric( V ) && size( V ,2) == 3
      [ov,d] = knnsearch( Moriginal.xyz , V );                 % validate against the ORIGINAL mesh
      if any( d > 1e-10 ), error('Constrained vertices coordinates seem to be too far. Try with MeshForceNode.'); end
      tv = orig2tidy( ov );
      if any( tv == 0 ), error('g4Remesh:constraint','%d constrained vertices were removed by MeshTidy.', sum(tv==0) ); end
      if isempty( R.Constraints ), R.SetExternalConstraints( g4.MeshConstraints() ); end
      FVcoords = [ FVcoords ; Moriginal.xyz( ov ,:) ];
      for v = tv(:).'
        constraintID = constraintID + 1;
        R.Constraints.SetOrUpdateVertexConstraint( int32( v - 1 ) , g4.VertexConstraint( true , constraintID ) );
      end

    elseif ischar( V ) && ( strcmpi( V , 'MeshBoundary' ) || strcmpi( V , 'SoftBoundary' ) || strcmpi( V , 'SlidingBoundary' ) )
      LS = g4.MeshBoundaryLoops( G );
      for l = 0:( LS.Count-1 )
        L = LS.Item(l);
        constraintID = constraintID + 1;
        g4.MeshConstraintUtil.ConstrainVtxLoopTo( R , L.Vertices , g4.DCurveProjectionTarget( L.ToCurve ) , constraintID );
      end
    
    elseif ischar( V ) && ( strcmpi( V , 'FixBoundaryEdges' ) || strcmpi( V , 'FixBoundary' ) || strcmpi( V , 'HardBoundary' ) )
      FVcoords = [ FVcoords ; M.xyz( meshBoundaryNodes(M) , :) ];
      g4.MeshConstraintUtil.FixAllBoundaryEdges(R);

    elseif ischar( V ) && strcmpi( V ,'PreserveBoundaryLoops' )
      FVcoords = [ FVcoords ; M.xyz( meshBoundaryNodes(M) , :) ];
      g4.MeshConstraintUtil.PreserveBoundaryLoops( R );
% Intenta rastrear la linea de frontera para identificar "Bucles" cerrados. Una vez identificado el bucle, fija sus vertices.
% Topologica global. Intenta entender la forma del agujero.
% Diferencia con la anterior: Si tu malla tiene un borde "sucio" (ej. un vertice compartido por dos agujeros, tipo "8"), esta funcion puede fallar al intentar seguir el camino, mientras que FixBoundaryEdges funcionaria igual.
% Hace lo mismo que la anterior, pero es menos robusta ante geometria mala.
      
    elseif ischar( V ) && strcmpi( V , 'FixBoundaryNodes' )
      V = find( meshBoundaryNodes(M) );
      if isempty( R.Constraints ), R.SetExternalConstraints( g4.MeshConstraints() ); end
      FVcoords = [ FVcoords ; M.xyz( V ,:) ];
      for v = V(:).'
        constraintID = constraintID + 1;
        R.Constraints.SetOrUpdateVertexConstraint( int32( v - 1 ) , g4.VertexConstraint( true , constraintID ) );
      end

    elseif isnumeric( V ) && size( V ,2) == 2
      if isempty( R.Constraints ), R.SetExternalConstraints( g4.MeshConstraints() ); end
      uv  = unique( V(:) ).';
      tuv = orig2tidy( uv );                                   % ORIGINAL indices -> tidied
      if any( tuv == 0 ), error('g4Remesh:constraint','%d constrained vertices were removed by MeshTidy.', sum(tuv==0) ); end
      FVcoords = [ FVcoords ; Moriginal.xyz( uv ,:) ];
      for v = tuv(:).'
        constraintID = constraintID + 1;
        R.Constraints.SetOrUpdateVertexConstraint( int32( v - 1 ) , g4.VertexConstraint( true , constraintID ) );
      end
      for f = 1:size(V,1)
        eid = G.FindEdge( int32( orig2tidy( V(f,1) ) - 1 ) , int32( orig2tidy( V(f,2) ) - 1 ) );
        if eid >= 0
          R.Constraints.SetOrUpdateEdgeConstraint( eid , g4.EdgeConstraint( g4.EdgeRefineFlags.FullyConstrained ) );
        end
      end
    
    elseif ischar( V ) && numel( V ) > 3 && strncmp( V , 'tri',3) && isfield( M , V )
      F = M.(V);
      [~,~,F] = unique( F );
      for p = unique( F(:).' )
        M = MeshGenerateIDs( M ,'xyz');
        REGION = MeshRemoveFaces( M , F ~= p );
        %disp( unique( REGION.(V) ) )
        REGION = meshSeparate( REGION );
        for r = 1:numel(REGION)
          RG = REGION{r};
          RG.xyz = M.xyz;
          RG.tri = reshape( RG.xyzID( RG.tri ) ,[],3);

          G4REG        = G4( RG );
          LS = g4.MeshBoundaryLoops( G4REG );
          for l = 0:( LS.Count-1 ), constraintID = constraintID + 1;
            L = LS.Item(l);
            g4.MeshConstraintUtil.ConstrainVtxLoopTo( R , L.Vertices  , g4.DCurveProjectionTarget( L.ToCurve ) , constraintID );
            %hplot3d( RG.xyz( double(L.Vertices) + 1 ,:) , '.-r');
          end

%           % --- Open spans (región toca el borde de la malla) ---      <-- NUEVO
%           for s = 0:( LS.SpanCount-1 ), constraintID = constraintID + 1;
%             S = LS.Spans.Item(s);
%             g4.MeshConstraintUtil.ConstrainVtxLoopTo( R , S.Vertices , g4.DCurveProjectionTarget( S.ToCurve ) , constraintID );
%           end

        end
      end

    elseif ischar( V ) && strcmpi( V ,'cotan' )
      R.SmoothType = g4.('Remesher+SmoothTypes').Cotan;

    elseif ischar( V ) && strcmpi( V ,'AfterRefinement' )
      R.ProjectionMode = g4.('Remesher+TargetProjectionMode').AfterRefinement;

    else, error('Invalid argument (%s).' , V );
    end
  
  end

  if isempty( EnableParallelSmooth ) || ~EnableParallelSmooth
    % DETERMINISTIC (also the default): disable EVERY parallel path. The .NET
    % Remesher defaults to EnableParallelSmooth=true AND, crucially,
    % EnableParallelProjection=true (probed on the shipped DLL) -- the explicit
    % 'deterministic' flag used to switch off only the smoother, leaving the
    % projection racing: with a projection target ('auto'!) two runs could
    % differ despite the flag's bit-exact promise.
    R.EnableParallelSmooth     = false;
    R.EnableParallelProjection = false;
    R.EnableSmoothInPlace      = false;
  else
    R.EnableParallelSmooth     = true;   %explicit 'EnableParallelSmooth': speed
  end                                    %over reproducibility, projection stays
                                         %at the .NET default (parallel too)
  R.SmoothSpeedT = SmoothSpeedT;

  R.Precompute(); 
  
  last_el = NaN;
  if nargout > 1, Q = cell(Nits,2); end
  for k = 1:Nits
    this_el = el( min( k , end ) );
    if VERBOSE && ( last_el ~= this_el || ( Nits >= 50 && ~rem( k , 5 ) ) ), fprintf('Iteration  %3d of %d (edgeLength: %g).\n', k , Nits , this_el ); end
    if last_el ~= this_el
      R.SetTargetEdgeLength( this_el ); last_el = this_el;
    end
    R.BasicRemeshPass();
    if nargout > 1
      Q{k,1} = PostProcess( R , M0 , M , FVcoords );
      Q{k,2} = this_el;
    end
  end

  R = PostProcess( R , M0 , M , FVcoords );
end

function P = PostProcess( R , M0 , M , FVcoords )
  
  P = g4ToMesh( R.Mesh );
  if ~isempty( FVcoords )
    % Buscamos dónde han acabado los vértices fijos en la malla final
    [ids, d] = knnsearch( P.xyz , FVcoords );
    w = d < 1e-10;
  
    % Tolerancia de seguridad (1e-10 es razonable para doubles)
    if any( ~w )
      warning('g4Remesh:ConstraintsDrift', 'WARNING! %d fixed vertices have drifted. Max Drift: %e. Check the order of MeshBoundary vs FixNodes.' , sum(~w) , max(d(~w)) );
    end

    P.xyz( ids(w) , : ) = FVcoords(w,:);
  end

  if ~~numel( fieldnames( M , '^tri.+') )
    e = vtkClosestElement( M0 , meshFacesCenter( P ) );
  end
  for f = fieldnames( M , '^tri.+').', f = f{1};
    P.(f) = M.(f)(e,:,:,:,:,:);
  end
  
end


function G = G4( M )
  try,      evalc( 'g4.DMesh3();' );
  catch,    g4DLL = fullfile( fileparts(mfilename('fullpath')) , 'geometry4Sharp.dll' );
            asm = NET.addAssembly( g4DLL );
  end

  %try,   G = g4.AcorysMatlabExtensions.CreateMesh( doble( vec( M.xyz.' ) ) , vec( int32( M.tri - 1).' ) ); return; end

  G = g4.DMesh3();

  V = double( M.xyz );
  for i= 1:size( M.xyz ,1)
      G.AppendVertex( g4.Vector3d( V(i,1) , V(i,2) , V(i,3) ) );
  end
  
  F = int32( M.tri - 1);
  for i= 1:size( M.tri ,1)
    G.AppendTriangle( F(i,1) , F(i,2) , F(i,3) );
  end

end
function M = g4ToMesh( M )
  G = g4.DMesh3( M );
  G.CompactInPlace();
  
  V = double( G.VerticesBuffer.GetBuffer() );
  V = V( 1:( G.MaxVertexID * 3 ) );
  V = reshape( V ,3,[]).';
  
  F = double( G.TrianglesBuffer.GetBuffer() );
  F = F( 1:( G.MaxTriangleID  * 3 ) );
  F = reshape( F ,3,[]).';
  F = F+1;
  
  M = struct( 'xyz' , double(V) , 'tri' , double(F) );
end
