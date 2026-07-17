function [M,B] = MeshFillSmallestHole( M , el , maxL )

  if nargin < 2, el = []; end
  if nargin < 3, maxL = Inf; end

  makeAll = false;
  if islogical( el ) && ~~el
    makeAll = true;
    el = [];
  end
  if isempty( el )
    [~,el] = meshEdges( M );
    el = median( el );
  end

  if makeAll
    while 1, try, M = MeshFillSmallestHole( M , el ); catch, break; end; end
    return;
  end

  B = MeshBoundary( M );
  if ~isfield( B ,'tri') || isempty( B.tri )
    error('No holes.');
  end
  B = meshSeparate( B );
  Ls = cellfun( @(b)sum(meshQuality(b,'length')) , B );
  if all( Ls > maxL ), error('No small holes.'); end
  B = B{ argmin( Ls ) };

  if size( B.xyz ,1) == 3 && size( B.tri ,1) == 3
    B.tri = [1,2,3];
  else
    H = B;
    B = Mesh( fillContoursMesh( B , el ) ,0);
    for f = fieldnames( H ).', f = f{1};
      if strcmp( f , 'xyz' ), continue; end
      if ~strncmp( f , 'xyz',3), continue; end
      if all( H.(f) == H.(f)(1,:) )
        B = MeshAddField( B , f , H.(f)(1,:) );
      end
    end
  end
  M = MeshWeld( M , B );

end
