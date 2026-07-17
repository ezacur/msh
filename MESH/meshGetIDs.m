function id = meshGetIDs( M ,W)

  if nargin < 2
    if isfield( M ,'triID') && isfield( M ,'xyzID'), error('XYZ or TRI must be specified.'); end
    if isfield( M ,'triID'), W = 'tri'; end
    if isfield( M ,'xyzID'), W = 'xyz'; end
  elseif ischar( W )
    switch lower( W )
      case {'t','tri'}, W = 'tri';
      case {'x','xyz'}, W = 'xyz';
    end
  else
    error( '''xyz'', ''x'', ''tri'' or ''t'' are the only valid keywords.');
  end

  if ~isfield( M , [W,'ID'] )
    error( 'No ''%sID'' field exists in M.' , W );
  end

  id = M.([W,'ID']);

end
