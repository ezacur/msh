function M = MeshColor( M , varargin )
%MESHCOLOR  Attach display options to a mesh struct (read later by plotMESH).
%
%   M = MeshColor( M , ... )  writes appearance options into FIELDS of the
%   mesh struct M (it does NOT draw anything). plotMESH / hplotMESH scan those
%   fields and turn them into patch properties (M.color -> FaceColor, and any
%   field named like a patch property -> that property). So the mesh carries
%   its own default look:
%       M = MeshColor( sphereMesh() , 'r' , 0.3 , "outer shell" );
%       plotMESH( M );                       % drawn red, alpha 0.3, that legend
%
%   Arguments are consumed left to right; each is dispatched by TYPE:
%       [r g b]      (numeric, 3 elems)  -> M.color      (the face color)
%       a            (numeric scalar)    -> M.FaceAlpha  (opacity 0..1)
%       "text"       (string, dbl-quote) -> M.DisplayName (legend label)
%       'r'|'g'|'b'  (char)              -> M.color = red/green/blue
%       'Name', val  (char + next arg)   -> M.(Name) = val  (any patch prop,
%                                           e.g. 'EdgeColor',[0 0 0])
%
%   NOTES / gotchas:
%     * color shortcuts are ONLY r/g/b (no k/w/c/m/y); for others pass a
%       triple ([1 1 0]) or 'FaceColor',spec.
%     * char vs string matters: 'r' (single quotes) = red, but "r" (double
%       quotes, a string) becomes the DisplayName "r".
%     * the positional numeric shortcuts are rigid: a 3-elem vector is ALWAYS
%       read as color and a scalar ALWAYS as FaceAlpha. Numerics of any other
%       size (e.g. an N x 3 per-vertex color matrix) are SILENTLY IGNORED --
%       pass those as an explicit 'FaceVertexCData',C pair instead.
%     * a trailing 'Name' with no following value errors (it consumes the
%       next argument).
%
% See also plotMESH, hplotMESH, Mesh.

  while ~isempty( varargin )
    V = varargin{1}; varargin(1) = [];

    if isnumeric(V) && numel(V) == 3
      M.color = V; continue; end
    
    if isnumeric(V) && numel(V) == 1
      M.FaceAlpha = V; continue; end
    
    if isstring(V)
      M.DisplayName = char(V); continue; end

    if ischar(V)
      switch V
        case 'r', M.color = [1,0,0]; continue;
        case 'g', M.color = [0,1,0]; continue;
        case 'b', M.color = [0,0,1]; continue;
      end; end

    if ischar(V)
      P = V; V = varargin{1}; varargin(1) = [];
      M.(P) = V; continue; end

  end

end
