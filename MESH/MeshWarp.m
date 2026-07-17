function M = MeshWarp( M , d )

  if size( d ,2) == 1
    if isfield( M ,'xyzNORMALS')
      d = M.xyzNORMALS .* d;
    else
      d = meshNormals( M ,'u') .* d;
    end
  end

  M.xyz = M.xyz + d;

end
