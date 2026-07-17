function M = MeshFlipFaces( M )

  M.tri = M.tri( : ,[1:end-2,end,end-1]);

end
