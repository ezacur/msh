function M = MeshCommonFaces( A , B )

  M = A;
  a = A.tri; a = sort( a ,2);
  b = B.tri; b = sort( b ,2);

  w = ismember( a , b ,'rows');
  M = MeshRemoveFaces( M , ~w ,true);


  


end
