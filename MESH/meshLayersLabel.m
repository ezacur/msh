function M = meshLayersLabel( M , f_or_v )

  if ~isstruct( M )
    M = struct( 'tri' , M );
  end
  if nargin < 2
    f_or_v = 'F';
  end

  switch upper( f_or_v )
    case {'F','FACES','T','TRI'}
      M = MeshAddField( M ,'triLayerID',0);

      L = M;
      L = struct( 'tri' , L.tri ,'triID' , (1:size(L.tri,1)).' );
      Lid = 0;
      while 1
        w = meshBoundaryElements( L );
        if ~any( w ), break; end
        Lid = Lid + 1;
        M.triLayerID( L.triID(w) ) = Lid;
        L.tri(   w ,:) = [];
        L.triID( w ,:) = [];
      end
      
      M = M.triLayerID;
    
    otherwise
      error('TRI or XYZ?');
  end




end
