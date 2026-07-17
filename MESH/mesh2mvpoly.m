function [ MP , L ] = mesh2mvpoly( M , type )


  if ~isfield( M , 'xyz' ), error('Coordinates must be specified in M.xyz'); end
  if ~isfield( M , 'tri' ), error('Connectivity must be specified in M.tri'); end

  NSD = size( M.xyz , 2 );
  switch lower(type)
    case 'linear'
      if size( M.tri , 2 ) ~= 2, error('Connectivity must have 2 columns for the LINEAR case.'); end
    
      MAT = [1 0;-1 1];
    
      f = 1; MP =   mvpoly( MAT * M.xyz( M.tri(f,:) , : ) , 1 , 1 );
      for f = 2:size(M.tri,1)
        MP = [ MP ; mvpoly( MAT * M.xyz( M.tri(f,:) , : ) , 1 , 1 ) ];
      end
      
      if nargout > 1
        Nxyz = size( M.xyz , 1 );
        L = [];
        for f = 1:size(M.tri,1)
          L = [ L ; kron(  eye(NSD) , MAT ) * kron( eye(NSD) , sparse( 1:2 , M.tri(f,:) , 1 , 2 , Nxyz ) ) ];
        end
      end
    
    case 'hermite'
      if size( M.tri , 2 ) ~= 2, error('Connectivity must have 2 columns for the LINEAR case.'); end
      if ~isfield( M , 'xyzd'  ), error('Derivatives must be specified in M.xyzd'); end
    
      MAT = [1,0,-0,-0;-0,0,1,-0;-3,3,-2,-1;2,-2,1,1];
    
      f = 1; MP =   mvpoly( MAT * [ M.xyz( M.tri(f,:) , : ) ; M.xyzd( M.tri(f,:) , : ) ] , 3 , 1 );
      for f = 2:size(M.tri,1)
        MP = [ MP ; mvpoly( MAT * [ M.xyz( M.tri(f,:) , : ) ; M.xyzd( M.tri(f,:) , : ) ] , 3 , 1 ) ];
      end
      
      if nargout > 1
        %L * vec([ M.xyz ; M.xyzd ]) - MP.coeffs(:)
        Nxyz = size( M.xyz , 1 );
        L = sparse([]);
        for f = 1:size(M.tri,1)
          tL = [];
          tL = [ tL ; sparse( 1:2 , M.tri(f,:)          , 1 , 2 , 2*Nxyz ) ];
          tL = [ tL ; sparse( 1:2 , M.tri(f,:) +   Nxyz , 1 , 2 , 2*Nxyz ) ];

          L = [ L ; kron( speye(NSD) , MAT ) * kron( speye(NSD) , tL ) ];
        end
      end
      
    case 'bilinear'
      if size( M.tri , 2 ) ~= 4, error('Connectivity must have 4 columns for the BILINEAR case.'); end

      MAT = [1 0 0 0;-1 1 0 0;-1 0 0 1;1 -1 1 -1];
      
      f = 1; MP =   mvpoly( MAT * M.xyz( M.tri(f,:) , : ) , [1 1] , 2 );
      for f = 2:size(M.tri,1)
        MP = [ MP ; mvpoly( MAT * M.xyz( M.tri(f,:) , : ) , [1 1] , 2 ) ];
      end

      if nargout > 1
        Nxyz = size( M.xyz , 1 );
        L = [];
        for f = 1:size(M.tri,1)
          L = [ L ; kron(  eye(NSD) , MAT ) * kron( eye(NSD) , sparse( 1:4 , M.tri(f,:) , 1 , 4 , Nxyz ) ) ];
        end
      end


    case 'bihermite'
      if size( M.tri , 2 ) ~= 4, error('Connectivity must have 4 columns for the BILINEAR case.'); end
      if ~isfield( M , 'xyzd1'  ), error('Derivatives with respect to 1 must be specified in M.xyzd1'); end
      if ~isfield( M , 'xyzd2'  ), error('Derivatives with respect to 2 must be specified in M.xyzd2'); end
      if ~isfield( M , 'xyzd12' ), error('Cross derivatives must be specified in M.xyzd12'); end

      H2P = [1,0,-0,-0;-0,0,1,-0;-3,3,-2,-1;2,-2,1,1];
      
      MAT = kron(H2P,H2P);
      %MAT = MAT*sparse(1:16,[1,2,5,6,4,3,8,7,9,10,13,14,12,11,16,15],1);
      MAT = MAT(:,[1,2,6,5,3,4,8,7,9,10,14,13,11,12,16,15]);
%       MAT = MAT(:,[1,5,6,2,3,7,8,4,9,13,14,10,11,15,16,12]);
      
      f = 1; MP =   mvpoly( MAT * [ M.xyz( M.tri(f,:) , : ) ; M.xyzd1( M.tri(f,:) , : ) ; M.xyzd2( M.tri(f,:) , : ) ; M.xyzd12( M.tri(f,:) , : ) ] , [3 3] , 2 );
      for f = 2:size(M.tri,1)
        MP = [ MP ; mvpoly( MAT * [ M.xyz( M.tri(f,:) , : ) ; M.xyzd1( M.tri(f,:) , : ) ; M.xyzd2( M.tri(f,:) , : ) ; M.xyzd12( M.tri(f,:) , : ) ] , [3 3] , 2 ) ];
      end

      if nargout > 1
        Nxyz = size( M.xyz , 1 );
        L = sparse([]);
        for f = 1:size(M.tri,1)
          tL = [];
          tL = [ tL ; sparse( 1:4 , M.tri(f,:)          , 1 , 4 , 4*Nxyz ) ];
          tL = [ tL ; sparse( 1:4 , M.tri(f,:) +   Nxyz , 1 , 4 , 4*Nxyz ) ];
          tL = [ tL ; sparse( 1:4 , M.tri(f,:) + 2*Nxyz , 1 , 4 , 4*Nxyz ) ];
          tL = [ tL ; sparse( 1:4 , M.tri(f,:) + 3*Nxyz , 1 , 4 , 4*Nxyz ) ];

          L = [ L ; kron( speye(NSD) , MAT ) * kron( speye(NSD) , tL ) ];
%           L = [ L ; kron(  eye(NSD) , MAT ) * kron( eye(NSD) , sparse( 1:16 , bsxfun(@plus , M.tri(f,:).' , (0:3)*Nxyz )  , 1 , 16 , 4*Nxyz ) ) ];
        end
      end


  end




end
