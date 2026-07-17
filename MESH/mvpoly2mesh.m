function M = mvpoly2mesh( MP , varargin )

  nD = numel( MP.g );
  nC = size( MP.C , numel(MP.g)+1 );
  plotTYPE = sprintf( '%d-%d' , nD , nC );
  sz = size( MP );
  dots = repmat({':'},1,numel( MP.g ));
  
  RES = ceil( realpow( prod( MP.g + 1 ) , 1/nD ) );
  [varargin,i,RES] = parseargs(varargin,'RESolution','$DEFS$',RES);
  RES(end+1:nD) = RES(end);
  
  M = struct('xyz',[],'tri',[]);
  switch plotTYPE
    case {'1-1'},
    case {'1-2','1-3'},
      M = struct('xyz',[],'tri',[]);

      faces = [ (1:RES-1).' , (2:RES).' ];

      c = linspace(0,1,RES(1)).';
      
      [v,fcn] = evaluate( MP([]) , c );
      for i = 1:prod(sz)
        M.tri = [ M.tri ; faces + size(M.xyz,1) ];
        M.xyz = [ M.xyz ; fcn( MP.C( dots{:} , : , i ) ) ];
      end

      [varargin,D] = parseargs(varargin,'d','$FORCE$',{true,false});
      if D
        M.xyzd = [];
        D = derivative( MP , 1 );
        [v,fcn] = evaluate( D([]) , c );
        for i = 1:prod(sz)
          M.xyzd = [ M.xyzd ; fcn( D.C( dots{:} , : , i ) ) ];
        end
      end
      
      if ~isempty( varargin ), warning('there are unused options in varargin'); end
      
    case {'2-1'},
    case {'2-2'},
    case {'2-3'},
      M = struct('xyz',[],'tri',[]);

      faces = reshape( bsxfun(@plus, ...
                              vec([         1 :    RES(1)-1    ;...
                                            2 :    RES(1)      ;...
                                     RES(1)+2 :  2*RES(1)      ;...
                                     RES(1)+1 :( 2*RES(1)-1 )  ]) ,...
                              (0:RES(2)-2)*RES(1) ) ,...
                       4 , [] ).';

      [varargin,TRIANGULATE] = parseargs(varargin,'TRIangulate','$FORCE$',{true,false});
      [varargin,symTRIANGULATE] = parseargs(varargin,'SYMmetrictriangulation','$FORCE$',{true,false});
      if symTRIANGULATE
        cNODE = ( 1:(RES(1)-1)*(RES(2)-1) )' + RES(1)*RES(2);

        faces = cat( 3 , [ faces(:,[1 2] ) cNODE ] ,...
                         [ faces(:,[2 3] ) cNODE ] ,...
                         [ faces(:,[3 4] ) cNODE ] ,...
                         [ faces(:,[4 1] ) cNODE ] );
        faces = permute(faces,[2 3 1]);
        faces = reshape(faces ,3,[] ).';
      elseif TRIANGULATE
        faces = cat(3,faces(:,[1 2 3]),faces(:,[3 4 1]));
        faces = permute(faces,[2 3 1]);
        faces = reshape(faces ,3,[] ).';
      end
     
      X = linspace(0,1,RES(1));
      Y = linspace(0,1,RES(2));
      c = ndmat( X , Y );
      if symTRIANGULATE
        cp = @(x) ( x(1:end-1) + x(2:end) )/2;
        c = [ c ; ndmat( cp(X) , cp(Y) ) ];
      end
      
      [v,fcn] = evaluate( MP([]) , c );
      for i = 1:prod(sz)
        M.tri = [ M.tri ; faces + size(M.xyz,1) ];
        M.xyz = [ M.xyz ; fcn( MP.C( dots{:} , : , i ) ) ];
      end

      [varargin,D] = parseargs(varargin,'d1','$FORCE$',{true,false});
      if D
        M.xyzd1 = [];
        D = derivative( MP , 1 );
        [v,fcn] = evaluate( D([]) , c );
        for i = 1:prod(sz)
          M.xyzd1 = [ M.xyzd1 ; fcn( D.C( dots{:} , : , i ) ) ];
        end
      end
      
      [varargin,D] = parseargs(varargin,'d2','$FORCE$',{true,false});
      if D
        M.xyzd2 = [];
        D = derivative( MP , 2 );
        [v,fcn] = evaluate( D([]) , c );
        for i = 1:prod(sz)
          M.xyzd2 = [ M.xyzd2 ; fcn( D.C( dots{:} , : , i ) ) ];
        end
      end
      
      [varargin,D] = parseargs(varargin,'d12','$FORCE$',{true,false});
      if D
        M.xyzd12 = [];
        D = derivative( MP , [ 1 , 2 ] );
        [v,fcn] = evaluate( D([]) , c );
        for i = 1:prod(sz)
          M.xyzd12 = [ M.xyzd12 ; fcn( D.C( dots{:} , : , i ) ) ];
        end
      end
      
      if ~isempty( varargin ), warning('there are unused options in varargin'); end
    case {'3-1'},
    case {'3-2'},
    case {'3-3'},
  end

  
  
  
end
