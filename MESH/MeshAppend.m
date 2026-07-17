function A = MeshAppend( varargin )
%
% m= AppendMeshes( m1 , m2 )
%

%%correct the celltype

  keepPARTS = false;
  if numel(varargin) && ( islogical( varargin{end} ) || ischar( varargin{end} ) )
    keepPARTS = varargin{end}; varargin(end) = [];
    if islogical( keepPARTS ) && numel( keepPARTS ) == 1
    elseif ischar( keepPARTS )
      switch lower( keepPARTS )
        case {'keepparts','kp','keep'}, keepPARTS = true;
        case {'removeparts','rm','remove'}, keepPARTS = false;
        otherwise, error('invalid specification of what to do with PARTS.');
      end
    else
      error('unknown specification of what to do with PARTS.');
    end
  end
    

  if numel( varargin ) && iscell( varargin{1} )
    A = MeshAppend( varargin{1}{:} , keepPARTS );
    return;
  end

  A = struct('xyz',[],'tri',[],'celltype',[],'xyzPART___',[],'triPART___',[]);
  %A = struct('xyz',[],'tri',[],'celltype',[]);

  for v = 1:numel(varargin)
    
    B = varargin{v};
    if isempty( B ), continue; end
    if ~isstruct( B ), error('only meshes are accepted as input'); end
    
    B.celltype = meshCelltype( B );
    if numel( B.celltype ) == 1
      B.celltype( 1:size(B.tri,1) , 1 ) = B.celltype;
    end
    
    B.xyzPART___( 1:size(B.xyz,1) ,1) = v;
    B.triPART___( 1:size(B.tri,1) ,1) = v;
    
    if isfield( B , 'xyzUV' ) && isfield( B , 'texture' )
      if ~isfield( A , 'xyzUV' )
        A.xyzUV = zeros( size(A.xyz,1) , size(B.xyzUV,2) ) + 0.5;
      end
      if ~isfield( A , 'texture' )
        if isempty( A.xyzUV ), A.texture = uint8(zeros(0,0,3));
        else,                  A.texture = uint8(cat(3,255,0,0));
        end
      end
    end
    if isfield( A , 'xyzUV' ) && isfield( A , 'texture' )
      if ~isfield( B , 'xyzUV' ),   B.xyzUV = zeros( size(B.xyz,1) , size(A.xyzUV,2) ) + 0.5; end
      if ~isfield( B , 'texture' ), B.texture = uint8(cat(3,255,0,0)); end
    end
    
    
    for f = fieldnames( B ).', f = f{1};
      if strcmp( f , 'celltype' )
        A.celltype = [ A.celltype ; B.celltype ];
        continue;
      end
      
      if strcmp( f , 'xyz' ), continue; end
      if strncmp( f , 'xyz' , 3 ) && ~strcmp( f , 'xyzUV' )
        if ~isfield( A , f )
          sz = size( B.(f) ); sz(1) = size( A.xyz , 1 );
          A.(f) = NaN( sz );
        end
        A.(f) = [ A.(f) ; B.(f) ];
        continue;
      end
      
      if strcmp( f , 'tri' ), continue; end
      if strncmp( f , 'tri' , 3 )
        if ~isfield( A , f )
          sz = size( B.(f) ); sz(1) = size( A.tri , 1 );
          if 0
          elseif iscell( B.(f) ),  val = {''};
          elseif isfloat( B.(f) ), val = NaN;
          elseif isnumeric( B.(f) ), val = 0;
          elseif islogical( B.(f) ), val = false;
          end
          A.(f) = repmat( val ,sz);
        end
        A.(f) = [ A.(f) ; B.(f) ];
        continue;
      end

      if strcmp( f , 'xyzUV' ) && isfield( B , 'xyzUV' ) && isfield( A , 'xyzUV' )
        if isfield( A , 'texture' ) && isfield( B , 'texture' ) && ( ~isequal( size(A.texture) , size(B.texture) ) || ~isequal( A.texture , B.texture ) )
          uv2ji = @(uv,T) [ uv(:,1) * size(T,2) + 0.5     , ( 1-uv(:,2) ) * size(T,1) + 0.5 ];
          ji2uv = @(ji,T) [ ( ji(:,1) - 0.5 ) / size(T,2) , 1 - ( ji(:,2) - 0.5 )/size(T,1) ];
          
          Aji = uv2ji( A.xyzUV , A.texture );
          Bji = uv2ji( B.xyzUV , B.texture ); Bji(:,2) = Bji(:,2) + size( A.texture ,1);
  
          A.texture = safecat( {1,'l',zeros([1,1],class(A.texture))} , A.texture , B.texture );
          A.xyzUV = ji2uv( [ Aji ; Bji ] , A.texture );
        else
          A.xyzUV = [ A.xyzUV ; B.xyzUV ];
        end
        continue;
      end
      
    end
    
    A.tri( 1:end , end+1:size(B.tri,2) ) = 0;
    B.tri( 1:end , end+1:size(A.tri,2) ) = 0;
    w = B.tri == 0;
    w = numel( A.tri ) + find(w);
    
    A.tri = [ A.tri ; B.tri + size( A.xyz , 1 ) ];
    A.tri( w ) = 0;
    
    nTRI = size( A.tri ,1);
    
    
    A.xyz = [ A.xyz ; B.xyz ];                     nXYZ = size( A.xyz ,1);
    

    
    for f = fieldnames( A ).', f = f{1};
      if strcmp( f , 'tri' ), continue; end
      if strcmp( f , 'xyz' ), continue; end
      if strncmp( f , 'xyz' , 3 ) && size( A.(f) ,1 ) < nXYZ
        if iscell( A.(f) )
          [ A.(f){ end+1:nXYZ ,:,:,:,:} ] = deal(NaN);
        else
          A.(f)(end+1:nXYZ,:,:,:,:) = NaN;
        end
        continue;
      end
      if strncmp( f , 'tri' , 3 ) && size( A.(f) ,1 ) < nTRI
        if iscell( A.(f) )
%           A.(f){nTRI,1} = [];
          [ A.(f){ end+1:nTRI ,:,:,:,:} ] = deal('');
        else
          A.(f)(end+1:nTRI,:,:,:,:) = NaN;
        end
        continue;
      end
    end        
    
    
  end

  if ~isempty( A.xyz )
    if ~all( A.celltype == A.celltype(1) )
      warning( 'celltypes look different');
    end
    A.celltype = A.celltype(1);
  end
  
%   if isfield( A , 'xyzUV' ) && isfield( A , 'texture' )
%     A.xyzUV = A.xyzUV - 1;
%     A.xyzUV = bsxfun( @rdivide , A.xyzUV , [ size( A.texture ,2) , size( A.texture ,1) ]-1 );
%     A.xyzUV(:,2) = 1 - A.xyzUV(:,2);
%   end



  if keepPARTS
    try, A.xyzPART = A.xyzPART___; end
    try, A.triPART = A.triPART___; end
  end
  try, A = rmfield( A , 'xyzPART___'); end
  try, A = rmfield( A , 'triPART___'); end
  
end
