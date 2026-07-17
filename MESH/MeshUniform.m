function M = MeshUniform( M , el )

  if numel( el ) < 2, el(2) = el(1)/2; end

  M0 = M; %struct('xyz',double(M.xyz),'tri',double(M.tri));

  nV = size( M.xyz ,1);
  for it = 1:100
    nVp = nV;

    EL = [  sum( ( M.xyz( M.tri(:,1) ,:) - M.xyz( M.tri(:,2) ,:) ).^2 ,2) ,...
            sum( ( M.xyz( M.tri(:,1) ,:) - M.xyz( M.tri(:,3) ,:) ).^2 ,2) ,...
            sum( ( M.xyz( M.tri(:,3) ,:) - M.xyz( M.tri(:,2) ,:) ).^2 ,2) ];

    M = MeshSubdivide( M , any( EL > el(2)^2 ,2) );
    nV = size( M.xyz ,1);
    if nV == nVp, break; end
  end

  %[P,IDS] = FarthestPointSampling( P , IDS , minD , maxN , D_fcn , VERBOSE )
  if el(1) > 0
    [~,nodeIDS] = FarthestPointSampling( M.xyz ,[],el(1));
  else
    [~,nodeIDS] = FarthestPointSampling( M.xyz ,[],0,-el(1));
  end

%   M = MeshSubdivide( M );
% 
%   try
% %     error(1);
%     D = ipd( M.xyz , M.xyz( nodeIDS ,:) );
%     [~,B] = min( D , [] ,2);
%   catch
%     B = zeros( size(M.xyz,1) , 1);
%     E = Inf( size(B) );
%     for c = 1:numel( nodeIDS )
%       D = sum( bsxfun( @minus , M.xyz , M.xyz( nodeIDS(c) ,:) ).^2 ,2);
%       w = D < E;
%       B(w) = c;
%       E(w) = D(w);
%     end
%   end

  nV = size( M.xyz ,1);
  [ES,EL] = meshEdges( M );

  G = graph( ES(:,1) , ES(:,2) , EL , nV );
%   W = sparse( double( ES(:,2) ) , double( ES(:,1) ) , EL , nV , nV );
%   G = bioinfo.internal.biograph2matlab( W ,'Directed',false);

  D = distances( G , nodeIDS , 'Method' , 'positive' );

  [~,B] = min( D , [] ,1);

  
  B = nodeIDS( B );

  B = B( M.tri );
  for f = fieldnames( M ).', f = f{1};
    if ~strncmp( f , 'tri' , 3 ), continue; end
    M = rmfield( M , f );
  end

  B( B(:,1) == B(:,2) ,:) = [];
  B( B(:,1) == B(:,3) ,:) = [];
  B( B(:,2) == B(:,3) ,:) = [];
  B = unique( B ,'rows');

  M.tri = B;


  b = zeros( size( M.xyz,1) ,1);
  b( M.tri(:) ) = 1;
  b( ~~b ) = 1:nnz(b);
  M.tri = b( M.tri );
  b = ~~b;
  for f = fieldnames( M ).', f = f{1};
    if ~strncmp( f , 'xyz' , 3 ), continue; end
    M.(f) = M.(f)( b ,:);
  end

  e = vtkClosestElement( struct('xyz',double(M0.xyz),'tri',double(M0.tri))  , double( meshFacesCenter( M ) ) );
  for f = fieldnames( M0 ).', f = f{1};
    if strcmp( f , 'tri' ), continue; end
    if strcmp( f , 'xyz' ), continue; end
    if strncmp( f , 'tri' ,3)
      M.(f) = M0.(f)(e,:,:,:,:,:);
    end
  end


end
