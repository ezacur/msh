function F = meshDilateField( M , F )

  if isstruct( M ) && isfield( M ,'tri') && ischar( F ) && isfield( M , F )

    M.(F) = meshDilateField( M , M.(F) );
    return;

  end


  if isnumeric( F ) && all( ismember(F,[0,1]) ), F = ~~F; end
  if ~islogical( F ), error( 'The field must be a logical.'); end
  if ~isvector( F ), error( 'The field must be a vector.' ); end
  if size( F ,2) ~= 1, error( 'The field must be a column vector.' ); end
  if numel(F) ~= size( M.tri ,1), error( 'The field must have size(M.tri,1) elements.' ); end


  w = ~~F;

  E = meshEsuE( M );

  F( any( E(w,:) ,1) ) = 1;

end
