function A = MeshRegistrationSimple( A , B , varargin )
%
% deforms A towards B
%
%
%


  B = Mesh(B,0);
  A = Mesh(A);

  LAMBDAS     = geospace( 1e-6 , 1e10 , 15 );
  PERCENTILES = 0.1;

  try, [varargin,~,LAMBDAS]     = parseargs(varargin,'Lambdas'            ,'$DEFS$',LAMBDAS     ); end
  try, [varargin,~,PERCENTILES] = parseargs(varargin,'PERCENTILES','PRCT' ,'$DEFS$',PERCENTILES ); end
  
  vprintf = @(varargin)fprintf( varargin{:} );
  
  LAMBDAS     = LAMBDAS(:);
  PERCENTILES = PERCENTILES(:);

  for it = 1:numel(LAMBDAS)
    %this step parameters
    tLAMBDA       = LAMBDAS( it );
    tPRCT         = PERCENTILES( min(it,end) );

    vprintf( '%4d - LAMBDA: %g , PRCT: %g\n' , it , tLAMBDA , tPRCT );
    
    onA = A.xyz;
    [~,onB] = vtkClosestElement( B , onA );

    d = fro( onA - onB ,2);
    vprintf('         connectors: mean( %g )  max( %g )\n', mean(d) , max(d) );
    TARGET = onA + ( onB - onA ) * tPRCT;

    if isinf( tLAMBDA )
      vprintf('         sticking ...');
      A.xyz  = TARGET;
    else
      vprintf('         deforming ...');
      A.xyz  = InterpolatingSplines( onA , TARGET , onA , 'r' , 'LAMBDA' , tLAMBDA , varargin{:} );
    end
    vprintf(' done\n');
  end

end
