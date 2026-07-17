function H = meshHistogram( M , V )
if 0


M = Mesh( ndmat(linspace(0,2,25),linspace(-1,1.5,35)) , 'delaunay');

V = fro( M.xyz ,2);

hist( V , 100 );
hplot( meshHistogram(M,V) , '-2r')


%%
end

  u = double( unique( single( V ) ) );
  m = ( u(1:end-1) + u(2:end) )/2;

  u = meshPrctile( M , V , u );  s = u(end);
  m = meshPrctile( M , V , m );

  x0 = u( 1:end-1 ,1);   y0 = u( 1:end-1 ,2)/s;
  xm = m(  :      ,1);   ym = m(  :      ,2)/s;
  x1 = u( 2:end   ,1);   y1 = u( 2:end   ,2)/s;
  
  d  = (x0 - x1) .* (x0 - xm) .* (x1 - xm);
  
  a  = ( xm    - x1    ).*y0 + ( x1    - x0    ).*ym + ( x0    - xm    ).*y1;   a = -a./d;
  b  = ( xm.^2 - x1.^2 ).*y0 + ( x1.^2 - x0.^2 ).*ym + ( x0.^2 - xm.^2 ).*y1;   b =  b./d;

  H = [ x0 , 2*a.*x0 + b ];
  H(end+1,:) = [ x1(end) , 0];
  H(1,2) = 0;

  w = find( ~isfinite( H(:,2) ) );
  H(w,2) = H( w-1 ,2);


end