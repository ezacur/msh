function [M,E] = MeshSymmetrize( M00 )

if 0

  y = Optimize1D( @(y)MeshSymmetrize__o2( transform( M ,'ry',y) ) , [-20 20] )
  M = MeshSymmetrize( transform( M ,'ry',y) );
  
  %%
  
end

  M0 = M00;

  M0 = MeshTidy( Mesh( M0 ,0) ,0,1);
  M  = M0;

  vtkClosestElement([],[]); CLEANUP = onCleanup( @()vtkClosestElement([],[]) );
  
  z = Optimize1D( @(z)ENER(z) , [ -90 , 90 ] ,...
    'METHODS','exhaustive','EXHAUSTIVE.MAX_ITS',12,'EXHAUSTIVE.N',7,'EXHAUSTIVE.ALLOWLOG',false , 'MIN' , -95 , 'MAX', 95 );

  
  [E,t] = ENER( z );
  R = [ R , [ t ; 0 ; 0 ] ; 0 , 0 , 0 , 1 ];
  M = transform( M00 , R );
  
  function [E,t] = ENER( z )
    R = [  cosd(z) , -sind(z) , 0 ;
           sind(z) ,  cosd(z) , 0 ;
                0  ,       0  , 1 ];

    M.xyz = M0.xyz * R.';
    vtkClosestElement([],[]);
    vtkClosestElement( M );

    xyz = M.xyz;
    a = - max( xyz(:,1) );
    b = - min( xyz(:,1) );
    
    xyz(:,1) = - xyz(:,1);
    
    [t,E] = Optimize1D( @(t)distance_t( xyz , t ) , [ a , b ] ,...
      'METHODS','golden','GOLDEN.MAX_ITS',12 , 'MIN' , a , 'MAX' , b );
      %'METHODS','exhaustive','EXHAUSTIVE.MAX_ITS',7,'EXHAUSTIVE.N',5,'EXHAUSTIVE.ALLOWLOG',false , 'MIN' , -b , 'MAX', -a );
      %'METHODS','golden','MAX_ITS',5 );
      %'METHODS','exhaustive','MAX_ITS',5,'N',5,'ALLOWLOG',false);
%     disp( [ E , t ] );
  end
  function E = distance_t( x , t )
    x(:,1) = x(:,1) - 2*t;
    [~,~,d] = vtkClosestElement( x );
    E = sum( d.^2 );
  end

end
