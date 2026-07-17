function R = g4BestRemesh( I , varargin )

  fcn = [];
  if isa( varargin{end} ,'function_handle')
    fcn = varargin{end};
    varargin(end) = [];
  end

  I = MeshTidy( I ,0,true);
  I = MeshFixCellOrientation( I );
  [~,Q] = g4Remesh( I , varargin{:} );
  
  Q = Q( [ Q{:,2} ] == Q{end,2} ,:);
  for q = 1:size(Q,1), %disp(q);
    M = Q{q,1};
    [rr,er,Ma,ma] = meshQuality( M ,'rr','edgeratio','maxangle','minangle');
    Q{q,3} =  prctile(  rr , 99);
    Q{q,4} =  prctile(  er , 99);
    Q{q,5} =  prctile(  Ma , 99);
    Q{q,6} = -prctile( -ma , 99);
    Q{q,7} = sum( Ma > 170 );
    Q{q,8} = sum( ma <   7 );
  end

  K = cell2mat( Q(:,[3,4,5,6,7,8]) );
  if 0
   plot( log10( K(:,1) ) , '.-r');
  hplot( log10( K(:,2) ) , '.-b');
  hplot(      ( K(:,3) ) , '.-m');
  hplot(      ( K(:,4) ) , '.-g');
  try,hplot(      ( K(:,5) ) , '.-y');end
  try,hplot(      ( K(:,6) ) , '.-c');end
  %%
  end
  [best,info] = consensus_select( K , [1,1,1,-1,1,1] ,'plot',false,'Uncertainty',[0.1,0.1,4,4,5,5],'AnalysisMode','soft','Weights',[1,1,1,1,2,2]);
  if ~info.is_pareto_optimal, best = info.dominated_by(1); end
  R = Q{best,1};

  if ~isempty( fcn ), R = fcn(R); end
  for it = 1:1
    R = MeshImproveObtuses( R , 160 , 5000 );
    R = MeshImproveAcutes(  R ,  10 , 50 );
  end
  if ~isempty( fcn ), R = fcn(R); end

  R = MeshTidy( R ,0,true);
  R = MeshFixCellOrientation( R );

%   R = g4Remesh( R , Q{end,2}*1.2 , 50 , varargin{3:end} );
%   R = g4Remesh( R , Q{end,2}*1.1 , 50 , varargin{3:end} );

  [~,Q] = g4Remesh( R , Q{end,2} , varargin{2:end} );
  for q = 1:size(Q,1)
    M = Q{q,1};
    [rr,er,Ma,ma] = meshQuality( M ,'rr','edgeratio','maxangle','minangle');
    Q{q,3} =  max(  rr );
    Q{q,4} =  max(  er );
    Q{q,5} =  max(  Ma );
    Q{q,6} = -max( -ma );
  end

  K = cell2mat( Q(:,[3,4,5,6]) );
  [best,info] = consensus_select( K , [1,1,1,-1] ,'plot',false,'Uncertainty',[0.1,0.1,4,4],'AnalysisMode','soft','Weight',[1,1,2,2]);
  if ~info.is_pareto_optimal, best = info.dominated_by(1); end
  R = Q{best,1};
  if ~isempty( fcn ), R = fcn(R); end

  R = MeshFixCellOrientation( R );
  

end