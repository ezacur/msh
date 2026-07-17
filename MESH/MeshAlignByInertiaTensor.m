function MA = MeshAlignByInertiaTensor( M0 , CENTER )
%MESHALIGNBYINERTIATENSOR  Reorient a mesh onto its principal axes of inertia.
%
%   MA = MeshAlignByInertiaTensor( M )
%   MA = MeshAlignByInertiaTensor( M , CENTER )
%   MA = MeshAlignByInertiaTensor( M , T )
%
%   Rotates the CLOSED triangular mesh M so its principal axes of inertia
%   (eigenvectors of the volume inertia tensor) align with the coordinate
%   axes. The inertia tensor is invariant under the 24 proper rotations of the
%   octahedral group (axis permutations and pairs of sign flips), so that
%   discrete ambiguity is resolved explicitly (see below).
%
%   CENTER (logical, default false):
%     false : keep the mesh at its original centroid (reorient in place).
%     true  : also translate the centroid to the origin.
%
%   NOT rotation-invariant (by design): among the 24 equivalent frames it
%   picks the one CLOSEST to M's current orientation (least reorientation), so
%   two rotations of the same object give two different results. It normalizes
%   the axes with minimal disruption; it does NOT produce a single canonical
%   pose.
%
%   TEMPLATE / registration mode:  MA = MeshAlignByInertiaTensor( M , T )  with
%   T a mesh struct brings M into T's frame (aligns M onto T through their
%   common inertial frame). Here the octahedral freedom is resolved
%   CONSISTENTLY: the 24 candidates are tried and the one best overlaying T is
%   kept, so it works for ARBITRARY relative pose (by exact vertex overlay when
%   M and T share the triangulation, else by a correspondence-free surface
%   distance via distanceFrom).
%
%   LIMITATION: for (near-)axisymmetric meshes two inertia eigenvalues
%   coincide, leaving a CONTINUOUS roll ambiguity about the symmetry axis that
%   these 24 discrete rotations cannot resolve; that roll is left arbitrary.
%
%   Requires a closed mesh (volume inertia via meshVolume). The template mode's
%   correspondence-free fallback additionally needs VTK (distanceFrom).
%
% See also meshVolume, MatchPoints, transform, logmrot.

  if nargin < 2, CENTER = false; end

  if isstruct( CENTER )

    T0 = CENTER;
    TA = MeshAlignByInertiaTensor( T0 ,true); rT = MatchPoints( double(TA.xyz) , double(T0.xyz) ,'Rt' );
    MA = MeshAlignByInertiaTensor( M0 ,true); rM = MatchPoints( double(MA.xyz) , double(M0.xyz) ,'Rt' );

    MA = transform( M0 , rM , minv( rT ) );

    % The octahedral sign/permutation freedom was resolved INDEPENDENTLY for M0
    % and T0 (each closest to its OWN orientation), so when they differ by more
    % than ~45 deg the composition lands in a different octahedral cell than T0
    % (a 90/180 deg residual). Snap MA into T0's cell: try the 24 rotations
    % about T0's centroid and keep the one that best overlays T0.
    [~,cT] = meshVolume( T0 ,'center' );
    RS = octahedralRotations();
    sameN = size( M0.xyz ,1) == size( T0.xyz ,1);
    bestD = Inf; best = MA;
    for r = 1:size(RS,3)
      MRr = transform( MA , 't',-cT , blkdiag( RS(:,:,r) ,1) , 't',cT );
      if sameN
        d = fro2( MRr.xyz - T0.xyz );                        %exact overlay (shared triangulation)
      else
        d = sum( distanceFrom( double(MRr.xyz) , T0 ).^2 );  %correspondence-free (needs VTK)
      end
      if d < bestD, bestD = d; best = MRr; end
    end
    MA = best;
    return;

  end



  M = M0;
  c0 = [0,0,0];
  for i = 1:4
    [~,c,it] = meshVolume( M ,'center');
    if i == 1, c0 = c; end
    M = transform( M , 't' , -c );
    [R,~] = eig(it); %disp( round(R*10)/10 ); disp(' ')
    M = transform( M , blkdiag( R.' , 1 ) );
    if det(R) < 0
      M = transform( M , diag([-1,1,1,1]) );
    end
  end

  if ~CENTER && any(c0)
    M = transform( M , 't' , c0 );
  end

  RS = octahedralRotations();
  %M = transform( M , blkdiag( RS(:,:,round( rand(1)*24 ) ) ,1) );
  D = Inf;
  for r = 1:size(RS,3)
    MR = transform( M , blkdiag( RS(:,:,r) ,1) );
    R = MatchPoints( double(MR.xyz) , double(M0.xyz) ,'Rt'); R = R(1:3,1:3);
    d = -trace( R );  %rotation size: trace(R)=1+2cos(angle), so max trace <=>
                      %min angle (closest to original). Cheaper than logmrot and
                      %free of its branch cut at +-pi.
    %d = fun( @(c)max( abs( invmaketransform( R , c ) ) ) , {'rxyz';'rxzy';'ryxz';'ryzx';'rzyx';'rzxy';'rxyx';'rxzx';'ryxy';'ryzy';'rzxz';'rzyz'} );
    if d < D
      D  = d;
      MA = MR;
    end
  end

%   MR = {};
%   Rs = ndmat( [ 0 90 -90 180 ] , [ 0 90 -90 180 ] , [ 0 90 -90 180 ] );
%   for r = 1:size(Rs,1)
%     MR{r,2} = Rs(r,:);
%     MR{r,3} = maketransform( 'rxyz', MR{r,2} );
%     MR{r,3} = round( MR{r,3} * 1000 ) / 1000;
%     MR{r,1} = transform( M , MR{r,3} );
%   end
%   try, for p = 1:350, MR( setdiff( find( fun( @(m)isequal(m,MR{p,3}) , MR(:,3) ) ) ,p) ,:) = []; end; end
% 
%   for r = 1:size(MR,1)
%     R = MatchPoints( MR{r,1}.xyz , M0.xyz ,'Rt' );
%     MR{r,4} = fro2( logm( R(1:3,1:3) ) );
%   end
%   MR = MR( order( [ MR{:,4} ] ) ,:);
%

end


function RS = octahedralRotations()
% The 24 proper (det=+1) rotations of the octahedral group = the axis
% permutations with an even number of sign flips. They are exactly the
% symmetries under which an inertia tensor is invariant, i.e. the discrete
% sign/permutation ambiguity of a principal-axis alignment.
  RS = reshape([1,0,0,0,-1,0,0,0,-1,1,0,0,1,0,0,0,0,1,0,-1,0,0,0,-1,0,0,1,0,-1,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,1,0,0,0,1,-1,0,0,0,0,1,0,0,-1,-1,0,0,0,-1,0,-1,0,0,1,0,0,0,1,0;0,1,0,1,0,0,0,1,0,0,0,-1,0,0,1,0,1,0,0,0,-1,1,0,0,1,0,0,0,0,1,-1,0,0,-1,0,0,0,0,1,0,-1,0,0,0,-1,-1,0,0,0,0,-1,0,-1,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,-1,0,1,0,0;0,0,1,0,0,1,1,0,0,0,1,0,0,-1,0,-1,0,0,1,0,0,0,-1,0,0,1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,-1,0,1,0,0,-1,0,0,0,1,0,0,0,-1,0,0,-1,0,0,-1,0,0,-1],[3,3,24]);
end
