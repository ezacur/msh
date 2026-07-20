function [ varargout ] = parseargs( in , varargin )
%PARSEARGS  Pull a keyword option (and its values) out of an argument list.
%
%   Scans the cell IN (typically a function's VARARGIN) for any of the given
%   KEYs and reports where it was found, the value(s) that follow it, and IN
%   with that key+values removed. Matching is CASE-INSENSITIVE; a CamelCase key
%   also matches its capitals-only abbreviation (e.g. 'FaceColor' matches 'FC').
%
%   pos = parseargs( IN , KEY1 , KEY2 , ... )
%       POS = 1-based position of the LAST matching key in IN, or 0 if none of
%       the keys is present. Several KEYs act as aliases for one option.
%
%   [out, pos, v1, ..., vN] = parseargs( IN , KEYS... )
%       Requests N value outputs (N = nargout-2). v1..vN are the N entries that
%       follow the key in IN; OUT is IN with the key and those N entries removed
%       ("consumed"); POS is the position (0 if absent). If the key appears more
%       than once, the LAST occurrence wins and every occurrence is removed.
%
%   ... = parseargs( ... , '$DEFS$' , d1 , ... , dN )
%       Default values for v1..vN when the key is absent (exactly N required).
%
%   ... = parseargs( ... , '$FORCE$' , {A,B} )
%       Override POS: return A when the key is present, B when absent (instead
%       of the numeric position / 0).
%   ... = parseargs( ... , '$FORCE$' )
%       Bare form: POS becomes TRUE when present and FALSE when absent (a handy
%       "is this flag set?" test). When both are used '$FORCE$' precedes '$DEFS$'.
%
%   Keys must be non-empty char row vectors. parseargs({}) is valid. Values are
%   read verbatim, so a key may take a value of any type.
%
%   Examples (inside a function that has VARARGIN):
%       [varargin,pos]  = parseargs( varargin , 'verbose' );
%       [varargin,~,n]  = parseargs( varargin , 'N' , 'num' , '$DEFS$' , 10 );
%       [varargin,tf]   = parseargs( varargin , 'Debug' , '$FORCE$' );  % tf=true/false
%
%   See also inputParser, nargin.

if 0
  %     1  2      3  4  5  6      7 8    9       10   11    12 13
  in = {0,'first',1,'s',2,'third',3,3.1,'second',2.2,'yes','N',Inf};
  [var,pos,val] = parseargs(in,'First','$DEFS$',-1);
  [var,pos,val] = parseargs(in,'Second','sec','$DEFS$',-2);
  [var,pos,val1,val2] = parseargs(in,'third','$DEFS$',-3,-3.1);
  [var,pos] = parseargs(in,'Yes','$FORCE$',{1,2})
  [var,pos] = parseargs(in,'No','$FORCE$',{1,2})
  [var,pos] = parseargs(in,'No')
  [var,pos] = parseargs(in,'maybe')
  [var,pos] = parseargs(in,'maybe','$FORCE$',{1,2})
end
if 0
in = {'-sergio',3,'hola',1,'ho',4,'chau',3,'mig',1,'ernesto',2,'salvador',7,8,'sergio','n',100,'n:200','miguel',6, '-sergio', 79, 'carve', 'bool'};
k1 = rand(3);
%%
parseargs({},'a')       %mal, no debe dar error, corregido
i = parseargs(in,'david','$FORCE$',{1500;-1500})  %que no error, corregido

[out i] = parseargs({}, 'a', '$FORCE$', {'jose', 4}, '$DEFS$', pi)   %cambiar el mensaje
[out i v,b] = parseargs({}, 'a', '$FORCE$', {'jose', 4}, '$DEFS$', pi)  %cambiar el mensaje
parseargs(in,'sergio','$FORCE$',1500)             %no dar error, corregido
i = parseargs(in,'sergio','$FORCE$',uint8(15))   %warning
i = parseargs(in,'sergio','$FORCE$',{uint8(15),uint8(3)})   %warning
i = parseargs(in,'sergio','$FORCE$',{uint8(15)})   %warning
i = parseargs(in,'sergio','$FORCE$',{-1,uint8(15)})   %warning
i = parseargs(in,'-sergio','$FORCE$',{1,uint8(0)})   %warning
[out i v] = parseargs(in,'sergio','$FORCE$',{1,uint8(0)})   %warning


i = parseargs({'-sergio'},'SERGIO')
i = parseargs({'so'},'SergiO')



parseargs({})
i = parseargs({})
[out,i,v] = parseargs({})                                               
[var i v] = parseargs({}, '$FORCE$', {'jose', 4}, '$DEFS$', pi)

 i = parseargs({},'a')                       
[out,i] = parseargs({}, 'a')
[out,i,v] = parseargs({}, 'a')
[out i v] = parseargs({}, 'a', '$FORCE$', {'jose', 4}, '$DEFS$', pi)
[out,i,v,b] = parseargs({}, 'a')
[out i v,b] = parseargs({}, 'a', '$FORCE$', {'jose', 4}, '$DEFS$', pi, k1)

i = parseargs(in,'hola')
[out,i] = parseargs(in,'hola')
[out,i] = parseargs(in,'hola','ho')       %% aca da mal, i deberia ser 3!!!
[out,i,v] = parseargs(in,'HOla','hol')  %% aca da mal, i deberia ser 3!!!
[out,i,v,b] = parseargs(in,'HOla','hol')  %% aca da mal, i deberia ser 3!!!

[out,i,v] = parseargs(in,'m','mi','migu','MIGuel')

[out,i,v1,v2] = parseargs(in,'SALvador')

[out,i,v1,v2] = parseargs(in,'m','mi','migu','MIGuel')  %%error... despues de miguel se esperan 2 valores

[out,i,v1] = parseargs(in,'matias')
[out,i,v1] = parseargs(in,'matias','$DEFS$',300)
[out,i,v1,v2] = parseargs(in,'matias','$DEFS$',300)   %%error, como se piden 2 salidas, deberia haber 2 $DEFS$
[out,i,v1,v2] = parseargs(in,'matias','$DEFS$',300,rand(8))





i = parseargs(in,'sergio','$FORCE$')
i = parseargs(in,'sergio','$FORCE$','carlitos')
i = parseargs(in,'david','$FORCE$','carlitos')
i = parseargs(in,'david','$FORCE$',1500)

[out,i,v] = parseargs(in,'david','$FORCE$',{1500,-1500},'$DEFS$',333)
[out,i,v] = parseargs(in,'MIGuel','$DEFS$',333,'$FORCE$',{'estaMIGUEL','NONONO_estaMIGUEL'},'$DEFS$',333)
[out,i,v] = parseargs(in,'MATias','$FORCE$',{'estaMATIAS','NONONO_estaMATIAS'},'$DEFS$',666)

[out,i,v] = parseargs(in,'n','$DEFS$',pi)   %devuelve 100, pero estaria perfecto si devuelve 200!!


end


  % A compiled parseargs_mex existed ~10 years ago; it is gone and will not
  % come back. The old code TRIED it on every call -- but invoking a missing
  % function throws + catches an exception (~75 us), which was ~74% of the whole
  % parseargs cost. Removed. If a MEX is ever reintroduced, DON'T pay the
  % exception per call: cache the existence check and guard the call, e.g.
  %   persistent USE_MEX
  %   if isempty( USE_MEX ), USE_MEX = exist('parseargs_mex','file') == 3; end
  %   if USE_MEX
  %     try, [ varargout{1:nargout} ] = parseargs_mex( in , varargin{:} ); return; end
  %   end

  if ~iscell( in ), error( 'IN should be a cell'); end
  
  Nouts           = max( nargout - 2 , 0 );
  POSITION_output = [];
  POSITION_set    = false;
  DEFAULTS_output = cell(1,Nouts);

  OPTSid = find( strcmp( varargin , '$DEFS$' ) | strcmp( varargin , '$FORCE$' ) ,1);
  if isempty( OPTSid ), OPTSid = numel( varargin ) + 1; end
  OPTS = varargin( OPTSid:end );

  o = 1;
  while o <= numel( OPTS )
    if strcmp( OPTS{o} , '$FORCE$' )
      POSITION_set = true;
      if numel( OPTS ) < o + 1
          POSITION_output = true;
      else
          POSITION_output = OPTS{o+1};
      end
          
      o = o + 1 + 1;
      continue;
    end
    if strcmp( OPTS{o} , '$DEFS$' )
      if numel( OPTS ) ~= o + Nouts
        error('specification of %d defaults is expected',Nouts);
      end
      DEFAULTS_output = OPTS(o+(1:Nouts));
      o = o + Nouts + 1;
      continue;
    end
    error('it is not expected to be here!!');
  end
  if POSITION_set
    if ~iscell( POSITION_output )
      POSITION_output = { POSITION_output };
    end
    if numel( POSITION_output ) == 1
      POSITION_output{2} = false;
    end
    if numel( POSITION_output ) > 2
      error('invalid $FORCE$ specification');
    end
  end
  
  
  KEYS = cell( 1 , 2*(OPTSid-1) );
  lk = 0;
  for k = 1:( OPTSid-1 )
    K = varargin{k};
    if ~ischarrow(K)
      error('only strings are allowed as keys');
    end
    lk = lk + 1; KEYS{ lk } = lower( K );
    
    K( K >= 'a' & K <= 'z' ) = [];
    if isempty( K ), continue; end
    lk = lk + 1; KEYS{ lk } = lower( K );
  end
  KEYS( (lk+1):end ) = [];
  

  toDELETE = [];
  POSITION = 0;
  i = 1;
  while i <= numel( in )
    if ~ischarrow( in{i} ), i = i+1; continue; end
    if ~any( strcmpi(  in{i}  , KEYS ) ), i = i+1; continue; end
  
    if numel( in ) < i+Nouts
      error('after key ''%s'', %d inputs are expected.',in{i},Nouts);
    end
    
    toDELETE = [ toDELETE , i:(i+Nouts) ];
    POSITION = i;
    i = i + Nouts + 1;
  end
  
  if POSITION == 0 && nargout < 2

    if POSITION_set
      varargout{1} = POSITION_output{2};
    else
      varargout{1} = 0;
    end
  
  elseif POSITION ~= 0 && nargout < 2
    
    if POSITION_set
      varargout{1} = POSITION_output{1};
    else
      varargout{1} = POSITION;
    end
    
  elseif POSITION == 0

    for i = 1:Nouts
      varargout{2+i} = DEFAULTS_output{ i }; %#ok<AGROW> 
    end
    if POSITION_set
      varargout{2} = POSITION_output{2};
    else
      varargout{2} = 0;
    end
    varargout{1} = in;

  else
    
    for i = 1:Nouts
      varargout{2+i} = in{ POSITION + i }; %#ok<AGROW> 
    end
    if POSITION_set
      varargout{2} = POSITION_output{1};
    else
      varargout{2} = POSITION;
    end
    in( toDELETE ) = [];
    varargout{1} = in;

  end
  
end
function s = ischarrow( x )
  s = ischar( x ) && ~isempty( x ) && ndims( x ) <= 2 && size( x , 1 ) == 1;
end
