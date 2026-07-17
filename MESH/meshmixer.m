function N = meshmixer( M , GROUPfield )

  MM_executable = '';
  if ~isfile( MM_executable )
    MM_executable = fullfile( fileparts( mfilename('fullpath') ) , 'meshmixer' , 'meshmixer.exe' );
  end
  if ~isfile( MM_executable )
    MM_executable = findfilename( fullfile( 'C:\' , 'Win' , 'Meshmixer' ) , { 'meshmixer.exe' } , false );
  end
  if ~isfile( MM_executable )
    try,
      MM_executable = getoption( 'MESHMIXER' , 'executable' );
      if isempty( MM_executable )
        setoption( 'MESHMIXER' , 'executable' , '"C:\Win\Meshmixer\meshmixer.exe"' );
      end
    end
    if isempty( MM_executable ) || ~isfile( MM_executable )
      try, setoption(); end
      error('Check the [MESHMIXER] executable option' );
    end
  end
  if isempty( MM_executable ) || ~isfile( MM_executable )
    error('Set the MESHMIXER_executable variable.' );
  end

  fn_in = inputname(1);
  if isempty( fn_in ), fn_in = 'mesh'; end

  [dn,CLEANOUT] = sanddir( 'MESHMIXER_????' );
  fn_in  = fullfile( dn , [ fn_in , '.obj' ] );
  fn_out = fullfile( dn , [ 'out' , '.obj' ] );

  try, M.triGROUP = M.(GROUPfield); end
  write_MMOBJ( M , fn_in );


  APPDATA = getenv('APPDATA');
  %APPDATA = fileparts( MM_executable );

  MM_ini_file = fullfile( APPDATA , 'Autodesk' , 'meshmixer.ini' );
  if isfile( MM_ini_file )  
    MM_opts = readFile( MM_ini_file );
    REWRITE_ini = onCleanup( @()fwrite( +MM_ini_file , strjoin( MM_opts , '\n' ) ) );

    MM_opts{end+1} = '[MainWindow]';
    MM_opts{end+1} = [ 'SaveFilePath=' , strrep( fn_out ,'\','/') ];

    MM_opts{end+1} = '[Options]';
    MM_opts{end+1} = 'ExportFileFormat=OBJ Format (*.obj)';
    MM_opts{end+1} = 'OBJGroupTag=g';
    MM_opts{end+1} = 'ReadOBJGroupsAsFaceGroups=true';
    MM_opts{end+1} = 'WriteFaceGroupsAsOBJGroups=true';
    MM_opts{end+1} = 'EnableOrthoCamera=true';
    MM_opts{end+1} = 'DefaultUpDirection=1';
    MM_opts{end+1} = 'DefaultNavigationMode=1';
    MM_opts{end+1} = 'EnableTransparentTarget=true';
    MM_opts{end+1} = 'WireframeLineWidth=1.000000';
    MM_opts{end+1} = 'ShowGroundPlaneGrid=false';
    MM_opts{end+1} = 'ShowPrintbed_MMView=false';
    MM_opts{end+1} = 'EnableUniformScale=false';
    MM_opts{end+1} = 'SaveToolSettings=false';
    MM_opts{end+1} = '';

    fwrite( +MM_ini_file , strjoin( MM_opts , '\n' ) );
  end

  [~,allPIDs] = system( sprintf('tasklist /fo list /fi "IMAGENAME eq %s"', filename( MM_executable ) ) );
  allPIDs = regexp( allPIDs , 'PID:\s*(\d*)' , 'tokens' );
  allPIDs = cellfun( @(t)str2double(t{1}) , allPIDs );
  allPIDs0 = allPIDs; 

  [a,b] = system( [ sprintf( 'cd /d "%s"  &' , dn ) ,...
                     '"' , MM_executable , '"' ,' ' ,...
                    '"' , fn_in , '"  &' ] );
  
  [~,allPIDs] = system( sprintf('tasklist /fo list /fi "IMAGENAME eq %s"', filename( MM_executable ) ) );
  allPIDs = regexp( allPIDs , 'PID:\s*(\d*)' , 'tokens' );
  allPIDs = cellfun( @(t)str2double(t{1}) , allPIDs );
  PID = setdiff( allPIDs , allPIDs0 );
  KILL = onCleanup( @()system( sprintf( 'taskkill /f /t /pid %d' , PID ) ) );
  
  %pause(1);

  wait4file( fn_out , Inf );

  clearvars( 'KILL' );
  clearvars( 'REWRITE_ini' );
  
  N = read_MMOBJ( fn_out );
  delete( fn_out );
  N = stickToPoints( M , N , 1e-4 );
  try
    N.(GROUPfield) = N.triGROUP;
    N = rmfields( N , 'triGROUP' );
  end

end
