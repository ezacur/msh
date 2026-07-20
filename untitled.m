cd c:\repos\msh\
addpath( 'c:\repos\msh\MESH\');
addpath( 'c:\repos\msh\tools\');
addpath( 'c:\repos\msh\uiTools\');

%%

H = load('H','H').H;
hP = plotMESH( H ,'ne','nf');
silhouette( hP ,'EdgeColor','r','LineWidth',3)
