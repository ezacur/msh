cd c:\repos\msh\
addpath( 'c:\repos\msh\MESH\');
addpath( 'c:\repos\msh\tools\');
addpath( 'c:\repos\msh\uiTools\');
addpath( 'c:\repos\msh\BVH\');     %motor de rayos de silhouette (bvhIntersectRay_mx)

%%

H = load('H','H').H;
hP = plotMESH( H ,'ne','nf');
silhouette( hP ,'EdgeColor','r','LineWidth',3,'UserData',[0.2,1])
