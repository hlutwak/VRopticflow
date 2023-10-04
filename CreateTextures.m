function [ds,pa] = CreateTextures(ds,pa)

InitializeMatlabOpenGL(1);
    bv = zeros(32);
    wv = ones(32);
    myimg = double(repmat([bv wv; wv bv],32,32) > 0.5);
    mytex = Screen('MakeTexture', ds.w, myimg, [], 1);
    [gltex, gltextarget] = Screen('GetOpenGLTexture', ds.w, mytex);

Screen('BeginOpenGL', ds.w); % Setup the OpenGL rendering context
glEnable(gltextarget); % Enable 2D texture mapping

mag_filter = GL_LINEAR; %GL_NEAREST; %_LINEAR;
min_filter = GL_LINEAR; %GL_NEAREST; %GL_LINEAR_MIPMAP_NEAREST; %_LINEAR


% %% create sphere with noise pattern
% target_texid = glGenTextures(1);
% 
% % Apply regular checkerboard pattern as texture:
% % % 1/f floor
% floorSize = 1024; % power of 2
% 
% [x,~] = meshgrid(-floorSize+1:floorSize,-floorSize+1:floorSize);
% % temp = zeros(floorSize*floorSize*3, 1);
% col = 0;
% 
% noysSlope = 1; %1.5;
% noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
% noys=repmat(noys,[ 1 1 3 ]);
% noys=permute(uint8(noys),[ 3 2 1 ]);
% 
% nSlices = 64;
% nStacks = 32;
% 
%     glBindTexture(gltextarget, gltex);
%     glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);
% %     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
% %     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
% %   glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.REPEAT);
%     glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.REPEAT);
%     glTexParameteri(gltextarget, GL.TEXTURE_WRAP_T, GL.REPEAT);
%     glTexParameteri(gltextarget, GL.TEXTURE_MIN_FILTER, GL.LINEAR_MIPMAP_LINEAR);
% 
% %     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
% %     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, floorSize, floorSize, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
% 
%         % Need mipmapping for trilinear filtering --> Create mipmaps:
%     if ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_EXT_framebuffer_object'))
%         % Ask the hardware to generate all depth levels automatically:
%         glGenerateMipmapEXT(GL.TEXTURE_2D);
%     else
%         % No hardware support for auto-mipmap-generation. Do it "manually":
% 
%         % Use GLU to compute the image resolution mipmap pyramid and create
%         % OpenGL textures ouf of it: This is slow, compared to glGenerateMipmapEXT:
%         r = gluBuild2DMipmaps(gltextarget, GL.LUMINANCE, size(myimg,1), size(myimg,2), GL.LUMINANCE, GL.UNSIGNED_BYTE, uint8(myimg));
%         if r>0
%             error('gluBuild2DMipmaps failed for some reason.');
%         end
%     end
%     glTexParameteri(gltextarget, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
%     
% 
%     quadratic=gluNewQuadric();
%     gluQuadricDrawStyle(quadratic, GLU_FILL);
%     gluQuadricNormals(quadratic, GLU_SMOOTH);
%     gluQuadricTexture(quadratic, GL_TRUE);
%     
%     ds.fixation = glGenLists(1);
%     glNewList(ds.fixation, GL.COMPILE);
%     gluSphere(quadratic, pa.fixationSize,nSlices,nStacks);
%     glEndList();
%     
% 
% 
% % ds.fixation = glGenLists(1);
% % glNewList(ds.fixation, GL.COMPILE);
% % gluQuadricTexture(quadratic, GL.TRUE); 
% % gluSphere(quadratic,pa.fixationSize,nSlices,nStacks);  
% % glEndList();


%% rest of the textures
glEnable(GL_TEXTURE_2D); % Enable 2D texture mapping
ds.wall_texid = glGenTextures(1); % this will be the surround texture with the fixation disk and fixation lines embedded
cubeside1_texid = glGenTextures(1); % this will be the texture for the large paddle faces
cubeside2_texid = glGenTextures(1); % this will be the texture for the small paddle faces
cubeside3_texid = glGenTextures(1); % this will be the texture for the large paddle faces
cubeside4_texid = glGenTextures(1); % this will be the texture for the small paddle faces
cubeside5_texid = glGenTextures(1); % this will be the texture for the large paddle faces
cubeside6_texid = glGenTextures(1); % this will be the texture for the small paddle faces

ds.floor_texid = glGenTextures(1); % this will be the floor texture
ds.ceiling_texid = glGenTextures(1); % this will be the ceiling texture
ds.roomwall_texid = glGenTextures(1); % this will be the wall texture
ds.correct = glGenLists(1);
ds.incorrect = glGenLists(1);

% %% Suround Texture - the fixation plane surround with fixation disk and fixation lines embedded
% 
% halfTextureSize = 1024;%ds.xc;  % needs to be in pixels for meshgrid
% 
% [x,y] = meshgrid(-halfTextureSize+1:halfTextureSize,-halfTextureSize+1:halfTextureSize);
% 
% noysSlope = 1.0; %1.5;
% noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
% noys=repmat(noys,[ 1 1 3 ]);
% noys=permute(uint8(noys),[ 3 2 1 ]);
% 
% xoffset = 0;
% yoffset = 0;
% pa.rmin_bg = 45.6874;% pixels 
% pa.rmax_bg = 137.7631;% pixels  
% pa.rstrip = 11.6268;% this cuts a strip into the fixation disk that has a height the size of the paddle height
% 
% % this code pokes out the transparent aperture
% opaque = ones(size(x'));
% for i = 1:length(xoffset)
%     opaque = min(opaque, ((sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) > pa.rmax_bg)  | ((abs(y'+yoffset(i)) > pa.rstrip) & sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) < pa.rmin_bg)));
% end
% noys(4,:,:) = shiftdim(255 .* opaque, -1); 
% 
% glBindTexture(GL_TEXTURE_2D, ds.wall_texid);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
% % glTexParameterfv(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
% % glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, halfTextureSize*2, halfTextureSize*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
% glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
% 
% corners = [0 0;
%     1 0;
%     1 1;
%     0 1];
% 
% canvas = max(ds.halfHeight, ds.halfWidth); % pick the larger of the two dimensions, since they are not equal
% 
% v=[-canvas -canvas 0 ;... 
%     canvas -canvas 0 ;...
%     canvas canvas 0 ;...
%     -canvas canvas 0]';
% 
% ds.surroundTexture = glGenLists(1);
% glNewList(ds.surroundTexture, GL.COMPILE);
% 
% glBegin(GL.POLYGON);
% for i = 1:4
%     glTexCoord2dv(corners(i,:));
%     glVertex3dv(v(:,i));
% end
% glEnd;
% 
% 
% % Add the fixation lines to the list as well using glLines instead of PTB's DrawLines - more efficient coding
% 
% glLineWidth(2);
% 
% glBegin(GL.LINES);
% glColor3f(0, 0, 0); % black diagonal line left
% glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
% glVertex3f(-pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);
% 
% glColor3f(1, 0, 0); % red diagonal line left
% glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
% glVertex3f(-pa.fixationHalfSquare, pa.fixationHalfSquare, 0);
% 
% glColor3f(1, 0, 0); % red diagonal line right
% glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
% glVertex3f(pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);
% 
% glColor3f(0, 0, 0); % black diagonal line right
% glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
% glVertex3f(pa.fixationHalfSquare, pa.fixationHalfSquare, 0);
% glEnd();
% 
% glEndList(); % 1/f noise texture is complete


%% Paddle textures

% Cube side 1
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoysSlope = 1.1; 
paddleNoys = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys=repmat(paddleNoys,[ 1 1 3 ]);
paddleNoys=permute(uint8(paddleNoys),[ 3 2 1 ]);
paddleNoys(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside1_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,size(x,2),size(x,1),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


% Cube side 2
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoys2 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys2=repmat(paddleNoys2,[ 1 1 3 ]);
paddleNoys2=permute(uint8(paddleNoys2),[ 3 2 1 ]);
paddleNoys2(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside2_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys2);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Cube side 3
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoys3 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys3=repmat(paddleNoys3,[ 1 1 3 ]);
paddleNoys3=permute(uint8(paddleNoys3),[ 3 2 1 ]);
paddleNoys3(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside3_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys3);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Cube side 4
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoys4 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys4=repmat(paddleNoys4,[ 1 1 3 ]);
paddleNoys4=permute(uint8(paddleNoys4),[ 3 2 1 ]);
paddleNoys4(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside4_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys4);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Cube side 5
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoys5 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys5=repmat(paddleNoys5,[ 1 1 3 ]);
paddleNoys5=permute(uint8(paddleNoys5),[ 3 2 1 ]);
paddleNoys5(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside5_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys5);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Cube side 6
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoys6 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys6=repmat(paddleNoys6,[ 1 1 3 ]);
paddleNoys6=permute(uint8(paddleNoys6),[ 3 2 1 ]);
paddleNoys6(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, cubeside6_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys6);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Create the full paddle list and compile it
ds.paddleList = glGenLists(1);
glNewList(ds.paddleList, GL.COMPILE);

% First add the large paddle faces
glBindTexture(GL.TEXTURE_2D,cubeside1_texid);
%Begin drawing of a new polygon:
glBegin(GL.QUADS); 
   % Right face
   glNormal3f(1.0, 0.0, 0.0);
    glTexCoord2f(1.0, 0.0); 
   glVertex3f(pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
    glTexCoord2f(1.0, 1.0); 
   glVertex3f(pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
    glTexCoord2f(0.0, 1.0); 
   glVertex3f(pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
    glTexCoord2f(0.0, 0.0); 
   glVertex3f(pa.paddleHalfWidth, -pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
glEnd();

glBindTexture(GL.TEXTURE_2D,cubeside2_texid);
glBegin(GL.QUADS); 
   % Left Face
    glNormal3f(-1.0, 0.0, 0.0);
    glTexCoord2f(0.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
    glTexCoord2f(1.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight, pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
    glTexCoord2f(1.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
    glTexCoord2f(0.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
glEnd();


glBindTexture(GL.TEXTURE_2D,cubeside3_texid);
glBegin(GL.QUADS); 
   % Top Face
   glNormal3f(0.0, 1.0, 0.0);
   glTexCoord2f(0.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
   glTexCoord2f(0.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
   glTexCoord2f(1.0, 0.0); 
   glVertex3f( pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
   glTexCoord2f(1.0, 1.0); 
   glVertex3f( pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
glEnd();

glBindTexture(GL.TEXTURE_2D,cubeside4_texid);
glBegin(GL.QUADS);
   % Bottom Face
   glNormal3f(0.0, -1.0, 0.0);
   glTexCoord2f(1.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
   glTexCoord2f(0.0, 1.0);
   glVertex3f( pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
    glTexCoord2f(0.0, 0.0); 
   glVertex3f( pa.paddleHalfWidth, -pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
    glTexCoord2f(1.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
glEnd();

% Now add the smaller paddle faces
glBindTexture(GL.TEXTURE_2D,cubeside5_texid);
%Begin drawing of a new polygon:
glBegin(GL.QUADS);
    % Front Face
   glNormal3f(0.0, 0.0, 1.0);
   glTexCoord2f(0.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
   glTexCoord2f(1.0, 0.0); 
   glVertex3f(pa.paddleHalfWidth, -pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
   glTexCoord2f(1.0, 1.0); 
   glVertex3f(pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
   glTexCoord2f(0.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight,  pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
glEnd();
  
glBindTexture(GL.TEXTURE_2D,cubeside6_texid);
glBegin(GL.QUADS);
   % Back Face
   glNormal3f(0.0, 0.0, -1.0);
   glTexCoord2f(1.0, 0.0); 
   glVertex3f(-pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Bottom Right Of The Texture and Quad
   glTexCoord2f(1.0, 1.0); 
   glVertex3f(-pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Right Of The Texture and Quad
   glTexCoord2f(0.0, 1.0); 
   glVertex3f(pa.paddleHalfWidth,  pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Top Left Of The Texture and Quad
   glTexCoord2f(0.0, 0.0); 
   glVertex3f(pa.paddleHalfWidth, -pa.paddleHalfHeight, -pa.paddleHalfDepth);  %% Bottom Left Of The Texture and Quad
glEnd();

glEndList();  % The paddle list is now complete


%% Floor (noise) texture - we will include a floor to help 'ground' participants 
% 

% % 1/f floor
floorSize = 4096; % power of 2

[x,y] = meshgrid(-floorSize+1:floorSize,-floorSize+1:floorSize);
% temp = zeros(floorSize*floorSize*3, 1);
col = 0;

noysSlope = 1; %1.5;
noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
noys=repmat(noys,[ 1 1 3 ]);
noys=permute(uint8(noys),[ 3 2 1 ]);
noys(4,:,:) = 255; 



glBindTexture(GL_TEXTURE_2D, ds.floor_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, floorSize, floorSize, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
% glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, floor_tex_data);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

mul = 1.0;


ds.floorTexture = glGenLists(1); 
glNewList(ds.floorTexture, GL.COMPILE);

glBegin(GL_QUADS);
glTexCoord2f(0, 0);
glVertex3f(-pa.floorWidth/2,pa.floorHeight,ds.viewingDistance);

glTexCoord2f(0 ,.5);
glVertex3f(-pa.floorWidth/2,pa.floorHeight,-pa.floorWidth);

glTexCoord2f(.5,.5);
glVertex3f(pa.floorWidth/2,pa.floorHeight,-pa.floorWidth);

glTexCoord2f(.5,0);
glVertex3f(pa.floorWidth/2,pa.floorHeight,ds.viewingDistance);
glEnd;

% glBegin(GL_QUADS);
% glTexCoord2f(0, 0);
% glVertex3f(-pa.floorWidth/2,pa.floorHeight,ds.viewingDistance);
% 
% glTexCoord2f(0 ,1);
% glVertex3f(-pa.floorWidth/2,pa.floorHeight,-pa.floorWidth);
% 
% glTexCoord2f(1,1);
% glVertex3f(pa.floorWidth/2,pa.floorHeight,-pa.floorWidth);
% 
% glTexCoord2f(1,0);
% glVertex3f(pa.floorWidth/2,pa.floorHeight,ds.viewingDistance);
% glEnd;

glEndList(); % done with the floor

%% Create Target Spheres

highcontrast_texid = glGenTextures(1); % this will be the high contrast target
midcontrast_texid = glGenTextures(1); % this will be the mid contrast target
lowcontrast_texid = glGenTextures(1); % this will be the low contrast target

nSlices = 64;
nStacks = 32;

quadratic=gluNewQuadric();
gluQuadricDrawStyle(quadratic, GLU_FILL);
gluQuadricNormals(quadratic, GLU_SMOOTH);
gluQuadricTexture(quadratic, GL_FALSE);

ds.highcontrastTarget = glGenLists(1);
glNewList(ds.highcontrastTarget, GL.COMPILE);
glColor4f(1,1,1,1);
gluSphere(quadratic,pa.fixationSize,nSlices,nStacks);  
glEndList();

ds.incorrect = glGenLists(1);
glNewList(ds.incorrect, GL.COMPILE);
glColor4f(1,0,0,pa.targetContrast(1));
gluSphere(quadratic,pa.targetSize,nSlices,nStacks);
glEndList();

ds.correct = glGenLists(1);
glNewList(ds.correct, GL.COMPILE);
glColor4f(0,1,0,pa.targetContrast(1));
gluSphere(quadratic,pa.targetSize,nSlices,nStacks);
glEndList();

% ds.dot = glGenLists(1);
% glNewList(ds.fixaiton, GL.COMPILE);
% glBegin(GL_POINTS);
% glColor4f(1.0, 1.0, 1.0, 1.0);
% glVertex3f(0.0, 0.0, 0.0);
% glEnd;
% glEndList();

ds.fixation = glGenLists(1);
glNewList(ds.fixation, GL.COMPILE);
glColor4f(1,1,1,1);
gluSphere(quadratic,pa.fixationSize,nSlices,nStacks);  
glEndList();


% Close the OpenGL rendering context
Screen('EndOpenGL', ds.w);


end