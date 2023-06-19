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


%% create sphere with noise pattern
target_texid = glGenTextures(1);

% Apply regular checkerboard pattern as texture:
% % 1/f floor
floorSize = 1024; % power of 2

[x,~] = meshgrid(-floorSize+1:floorSize,-floorSize+1:floorSize);
% temp = zeros(floorSize*floorSize*3, 1);
col = 0;

noysSlope = 1; %1.5;
noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
noys=repmat(noys,[ 1 1 3 ]);
noys=permute(uint8(noys),[ 3 2 1 ]);

nSlices = 64;
nStacks = 32;

    glBindTexture(gltextarget, gltex);
    glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);
%     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
%     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
%   glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.REPEAT);
    glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.REPEAT);
    glTexParameteri(gltextarget, GL.TEXTURE_WRAP_T, GL.REPEAT);
    glTexParameteri(gltextarget, GL.TEXTURE_MIN_FILTER, GL.LINEAR_MIPMAP_LINEAR);

%     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
%     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, floorSize, floorSize, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);

        % Need mipmapping for trilinear filtering --> Create mipmaps:
    if ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_EXT_framebuffer_object'))
        % Ask the hardware to generate all depth levels automatically:
        glGenerateMipmapEXT(GL.TEXTURE_2D);
    else
        % No hardware support for auto-mipmap-generation. Do it "manually":

        % Use GLU to compute the image resolution mipmap pyramid and create
        % OpenGL textures ouf of it: This is slow, compared to glGenerateMipmapEXT:
        r = gluBuild2DMipmaps(gltextarget, GL.LUMINANCE, size(myimg,1), size(myimg,2), GL.LUMINANCE, GL.UNSIGNED_BYTE, uint8(myimg));
        if r>0
            error('gluBuild2DMipmaps failed for some reason.');
        end
    end
    glTexParameteri(gltextarget, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
    

    quadratic=gluNewQuadric();
    gluQuadricDrawStyle(quadratic, GLU_FILL);
    gluQuadricNormals(quadratic, GLU_SMOOTH);
    gluQuadricTexture(quadratic, GL_TRUE);
    
    ds.fixation = glGenLists(1);
    glNewList(ds.fixation, GL.COMPILE);
    gluSphere(quadratic, pa.fixationSize,nSlices,nStacks);
    glEndList();
    


% ds.fixation = glGenLists(1);
% glNewList(ds.fixation, GL.COMPILE);
% gluQuadricTexture(quadratic, GL.TRUE); 
% gluSphere(quadratic,pa.fixationSize,nSlices,nStacks);  
% glEndList();


%% rest of the textures
glEnable(GL_TEXTURE_2D); % Enable 2D texture mapping
ds.wall_texid = glGenTextures(1); % this will be the surround texture with the fixation disk and fixation lines embedded
largepaddle_texid = glGenTextures(1); % this will be the texture for the large paddle faces
smallpaddle_texid = glGenTextures(1); % this will be the texture for the small paddle faces
ds.floor_texid = glGenTextures(1); % this will be the floor texture
ds.ceiling_texid = glGenTextures(1); % this will be the ceiling texture
ds.roomwall_texid = glGenTextures(1); % this will be the wall texture
ds.correct = glGenLists(1);
ds.incorrect = glGenLists(1);
%% Suround Texture - the fixation plane surround with fixation disk and fixation lines embedded

halfTextureSize = 1024;%ds.xc;  % needs to be in pixels for meshgrid

[x,y] = meshgrid(-halfTextureSize+1:halfTextureSize,-halfTextureSize+1:halfTextureSize);

noysSlope = 1.0; %1.5;
noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
noys=repmat(noys,[ 1 1 3 ]);
noys=permute(uint8(noys),[ 3 2 1 ]);

xoffset = 0;
yoffset = 0;
pa.rmin_bg = 45.6874;% pixels 
pa.rmax_bg = 137.7631;% pixels  
pa.rstrip = 11.6268;% this cuts a strip into the fixation disk that has a height the size of the paddle height

% this code pokes out the transparent aperture
opaque = ones(size(x'));
for i = 1:length(xoffset)
    opaque = min(opaque, ((sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) > pa.rmax_bg)  | ((abs(y'+yoffset(i)) > pa.rstrip) & sqrt((x'+xoffset(i)).^2+(y'+yoffset(i)).^2) < pa.rmin_bg)));
end
noys(4,:,:) = shiftdim(255 .* opaque, -1); 

glBindTexture(GL_TEXTURE_2D, ds.wall_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
% glTexParameterfv(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
% glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, halfTextureSize*2, halfTextureSize*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

corners = [0 0;
    1 0;
    1 1;
    0 1];

canvas = max(ds.halfHeight, ds.halfWidth); % pick the larger of the two dimensions, since they are not equal

v=[-canvas -canvas 0 ;... 
    canvas -canvas 0 ;...
    canvas canvas 0 ;...
    -canvas canvas 0]';

ds.surroundTexture = glGenLists(1);
glNewList(ds.surroundTexture, GL.COMPILE);

glBegin(GL.POLYGON);
for i = 1:4
    glTexCoord2dv(corners(i,:));
    glVertex3dv(v(:,i));
end
glEnd;


% Add the fixation lines to the list as well using glLines instead of PTB's DrawLines - more efficient coding

glLineWidth(2);

glBegin(GL.LINES);
glColor3f(0, 0, 0); % black diagonal line left
glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
glVertex3f(-pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);

glColor3f(1, 0, 0); % red diagonal line left
glVertex3f(-pa.fixationHalfSquare-pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
glVertex3f(-pa.fixationHalfSquare, pa.fixationHalfSquare, 0);

glColor3f(1, 0, 0); % red diagonal line right
glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, -pa.fixationHalfSquare-pa.fixationLineLength, 0);
glVertex3f(pa.fixationHalfSquare, -pa.fixationHalfSquare, 0);

glColor3f(0, 0, 0); % black diagonal line right
glVertex3f(pa.fixationHalfSquare+pa.fixationLineLength, pa.fixationHalfSquare+pa.fixationLineLength, 0);
glVertex3f(pa.fixationHalfSquare, pa.fixationHalfSquare, 0);
glEnd();

glEndList(); % 1/f noise texture is complete


%% Paddle textures

% Large Paddle Faces
% pmsize = 42.7350; 
% pmsize2 = 21.3675; 
pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoysSlope = 1.5; 
paddleNoys = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys=repmat(paddleNoys,[ 1 1 3 ]);
paddleNoys=permute(uint8(paddleNoys),[ 3 2 1 ]);
paddleNoys(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, largepaddle_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,size(x,2),size(x,1),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


% Small Paddle Faces

pmsize = 512;
pmsize2 = 512;

[x,y] = meshgrid(-pmsize:pmsize,-pmsize2:pmsize2);

paddleNoysSlope = 1.5;
paddleNoys2 = 255.*oneoverf(paddleNoysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
paddleNoys2=repmat(paddleNoys2,[ 1 1 3 ]);
paddleNoys2=permute(uint8(paddleNoys2),[ 3 2 1 ]);
paddleNoys2(4,:,:) = 255; 

glBindTexture(GL_TEXTURE_2D, smallpaddle_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size(x,1),size(x,2),0, GL_RGBA, GL_UNSIGNED_BYTE, paddleNoys2);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

% Create the full paddle list and compile it
ds.paddleList = glGenLists(1);
glNewList(ds.paddleList, GL.COMPILE);

% First add the large paddle faces
glBindTexture(GL.TEXTURE_2D,largepaddle_texid);
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
glBindTexture(GL.TEXTURE_2D,smallpaddle_texid);
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




%% Room textures (floor, ceiling, and walls)

roomSize = 0.45*4;
%% Floor (checkerboard) texture - we will include a floor to help 'ground' participants 
% 

% % 1/f floor
floorSize = 512; % power of 2

[x,y] = meshgrid(-floorSize+1:floorSize,-floorSize+1:floorSize);
% temp = zeros(floorSize*floorSize*3, 1);
col = 0;

noysSlope = 1.0; %1.5;
noys = 255.*oneoverf(noysSlope, size(x,1), size(x,2)); % oneoverf -> [0:1]
noys=repmat(noys,[ 1 1 3 ]);
noys=permute(uint8(noys),[ 3 2 1 ]);
noys(4,:,:) = 255; 



glBindTexture(GL_TEXTURE_2D, ds.floor_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, floorSize, floorSize, 0, GL_RGBA, GL_UNSIGNED_BYTE, noys);
% glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, floor_tex_data);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

mul = 1.0;

% ds.floorTexture = glGenLists(1); 
% glNewList(ds.floorTexture, GL.COMPILE);

% glBegin(GL_QUADS);
% glTexCoord2f(0.0 ,0.0);
% glVertex3f(ds.halfWidth*4,pa.floorHeight,0);
% 
% glTexCoord2f(0.0 ,mul*roomSize);
% glVertex3f(ds.halfWidth*4,pa.floorHeight,-roomSize);
% 
% glTexCoord2f(mul*ds.halfWidth*4,mul*roomSize);
% glVertex3f(-ds.halfWidth*4,pa.floorHeight,-roomSize);
% 
% glTexCoord2f(mul*ds.halfWidth*4,0.0);
% glVertex3f(-ds.halfWidth*4,pa.floorHeight,0);
% glEnd;
% 

ds.floorTexture = glGenLists(1); 
glNewList(ds.floorTexture, GL.COMPILE);

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

glBegin(GL_QUADS);
glTexCoord2f(0, 0);
glVertex3f(-pa.floorWidth/2,-pa.floorWidth/2,-5);

glTexCoord2f(0 ,1);
glVertex3f(-pa.floorWidth/2,pa.floorWidth/2,-5);

glTexCoord2f(1,1);
glVertex3f(pa.floorWidth/2,pa.floorWidth/2,-5);

glTexCoord2f(1,0);
glVertex3f(pa.floorWidth/2,-pa.floorWidth/2,-5);
glEnd;

glEndList(); % done with the floor

% OLD FLOOR
% floorSize = 1024; % power of 2
% temp = zeros(floorSize*floorSize*3, 1);
% col = 0;
% 
% for i=1:floorSize
%     for j=1:floorSize
%         if (i<floorSize/2 && j<floorSize/2) || (i>=floorSize/2 && j>=floorSize/2)
%             % white
%             col = 180;
%         else
%             %gray
%             col = 80;
%         end
%         temp((i-1)*floorSize*3 + (j-1)*3+1, 1) = col;
%         temp((i-1)*floorSize*3 + (j-1)*3+2, 1) = col;
%         temp((i-1)*floorSize*3 + (j-1)*3+3, 1) = col;
%     end
% end
% 
% floor_tex_data = uint8(temp);
% 
% glBindTexture(GL_TEXTURE_2D, ds.floor_texid);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
% glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
% glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
% glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, floor_tex_data);
% glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
% 
% mul = 1.0;
% 
% ds.floorTexture = glGenLists(1); 
% glNewList(ds.floorTexture, GL.COMPILE);
% 
% glBegin(GL_QUADS);
% glTexCoord2f(0.0 ,0.0);
% glVertex3f(ds.halfWidth*4,pa.floorHeight,roomSize);
% 
% glTexCoord2f(0.0 ,mul*roomSize);
% glVertex3f(ds.halfWidth*4,pa.floorHeight,-roomSize);
% 
% glTexCoord2f(mul*ds.halfWidth*4,mul*roomSize);
% glVertex3f(-ds.halfWidth*4,pa.floorHeight,-roomSize);
% 
% glTexCoord2f(mul*ds.halfWidth*4,0.0);
% glVertex3f(-ds.halfWidth*4,pa.floorHeight,roomSize);
% glEnd;
% 
% glEndList(); % done with the floor

 
%% Ceiling texture

ceilingSize = 1024;
temp = zeros(ceilingSize*ceilingSize*3, 1);
col = 0;

    for i=1:ceilingSize
        for j=1:ceilingSize
            if ((abs(ceilingSize/2-i) <= (ceilingSize/2)-3) && (abs(ceilingSize/2-j) <= (ceilingSize/2)-3))
                % white
                col = 80;%180;
            else
                %gray
                col = 180;%80;
            end
            temp((i-1)*ceilingSize*3 + (j-1)*3+1, 1) = col;
            temp((i-1)*ceilingSize*3 + (j-1)*3+2, 1) = col;
            temp((i-1)*ceilingSize*3 + (j-1)*3+3, 1) = col;
        end
    end

ceiling_tex_data = uint8(temp);

glBindTexture(GL_TEXTURE_2D, ds.ceiling_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, ceiling_tex_data);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

ds.ceilingTexture = glGenLists(1); 
glNewList(ds.ceilingTexture, GL.COMPILE);

glBegin(GL_QUADS);
glTexCoord2f(0.0 ,0.0);
glVertex3f(ds.halfWidth*4,pa.ceilingHeight,roomSize);

glTexCoord2f(0.0 ,mul*roomSize);
glVertex3f(ds.halfWidth*4,pa.ceilingHeight,-roomSize);

glTexCoord2f(mul*ds.halfWidth*4,mul*roomSize);
glVertex3f(-ds.halfWidth*4,pa.ceilingHeight,-roomSize);

glTexCoord2f(mul*ds.halfWidth*4,0.0);
glVertex3f(-ds.halfWidth*4,pa.ceilingHeight,roomSize);
glEnd;

glEndList(); % done with the ceiling


%% Wall textures

wallSize = 1024;
temp = zeros(wallSize*wallSize*3, 1);
col = 0;

for i=1:wallSize
    for j=1:wallSize
        if (i <= wallSize/2)
            if ((j >= 3) && (j <= wallSize-3) && (i >= 2) && (i <= wallSize-2))
                col = 180;
            else
                col = 60;
            end
        else
            if (abs(wallSize/2-j) <= 2)
                col = 60;
            elseif ((i <= wallSize/2) || (i >= wallSize-2))
                col = 60;
            else
                col = 180;
            end
        end
        
        temp((i-1)*wallSize*3 + (j-1)*3+1, 1) = col;
        temp((i-1)*wallSize*3 + (j-1)*3+2, 1) = col;
        temp((i-1)*wallSize*3 + (j-1)*3+3, 1) = col;
    end
end

wall_tex_data = uint8(temp);
glBindTexture(GL_TEXTURE_2D, ds.roomwall_texid);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, wall_tex_data);
glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

mul = 2.0;

ds.wallTexture = glGenLists(1); 
glNewList(ds.wallTexture, GL.COMPILE);


% right wall
glBegin(GL_QUADS);
glTexCoord2f(mul*roomSize,mul*ds.halfHeight*2);
glVertex3f(ds.halfWidth*4,pa.floorHeight,roomSize);

glTexCoord2f(mul*roomSize,0.0);
glVertex3f(ds.halfWidth*4,ds.halfHeight*2-pa.floorHeight,roomSize);

glTexCoord2f(0.0,0.0);
glVertex3f(ds.halfWidth*4,ds.halfHeight*2-pa.floorHeight,-roomSize);

glTexCoord2f(0.0,mul*ds.halfHeight*2);
glVertex3f(ds.halfWidth*4,pa.floorHeight,-roomSize);

% left wall

glTexCoord2f(mul*roomSize ,mul*ds.halfHeight*2);
glVertex3f(-ds.halfWidth*4,pa.floorHeight,-roomSize);

glTexCoord2f(mul*roomSize,0.0);
glVertex3f(-ds.halfWidth*4,ds.halfHeight*2-pa.floorHeight,-roomSize);

glTexCoord2f(0.0,0.0);
glVertex3f(-ds.halfWidth*4,ds.halfHeight*2-pa.floorHeight,roomSize);

glTexCoord2f(0.0,mul*ds.halfHeight*2);
glVertex3f(-ds.halfWidth*4,pa.floorHeight,roomSize);
glEnd;

glEndList(); % done with the walls



%% Create Target Spheres

highcontrast_texid = glGenTextures(1); % this will be the high contrast target
midcontrast_texid = glGenTextures(1); % this will be the mid contrast target
lowcontrast_texid = glGenTextures(1); % this will be the low contrast target

nSlices = 64;
nStacks = 32;

% quadratic=gluNewQuadric();
% gluQuadricDrawStyle(quadratic, GLU_FILL);
% gluQuadricNormals(quadratic, GLU_SMOOTH);
% gluQuadricTexture(quadratic, GL_FALSE);

ds.highcontrastTarget = glGenLists(1);
glNewList(ds.highcontrastTarget, GL.COMPILE);
glColor4f(1,1,1,pa.targetContrast(1));
gluSphere(quadratic,pa.fixationSize,nSlices,nStacks);  
glEndList();

ds.midcontrastTarget = glGenLists(1);
glNewList(ds.midcontrastTarget, GL.COMPILE);
glColor4f(1,1,1,pa.targetContrast(2));
gluSphere(quadratic,pa.targetSize,nSlices,nStacks);  
glEndList();

ds.lowcontrastTarget = glGenLists(1);
glNewList(ds.lowcontrastTarget, GL.COMPILE);
glColor4f(1,1,1,pa.targetContrast(3));
gluSphere(quadratic,pa.targetSize,nSlices,nStacks);  
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

ds.dot = glGenLists(1);
glNewList(ds.dot, GL.COMPILE);
glBegin(GL_POINTS);
glColor4f(1.0, 1.0, 1.0, 1.0);
glVertex3f(0.0, 0.0, 0.0);
glEnd;
glEndList();


% Close the OpenGL rendering context
Screen('EndOpenGL', ds.w);


end