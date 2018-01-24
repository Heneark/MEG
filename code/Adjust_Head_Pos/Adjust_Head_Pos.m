function Adjust_Head_Pos()
% 
% cd Y:\\Epilepto\CAT_Aurelie\MEG\
% 
% % necessaire pour les fonction d'exportation CTF
% addpath(genpath('C:\\Data\Sauvegarde\Programmes_Matlab\spm12\'))

%% Init

separator = '\';

global delta_mov

Header             = {};
Header.blockoffset = [];
GUI = {};
DSpath = {};

pos_ref = 1;   % sample de reference
delta_mov = .05;  % mouvement acceptable en cm


%% Select data
[Path] = uigetdir('.', 'Select DS directory');
[a, b, c] = fileparts(Path);
clear Path

while strcmp(c, '.ds')
    DSpath{end+1} = [a separator b c];
    [Path] = uigetdir(a, 'Select DS directory');
    [a, b, c] = fileparts(Path);
    clear Path
end


%% Import data

data = [];
for xi_dsfile = 1 : length(DSpath)
    
    DS = readCTFds(DSpath{xi_dsfile});
    CharLabel = DS.res4.chanNames;
    CharLabel = deblank(CharLabel);
    
    IndexHC = strmatch('HLC', CharLabel);
    % 27 mesures :
    % standard nasion coil position relative to dewar (m):
    % x
    % y
    % z
    % standard left ear coil position relative to dewar (m):
    % x
    % y
    % z
    % standard right ear coil position relative to dewar (m)
    % x
    % y
    % z
    %
    % Les suivantes sont :
    % measured nasion coil position relative to dewar (cm): x,y,z
    % measured left ear coil position relative to dewar (cm): x,y,z
    % measured right ear coil position relative to dewar (cm): x,y,z
    % measured nasion coil position relative to head (cm): x,y,z
    % measured left ear coil position relative to head (cm): x,y,z
    % measured right ear coil position relative to head (cm): x,y,z
    
    % Donc on ne concerve que les 9 premiers channles
    data = [data; getCTFdata(DS, [], IndexHC(1:9), 'ft', 'double')];
    Header.blockoffset(end+1) = length(data);
    Header.sample_rate = DS.res4.sample_rate;
    clear DS CharLabel IndexHC
    
end

% mesure sont en mètre, transformation en cm
data = data*100;
Header.unit = 'cm';
Header.HClabel = {'Nasion_X'; 'Nasion_Y'; 'Nasion_Z'; 'LeftEar_X'; 'LeftEar_Y'; 'LeftEar_Z'; 'RightEar_X'; 'RightEar_Y'; 'RightEar_Z'};


% position par raport au sample de ref
data = data-repmat(data(1,:),size(data,1),1);

% distance eucclidienne
D_eucl_nas = sqrt((0-data(:,1)).^2 + (0-data(:,2)).^2 + (0-data(:,3)).^2);
D_eucl_ELeft = sqrt((0-data(:,4)).^2 + (0-data(:,5)).^2 + (0-data(:,6)).^2);
D_eucl_ERight = sqrt((0-data(:,7)).^2 + (0-data(:,8)).^2 + (0-data(:,9)).^2);

%% Display

y_mean = 0;

GUI.f = figure;

subplot(3,1,1)
hold on
GUI.nasion.goodarea = fill([1 size(data,1) size(data,1) 1], [y_mean+0.05 y_mean+0.05 y_mean-0.05 y_mean-0.05], 'r');
GUI.nasion.goodarea.Tag = 'goodarea_nas';
GUI.nasion.goodarea.EdgeAlpha = 0;
GUI.nasion.goodarea.FaceAlpha = .05;
GUI.nasion.plot_nasion = plot(D_eucl_nas);
GUI.nasion.meanref = plot([1 size(data,1)], [y_mean y_mean], 'k:');
GUI.nasion.meanref.Tag = 'meanref_nas';
grille = GUI.nasion.plot_nasion.Parent;
GUI.nasion.goodarea.ButtonDownFcn = @clicy;
grille.ButtonDownFcn = @clicy;
grille.YLim = [min(D_eucl_nas)-delta_mov max(D_eucl_nas)+delta_mov];
grille.XLim = [0 length(D_eucl_nas)];
title('Distance du Nasion     ref sample #1')

subplot(3,1,2)
hold on
GUI.ELeft.goodarea = fill([1 size(data,1) size(data,1) 1], [y_mean+0.05 y_mean+0.05 y_mean-0.05 y_mean-0.05], 'r');
GUI.ELeft.goodarea.Tag = 'goodarea_ELeft';
GUI.ELeft.goodarea.EdgeAlpha = 0;
GUI.ELeft.goodarea.FaceAlpha = .05;
GUI.ELeft.plot_ELeft = plot(D_eucl_ELeft);
GUI.ELeft.meanref = plot([1 size(data,1)], [y_mean y_mean], 'k:');
GUI.ELeft.meanref.Tag = 'meanref_ELeft';
grille = GUI.ELeft.plot_ELeft.Parent;
GUI.ELeft.goodarea.ButtonDownFcn = @clicy;
grille.ButtonDownFcn = @clicy;
grille.YLim = [min(D_eucl_ELeft)-delta_mov max(D_eucl_ELeft)+delta_mov];
grille.XLim = [0 length(D_eucl_ELeft)];
title('Distance du ELeft     ref sample #1')

subplot(3,1,3)
hold on
GUI.ERight.goodarea = fill([1 size(data,1) size(data,1) 1], [y_mean+0.05 y_mean+0.05 y_mean-0.05 y_mean-0.05], 'r');
GUI.ERight.goodarea.Tag = 'goodarea_ERight';
GUI.ERight.goodarea.EdgeAlpha = 0;
GUI.ERight.goodarea.FaceAlpha = .05;
GUI.ERight.plot_ERight = plot(D_eucl_ERight);
GUI.ERight.meanref = plot([1 size(data,1)], [y_mean y_mean], 'k:');
GUI.ERight.meanref.Tag = 'meanref_ERight';
grille = GUI.ERight.plot_ERight.Parent;
GUI.ERight.goodarea.ButtonDownFcn = @clicy;
grille.ButtonDownFcn = @clicy;
grille.YLim = [min(D_eucl_ERight)-delta_mov max(D_eucl_ERight)+delta_mov];
grille.XLim = [0 length(D_eucl_ERight)];
title('Distance du ERight     ref sample #1')

function clicy(gcbo, eventdata, handles)

global delta_mov

h = gca;
set(findobj('Tag', 'goodarea_nas'), 'YData', [h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)-delta_mov h.CurrentPoint(1,2)-delta_mov])
set(findobj('Tag', 'meanref_nas'), 'YData', [h.CurrentPoint(1,2) h.CurrentPoint(1,2)])

set(findobj('Tag', 'goodarea_ELeft'), 'YData', [h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)-delta_mov h.CurrentPoint(1,2)-delta_mov])
set(findobj('Tag', 'meanref_ELeft'), 'YData', [h.CurrentPoint(1,2) h.CurrentPoint(1,2)])

set(findobj('Tag', 'goodarea_ERight'), 'YData', [h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)+delta_mov h.CurrentPoint(1,2)-delta_mov h.CurrentPoint(1,2)-delta_mov])
set(findobj('Tag', 'meanref_ERight'), 'YData', [h.CurrentPoint(1,2) h.CurrentPoint(1,2)])

get(gcf,'SelectionType')

% assignin('caller', 'new_x',  h.CurrentPoint(1,1));
% assignin('caller', 'new_y',  h.CurrentPoint(1,2));

