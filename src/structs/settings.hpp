#pragma once

enum CAM_MODE { CAM_ARC, CAM_FREE };
enum VOX_MODE { VOX_NONE, VOX_ABS, VOX_EMIT };
enum TEXT_MODE { TEXT_NONE, TEXT_HUD, TEXT_ALL };
enum GRID_MODE { GRID_NONE, GRID_PARTIAL, GRID_ALL };
enum GEOM_MODE { GEOM_NONE, GEOM_WHITE, GEOM_COLOR };

enum MENU_ITEM {
	MENU_DISPLAYMODE,
	MENU_CAMERAMODE,
	MENU_GRIDMODE,
	MENU_TEXTMODE,
	MENU_GEOMMODE,
	MENU_TOGGLE_FRAME,
	MENU_TOGGLE_PATHS,
	MENU_TOGGLE_BOUNDS
};

struct Settings {
	bool drawpaths;
	bool drawframe;
	bool drawbounds;

	CAM_MODE cammode;
	VOX_MODE voxmode;
	TEXT_MODE textmode;
	GRID_MODE gridmode;
	GEOM_MODE geommode;

	Settings() {
		drawpaths = true;
		drawframe = false;
		drawbounds = true;

		cammode = CAM_ARC;
		voxmode = VOX_NONE;
		textmode = TEXT_ALL;
		gridmode = GRID_NONE;
		geommode = GEOM_WHITE;
	}
};
