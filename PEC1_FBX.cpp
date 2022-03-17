#include <SDL.h>
#include <SDL_image.h>
#include <SDL_mixer.h>
#include <iostream>
#include <cmath>

#include "vector.h"
#include "matrix.h"

//Screen dimension constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

//The window we'll be rendering to
SDL_Window* window = NULL;
//The surface contained by the window
SDL_Surface* screenSurface = NULL;

#define FPS 60
int lastTime = 0, currentTime, deltaTime;
float msFrame = 1 / (FPS / 1000.0f);

// size of the spot light
#define LIGHTSIZE 2.4f
#define LIGHT_PIXEL_RES 256
// contains the precalculated spotlight
unsigned char *light; 
// image color
SDL_Surface *image;
// imagen bump
SDL_Surface *bump;

// define the distortion buffer movement
int windowx1, windowy1, windowx2, windowy2, windowZ, windowx3, windowy3, windowx4, windowy4, windowx5, windowy5;

bool initSDL();
void update();
void render();

void close();
void waitTime();

typedef enum { NIEVE = 0, LAVA, FUEGO, PARTICULAS, DISTORSION, BUMP_MAP } Efecto;

Efecto efectoActual;

// bump map
void initBumpMap();
void updateBumpMap();
void renderBumpMap();
void closeBumpMap();

void Compute_Light();
void Bump();

// nieve
#define MAXSNOW 256

struct TSnow
{
	float x, y;             // posicion
	unsigned char plane;    // plano
};
TSnow* snow;

void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel);

void initSnow();
void updateSnow();
void renderSnow();
void closeSnow();

// lava
unsigned char* plasma1;
unsigned char* plasma2;
long src1, src2;

struct RGBColor { unsigned char R, G, B; };
RGBColor palette[256];

void initPlasma();
void updatePlasma();
void renderPlasma();
void closePlasma();

void buildPalettePlasma();

// fuego
unsigned char* fire1;
unsigned char* fire2;
unsigned char* tmp;

void initFire();
void updateFire();
void renderFire();
void closeFire();

void buildPaletteFuego();
void Shade_Pal(int s, int e, int r1, int g1, int b1, int r2, int g2, int b2);

void Blur_Up(unsigned char* src, unsigned char* dst);
void NuevosPuntosFuego(unsigned char* dst);

// particulas
#define MAXPTS 1500

MATRIX obj;
float base_dist;
VECTOR* pts;

int scaleX[SCREEN_WIDTH];
int scaleY[SCREEN_HEIGHT];
int index;

SDL_Surface* secondScreen;

void initParticles();
void updateParticles();
void renderParticles();
void closeParticles();

void RescaleParticula(SDL_Surface* src, SDL_Surface* dst);
void BlurParticula(SDL_Surface* src, SDL_Surface* dst);
void DrawParticula(SDL_Surface* where, VECTOR v);

// distorsion
char* dispX, * dispY;

void initDistortion();
void updateDistortion();
void renderDistortion();
void closeDistorsion();

void precalculate();
void Distort();
void Distort_Bili();

//musica
Mix_Music* mySong;
#define BPM_MUSIC 128
#define MSEG_BPM (60000 / BPM_MUSIC)
#define FLASH_MAX_TIME 300
int flashtime;
int MusicCurrentTime;
int MusicCurrentTimeBeat;
int MusicCurrentBeat;
int MusicPreviousBeat;

void initMusic();
void updateMusic();
void renderMusic();

// transiciones
Uint32 Backgroundcolor;
bool transicion;
SDL_Rect transRec;
int transVelocidad;
void renderTransicion();
void updateTransicion();

int main(int argc, char* args[])
{
	//Start up SDL and create window
	if (!initSDL())
	{
		std::cout << "Failed to initialize!\n";
		return 1;
	}
	else
	{
		IMG_Init(IMG_INIT_PNG);
		initMusic();
		transicion = false;
		transVelocidad = 5;
		Backgroundcolor = SDL_MapRGB(screenSurface->format, 255, 255, 255);

		efectoActual = PARTICULAS;
		switch (efectoActual)
		{
		case NIEVE:
			initSnow();
			break;

		case LAVA:
			initPlasma();
			break;

		case FUEGO:
			initFire();
			break;

		case PARTICULAS:
			initParticles();
			break;

		case DISTORSION:
			initDistortion();
			break;

		case BUMP_MAP:
			initBumpMap();
			break;
		}

		//Main loop flag
		bool quit = false;

		//Event handler
		SDL_Event e;

		//While application is running
		while (!quit)
		{
			//Handle events on queue
			while (SDL_PollEvent(&e) != 0)
			{
				if (e.type == SDL_KEYDOWN) {
					if (e.key.keysym.scancode == SDL_SCANCODE_ESCAPE) {
						quit = true;
					}
				}
				//User requests quit
				if (e.type == SDL_QUIT)
				{
					quit = true;
				}
			}

			// updates all
			update();

			//Render
			render();

			//Update the surface
			SDL_UpdateWindowSurface(window);
			waitTime();
		}
	}

	//Free resources and close SDL
	close();

	return 0;
}

bool initSDL() {
	//Initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
	{
		std::cout << "SDL could not initialize! SDL_Error: %s\n" << SDL_GetError();
		return false;
	}
	//Create window
	window = SDL_CreateWindow("PEC1 - Efectos visuales y sonoros", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);

	if (window == NULL)
	{
		std::cout << "Window could not be created! SDL_Error: %s\n" << SDL_GetError();
		return false;
	}
	//Get window surface
	screenSurface = SDL_GetWindowSurface(window);
	return true;
}

void update() {
	updateMusic();

	switch (efectoActual)
	{
	case NIEVE:
		updateSnow();
		break;

	case LAVA:
		updatePlasma();
		break;

	case FUEGO:
		updateFire();
		break;

	case PARTICULAS:
		updateParticles();

	case DISTORSION:
		updateDistortion();
		break;

	case BUMP_MAP:
		updateBumpMap();
		break;
	}

	if (transicion) updateTransicion();
}

void render() {
	renderMusic();

	switch (efectoActual)
	{
	case NIEVE:
		SDL_FillRect(screenSurface, NULL, SDL_MapRGB(screenSurface->format, 0, 0, 128));
		renderSnow();
		break;

	case LAVA:
		renderPlasma();
		break;

	case FUEGO:
		renderFire();
		break;

	case PARTICULAS:
		renderParticles();
		break;

	case DISTORSION:
		renderDistortion();
		break;

	case BUMP_MAP:
		renderBumpMap();
		break;
	}

	if (transicion) renderTransicion();
}

void close() {
	switch (efectoActual)
	{
	case NIEVE:
		closeSnow();
		break;

	case LAVA:
		closePlasma();
		break;

	case FUEGO:
		closeFire();
		break;

	case PARTICULAS:
		closeParticles();
		break;

	case DISTORSION:
		closeDistorsion();
		break;

	case BUMP_MAP:
		closeBumpMap();
		break;
	}

	Mix_HaltMusic();
	Mix_FreeMusic(mySong);
	Mix_Quit();
	Mix_CloseAudio();

	//Destroy window
	SDL_DestroyWindow(window);
	//Quit SDL subsystems
	SDL_Quit();
}

void waitTime() {
	currentTime = SDL_GetTicks();
	deltaTime = currentTime - lastTime;
	if (deltaTime < (int)msFrame) {
		SDL_Delay((int)msFrame - deltaTime);
	}
	lastTime = currentTime;
}

void initBumpMap() {
	// contains the image of the spotlight
	light = (unsigned char*)malloc(LIGHT_PIXEL_RES * LIGHT_PIXEL_RES);
	// generate the light pattern
	Compute_Light();
	// load the color image
	SDL_Surface* temp = IMG_Load("imagen_bump.png");
	if (temp == NULL) {
		std::cout << "Image can be loaded! " << IMG_GetError();
		close();
		exit(1);
	}
	image = SDL_ConvertSurfaceFormat(temp, SDL_PIXELFORMAT_ARGB8888, 0);
	// load the bump image
	temp = IMG_Load("bump_bump.png");
	if (temp == NULL) {
		std::cout << "Image can be loaded! " << IMG_GetError();
		close();
		exit(1);
	}
	bump = SDL_ConvertSurfaceFormat(temp, SDL_PIXELFORMAT_ARGB8888, 0);
}

void updateBumpMap() {
	// move the lights.... more sines :)
	windowx1 = (int)((LIGHT_PIXEL_RES / 2) * cos((double)currentTime / 640)) - 20;
	windowy1 = (int)((LIGHT_PIXEL_RES / 2) * sin((double)-currentTime / 450)) + 20;

	windowx2 = (int)((LIGHT_PIXEL_RES / 2) * cos((double)-currentTime / 510)) - 20;
	windowy2 = (int)((LIGHT_PIXEL_RES / 2) * sin((double)currentTime / 710)) + 20;

	windowx3 = (int)((LIGHT_PIXEL_RES / 2) * cos((double)-currentTime / 320)) - 20;
	windowy3 = (int)((LIGHT_PIXEL_RES / 2) * sin((double)currentTime / 560)) + 20;

	windowx4 = (int)((LIGHT_PIXEL_RES / 2) * cos((double)-currentTime / 180)) - 20;
	windowy4 = (int)((LIGHT_PIXEL_RES / 2) * sin((double)currentTime / 460)) + 20;

	windowx5 = (int)((LIGHT_PIXEL_RES / 2) * cos((double)-currentTime / 820)) - 20;
	windowy5 = (int)((LIGHT_PIXEL_RES / 2) * sin((double)currentTime / 480)) + 20;

	windowZ = 192 + (int)(((LIGHT_PIXEL_RES / 2) - 1) * sin((double)currentTime / 1120));
}

void renderBumpMap() {
	// draw the bumped image
	Bump();
}

/*
* generate a "spot light" pattern
*/
void Compute_Light()
{
	for (int j = 0; j < LIGHT_PIXEL_RES; j++)
		for (int i = 0; i < LIGHT_PIXEL_RES; i++)
		{
			// get the distance from the centre
			float dist = (float)((LIGHT_PIXEL_RES / 2) - i) * ((LIGHT_PIXEL_RES / 2) - i) + ((LIGHT_PIXEL_RES / 2) - j) * ((LIGHT_PIXEL_RES / 2) - j);
			if (fabs(dist) > 1) dist = sqrt(dist);
			// then fade if according to the distance, and a random coefficient
			int c = (int)(LIGHTSIZE * dist) + (rand() & 7) - 3;
			// clip it
			if (c < 0) c = 0;
			if (c > 255) c = 255;
			// and store it
			light[(j * LIGHT_PIXEL_RES) + i] = 255 - c;
		}
}

/*
* this needs a bump map and a colour map, and 5 light coordinates
* it computes the output colour with the look up table
*/
void Bump()
{
	int i, j, px, py, x, y, c;
	// setup the offsets in the buffers
	Uint8* dst;
	Uint8* initbuffer = (Uint8*)screenSurface->pixels;
	int bpp = screenSurface->format->BytesPerPixel;
	Uint32* imagebuffer = (Uint32*)image->pixels;
	int bppImage = image->format->BytesPerPixel;
	Uint32* bumpbuffer = (Uint32*)bump->pixels;
	int bppBump = bump->format->BytesPerPixel;

	// we skip the first line since there are no pixels above
	// to calculate the slope with
	// loop for all the other lines
	SDL_LockSurface(screenSurface);
	for (j = 1; j < SCREEN_HEIGHT; j++)
	{
		// likewise, skip first pixel since there are no pixels on the left
		dst = initbuffer + j * screenSurface->pitch;
		*(Uint32*)dst = 0;
		for (i = 1; i < SCREEN_WIDTH; i++)
		{
			// calculate coordinates of the pixel we need in light map
			// given the slope at this point, and the zoom coefficient
			SDL_Color Colors[3]; // 0=left pixel 1=center pixel 2=up pixel
			SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + j * image->pitch + (i - 1) * bppImage), image->format, &Colors[0].r, &Colors[0].g, &Colors[0].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + j * image->pitch + i * bppImage), image->format, &Colors[1].r, &Colors[1].g, &Colors[1].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + (j - 1) * image->pitch + i * bppImage), image->format, &Colors[2].r, &Colors[2].g, &Colors[2].b);
			px = (i * windowZ >> 8) + Colors[0].r - Colors[1].r;
			py = (j * windowZ >> 8) + Colors[2].r - Colors[1].r;

			// add the movement of the light 1
			x = px + windowx1;
			y = py + windowy1;
			// check if the coordinates are inside the light buffer
			if ((y >= 0) && (y < LIGHT_PIXEL_RES) && (x >= 0) && (x < LIGHT_PIXEL_RES))
				// if so get the pixel
				c = light[(y * LIGHT_PIXEL_RES) + x];
			// otherwise assume intensity 0
			else c = 0;

			// now do the same for the light 2
			x = px + windowx2;
			y = py + windowy2;
			// this time we add the light's intensity to the first value
			if ((y >= 0) && (y < LIGHT_PIXEL_RES) && (x >= 0) && (x < LIGHT_PIXEL_RES))
				c += light[(y * LIGHT_PIXEL_RES) + x];

			// now do the same for the light 3
			x = px + windowx3;
			y = py + windowy3;
			// this time we add the light's intensity to the first value
			if ((y >= 0) && (y < LIGHT_PIXEL_RES) && (x >= 0) && (x < LIGHT_PIXEL_RES))
				c += light[(y * LIGHT_PIXEL_RES) + x];

			// now do the same for the light 4
			x = px + windowx4;
			y = py + windowy4;
			// this time we add the light's intensity to the first value
			if ((y >= 0) && (y < LIGHT_PIXEL_RES) && (x >= 0) && (x < LIGHT_PIXEL_RES))
				c += light[(y * LIGHT_PIXEL_RES) + x];

			// now do the same for the light 5
			x = px + windowx5;
			y = py + windowy5;
			// this time we add the light's intensity to the first value
			if ((y >= 0) && (y < LIGHT_PIXEL_RES) && (x >= 0) && (x < LIGHT_PIXEL_RES))
				c += light[(y * LIGHT_PIXEL_RES) + x];
			// make sure it's not too big
			if (c > 255) c = 255;

			// look up the colour multiplied by the light coeficient
			SDL_Color ColorBump; // 0=left pixel 1=center pixel 2=up pixel
			SDL_GetRGB(*(Uint32*)((Uint8*)bump->pixels + j * bump->pitch + i * bppBump), bump->format, &ColorBump.r, &ColorBump.g, &ColorBump.b);

			Uint32 Color[3]; // 0=R  1=G  2=B
			Color[0] = (Uint32)((((Colors[0].r * ColorBump.r) / 255)) * (c / 255.0f));
			Color[1] = (Uint32)((((Colors[0].g * ColorBump.g) / 255)) * (c / 255.0f));
			Color[2] = (Uint32)((((Colors[0].b * ColorBump.b) / 255)) * (c / 255.0f));
			if (Color[0] > 255) Color[0] = 255;
			if (Color[1] > 255) Color[1] = 255;
			if (Color[2] > 255) Color[2] = 255;
			Uint32 resultColor = SDL_MapRGB(image->format, Color[0], Color[1], Color[2]);
			*(Uint32*)dst = resultColor;
			dst += bpp;
		}
	}
	SDL_UnlockSurface(screenSurface);
}

void closeBumpMap() {
	free(light);
	SDL_FreeSurface(image);
	SDL_FreeSurface(bump);
}

void initSnow() {
	// allocate memory 
	snow = new TSnow[MAXSNOW];
	// generación random de copos
	for (int i = 0; i < MAXSNOW; i++)
	{
		snow[i].x = (float)(rand() % SCREEN_WIDTH);
		snow[i].y = (float)(rand() % SCREEN_HEIGHT);
		snow[i].plane = rand() % 3; 
	}

}

void updateSnow() {
	for (int i = 0; i < MAXSNOW; i++)
	{
		// movimiento regular descendente, pequeño y aletorio movimiento lateral
		snow[i].y += (deltaTime + (float)snow[i].plane) * 0.1f;
		snow[i].x += ((float)(rand() % 21) - 10) * 0.15f;

		// si el copo se sale de la pantalla verticalmente
		if (snow[i].y > SCREEN_HEIGHT)
		{
			// vuelve a generarse arriba
			snow[i].y = 0;
			// en posicion aletoria horizontal
			snow[i].x = (float)(rand() % SCREEN_WIDTH);
		}
	}
}

void renderSnow() {
	for (int i = 0; i < MAXSNOW; i++)
	{
		// el color depende del plano
		unsigned int color = 0;
		switch (1 + snow[i].plane) {
		case 1:
			color = 0xFF606060; // dark grey
			break;
		case 2:
			color = 0xFFC2C2C2; // light grey
			break;
		case 3:
			color = 0xFFFFFFFF; // white
			break;
		}

		// el tamaño tambien depende del plano
		// (como no hay breaks cuando se ejecuta uno se ejecutan todos los siguientes, de esta forma a más cercano más grande el copo)
		switch (1 + snow[i].plane) {
		case 3:
			putpixel(screenSurface, (int)snow[i].x + 2, (int)snow[i].y, color);
			putpixel(screenSurface, (int)snow[i].x + 2, (int)snow[i].y + 1, color);
			putpixel(screenSurface, (int)snow[i].x, (int)snow[i].y + 2, color);
			putpixel(screenSurface, (int)snow[i].x + 1, (int)snow[i].y + 2, color);
			putpixel(screenSurface, (int)snow[i].x + 2, (int)snow[i].y + 2, color);


		case 2:
			putpixel(screenSurface, (int)snow[i].x + 1, (int)snow[i].y, color);
			putpixel(screenSurface, (int)snow[i].x, (int)snow[i].y + 1, color);
			putpixel(screenSurface, (int)snow[i].x + 1, (int)snow[i].y + 1, color);

		case 1:
			putpixel(screenSurface, (int)snow[i].x, (int)snow[i].y, color);
		}
	}
}

/*
* Set the pixel at (x, y) to the given value
* NOTE: The surface must be locked before calling this!
*/
void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
	// Clipping
	if ((x < 0) || (x >= SCREEN_WIDTH) || (y < 0) || (y >= SCREEN_HEIGHT))
		return;

	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to set */
	Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

	switch (bpp) {
	case 1:
		*p = pixel;
		break;

	case 2:
		*(Uint16*)p = pixel;
		break;

	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
			p[0] = (pixel >> 16) & 0xff;
			p[1] = (pixel >> 8) & 0xff;
			p[2] = pixel & 0xff;
		}
		else {
			p[0] = pixel & 0xff;
			p[1] = (pixel >> 8) & 0xff;
			p[2] = (pixel >> 16) & 0xff;
		}
		break;

	case 4:
		*(Uint32*)p = pixel;
		break;
	}
}

void closeSnow() {
	delete[](snow);
}

void initPlasma() {
	plasma1 = (unsigned char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);
	plasma2 = (unsigned char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);

	int i, j, dst = 0;
	for (j = 0; j < (SCREEN_HEIGHT * 2); j++)
	{
		for (i = 0; i < (SCREEN_WIDTH * 2); i++)
		{
			// creación de los dos plasmas
			plasma1[dst] = (unsigned char)(75 * (sin((double)hypot(SCREEN_HEIGHT - j, SCREEN_WIDTH - i) / 16)));
			plasma2[dst] = (unsigned char)(89 + 150 * sin((float)i / (37 + 15 * cos((float)j / 74))) * cos((float)j / (31 + 11 * sin((float)i / 57))));
			dst++;
		}
	}
}

void updatePlasma() {
	buildPalettePlasma();

	// movimiento de los dos "plasmas"
	windowx1 = (SCREEN_WIDTH / 2) + (int)(((SCREEN_WIDTH / 2) - 1) * cos((double)currentTime / 1356));
	windowx2 = (SCREEN_WIDTH / 2) + (int)(((SCREEN_WIDTH / 2) - 1) * sin((double)-currentTime / 2587));
	windowy1 = (SCREEN_HEIGHT / 2) + (int)(((SCREEN_HEIGHT / 2) - 1) * sin((double)currentTime / 2563));
	windowy2 = (SCREEN_HEIGHT / 2) + (int)(((SCREEN_HEIGHT / 2) - 1) * cos((double)-currentTime / 1254));

	src1 = windowy1 * (SCREEN_WIDTH * 2) + windowx1;
	src2 = windowy2 * (SCREEN_WIDTH * 2) + windowx2;
}

void renderPlasma() {
	Uint8* dst;
	long i, j;
	Uint8* initbuffer = (Uint8*)screenSurface->pixels;
	int bpp = screenSurface->format->BytesPerPixel;

	SDL_LockSurface(screenSurface);

	dst = initbuffer;
	for (j = 0; j < SCREEN_HEIGHT; j++)
	{
		dst = initbuffer + j * screenSurface->pitch;
		for (i = 0; i < SCREEN_WIDTH; i++)
		{
			// el valor del pixel es la suma de los dos plasmas
			unsigned int Color = 0;
			int indexColor = (plasma1[src1] + plasma2[src2]) % 256;
			Color = 0xFF000000 + (palette[indexColor].R << 16) + (palette[indexColor].G << 8) + palette[indexColor].B;
			*(Uint32*)dst = Color;

			dst += bpp;
			src1++; src2++;
		}
		src1 += SCREEN_WIDTH; src2 += SCREEN_WIDTH;
	}
	SDL_UnlockSurface(screenSurface);
}

void buildPalettePlasma() {
	for (int i = 0; i < 256; i++)
	{
		palette[i].R = (unsigned char)(128 + 127 * cos(i * M_PI / 128 + (double)currentTime / 740)) * 0.25;
		palette[i].G = (unsigned char)(128 + 127 * sin(i * M_PI / 128 + (double)currentTime / 630));
		palette[i].B = (unsigned char)(128 - 127 * cos(i * M_PI / 128 + (double)currentTime / 810)) * 0.25;
	}
}

void closePlasma() {
	free(plasma1);
	free(plasma2);
}

void initFire() {
	buildPaletteFuego();

	// two fire buffers
	fire1 = (unsigned char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT);
	fire2 = (unsigned char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT);

	// clear the buffers
	memset(fire1, 0, SCREEN_WIDTH * SCREEN_HEIGHT);
	memset(fire2, 0, SCREEN_WIDTH * SCREEN_HEIGHT);
}

void updateFire() {
	// swap our two fire buffers
	tmp = fire1;
	fire1 = fire2;
	fire2 = tmp;

	// heat the fire
	NuevosPuntosFuego(fire1);
	// apply the filter
	Blur_Up(fire1, fire2);
}

void renderFire() {
	Uint8* dst;
	int src = 0;
	long i, j;
	Uint8* initbuffer = (Uint8*)screenSurface->pixels;
	int bpp = screenSurface->format->BytesPerPixel;

	SDL_LockSurface(screenSurface);
	dst = initbuffer;

	for (j = 0; j < (SCREEN_HEIGHT - 3); j++) // Menos las 3 ultimas lineas
	{
		dst = initbuffer + j * screenSurface->pitch;

		for (i = 0; i < SCREEN_WIDTH; i++)
		{
			unsigned int Color = 0;
			int indexColor = fire2[src];
			Color = 0xFF000000 + (palette[indexColor].R << 16) + (palette[indexColor].G << 8) + palette[indexColor].B;
			*(Uint32*)dst = Color;
			dst += bpp;
			src++;
		}
	}
	SDL_UnlockSurface(screenSurface);
}

void buildPaletteFuego() {
	Shade_Pal(128, 255, 0, 0, 0, 0, 32, 64);
	Shade_Pal(64, 127, 0, 32, 64, 0, 255, 0);
	Shade_Pal(48, 63, 0, 255, 0, 255, 255, 0);
	Shade_Pal(24, 47, 255, 255, 0, 255, 255, 255);
	Shade_Pal(0, 23, 255, 255, 255, 255, 255, 255);
}

void Shade_Pal(int s, int e, int r1, int g1, int b1, int r2, int g2, int b2)
{
	int i;
	float k;
	for (i = 0; i <= e - s; i++)
	{
		k = (float)i / (float)(e - s);
		palette[s + i].R = (int)(r1 + (r2 - r1) * k);
		palette[s + i].G = (int)(g1 + (g2 - g1) * k);
		palette[s + i].B = (int)(b1 + (b2 - b1) * k);
	}
}

/*
* Crea nuevos puntos en toda la pantalla
*/
void NuevosPuntosFuego(unsigned char* dst)
{
	int j = (rand() % 1500);

	// Genera puntos en un área centrada equivalente a 1/5 de la pantalla
	for (int i = 0; i < j; i++)
	{
		int h = (rand() % SCREEN_HEIGHT);
		int w = (rand() % SCREEN_WIDTH);

		dst[SCREEN_WIDTH * h + w] = 255;
	}
}

void Blur_Up(unsigned char* src, unsigned char* dst)
{
	int offs = 0;
	unsigned char b;
	// primeras dos filas con 0
	memset(&dst[offs], 0, SCREEN_WIDTH * 2);
	offs += SCREEN_WIDTH * 2;

	for (int j = 2; j < (SCREEN_HEIGHT - 2); j++)
	{
		// primeros dos pixeles de la fila como 0
		dst[offs] = 0; offs++;
		dst[offs] = 0; offs++;
		
		for (int i = 2; i < (SCREEN_WIDTH - 2); i++)
		{
			// calculo de la media de la mascara
			b = (int)(
				+src[offs - ((SCREEN_WIDTH * 2) + 2)] + src[offs - ((SCREEN_WIDTH * 2) + 1)] * 0.5
				+ src[offs - (SCREEN_WIDTH + 2)] * 0.5 + src[offs - (SCREEN_WIDTH + 1)] * 2.0f + src[offs - SCREEN_WIDTH] * 0.5
				+ src[offs + 1] * 0.5) / 5;

			// guardado del valor del pixel
			dst[offs] = b;
			offs++;
		}
		// ultimos dos pixeles de la fila como 0
		dst[offs] = 0; offs++;
		dst[offs] = 0; offs++;
	}

	// ultimas dos filas con 0
	memset(&dst[offs], 0, SCREEN_WIDTH * 2);
}

void closeFire() {
	free(fire1);
	free(fire2);
}

void initParticles() {
	// generacion de los puntos
	pts = new VECTOR[MAXPTS];
	for (int i = 0; i < MAXPTS; i++) {
		pts[i] = (rotX(2.0f * M_PI * sin((float)i / 100))
			* rotY(2.0f * M_PI * cos((float)i + 100))
			* rotZ(-2.0f * M_PI * cos((float)i + 100)))
			* VECTOR(64 + 16 * sin((float)i + 100), 0, 0);
	}

	secondScreen = SDL_CreateRGBSurface(0, SCREEN_WIDTH, SCREEN_HEIGHT, 32, 0, 0, 0, 0);
}

void updateParticles() {
	// recompute parameters for image rescaling
	int sx = (int)((SCREEN_WIDTH / 2) - (SCREEN_WIDTH / 4) * sin((float)currentTime / 5590)),
		sy = (int)((SCREEN_HEIGHT / 2) + (SCREEN_HEIGHT / 4) * sin((float)currentTime / 6110));
	for (int i = 0; i < SCREEN_WIDTH; i++) scaleX[i] = (int)(sx + (i - sx) * 0.85f);
	for (int i = 0; i < SCREEN_HEIGHT; i++) scaleY[i] = (int)(sy + (i - sy) * 0.85f);

	// setup the position of the object
	base_dist = 128 + 64 * sin((float)currentTime / 5632);
	obj = rotX(2.0f * M_PI * sin((float)currentTime / 4521))
		* rotY(2.0f * M_PI * cos((float)currentTime / 1593))
		* rotZ(-2.0f * M_PI * sin((float)currentTime / 2365));
}

void renderParticles() {
	// rescale the image
	RescaleParticula(screenSurface, secondScreen);

	// blur it
	BlurParticula(secondScreen, screenSurface);

	// draw the particles
	for (int i = 0; i < MAXPTS; i++) {
		DrawParticula(screenSurface, obj * pts[i]);
	}
}

void RescaleParticula(SDL_Surface* src, SDL_Surface* dst)
{
	Uint8* dstPixel;
	Uint8* initbuffer = (Uint8*)src->pixels;
	Uint8* finalbuffer = (Uint8*)dst->pixels;
	int initbufferbpp = src->format->BytesPerPixel;
	int finalbufferbpp = dst->format->BytesPerPixel;

	for (int j = 0; j < SCREEN_HEIGHT; j++)
	{
		dstPixel = finalbuffer + j * dst->pitch;

		for (int i = 0; i < SCREEN_WIDTH; i++)
		{
			// get value from pixel in scaled image, and store
			SDL_Color resultColor;
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + scaleY[j] * src->pitch + scaleX[i] * initbufferbpp), src->format, &resultColor.r, &resultColor.g, &resultColor.b);
			*(Uint32*)dstPixel = SDL_MapRGB(dst->format, resultColor.r, resultColor.g, resultColor.b);
			dstPixel += finalbufferbpp;
		}
	}
}

void BlurParticula(SDL_Surface* src, SDL_Surface* dst)
{
	Uint8* dstPixel;
	Uint8* initbuffer = (Uint8*)src->pixels;
	Uint8* finalbuffer = (Uint8*)dst->pixels;
	int initbufferbpp = src->format->BytesPerPixel;
	int finalbufferbpp = dst->format->BytesPerPixel;

	//first line black
	memset(finalbuffer, 0, SCREEN_WIDTH * finalbufferbpp);

	for (int j = 1; j < (SCREEN_HEIGHT - 1); j++)
	{
		dstPixel = finalbuffer + j * dst->pitch;

		// set first pixel of the line to 0
		*(Uint32*)dstPixel = 0;
		dstPixel += finalbufferbpp;

		// calculate the filter for all the other pixels
		for (int i = 1; i < (SCREEN_WIDTH - 1); i++)
		{
			// calculate the average
			SDL_Color resultColor[8];
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j - 1) * src->pitch + (i - 1) * initbufferbpp), src->format, &resultColor[0].r, &resultColor[0].g, &resultColor[0].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j - 1) * src->pitch + (i)*initbufferbpp), src->format, &resultColor[1].r, &resultColor[1].g, &resultColor[1].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j - 1) * src->pitch + (i + 1) * initbufferbpp), src->format, &resultColor[2].r, &resultColor[2].g, &resultColor[2].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j)*src->pitch + (i - 1) * initbufferbpp), src->format, &resultColor[3].r, &resultColor[3].g, &resultColor[3].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j)*src->pitch + (i + 1) * initbufferbpp), src->format, &resultColor[4].r, &resultColor[4].g, &resultColor[4].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j + 1) * src->pitch + (i - 1) * initbufferbpp), src->format, &resultColor[5].r, &resultColor[5].g, &resultColor[5].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j + 1) * src->pitch + (i)*initbufferbpp), src->format, &resultColor[6].r, &resultColor[6].g, &resultColor[6].b);
			SDL_GetRGB(*(Uint32*)((Uint8*)src->pixels + (j + 1) * src->pitch + (i + 1) * initbufferbpp), src->format, &resultColor[7].r, &resultColor[7].g, &resultColor[7].b);
			SDL_Color medianColor = { 0,0,0 };
			for (int c = 0; c < 8; c++) {
				medianColor.r += resultColor[c].r;
				medianColor.g += resultColor[c].g;
				medianColor.b += resultColor[c].b;
			}
			medianColor.r = medianColor.r / 8;
			medianColor.g = medianColor.g / 8;
			medianColor.b = medianColor.b / 8;

			// store the pixel
			Uint32 finalColor = SDL_MapRGB(dst->format, medianColor.r, medianColor.g, medianColor.b);
			*(Uint32*)dstPixel = finalColor;
			dstPixel += finalbufferbpp;
		}

		// set last pixel of the line to 0
		*(Uint32*)dstPixel = 0;
		dstPixel += finalbufferbpp;
	}
}

void DrawParticula(SDL_Surface* where, VECTOR v)
{
	// coordenadas
	float iz = 1 / (v[2] + base_dist),
		x = (SCREEN_WIDTH / 2) + (SCREEN_WIDTH / 2) * v[0] * iz,
		y = (SCREEN_HEIGHT / 2) + (SCREEN_HEIGHT / 2) * v[1] * iz;

	// clipping
	if ((x < 0) || (x > (SCREEN_WIDTH - 1)) || (y < 0) || (y > (SCREEN_HEIGHT - 1))) return;
	// convert to fixed point to help antialiasing
	int sx = (int)(x * 8.0f),
		sy = (int)(y * 8.0f);

	// calculo del color segun la Z, solo se usa el canal verde, y muy diluido
	Uint32 colorZ = (int)(iz * 30000);
	if (colorZ > 255) colorZ = 255;
	/*Uint32 color = 0xFF000000 + (colorZ << 8);*/
	Uint32 color = 0xFF000000 + (colorZ << index);
	index += 8;
	if (index > 16) index = 0;

	// draw particle
	*((Uint32*)((Uint8*)where->pixels + (sy / 8) * where->pitch + (sx / 8) * where->format->BytesPerPixel)) = color;
}

void closeParticles() {
	SDL_FreeSurface(secondScreen);
}

void initDistortion() {
	// two buffers
	dispX = (char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);
	dispY = (char*)malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);
	// create two distortion functions
	precalculate();
	// load the background image
	SDL_Surface* temp = IMG_Load("Titulo.png");
	if (temp == NULL) {
		std::cout << "Image can be loaded! " << IMG_GetError();
		close();
		exit(1);
	}
	image = SDL_ConvertSurfaceFormat(temp, SDL_PIXELFORMAT_ARGB8888, 0);
}

void updateDistortion() {
	// move distortion buffer
	windowx1 = (SCREEN_WIDTH / 2) + (int)(((SCREEN_WIDTH / 2) - 1) * cos((double)currentTime / 2050));
	windowx2 = (SCREEN_WIDTH / 2) + (int)(((SCREEN_WIDTH / 2) - 1) * sin((double)-currentTime / 1970));
	windowy1 = (SCREEN_HEIGHT / 2) + (int)(((SCREEN_HEIGHT / 2) - 1) * sin((double)currentTime / 2310));
	windowy2 = (SCREEN_HEIGHT / 2) + (int)(((SCREEN_HEIGHT / 2) - 1) * cos((double)-currentTime / 2240));
}

void renderDistortion() {
	// draw the effect showing without filter and with filter each 2 seconds
	if ((currentTime & 2048) < 1024) {
		Distort();
	}
	else {
		Distort_Bili();
	}
}

/*
* calculate a distorion function for X and Y in 5.3 fixed point
*/
void precalculate()
{
	int i, j, dst;
	dst = 0;
	for (j = 0; j < (SCREEN_HEIGHT * 2); j++)
	{
		for (i = 0; i < (SCREEN_WIDTH * 2); i++)
		{
			float x = (float)i;
			float y = (float)j;
			
			// como solo queremos distorsion vertical, ponemos la horizontal a 0
			dispX[dst] = (signed char)(0);
			dispY[dst] = (signed char)(15 * (cos(x / 300) + cos(x * y / 15000)));
			dst++;
		}
	}
}

/*
*   copy an image to the screen with added distortion.
*   no bilinear filtering.
*/
void Distort()
{
	// setup the offsets in the buffers
	Uint8* dst;
	int	src1 = windowy1 * (SCREEN_WIDTH * 2) + windowx1,
		src2 = windowy2 * (SCREEN_WIDTH * 2) + windowx2;
	int dX, dY;
	Uint8* initbuffer = (Uint8*)screenSurface->pixels;
	int bpp = screenSurface->format->BytesPerPixel;
	Uint8* imagebuffer = (Uint8*)image->pixels;
	int bppImage = image->format->BytesPerPixel;

	SDL_LockSurface(screenSurface);
	// loop for all lines
	for (int j = 0; j < SCREEN_HEIGHT; j++)
	{
		dst = initbuffer + j * screenSurface->pitch;
		// for all pixels
		for (int i = 0; i < SCREEN_WIDTH; i++)
		{
			// get distorted coordinates, use the integer part of the distortion
			// buffers and truncate to closest texel
			dY = j + (dispY[src1] >> 3);
			dX = i + (dispX[src2] >> 3);
			// check the texel is valid
			if ((dY >= 0) && (dY < (SCREEN_HEIGHT - 1)) && (dX >= 0) && (dX < (SCREEN_WIDTH - 1)))
			{
				// copy it to the screen
				Uint8* p = (Uint8*)imagebuffer + dY * image->pitch + dX * bppImage;
				*(Uint32*)dst = *(Uint32*)p;
			}
			// otherwise, just set it to black
			else *(Uint32*)dst = 0;
			// next pixel
			dst += bpp;
			src1++; src2++;
		}
		// next line
		src1 += SCREEN_WIDTH;
		src2 += SCREEN_WIDTH;
	}
	SDL_UnlockSurface(screenSurface);
}

/*
*   copy an image to the screen with added distortion.
*   with bilinear filtering.
*/
void Distort_Bili()
{
	// setup the offsets in the buffers
	Uint8* dst;
	int src1 = windowy1 * (SCREEN_WIDTH * 2) + windowx1,
		src2 = windowy2 * (SCREEN_WIDTH * 2) + windowx2;
	int dX, dY, cX, cY;
	Uint8* initbuffer = (Uint8*)screenSurface->pixels;
	int bpp = screenSurface->format->BytesPerPixel;
	Uint32* imagebuffer = (Uint32*)image->pixels;
	int bppImage = image->format->BytesPerPixel;

	SDL_LockSurface(screenSurface);
	// loop for all lines
	for (int j = 0; j < SCREEN_HEIGHT; j++)
	{
		dst = initbuffer + j * screenSurface->pitch;
		// for all pixels
		for (int i = 0; i < SCREEN_WIDTH; i++)
		{
			// get distorted coordinates, by using the truncated integer part
			// of the distortion coefficients
			dY = j + (dispY[src1] >> 3);
			dX = i + (dispX[src2] >> 3);
			// get the linear interpolation coefficiants by using the fractionnal
			// part of the distortion coefficients
			cY = dispY[src1] & 0x7;
			cX = dispX[src2] & 0x7;
			// check if the texel is valid
			if ((dY >= 0) && (dY < (SCREEN_HEIGHT - 1)) && (dX >= 0) && (dX < (SCREEN_WIDTH - 1)))
			{
				// load the 4 surrounding texels and multiply them by the
				// right bilinear coefficients, then get rid of the fractionnal
				// part by shifting right by 6
				SDL_Color Colorvalues[4];
				SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + dY * image->pitch + dX * bppImage), image->format, &Colorvalues[0].r, &Colorvalues[0].g, &Colorvalues[0].b);
				SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + dY * image->pitch + (dX + 1) * bppImage), image->format, &Colorvalues[1].r, &Colorvalues[1].g, &Colorvalues[1].b);
				SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + (dY + 1) * image->pitch + dX * bppImage), image->format, &Colorvalues[2].r, &Colorvalues[2].g, &Colorvalues[2].b);
				SDL_GetRGB(*(Uint32*)((Uint8*)image->pixels + (dY + 1) * image->pitch + (dX + 1) * bppImage), image->format, &Colorvalues[3].r, &Colorvalues[3].g, &Colorvalues[3].b);
				Uint32 Color[3]; // 0=R  1=G  2=B
				Color[0] = (Colorvalues[0].r * (0x8 - cX) * (0x8 - cY) +
					Colorvalues[1].r * cX * (0x8 - cY) +
					Colorvalues[2].r * (0x8 - cX) * cY +
					Colorvalues[3].r * cX * cY) >> 6;
				Color[1] = (Colorvalues[0].g * (0x8 - cX) * (0x8 - cY) +
					Colorvalues[1].g * cX * (0x8 - cY) +
					Colorvalues[2].g * (0x8 - cX) * cY +
					Colorvalues[3].g * cX * cY) >> 6;
				Color[2] = (Colorvalues[0].b * (0x8 - cX) * (0x8 - cY) +
					Colorvalues[1].b * cX * (0x8 - cY) +
					Colorvalues[2].b * (0x8 - cX) * cY +
					Colorvalues[3].b * cX * cY) >> 6;
				Uint32 resultColor = SDL_MapRGB(image->format, Color[0], Color[1], Color[2]);
				*(Uint32*)dst = resultColor;
			}
			// otherwise, just make it black
			else *(Uint32*)dst = 0;
			dst += bpp;
			src1++; src2++;
		}
		// next line
		src1 += SCREEN_WIDTH;
		src2 += SCREEN_WIDTH;
	}
	SDL_UnlockSurface(screenSurface);
}

void closeDistorsion() {
	free(dispX);
	free(dispY);
	SDL_FreeSurface(image);
}

void initMusic() {
	Mix_OpenAudio(44100, MIX_DEFAULT_FORMAT, 2, 4096);
	Mix_Init(MIX_INIT_OGG);
	mySong = Mix_LoadMUS("atmo.ogg");
	if (!mySong) {
		std::cout << "Error loading Music: " << Mix_GetError() << std::endl;
		close();
		exit(1);
	}
	Mix_PlayMusic(mySong, 0);
	flashtime = 0;
	MusicCurrentTime = 0;
	MusicCurrentTimeBeat = 0;
	MusicCurrentBeat = 0;
	MusicPreviousBeat = -1;
}

void updateMusic() {
	MusicCurrentTime += deltaTime;
	MusicCurrentTimeBeat += deltaTime;
	MusicPreviousBeat = MusicCurrentBeat;

	if (MusicCurrentTimeBeat >= MSEG_BPM) {
		MusicCurrentTimeBeat = 0;
		MusicCurrentBeat++;
		flashtime = FLASH_MAX_TIME;

		if (MusicCurrentBeat % 20 == 0) {
			transicion = true;

			switch (efectoActual)
			{
			case PARTICULAS:
				transRec.x = SCREEN_WIDTH;
				transRec.y = 0;
				transRec.w = 0;
				transRec.h = SCREEN_HEIGHT;
				Backgroundcolor = SDL_MapRGB(screenSurface->format, 0, 0, 0);
				break;

			case BUMP_MAP:
				transRec.x = 0;
				transRec.y = SCREEN_HEIGHT;
				transRec.w = SCREEN_WIDTH;
				transRec.h = 0;
				Backgroundcolor = SDL_MapRGB(screenSurface->format, 255, 255, 255);
				break;

			case FUEGO:
				transRec.x = SCREEN_WIDTH;
				transRec.y = 0;
				transRec.w = 0;
				transRec.h = SCREEN_HEIGHT;
				Backgroundcolor = SDL_MapRGB(screenSurface->format, 0, 0, 128);
				break;

			case NIEVE:
				transRec.x = 0;
				transRec.y = 0;
				transRec.w = 0;
				transRec.h = SCREEN_HEIGHT;
				Backgroundcolor = SDL_MapRGB(screenSurface->format, 0, 255, 0);
				break;

			case LAVA:
				transRec.x = 0;
				transRec.y = 0;
				transRec.w = SCREEN_WIDTH;
				transRec.h = 0;
				Backgroundcolor = SDL_MapRGB(screenSurface->format, 0, 0, 0);
				break;
			}
		}
	}

	if (!Mix_PlayingMusic()) {
		close();
		exit(0);
	}
}

void renderMusic() {
	std::cout << "Beat Time: " << MusicCurrentTimeBeat;
	std::cout << "\tBeat: " << MusicCurrentBeat;
	std::cout << "\tEffect Time: " << flashtime;
	std::cout << "\tMusic Time: " << MusicCurrentTime;
	std::cout << std::endl;
}

void renderTransicion() {
	SDL_FillRect(screenSurface, &transRec, Backgroundcolor);
}

void updateTransicion() {
	if (transRec.w < SCREEN_WIDTH) transRec.w += transVelocidad;
	else if (transRec.h < SCREEN_HEIGHT) transRec.h += transVelocidad;
	if (transRec.x > 0) {
		for (int i = 0; i < transVelocidad; i++) if (transRec.x > 0) transRec.x--;
	}
	else if (transRec.y > 0) {
		for (int i = 0; i < transVelocidad; i++) if (transRec.y > 0) transRec.y--;
	}

	if (transRec.w == SCREEN_WIDTH && transRec.h == SCREEN_HEIGHT == transRec.x == 0 && transRec.y == 0)
	{
		transicion = false;
		switch (efectoActual)
		{
		case PARTICULAS:
			closeParticles();
			initBumpMap();
			efectoActual = BUMP_MAP;
			break;

		case BUMP_MAP:
			closeBumpMap();
			initFire();
			efectoActual = FUEGO;
			break;

		case FUEGO:
			closeFire();
			initSnow();
			efectoActual = NIEVE;
			break;

		case NIEVE:
			closeSnow();
			initPlasma();
			efectoActual = LAVA;
			break;

		case LAVA:
			closePlasma();
			initDistortion();
			efectoActual = DISTORSION;
			break;
		}
	}
}
