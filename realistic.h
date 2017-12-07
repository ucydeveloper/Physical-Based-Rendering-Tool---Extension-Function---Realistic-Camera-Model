#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

// Example representation of an autofocus zone.
class AfZone {
	public:
	  // from the data file
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};

struct LensInterface {
    float radius, thickness, index_of_refraction, aperture;
};

class RealisticCamera : public Camera {
public:
    RealisticCamera(const AnimatedTransform &cam2world,
        float hither, float yon, float sopen,
        float sclose, float filmdistance, float aperture_diameter,
        const string &specfile, const string &autofocusfile,
        float filmdiag, Film *film, const string &bokeh, int n_blades);
    ~RealisticCamera();

    float GenerateRay(const CameraSample &sample, Ray *) const;
    //void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
    void  ParseAfZones(const string& filename);

private:
    bool  autofocus;
    vector<AfZone> afZones;
    float ShutterOpen;
    float ShutterClose;
    Film * film;

    float units_per_pixel;
    vector<LensInterface> lenses;
    bool  RayLensRefraction(Ray *r, const LensInterface & lens, float origin_z, float n1, float n2) const;
    float AutoFocusScore(Renderer * renderer, const Scene * scene, Sample * origSample, bool quick);

    int bokeh_type;
    int bokeh_n_blades;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif