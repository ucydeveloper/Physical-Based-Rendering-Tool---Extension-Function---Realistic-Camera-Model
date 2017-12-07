#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <deque>
#include <string>

using namespace std;

#define ABS(n) ((n) < 0 ? -(n) : (n))

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
            const AnimatedTransform &cam2world, Film *film) {
    // Extract common camera parameters from \use{ParamSet}
    float hither = params.FindOneFloat("hither", -1);
    float yon = params.FindOneFloat("yon", -1);
    float shutteropen = params.FindOneFloat("shutteropen", -1);
    float shutterclose = params.FindOneFloat("shutterclose", -1);
    
    // Realistic camera-specific parameters
    string specfile = params.FindOneString("specfile", "");
    float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
    float fstop = params.FindOneFloat("aperture_diameter", 1.0);
    float filmdiag = params.FindOneFloat("filmdiag", 35.0);
    string autofocusfile = params.FindOneString("af_zones", "");

    assert(hither != -1 && yon != -1 && shutteropen != -1 &&
            shutterclose != -1 && filmdistance!= -1);

    string bokeh = params.FindOneString("bokeh", "");
    int n_blades = params.FindOneInt("n_blades", 6);

    if (specfile == "") {
        Severe( "No lens spec file supplied!\n" );
    }
    return new RealisticCamera(cam2world, hither, yon,
        shutteropen, shutterclose, filmdistance, fstop,
        specfile, autofocusfile, filmdiag, film, bokeh, n_blades);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter,
                                 const string &specfile,
                                 const string &autofocusfile,
                                 float filmdiag,
                                 Film *f,
                                 const string &bokeh,
                                 int n_blades
                                 )
                                 : Camera(cam2world, sopen, sclose, f),
                                   ShutterOpen(sopen),
                                   ShutterClose(sclose),
                                   film(f)
{

    // YOUR CODE HERE -- build and store datastructures representing the given lens
    // and film placement.
        
    // parse the lens data file
    ifstream specfile_stream;
    specfile_stream.open(specfile.c_str());
    if (!specfile_stream) {
        fprintf(stderr, "Cannot open file %s\n", specfile.c_str());
        exit (-1);
    }
    char line[512];
    while (!specfile_stream.eof()) {
        specfile_stream.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
            lenses.resize(lenses.size()+1);
            LensInterface& li = lenses[lenses.size()-1];
            sscanf(line, "%f %f %f %f\n", &li.radius, &li.thickness, &li.index_of_refraction, &li.aperture);
            li.radius *= -1;
        }
    }
    lenses.back().thickness = filmdistance;

    float pixeldiag = f->xResolution * f->xResolution + f->yResolution * f->yResolution;
    pixeldiag = sqrtf(pixeldiag);
    units_per_pixel = filmdiag/pixeldiag;

    // If 'autofocusfile' is the empty string, then you should do
    // nothing in any subsequent call to AutoFocus()
    autofocus = false;

    if (autofocusfile.compare("") != 0)  {
        ParseAfZones(autofocusfile);
        autofocus = true;
    }

    bokeh_type = 0;
    bokeh_n_blades = n_blades;
    if (bokeh.compare("") != 0) {
        if (bokeh.compare("blades") == 0) {
            bokeh_type = 1;
        } else if (bokeh.compare("half_circle") == 0) {
            bokeh_type = 2;
        } else if (bokeh.compare("heart") == 0) {
            bokeh_type = 3;
        }
    }
    printf("bokeh type is: %i for %s\n", bokeh_type, bokeh.c_str());
    printf("n_blades is %i\n", bokeh_n_blades);
}


// parses the AF zone file
void RealisticCamera::ParseAfZones(const string& filename)
{
    ifstream specfile(filename.c_str());
    if (!specfile) {
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        exit (-1);
    }

    char line[512];

    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
            afZones.resize(afZones.size()+1);
            AfZone& zone = afZones[afZones.size()-1];
            sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
        }
    }

    printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

RealisticCamera::~RealisticCamera()
{

}

// origin_z is the origin of the sphere in camera coordinates
bool RealisticCamera::RayLensRefraction(Ray *r, const LensInterface & lens, float origin_z, float n1, float n2) const {
    float phi;
    Point phit;
    float radius = lens.radius;

    // Transform _Ray_ to "object space" of sphere by doing a translation
    Ray ray(*r);
    ray.o.z -= origin_z;
    //printf("radius + thickness = %f\n", radius + lens.thickness);

    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    float thit = radius > 0 ? t0 : t1;

    // Compute sphere hit position and $\phi$
    phit = ray(thit);
    if (phit.x * phit.x + phit.y * phit.y > (lens.aperture/2)*(lens.aperture/2)) {
        return false;
    }
    if (phit.x == 0.f && phit.y == 0.f) phit.x = 1e-5f * radius;
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;
    
    // use hit position to calculate normal

    Normal n = Normalize(Normal(phit.x, phit.y, phit.z));
    if (radius < 0) n = -n;

    // use snells law to update ray
    float costheta1 = Dot(ray.d, -n);
    float theta1 = acosf(costheta1);
    if (isnan(theta1)) {
        return false;
    }

    // Direction from the ray to the normal tangent to the surface
    Vector dtan = Normalize(ray.d + Vector(costheta1 * n));

    float theta2 = theta1 * n1 / n2;
    ray.o = phit;
    ray.d = sinf(theta2) * dtan - Vector(cosf(theta2) * n);

    *r = ray;
    r->o.z += origin_z;

    return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
  
  // YOUR CODE HERE -- make that ray!
    float time = Lerp(sample.time, shutterOpen, shutterClose);

    float lensU, lensV;
    float aperture_radius = lenses.back().aperture/2;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    lensU *= aperture_radius;
    lensV *= aperture_radius;

    Point p = Point((sample.imageX - film->xResolution/2) * units_per_pixel,(sample.imageY - film->yResolution/2) * units_per_pixel,0);
    Point d = Point (lensU,lensV,lenses.back().thickness);

        *ray = Ray(p, d - p, 0.f, INFINITY, time);
    float L = ray->d.Length();
    ray->d = ray->d/L;
    float A = M_PI * aperture_radius * aperture_radius;
    float Z = lenses.back().thickness;
    float costheta = Dot(ray->d, Vector(0,0,1));
    float costhetapow2 = costheta * costheta;
    float costhetapow4 = costhetapow2 * costhetapow2;

    float E = A/(Z * Z) * costhetapow4;

    // find intersection with first (last in queue) lens
                          
    float offset = 0;
    float n1 = 1.0, n2 = 1.0;
    for (int i = lenses.size() - 1; i >= 0; i--) {
        const LensInterface & lens = lenses[i];
        offset += lens.thickness;
        if (lens.radius == 0) {
            // intersect ray with offset plane
            float t = (offset - ray->o.z)/ray->d.z;
            Point phit = (*ray)(t);
            // see if the aperture stops the ray
            switch (bokeh_type) {
                case 1:
                    {
                        // convert to polar coordinates
                        float r = sqrtf(phit.x * phit.x + phit.y * phit.y);
                        float theta = atan(phit.y/phit.x);
                        // get it into the top slice
                        while (theta > M_PI/bokeh_n_blades) theta -= 2*M_PI/bokeh_n_blades; 
                        while (theta < -M_PI/bokeh_n_blades) theta += 2*M_PI/bokeh_n_blades; 
                        float y = r * cos(theta);
                        // cutoff based on y line
                        if (y > lens.aperture/2 * cosf(M_PI/bokeh_n_blades)) {
                            return 0.0f;
                        }
                    }
                    break;
                case 2: // half_circle
                    if (phit.y < 0 || phit.x * phit.x + phit.y * phit.y > (lens.aperture/2)*(lens.aperture/2)) {
                        return 0.0f;
                    }
                    break;
                case 3: // heart
                    {
                        // convert to polar coordinates
                        float r = sqrtf(phit.x * phit.x + phit.y * phit.y);
                        r /= lens.aperture/2;
                        r *= 4;
                        float theta = atan2(phit.y, phit.x);

                        if (r > 2 - 2*sinf(theta) + sinf(theta) * sqrtf(ABS(cosf(theta))) / (sinf(theta) + 1.4)) {
                             return 0.0f;
                        }
                        // heart equation taken from wolfram alpha
                    }
                    break;
                default:
                    if (phit.x * phit.x + phit.y * phit.y > (lens.aperture/2)*(lens.aperture/2)) {
                        return 0.0f;
                    }
            }
            // otherwise continue
            continue;
        }
        n1 = n2;
        n2 = i == 0 ? 1.0f : lenses[i-1].index_of_refraction;
        if (n2 == 0) {
            n2 = 1;
        }

        if (RayLensRefraction(ray, lens, offset + lens.radius, n1, n2) == false) {
            return 0.0f;
        }
    }
    
    // xaxis flip hack:
    ray->o.x *= -1;
    ray->d.x *= -1;
    // put the front lens at 0,0
    ray->o.z -= offset;

    CameraToWorld(*ray, ray);
    ray->d = Normalize(ray->d);

    // GenerateRay() should return the weight of the generated ray
    return E;
}

