#include "random.h"
#include "math.h"
#include <cmath>
#include "camera.h"

Random::Random(unsigned int seed) : seed(seed), generator(seed) {}

/**********************
 * Generates a random number from 0 to 1
 * @return number from (0,1)
 * ********************/
double Random::GenerateUniformFloat() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(this->generator);
}

/****************************
 * Generates a random point on a disc of radius 1
 * Sqrt is used to make more uniform distribution on the disc otherwise the randomly generated numbers will be more concentrated to center 
 * @return random point from a disc of radius 1 
 * ****************************/
Vec3 Random::GenerateUniformPointDisc() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    double r = std::sqrt(distribution(generator));
    double theta = distribution(generator) * 2 * M_PI;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    return Vec3(x, y, 0);
}


/*********************
 * Generates a random point on a sphere of radius 1
 * @return point on sphere of radius 1
 *******************/
Vec3 Random::GenerateUniformPointSphere() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float u = distribution(generator);
    float v = distribution(generator);

    float theta = u * 2 * M_PI;                         // Using Spherical polar coordinates 
    float cos_phi = 2 * v - 1;
    float sin_phi = std::sqrt(1 - cos_phi * cos_phi);

    Vec3 ret;

    ret.x = cos_phi;
    ret.y = std::sin(theta) * sin_phi;
    ret.z = std::cos(theta) * sin_phi;
    return ret;
}

Vec3 Random::GenerateCosineWeightedUniformPointSphere(const Vec3 & n){
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float u = distribution(generator);
    float v = distribution(generator);

    float theta = 2 * M_PI * u;
    float r = std::sqrt(v);

   
    Vec3 ret;
    ret.x = r * std::cos(theta);
    ret.y = r * std::sin(theta);
    ret.z = std::sqrt(1.0f - v);

    
    Vec3 nn = n.normalized();  
    float c1 = 1 / (1 + nn.z);
    float c2 = nn.y * c1;
    float c3 = nn.x * c1;

    Mat3 m;
    m[0][0] = 1 - nn.x * c3;
    m[0][1] = -nn.x * c2;
    m[0][2] = nn.x;
    m[1][0] = -nn.x * c2;
    m[1][1] = 1 - nn.y * c2;
    m[1][2] = nn.y;
    m[2][0] = -nn.x;
    m[2][1] = -nn.y;
    m[2][2] = nn.z;

    return m * ret;
}

/**
 * @brief Generates an importance-sampled half-vector using the GGX distribution.
 * 
 * This function samples a half-vector `h` from the GGX microfacet distribution. 
 * The distribution is used in physically-based rendering for sampling specular reflection 
 * based on roughness and the angle between the surface normal and the incoming light.
 * 
 * @param n The surface normal, which is used to align the sampled half-vector.
 * @param roughness The roughness of the material, affecting the width of the GGX distribution.
 * 
 * @return Vec3 The generated half-vector in world space, normalized.
 */

Vec3 Random::GenerateImportanceSampleGGX(const Vec3& n, float roughness) {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float u1 = distribution(generator);
    float u2 = distribution(generator);

    float a2 = roughness * roughness;
    float phi = 2.0f * M_PI * u1;
    float cos_theta = sqrt((1 - u2) / (1 + (a2 - 1) * u2));
    float sin_theta = sqrt(1 - cos_theta * cos_theta);

    Vec3 ret = Vec3(sin_theta * cosf(phi), sin_theta * sinf(phi), cos_theta);

    Vec3 nn = n.normalized();  
    float c1 = 1 / (1 + nn.z);
    float c2 = nn.y * c1;
    float c3 = nn.x * c1;

    Mat3 m;
    m[0][0] = 1 - nn.x * c3;
    m[0][1] = -nn.x * c2;
    m[0][2] = nn.x;
    m[1][0] = -nn.x * c2;
    m[1][1] = 1 - nn.y * c2;
    m[1][2] = nn.y;
    m[2][0] = -nn.x;
    m[2][1] = -nn.y;
    m[2][2] = nn.z;

    return m * ret;
}



/*****************************
 *  Generates a random point on a hemisphere around a normal vector
 * @par n ,A 3D vector representing a normal passing through hemisphere 
 * @return random point on the hemisphere 
 *  */
Vec3 Random::GenerateUniformPointHemisphere(const Vec3 &n) {
    Vec3 ret = GenerateUniformPointSphere();
    float factor = dot(ret, n) > 0 ? 1 : -1;

    return ret * factor;
}
