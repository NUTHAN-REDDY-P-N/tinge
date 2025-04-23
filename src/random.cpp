#include "random.h"
#include "math.h"
#include <cmath>
#include "camera.h"


Random::Random(unsigned int seed) : seed(seed), generator(seed) {}

double Random::GenerateUniformFloat() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(this->generator);
}

Vec3 Random::GenerateUniformPointDisc() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    double r = std::sqrt(distribution(generator));
    double theta = distribution(generator) * 2 * M_PI;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    return Vec3(x, y, 0);
}


Vec3 Random::GenerateUniformPointSphere() {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float u = distribution(generator);
    float v = distribution(generator);

    float theta = u * 2 * M_PI;
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

Ray Random::sample_wi(const Ray &wo, const Vec3 &at, const Vec3 &n, Random &random_gen, float refractive_index, float roughness) {
   
    auto clamp = [](float value, float min, float max) {
        return (value < min) ? min : (value > max) ? max : value;
    };

 
    float cos_theta_o = clamp(dot(wo.direction, n), -1.0f, 1.0f);

   
    float r0 = pow((1 - refractive_index) / (1 + refractive_index), 2);
    float F = r0 + (1 - r0) * pow(1 - fabs(cos_theta_o), 5);

    if (random_gen.GenerateUniformFloat() < F) {
    
        float phi = 2.0f * M_PI * random_gen.GenerateUniformFloat();
        float r = random_gen.GenerateUniformFloat();
        float a2 = roughness * roughness;
        float cos_theta = sqrtf((1 - r) / (1 + (a2 - 1) * r));
        float sin_theta = sqrtf(1 - cos_theta * cos_theta);

      
        Vec3 h_local(sin_theta * cosf(phi), sin_theta * sinf(phi), cos_theta);

       
        Vec3 nn = n.normalized();
        float c1 = 1.0f / (1.0f + nn.z);
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

        Vec3 h = m * h_local;
        Ray reflected_dir = reflect(wo, h,n);

      
        return reflected_dir; 
    } else {
     
        Vec3 dir = random_gen.GenerateCosineWeightedUniformPointSphere(n);

        return Ray(at, normalize(dir)); 
    }
}


// Generate a random point on a hemisphere around a normal vector
Vec3 Random::GenerateUniformPointHemisphere(const Vec3 &n) {
    Vec3 ret = GenerateUniformPointSphere();
    float factor = dot(ret, n) > 0 ? 1 : -1;

    return ret * factor;
}
